#!/usr/bin/env python
# encoding: utf-8

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

"""
lxcat2yaml.py: Convert lxcat-format cross-section data to YAML files

Usage:
    lxcat2yaml [--input=<filename>]
               [--name=<name>]

Example:
    lxcat2yaml --input=lxcat.txt

If the output file name is not given, an output file with the same name as the
input file, with the extension changed to '.yaml'.
"""

from collections import defaultdict, OrderedDict
import logging
import os.path
import sys
import numpy as np
import re
import itertools
import getopt
import textwrap
from email.utils import formatdate

try:
    import ruamel_yaml as yaml
except ImportError:
    from ruamel import yaml

BlockMap = yaml.comments.CommentedMap

class Parser:
    def __init__(self):
        self.processed_units = False
        self.energy_units = 'cal/mol'  # for the current REACTIONS section
        self.output_energy_units = 'cal/mol'  # for the output file
        self.quantity_units = 'mol'  # for the current REACTIONS section
        self.output_quantity_units = 'mol'  # for the output file
        self.motz_wise = None
        self.warning_as_error = True

        self.elements = []
        self.element_weights = {}  # for custom elements only
        self.species_list = []  # bulk species only
        self.species_dict = {}  # bulk and surface species
        self.surfaces = []
        self.reactions = []
        self.header_lines = []
        self.extra = {}  # for extra entries
        self.files = []  # input file names

    def warn(self, message):
        if self.warning_as_error:
            raise InputError(message)
        else:
            logging.warning(message)

    def write_yaml(self, name='gas', out_name='mech.yaml'):
        emitter = yaml.YAML()
        emitter.width = 70

        emitter.register_class(Species)
        emitter.register_class(Nasa7)
        emitter.register_class(Nasa9)
        emitter.register_class(TransportData)
        emitter.register_class(Reaction)

        with open(out_name, 'w') as dest:
            have_transport = True
            for s in self.species_list:
                if not s.transport:
                    have_transport = False

            surface_names = []
            n_reacting_phases = 0
            if self.reactions:
                n_reacting_phases += 1
            for surf in self.surfaces:
                surface_names.append(surf.name)
                if surf.reactions:
                    n_reacting_phases += 1

            # Write header lines
            desc = '\n'.join(line.rstrip() for line in self.header_lines)
            desc = desc.strip('\n')
            desc = textwrap.dedent(desc)
            if desc.strip():
                emitter.dump({'description': yaml.scalarstring.PreservedScalarString(desc)}, dest)

            # Additional information regarding conversion
            files = [os.path.basename(f) for f in self.files]
            metadata = BlockMap([
                ('generator', 'ck2yaml'),
                ('input-files', FlowList(files)),
                ('cantera-version', '2.5.0a4'),
                ('date', formatdate(localtime=True)),
            ])
            if desc.strip():
                metadata.yaml_set_comment_before_after_key('generator', before='\n')
            emitter.dump(metadata, dest)

            # Write extra entries
            if self.extra:
                extra = BlockMap(self.extra)
                key = list(self.extra.keys())[0]
                extra.yaml_set_comment_before_after_key(key, before='\n')
                emitter.dump(extra, dest)

            units = FlowMap([('length', 'cm'), ('time', 's')])
            units['quantity'] = self.output_quantity_units
            units['activation-energy'] = self.output_energy_units
            units_map = BlockMap([('units', units)])
            units_map.yaml_set_comment_before_after_key('units', before='\n')
            emitter.dump(units_map, dest)

            phases = []
            reactions = []
            if name is not None:
                phase = BlockMap()
                phase['name'] = name
                phase['thermo'] = 'ideal-gas'
                phase['elements'] = FlowList(self.elements)
                phase['species'] = FlowList(S.label for S in self.species_list)
                if self.reactions:
                    phase['kinetics'] = 'gas'
                    if n_reacting_phases == 1:
                        reactions.append(('reactions', self.reactions))
                    else:
                        rname = '{}-reactions'.format(name)
                        phase['reactions'] = [rname]
                        reactions.append((rname, self.reactions))
                if have_transport:
                    phase['transport'] = 'mixture-averaged'
                phase['state'] = FlowMap([('T', 300.0), ('P', '1 atm')])
                phases.append(phase)

            for surf in self.surfaces:
                # Write definitions for surface phases
                phase = BlockMap()
                phase['name'] = surf.name
                phase['thermo'] = 'ideal-surface'
                phase['elements'] = FlowList(self.elements)
                phase['species'] = FlowList(S.label for S in surf.species_list)
                phase['site-density'] = surf.site_density
                if self.motz_wise is not None:
                    phase['Motz-Wise'] = self.motz_wise
                if surf.reactions:
                    phase['kinetics'] = 'surface'
                    if n_reacting_phases == 1:
                        reactions.append(('reactions', surf.reactions))
                    else:
                        rname = '{}-reactions'.format(surf.name)
                        phase['reactions'] = [rname]
                        reactions.append((rname, surf.reactions))
                phase['state'] = FlowMap([('T', 300.0), ('P', '1 atm')])
                phases.append(phase)

            if phases:
                phases_map = BlockMap([('phases', phases)])
                phases_map.yaml_set_comment_before_after_key('phases', before='\n')
                emitter.dump(phases_map, dest)

            # Write data on custom elements
            if self.element_weights:
                elements = []
                for name, weight in sorted(self.element_weights.items()):
                    E = BlockMap([('symbol', name), ('atomic-weight', weight)])
                    elements.append(E)
                elementsMap = BlockMap([('elements', elements)])
                elementsMap.yaml_set_comment_before_after_key('elements', before='\n')
                emitter.dump(elementsMap, dest)

            # Write the individual species data
            all_species = list(self.species_list)
            for species in all_species:
                if species.composition is None:
                    raise InputError('No thermo data found for '
                                     'species {!r}'.format(species.label))

            for surf in self.surfaces:
                all_species.extend(surf.species_list)
            speciesMap = BlockMap([('species', all_species)])
            speciesMap.yaml_set_comment_before_after_key('species', before='\n')
            emitter.dump(speciesMap, dest)

            # Write the reactions section(s)
            for label, R in reactions:
                reactionsMap = BlockMap([(label, R)])
                reactionsMap.yaml_set_comment_before_after_key(label, before='\n')
                emitter.dump(reactionsMap, dest)

        # Names of surface phases need to be returned so they can be imported as
        # part of mechanism validation
        return surface_names

def main(argv):

    longOptions = ['input=', 'thermo=', 'transport=', 'surface=', 'name=',
                   'extra=', 'output=', 'permissive', 'help', 'debug', 'quiet',
                   'no-validate', 'id=']

    try:
        optlist, args = getopt.getopt(argv, 'dh', longOptions)
        options = dict()
        for o,a in optlist:
            options[o] = a

        if args:
            raise getopt.GetoptError('Unexpected command line option: ' +
                                     repr(' '.join(args)))

    except getopt.GetoptError as e:
        print('ck2yaml.py: Error parsing arguments:')
        print(e)
        print('Run "ck2yaml.py --help" to see usage help.')
        sys.exit(1)

    if not options or '-h' in options or '--help' in options:
        print(__doc__)
        sys.exit(0)

    input_file = options.get('--input')

    if '--output' in options:
        out_name = options['--output']
        if not out_name.endswith('.yaml') and not out_name.endswith('.yml'):
            out_name += '.yaml'
    elif input_file:
        out_name = os.path.splitext(input_file)[0] + '.yaml'
    else:
        out_name = os.path.splitext(thermo_file)[0] + '.yaml'

    try:
        import cantera as ct
    except ImportError:
        print('WARNING: Unable to import Cantera Python module. Output '
              'mechanism has not been validated')
        sys.exit(0)


def script_entry_point():
    main(sys.argv[1:])

if __name__ == '__main__':
    main(sys.argv[1:])
