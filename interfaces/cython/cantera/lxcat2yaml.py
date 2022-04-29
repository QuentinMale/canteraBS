#!/usr/bin/env python3
# encoding: utf-8

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

"""
lxcat2yaml.py: Convert the LXCat integral cross-section data in XML format to YAML format

Usage:
    lxcat2yaml [--input=<filename>]
               [--database=<database name>]
               [--mech=<filename>]
               [--phase=<phase name>]
               [--insert]
               [--output=<filename>]

Example:
    lxcat2yaml --input=mycs.xml --database=itikawa --mech=oxygen-plasma.yaml
               --phase=isotropic-electron-energy-plasma --insert
               --output=oxygen-itikawa-plasma.yaml

If the output file name is not given, an output file with the same name as the
input file, with the extension changed to '.yaml'.
"""

from pathlib import Path
import argparse
import xml.etree.ElementTree as etree
from typing import Union, Optional
import sys
try:
    import ruamel_yaml as yaml
except ImportError:
    from ruamel import yaml
try:
    import cantera as ct
except ImportError:
    print("Could not load Cantera, so the mechanism file cannot be used.")

BlockMap = yaml.comments.CommentedMap

# A class of YAML data for collision of a target species
class Process:
    def __init__(self, equation, energy_levels, cross_section):
        self.equation = equation
        self.energy_levels = energy_levels
        self.cross_section = cross_section

    @classmethod
    def to_yaml(cls, representer, node):
        out = BlockMap([('equation', node.equation),
                        ('type', 'electron-collision-plasma'),
                        ('energy-levels', node.energy_levels),
                        ('cross-section', node.cross_section),
                        ])
        return representer.represent_dict(out)

# Define YAML emitter
emitter = yaml.YAML()
emitter.register_class(Process)

# Return indices of a child name
def get_children(parent, child_name):
    return [child for child in parent if child.tag.find(child_name) != -1]

def FlowList(*args, **kwargs):
    """A YAML sequence that flows onto one line."""
    lst = yaml.comments.CommentedSeq(*args, **kwargs)
    lst.fa.set_flow_style()
    return lst

class IncorrectXMLNode(LookupError):
    def __init__(self, message: str = "", node: Optional[etree.Element] = None):
        """Error raised when a required node is incorrect in the XML tree.

        :param message:
            The error message to be displayed to the user.
        :param node:
            The XML node from which the requested node is incorrect.
        """
        if node is not None:
            # node str
            node_str = etree.tostring(node, encoding="unicode")

            # Print the XML node
            if message:
                message += "\n" + node_str
            else:
                message = "\n" + node_str

        super().__init__(message)

def convert(
        inpfile: Union[str, Path] = None,
        database: str = None,
        mechfile: Union[str, Path] = None,
        phase: str = None,
        insert: bool = True,
        outfile: Union[str, Path] = None,
    ) -> None:
    """Convert an LXCat XML file to a YAML file.

    :param inpfile:
        The input LXCat file name. Exclusive with ``text``, only one of the two can be
        specified.
    :param database:
        The name of the database. E.g. itikawa.
    :param mechfile:
        The reaction mechanism file. This option requires using the Cantera library.
    :param phase:
        The phase name of the mechanism file. This option requires using mechfile.
    :param insert:
        The flag of whether to insert the collision reactions or not.
    :param outfile:
        The output YAML file name.

    All files are assumed to be relative to the current working directory of the Python
    process running this script.
    """
    if inpfile is not None:
        inpfile = Path(inpfile)
        lxcat_text = inpfile.read_text().lstrip()
        if outfile is None:
            outfile = inpfile.with_suffix(".yaml")
    else:
        raise ValueError("'inpfile' must be specified")

    if insert == True and mechfile is None:
        raise ValueError("'mech' must be specified if 'insert' is used")

    gas = None
    if mechfile is not None:
        mechfile = Path(mechfile)
        try:
            gas = ct.Solution(mechfile, phase, transport_model=None)
        except ImportError:
            print("Could not load Cantera, so the mechanism file cannot be used.")
            sys.exit(0)
    elif phase is not None:
        raise ValueError("'mech' must be specified if 'phase' is used.")

    xml_tree = etree.fromstring(lxcat_text)

    # If insert key word is used, create a process list,
    # and append all processes together
    process_list = None
    if not insert:
        process_list = []

    for database_node in xml_tree:
        if database_node.attrib["id"] != database:
            continue

        # Get groups node
        groups_node = get_children(database_node, "groups")[0]

        for group in groups_node:
            for process in get_children(group, "processes")[0]:
                registerProcess(process, process_list, gas)

    if not insert:
        # Put process list in collision node
        collision_node = {"collisions": process_list}
        with Path(outfile).open("w") as output_file:
            emitter.dump(collision_node, output_file)
    else:
        gas.write_yaml(outfile)

def registerProcess(process, process_list, gas):
    # Parse the threshold
    threshold = 0.0
    parameters_node = get_children(process, "parameters")[0]
    if len(get_children(parameters_node, "parameter")) == 1:
        parameter = get_children(parameters_node, "parameter")[0]
        if parameter.attrib["name"] == 'E':
            threshold = float(parameter.text)

    # Parse the equation
    product_array=[]

    for product_node in get_children(process, "products")[0]:
        if product_node.tag.find("electron") != -1:
            product_array.append(gas.electron_species_name)
        if product_node.tag.find("molecule") != -1:
            product_name = product_node.text
            if "state" in product_node.attrib:
                state = product_node.attrib["state"].replace(" ","-")
                # State is appended in a parenthesis
                product_name += f"({state})"
            if "charge" in product_node.attrib:
                charge = int(product_node.attrib["charge"])
                if charge > 0:
                    product_name += charge*"+"
                else:
                    product_name += -charge*"-"

            # Filter the collision based on the existed species in the mechanism file
            if gas is None or product_name in gas.species_names:
                product_array.append(product_name)
            else:
                return

    for reactant_node in get_children(process, "reactants")[0]:
        if reactant_node.tag.find("molecule") != -1:
            reactant = reactant_node.text

    products = " + ".join(product_array)
    equation = f"{reactant} + {gas.electron_species_name} => {products}"

    # Parse the cross-section data
    data_x_node = process.find("data_x")
    if data_x_node is None:
        raise IncorrectXMLNode("The 'process' node requires the 'data_x' node.", process)

    data_y_node = process.find("data_y")
    if data_y_node is None:
        raise IncorrectXMLNode("The 'process' node requires the 'data_y' node.", process)

    energy_levels = FlowList(map(float, data_x_node.text.split()))
    cross_section = FlowList(map(float, data_y_node.text.split()))

    # Edit energy levels and cross section
    if len(energy_levels) != len(cross_section):
        raise IncorrectXMLNode("Energy levels (data_x) and cross section "
                                "(data_y) must have the same length.", process)

    if energy_levels[0] > threshold:
        # Use FlowList again to ensure correct YAML format
        energy_levels = FlowList([threshold, *energy_levels])
        cross_section = FlowList([0.0, *cross_section])
    else:
        cross_section[0] = 0.0

    # If insert mode is on, add the process as a reaction to the gas object.
    if gas is not None:
        R = ct.Reaction(
            equation=equation,
            rate=ct.ElectronCollisionPlasmaRate(energy_levels=energy_levels,
                                                cross_section=cross_section))
        gas.add_reaction(R)

    # If insert mode is off, process_list is used to store the data.
    if process_list is not None:
        process_list.append(Process(equation=equation,
                                    energy_levels=energy_levels,
                                    cross_section=cross_section))

def main():
    """Parse command line arguments and pass them to `convert`."""
    parser = argparse.ArgumentParser(
        description="Convert the LXCat integral cross-section data in XML format (LXCATML) to YAML format",
        epilog=(
            "The 'mech-file' argument is optional. It is used to filter the cross-section "
            "data. If it is given, the 'insert' argument can be used. If 'insert' is used, "
            "the output file will contain the mechanism file with the cross-section data "
            "inserted. The 'output' argument is optional. If it is not given, an output "
            "file with the same name as the input file is used, with the extension "
            "changed to '.yaml'."
        ),
    )
    parser.add_argument("--input", required=True, type=str, help="The input LXCatML filename. Must be specified.")
    parser.add_argument("--database", required=True, type=str, help="The name of the database. Must be specified.")
    parser.add_argument("--mech", nargs="?", help="The mechanism filename. Optional.")
    parser.add_argument("--phase", nargs="?", help="The phase name of the mechanism file. Optional.")
    parser.add_argument("--insert", action='store_true', help="The flag for whether the collisions are inserted into the mechanism file. Optional.")
    parser.add_argument("--output", nargs="?", help="The output YAML filename. Optional.")

    args = parser.parse_args()
    input_file = Path(args.input)
    output_file = args.output or input_file.with_suffix(".yaml")

    convert(input_file, args.database, args.mech, args.phase, args.insert, output_file)

    if args.insert:
        # Test mechanism can be loaded back into Cantera
        try:
            print("Validating mechanism...", end="")
            ct.Solution(output_file, args.phase, transport_model=None)
            print("PASSED.")
        except ImportError:
            print("Could not load Cantera, so the converted file could not be validated.")
            sys.exit(0)
        except RuntimeError as e:
            print("FAILED.")
            print(e)
            sys.exit(1)

if __name__ == "__main__":
    main()
