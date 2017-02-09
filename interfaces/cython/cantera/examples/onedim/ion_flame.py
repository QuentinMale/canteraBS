"""
A freely-propagating, premixed hydrogen flat flame with multicomponent
transport properties.
"""

import cantera as ct
import numpy as np

# Simulation parameters
p = ct.one_atm  # pressure [Pa]
Tin = 300.0  # unburned gas temperature [K]
reactants = 'CH4:0.216, O2:2'  # premixed gas composition
width = 0.05  # m
loglevel = 1  # amount of diagnostic output (0 to 8)
offset = 5    #the offset for the species

# IdealGasMix object used to compute mixture properties, set to the state of the
# upstream fuel-air mixture
gas = ct.Solution('gri30_ion.xml')
gas.TPX = Tin, p, reactants

# Set up flame object
f = ct.IonFlame(gas, width=width)
f.set_refine_criteria(ratio=3, slope=0.06, curve=0.12)
f.show_solution()
# Set tolerance for ions
f.flame.set_steady_tolerances(**{'HCO+':(1.0e-5, 1.0e-16)})
f.flame.set_transient_tolerances(**{'HCO+':(1.0e-5, 1.0e-18)})
f.flame.set_steady_tolerances(**{'H3O+':(1.0e-5, 1.0e-13)})
f.flame.set_transient_tolerances(**{'H3O+':(1.0e-5, 1.0e-15)})
f.flame.set_steady_tolerances(**{'E':(1.0e-5, 1.0e-16)})
f.flame.set_transient_tolerances(**{'E':(1.0e-5, 1.0e-18)})

# phase one
gas.set_multiplier(325,0.01)
f.solve(loglevel=loglevel, auto=True)

# phase two
f.flame.set_solvingPhase(2)
f.energy_enabled = False
f.velocity_enabled = False
gas.set_multiplier(325,1.0)
# Solve with Prager's ambipolar diffusion model
f.solve(loglevel)
f.energy_enabled = True
f.velocity_enabled = True
f.solve(loglevel)

# phase three
f.flame.set_solvingPhase(3)
f.poisson_enabled = True
# Solve with electric force
f.solve(loglevel)

f.save('CH4_adiabatic.xml', 'mix', 'solution with mixture-averaged transport')
f.show_solution()
print('mixture-averaged flamespeed = {0:7f} m/s'.format(f.u[0]))

# write the velocity, temperature, density, and mole fractions to a CSV file
f.write_csv('CH4_adiabatic.csv', quiet=False)
f.write_electric_csv('CH4_adiabatic_ion.csv', quiet=False)
