"""

Requires: cantera >= XX.

.. tags:: Python, plasma
"""


import cantera as ct
import numpy as np


# solution
gas = ct.Solution('air-plasma_Phelps.yaml')
gas.TPX = 300., 101325., 'N2:0.79, O2:0.21, N2+:1E-9, Electron:1E-9'
gas.EN = 200.0 * 1e-21 # Reduced electric field [V.m^2]
gas.update_EEDF()

# plasma reactor
r = ct.PlasmaReactor(gas)
r.dis_vol = 5e-3*np.pi*(1e-3)**2/4
print(r.dis_vol)
print(r.dis_power)

sim = ct.ReactorNet([r])
sim.verbose = True

