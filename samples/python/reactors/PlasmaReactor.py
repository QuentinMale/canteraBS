"""

Requires: cantera >= XX.

.. tags:: Python, plasma
"""


import cantera as ct


# solution
gas = ct.Solution('air-plasma_Phelps.yaml')
gas.TPX = 300., 101325., 'N2:0.79, O2:0.21, N2+:1E-10, Electron:1E-10'

# plasma reactor
r = ct.PlasmaReactor(gas)

sim = ct.ReactorNet([r])
sim.verbose = True

