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

dt_max = 1.e-10
t_end = 100 * dt_max
states = ct.SolutionArray(gas, extra=['t'])

print('{:10s} {:10s} {:10s} {:14s}'.format(
    't [s]', 'T [K]', 'P [Pa]', 'u [J/kg]'))
while sim.time < t_end:
    sim.advance(sim.time + dt_max)
    states.append(r.thermo.state, t=sim.time*1e3)
    print('{:10.3e} {:10.3f} {:10.3f} {:14.6f}'.format(
            sim.time, r.T, r.thermo.P, r.thermo.u))
