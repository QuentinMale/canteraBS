"""
EEDF calculation
==============
Compute EEDF with two term approximation solver at constant E/N.

Requires: cantera >= XX.

.. tags:: Python, plasma
"""


import cantera as ct

gas = ct.Solution('air-plasma_Phelps.yaml')
gas.TPX = 300., 101325., 'N2:0.79,O2:0.21,N2+:1E-10,Electron:1E-10'
gas.EN = 200.0 * 1e-21 # Reduced electric field [V.m^2]
gas.update_EEDF()

