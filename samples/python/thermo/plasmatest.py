"""
EEDF calculation
==============
Compute EEDF with two term approximation solver at constant E/N.
Compare with results from BOLOS.

Requires: cantera >= XX.

.. tags:: Python, plasma
"""


import matplotlib.pyplot as plt
import cantera as ct

gas = ct.Solution('air-plasma_Phelps.yaml')
gas.TPX = 300., 101325., 'N2:0.79, O2:0.21, N2+:1E-10, Electron:1E-10'
gas.EN = 200.0 * 1e-21 # Reduced electric field [V.m^2]
gas.update_EEDF()

grid = gas.electron_energy_levels
eedf = gas.electron_energy_distribution

# results from BOLOS
cgrid = [6.000e-02, 1.282e+00, 2.504e+00, 3.726e+00, 4.948e+00, 6.170e+00, 7.392e+00, 8.614e+00, 9.836e+00, 1.106e+01, 1.228e+01, 1.350e+01, 1.472e+01, 1.595e+01, 1.717e+01, 1.839e+01, 1.961e+01, 2.083e+01, 2.206e+01, 2.328e+01, 2.450e+01, 2.572e+01, 2.694e+01, 2.817e+01, 2.939e+01, 3.061e+01, 3.183e+01, 3.306e+01, 3.428e+01, 3.550e+01, 3.672e+01, 3.794e+01, 3.917e+01, 4.039e+01, 4.161e+01, 4.283e+01, 4.405e+01, 4.528e+01, 4.650e+01, 4.772e+01, 4.894e+01, 5.016e+01, 5.139e+01, 5.261e+01, 5.383e+01, 5.505e+01, 5.627e+01, 5.750e+01, 5.872e+01, 5.994e+01]
cf0 = [1.445e-01, 1.332e-01, 6.907e-02, 4.753e-02, 3.908e-02, 3.160e-02, 2.481e-02, 1.874e-02, 1.350e-02, 9.381e-03, 6.258e-03, 4.004e-03, 2.488e-03, 1.548e-03, 9.596e-04, 5.917e-04, 3.622e-04, 2.198e-04, 1.323e-04, 7.904e-05, 4.695e-05, 2.773e-05, 1.627e-05, 9.491e-06, 5.504e-06, 3.182e-06, 1.831e-06, 1.048e-06, 5.975e-07, 3.390e-07, 1.916e-07, 1.078e-07, 6.039e-08, 3.371e-08, 1.876e-08, 1.041e-08, 5.759e-09, 3.177e-09, 1.750e-09, 9.632e-10, 5.297e-10, 2.910e-10, 1.598e-10, 8.780e-11, 4.831e-11, 2.672e-11, 1.499e-11, 8.773e-12, 5.762e-12, 4.822e-12]

fig, ax = plt.subplots()

ax.plot(grid, eedf, label='CANTERA')
ax.plot(cgrid, cf0, label='BOLOS')

ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlim(1e-2, 1e2)
ax.set_ylim(1e-10, 1e4)

ax.legend()

plt.show()
