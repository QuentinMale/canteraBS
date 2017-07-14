"""
A freely-propagating, premixed methane-air flat flame with charged species.
"""

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

# Simulation parameters
p = ct.one_atm  # pressure [Pa]
Tin = 300.0  # unburned gas temperature [K]
reactants = 'CH4:1, O2:2, N2:7.52'  # premixed gas composition
width = 0.05  # m
loglevel = 1  # amount of diagnostic output (0 to 8)

# IdealGasMix object used to compute mixture properties, set to the state of the
# upstream fuel-air mixture
gas = ct.Solution('ch4_plasma.cti')
gas.TPX = Tin, p, reactants

# Set up flame object
f = ct.IonFlame(gas, width=width)
f.set_refine_criteria(ratio=3, slope=0.06, curve=0.12)
#f.show_solution() #this will cause error, because it gives a bad initial guess

def plotGenerator():
    mobility_e = []
    diffCoeff_e = []
    Te = []
    X_E = []
    X_O2n = []
    for n in range(f.flame.n_points):
        f.set_gas_state(n)
        mobility_e.append(f.flame.get_elecMobility(n))
        diffCoeff_e.append(f.flame.get_elecDiffCoeff(n))
        Te.append(f.flame.get_elecTemperature(n))
        X_E.append(f.gas.X[f.gas.species_index('E')])
        X_O2n.append(f.gas.X[f.gas.species_index('O2^-')])

    plt.figure(1)
    plt.plot(f.grid, Te)
    plt.ylabel('Temperature [K]')
    plt.title('electron temperature')
    plt.figure(2)
    plt.plot(f.grid, f.E)
    plt.ylabel('E [V/m]')
    plt.title('electric field')
    plt.figure(3)
    plt.plot(f.grid, mobility_e)
    plt.ylabel('mu [m2/(V-s)]')
    plt.title('electron mobility')
    plt.figure(4)
    plt.plot(f.grid, diffCoeff_e)
    plt.ylabel('D [m2/s]')
    plt.title('electron diffusion coefficient')
    plt.figure(5)
    plt.plot(f.grid, X_E)
    plt.ylabel('mole fraction')
    plt.title('electron concentration')
    plt.figure(6)
    plt.plot(f.grid, X_O2n)
    plt.ylabel('mole fraction')
    plt.title('O2^- concentration')
    plt.show()

gas.set_multiplier(104,0.0)
# for i in range(104,109): #this method will produce bad guess
#     gas.set_multiplier(i,0.0)
#------------------------------------------------------------------
# stage 1: Frozen-Ion Method
#
# In this stage, the diffusion of ions is turned off due to
# the fast diffusion rate of electron without internal 
# electric forces (ambi-polar diffusion effect).
f.ionSolve(loglevel=loglevel, auto=True)
#plotGenerator()

#-------------------------------------------------------------------
# trun on chemi-ionization
psm = np.logspace(-5,0,6)
gas.set_multiplier(104,psm[0])
# for i in range(104,109):
#     gas.set_multiplier(i,1e-5)
# turn on plasma
print('begin to run with zdplaskin')
f.plasma_enabled = True
f.flame.set_transversElecField(2.5e4)
#f.flame.set_elecFrequency(1e10)
f.flame.set_plasmaSourceMultiplier(psm[0])
f.ionSolve(loglevel=loglevel, stage=1, enable_energy=True)
#plotGenerator()
#--------------------------------------------------------------------
# stage 2: Charge-Neutrality Model
#
# The second stage uses charge neutrality model, which assume
# zero charge flux throughout the domain, to calculate drift flux.
# The drift flux is added to the total flux of ions.
f.ionSolve(loglevel=loglevel, stage=2, enable_energy=False)
f.ionSolve(loglevel=loglevel, stage=2, enable_energy=True)
#plotGenerator()
#--------------------------------------------------------------------
# stage 3: Poisson Equation Method
#
# The third stage evaluates drift flux from electric field calculated 
# from Poisson's equation, which is solved together with other equations.
# Poisson's equation is coupled because the total charge densities 
# depends on the species' concentration.
f.flame.enable_elecHeat(True)
for i in psm:
    f.flame.set_plasmaSourceMultiplier(i)
    gas.set_multiplier(104,i)
    f.ionSolve(loglevel=loglevel, stage=3, enable_energy=False)
    f.ionSolve(loglevel=loglevel, stage=3, enable_energy=True)
plotGenerator()

f.save('CH4_adiabatic.xml', 'mix', 'solution with mixture-averaged transport')
#f.show_solution()
print('mixture-averaged flamespeed = {0:7f} m/s'.format(f.u[0]))

# write the velocity, temperature, density, and mole fractions to a CSV file
f.write_csv('CH4_adiabatic.csv', quiet=False)
