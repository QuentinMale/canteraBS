import os
from bolos import parser, grid, solver

def eedf(gas_species, gas_molefraction, gas_temperature):
    #plugin = __import__(bolos)
    gr = grid.LinearGrid(0, 60., 200)

    boltzmann = solver.BoltzmannSolver(gr)

    dir_path = os.path.dirname(os.path.realpath(__file__))
    with open(dir_path + '/data/bolsig_lxcat.dat') as fp:
        processes = parser.parse(fp)
    boltzmann.load_collisions(processes)

    for i, name in enumerate(gas_species):
        list_name = list(name)
        for j in list_name:
            if j == '+' or j == '-':
                k = list_name.index(j)
                if list_name[k-1] != '^':
                    list_name.insert(list_name.index(j), '^')
        list_name = ''.join(list_name)
        gas_species[i] = list_name

    gas_list = dict(zip(gas_species, gas_molefraction))

    for key in gas_list:
        boltzmann.target[key].density = gas_list[key]

    boltzmann.kT = gas_temperature * solver.KB / solver.ELECTRONVOLT
    # not very important for now
    boltzmann.EN = 500 * solver.TOWNSEND

    boltzmann.init()
    fMaxwell = boltzmann.maxwell(2.0)
    f = boltzmann.converge(fMaxwell, maxn=100, rtol=1e-5)

    # Calculate the mean energy according to the first EEDF
    mean_energy = boltzmann.mean_energy(f)

    # Set a new grid extending up to 15 times the mean energy.
    # Now we use a quadritic grid instead of a linear one.
    newgrid = grid.QuadraticGrid(0, 15 * mean_energy, 200)

    # Set the new grid and update the internal
    boltzmann.grid = newgrid
    boltzmann.init()

    # Calculate an EEDF in the new grid by interpolating the old one
    finterp = boltzmann.grid.interpolate(f, gr)

    # Iterate until we have a new solution
    f1 = boltzmann.converge(finterp, maxn=200, rtol=1e-5)

    #mun = boltzmann.mobility(f1)
    #diffn = boltzmann.diffusion(f1)

    # Obtain the reaction rate for impact ionization of molecular nitrogen.
    #k = boltzmann.rate(f1, "N2 -> N2^+")