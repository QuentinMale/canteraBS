/*!
 * @file flamespeed.cpp
 * C++ demo program to compute flame speeds using GRI-Mech.
 */

#include "cantera/oneD/Sim1D.h"
#include "cantera/oneD/Inlet1D.h"
#include "cantera/oneD/StFlow.h"
#include "cantera/oneD/IonFlow.h"
#include "cantera/IdealGasMix.h"
#include "cantera/transport.h"
#include "cantera/numerics/funcs.h"
#include <fstream>

using namespace std;
using namespace Cantera;
using fmt::print;

int flamespeed(double phi)
{
    try {
        IdealGasMix gas("gri30_ion.xml","gri30_mix");

        doublereal temp = 300.0; // K
        doublereal pressure = 1.0*OneAtm; //atm
        doublereal uin = 0.6; //m/sec

        size_t nsp = gas.nSpecies();
        vector_fp x(nsp, 0.0);

        doublereal C_atoms = 1.0;
        doublereal H_atoms = 4.0;
        doublereal ax = C_atoms + H_atoms / 4.0;
        doublereal fa_stoic = 1.0 / ax;
        x[gas.speciesIndex("CH4")] = 1.0;
        x[gas.speciesIndex("O2")] = 1.0 / phi / fa_stoic;

        gas.setState_TPX(temp, pressure, x.data());
        doublereal rho_in = gas.density();

        vector_fp yin(nsp, 0.0);
        gas.getMassFractions(&yin[0]);

        gas.equilibrate("HP");
        vector_fp yout(nsp, 0.0);
        gas.getMassFractions(&yout[0]);

        doublereal rho_out = gas.density();
        doublereal Tad = gas.temperature();
        print("phi = {}, Tad = {}\n", phi, Tad);

        vector_fp mw = gas.molecularWeights();

        //=============  build each domain ========================


        //-------- step 1: create the flow -------------

        IonFlow flow(&gas);

        // create an initial grid
        int nz = 6;
        doublereal lz = 0.1;
        vector_fp z(nz);
        doublereal dz = lz/((doublereal)(nz-1));
        for (int iz = 0; iz < nz; iz++) {
            z[iz] = ((doublereal)iz)*dz;
        }

        flow.setupGrid(nz, &z[0]);

        // specify the objects to use to compute kinetic rates and
        // transport properties

        unique_ptr<Transport> trmix(newTransportMgr("Mix", &gas));

        flow.setTransport(*trmix);
        flow.setKinetics(gas);
        flow.setPressure(pressure);

        //------- step 2: create the inlet  -----------------------

        Inlet1D inlet;

        inlet.setMoleFractions(x.data());
        doublereal mdot=uin*rho_in;
        inlet.setMdot(mdot);
        inlet.setTemperature(temp);

        //------- step 3: create the outlet  ---------------------

        Outlet1D outlet;

        //=================== create the container and insert the domains =====

        vector<Domain1D*> domains { &inlet, &flow, &outlet };
        Sim1D flame(domains);

        //----------- Supply initial guess----------------------

        vector_fp locs{0.0, 0.3, 0.7, 1.0};
        vector_fp value;

        double uout = inlet.mdot()/rho_out;
        value = {uin, uin, uout, uout};
        flame.setInitialGuess("u",locs,value);
        value = {temp, temp, Tad, Tad};
        flame.setInitialGuess("T",locs,value);

        for (size_t k=0; k<nsp; k++) {
            value = {yin[k], yin[k], yout[k], yout[k]};
            flame.setInitialGuess(gas.speciesName(k),locs,value);
        }

        inlet.setMoleFractions(x.data());
        inlet.setMdot(mdot);
        inlet.setTemperature(temp);

        flame.showSolution();

        int flowdomain = 1;
        double ratio = 10.0;
        double slope = 0.08;
        double curve = 0.1;

        flame.setRefineCriteria(flowdomain,ratio,slope,curve);

        int loglevel=1;

        // Solve freely propagating flame

        // Linearly interpolate to find location where this temperature would
        // exist. The temperature at this location will then be fixed for
        // remainder of calculation.
        flame.setFixedTemperature(0.5 * (temp + Tad));
        flow.solveEnergyEqn();
        bool refine_grid = true;

        // set tolerances for ions
        flow.setSteadyTolerances(1.0e-5,1.0e-16,c_offset_Y + gas.speciesIndex("HCO+"));
        flow.setTransientTolerances(1.0e-5,1.0e-18,c_offset_Y + gas.speciesIndex("HCO+"));
        flow.setSteadyTolerances(1.0e-5,1.0e-13,c_offset_Y + gas.speciesIndex("H3O+"));
        flow.setTransientTolerances(1.0e-5,1.0e-15,c_offset_Y + gas.speciesIndex("H3O+"));
        flow.setSteadyTolerances(1.0e-5,1.0e-16,c_offset_Y + gas.speciesIndex("E"));
        flow.setTransientTolerances(1.0e-5,1.0e-18,c_offset_Y + gas.speciesIndex("E"));

        // phase one
        // turn down the reaction for production of electron
        for (size_t i = 325; i < 326; i++) {
            gas.setMultiplier(i,0.01);
        }
        
        flame.solve(loglevel,refine_grid);
        double flameSpeed_mix = flame.value(flowdomain,flow.componentIndex("u"),0);
        print("Flame speed with mixture-averaged transport: {} m/s\n",
              flameSpeed_mix);       

        // phase two
        flow.setSolvingPhase(2);
        flow.fixTemperature();
        flow.fixVelocity();

        // use full reaction rate
        for (size_t i = 325; i < 326; i++) {
            gas.setMultiplier(i,1.0);
        }
        flame.solve(loglevel, refine_grid);

        // turn on the energy equation
        flow.solveEnergyEqn();
        flow.solveVelocity();
        flame.solve(loglevel, refine_grid);

        // phase three
        flow.setSolvingPhase(3);
        flow.solvePoissonEqn();
        flame.solve(loglevel, refine_grid);

        // post simulation       
        vector_fp rho_vec,HCOxvec,H3Oxvec,Evec;

        for (size_t n = 0; n < flow.nPoints(); n++) {
            // save the concentrations for charged species in number density
            double rho = flow.density(n);
            size_t k = gas.speciesIndex("HCO+");
            rho_vec.push_back(rho);
            HCOxvec.push_back(2e-4*Avogadro*rho*flame.value(flowdomain,k+c_offset_Y,n)/mw[k]);
            k = gas.speciesIndex("H3O+");
            H3Oxvec.push_back(1e-6*Avogadro*rho*flame.value(flowdomain,k+c_offset_Y,n)/mw[k]);
            k = gas.speciesIndex("E");
            Evec.push_back(1e-6*Avogadro*rho*flame.value(flowdomain,k+c_offset_Y,n)/mw[k]);
        } 

        print("\n{:9s}\t{:8s}\t{:8s}\t{:8s}\t{:5s}\n",
              "z (m)", "T (K)", "U (m/s)", "eP[V]", "Y(E)");

        vector_fp Tvec,ePvec,eFvec,Uvec;

        for (size_t n = 0; n < flow.nPoints(); n++) {
            size_t k = gas.speciesIndex("E");
            Tvec.push_back(flame.value(flowdomain,flow.componentIndex("T"),n));
            ePvec.push_back(flame.value(flowdomain,flow.componentIndex("ePotential"),n));
            Uvec.push_back(flame.value(flowdomain,flow.componentIndex("u"),n));

            print("{:9.6f}\t{:8.3f}\t{:8.3f}\t{:8.3f}\t{:11.3e}\n",
                  flow.grid(n), Tvec[n], Uvec[n], ePvec[n],flame.value(flowdomain,k+c_offset_Y,n));
        }

        for (size_t n = 0; n < flow.nPoints()-1; n++) {
            // save electrical field strength in V/cm
            eFvec.push_back(1e-2*(ePvec[n]-ePvec[n+1])/(flow.grid(n+1)-flow.grid(n)));
        }      
        eFvec.push_back(eFvec[flow.nPoints()-2]);

        print("\nAdiabatic flame temperature from equilibrium is: {}\n", Tad);
        print("Flame speed for phi={} is {} m/s.\n", phi, Uvec[0]);

        ofstream outfile("flamespeed.csv", ios::trunc);
        outfile << "  Grid, Temperature, Uvec, density, HCO+(x200), H3O+, E, ePotential, efield\n";
        
        for (size_t n = 0; n < flow.nPoints(); n++) {
            print(outfile, " {:11.3e}, {:11.3e}, {:11.3e}, {:11.3e}, {:11.3e}, {:11.3e}, {:11.3e}, {:11.3e}, {:11.3e}\n",
                  flow.grid(n), Tvec[n], Uvec[n], rho_vec[n], HCOxvec[n], H3Oxvec[n], Evec[n], ePvec[n], eFvec[n]);
        }
    } catch (CanteraError& err) {
        cerr << err.what() << endl;
        cerr << "program terminating." << endl;
        return -1;
    }
    return 0;
}

int main()
{
    double phi;
    print("Enter phi: ");
    cin >> phi;
    return flamespeed(phi);
}
