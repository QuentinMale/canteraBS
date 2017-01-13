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
#include <fstream>

using namespace std;
using namespace Cantera;
using fmt::print;

int flamespeed(double phi)
{
    try {
        IdealGasMix gas("methane_ion.xml");
        IdealGasMix gasCopy("methane_ion.xml");

        doublereal temp = 300.0; // K
        doublereal pressure = 1.0*OneAtm; //atm
        doublereal uin = 0.3; //m/sec

        size_t nsp = gas.nSpecies();
        vector_fp x(nsp, 0.0);

        doublereal C_atoms = 1.0;
        doublereal H_atoms = 4.0;
        doublereal ax = C_atoms + H_atoms / 4.0;
        doublereal fa_stoic = 1.0 / (4.76 * ax);
        x[gas.speciesIndex("CH4")] = 1.0;
        x[gas.speciesIndex("O2")] = 0.21 / phi / fa_stoic;
        x[gas.speciesIndex("N2")] = 0.79 / phi/ fa_stoic;

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

        std::unique_ptr<Transport> trmix(newTransportMgr("Mix", &gas));

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

        std::vector<Domain1D*> domains { &inlet, &flow, &outlet };
        Sim1D flame(domains);

        //----------- Supply initial guess----------------------

        vector_fp locs{0.0, 0.3, 0.7, 1.0};
        vector_fp value;

        double uout = inlet.mdot()/rho_out;
        value = {uin, uin, uout, uout};
        flame.setInitialGuess("u",locs,value);
        value = {temp, temp, Tad, Tad};
        flame.setInitialGuess("T",locs,value);

        for (size_t i=0; i<nsp; i++) {
            value = {yin[i], yin[i], yout[i], yout[i]};
            flame.setInitialGuess(gas.speciesName(i),locs,value);
        }

        inlet.setMoleFractions(x.data());
        inlet.setMdot(mdot);
        inlet.setTemperature(temp);

        flame.showSolution();

        int flowdomain = 1;
        double ratio = 5.0;
        double slope = 0.04;
        double curve = 0.05;

        flame.setRefineCriteria(flowdomain,ratio,slope,curve);

        int loglevel=1;

        // Solve freely propagating flame

        // Linearly interpolate to find location where this temperature would
        // exist. The temperature at this location will then be fixed for
        // remainder of calculation.
        flame.setFixedTemperature(0.5 * (temp + Tad));
        flow.solveEnergyEqn();
        bool refine_grid = true;

        // not sure what to do!!
        for (size_t i = 184; i < 186; i++) {
            gas.setMultiplier(i,0.0);
        }
        
        for (size_t i = 186; i < 190; i++) {
            gas.setMultiplier(i,0);
        }
        
        flame.solve(loglevel,refine_grid);
        double flameSpeed_mix = flame.value(flowdomain,flow.componentIndex("u"),0);
        print("Flame speed with mixture-averaged transport: {} m/s\n",
              flameSpeed_mix);
        
        // finish standard Cantera calculation
        // generate initial guess for the charged species
        // the list for important species
        vector_fp list;
        list.push_back(gas.speciesIndex("H2O"));
        list.push_back(gas.speciesIndex("CO2"));
        list.push_back(gas.speciesIndex("N2"));
        list.push_back(gas.speciesIndex("O2"));
        list.push_back(gas.speciesIndex("H2"));
        list.push_back(gas.speciesIndex("CO"));
        list.push_back(gas.speciesIndex("H"));
        list.push_back(gas.speciesIndex("OH"));
        list.push_back(gas.speciesIndex("O"));
        list.push_back(gas.speciesIndex("CH"));
        list.push_back(gas.speciesIndex("CH3"));
        list.push_back(gas.speciesIndex("CH4"));
        list.push_back(gas.speciesIndex("C2H2"));
        list.push_back(gas.speciesIndex("C2H4"));

        for (size_t n = 0; n < flow.nPoints(); n++) {
            double locTemp = flame.value(flowdomain,flow.componentIndex("T"),n);
            vector_fp y(nsp, 0.0);
            for (size_t k : list) {
                y[k] = flame.value(flowdomain,k+c_offset_Y,n);
            }
            gasCopy.setState_TPY(locTemp,pressure,y.data());
            gasCopy.equilibrate("HP");
            gasCopy.getMassFractions(&y[0]);
            for (size_t k : flow.chargeList()) {
                flame.setValue(flowdomain, c_offset_Y + k, n, y[k]);
            }
        }
        
        //****************** end of phase 1*****************
        // set solving phase for IonFlow
        flow.setSolvingPhase(2);
        refine_grid = false;
        flow.fixTemperature();
        flow.fixVelocity();
        
        for (size_t j = 0; j < 11; j++) {
            for (size_t i = 184; i < 190; i++) {
                gas.setMultiplier(i,0.1*j);
            }
            flame.solve(loglevel, refine_grid);
        }
        
        // turn on energy equation
        flow.solveEnergyEqn();
        flow.solveVelocity();
        for (size_t j = 0; j < 1; j++) {
            flame.solve(loglevel, refine_grid);
        }
        
        /*
        //****************** end of phase 2 ****************
        // set solving phase for IonFlow
        flow.setSolvingPhase(3);
        locs = {0, 0.3, 0.5, 0.7, 1.0};
        value = {0, 1, 10, 1, 0};
        flame.setInitialGuess("ePotential",locs,value);

        refine_grid = true;

        for (size_t i = 0; i < 1; i++) {
            // now enable Poisson's equation but hold y and T constant
            flow.solvePoissonEqn();
            flow.enableElectric(false);
            flow.fixSpeciesMassFrac();
            flame.solve(loglevel, refine_grid);

            // now enable electric effect but hold E and V constant
            flow.enableElectric(true);

        }*/

        vector_fp zvec,Tvec,HCOxvec,H3Oxvec,Evec,ePvec,eFvec,Uvec;

        for (size_t n = 0; n < flow.nPoints(); n++) {
            double rho = flow.density(n);
            size_t k = gas.speciesIndex("HCO+");
            HCOxvec.push_back(1e-6*Avogadro*rho*flame.value(flowdomain,k+c_offset_Y,n)/mw[k]);
            k = gas.speciesIndex("H3O+");
            H3Oxvec.push_back(1e-6*Avogadro*rho*flame.value(flowdomain,k+c_offset_Y,n)/mw[k]);
            k = gas.speciesIndex("E");
            Evec.push_back(1e-6*Avogadro*rho*flame.value(flowdomain,k+c_offset_Y,n)/mw[k]);
        } 

        print("\n{:9s}\t{:8s}\t{:8s}\t{:5s}\n",
              "z (m)", "T (K)", "U (m/s)", "eP[V]");

        Tvec.resize(flow.nPoints(), 0.0);
        zvec.resize(flow.nPoints(), 0.0);

        for (size_t n = 0; n < flow.nPoints(); n++) {
            Tvec[n] = flame.value(flowdomain,flow.componentIndex("T"),n);
            ePvec.push_back(flame.value(flowdomain,flow.componentIndex("ePotential"),n));
            Uvec.push_back(flame.value(flowdomain,flow.componentIndex("u"),n));
            zvec[n] = flow.grid(n);
            print("{:9.6f}\t{:8.3f}\t{:8.3f}\t{:5.3f}\n",
                  flow.grid(n), Tvec[n], Uvec[n], ePvec[n]);
        }

        eFvec.push_back(0);
        for (size_t n = 1; n < flow.nPoints()-1; n++) {
            eFvec.push_back(-(ePvec[n+1]-ePvec[n-1])/(flow.grid(n+1)-flow.grid(n-1)));
        }      
        eFvec.push_back(eFvec[flow.nPoints()-2]);
        eFvec[0] = eFvec[1];

        print("\nAdiabatic flame temperature from equilibrium is: {}\n", Tad);
        print("Flame speed for phi={} is {} m/s.\n", phi, Uvec[0]);

        std::ofstream outfile("flamespeed.csv", std::ios::trunc);
        outfile << "  Grid, Temperature, Uvec, HCO+, H3O+, E, ePotential, efield\n";
        
        for (size_t n = 0; n < flow.nPoints(); n++) {
            print(outfile, " {:11.3e}, {:11.3e}, {:11.3e}, {:11.3e}, {:11.3e}, {:11.3e}, {:11.3e}, {:11.3e}\n",
                  flow.grid(n), Tvec[n], Uvec[n], HCOxvec[n], H3Oxvec[n], Evec[n], ePvec[n], eFvec[n]);
        }
    } catch (CanteraError& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << "program terminating." << std::endl;
        return -1;
    }
    return 0;
}

int main()
{
    double phi;
    print("Enter phi: ");
    std::cin >> phi;
    return flamespeed(phi);
}
