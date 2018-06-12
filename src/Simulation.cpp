/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "Simulation.h"
#include "Setup.h"          //For setup object
// GJS
#include "barebones_Replica.cpp"
// GJS

#include "EnergyTypes.h"
#include "PSFOutput.h"
#include <iostream>
#include <iomanip>
#include <omp.h>

Simulation::Simulation(char const*const configFileName)
{
  //NOTE:
  //IMPORTANT! Keep this order...
  //as system depends on staticValues, and cpu sometimes depends on both.
  barebones_Replica* re;

  Setup set;
  set.Init(configFileName);
// GJS
  usingRE = set.config.sys.usingRE;
  replica_temps = set.config.sys.T.replica_temps;
  num_replicas = replica_temps.size(); 
// GJS
  totalSteps = set.config.sys.step.total;
  staticValues = new StaticVals(set);
  system = new System(*staticValues);
  staticValues->Init(set, *system);
  system->Init(set);
  //recall Init for static value for initializing ewald since ewald is
  //initialized in system
  staticValues->InitOver(set, *system);
  cpu = new CPUSide(*system, *staticValues);
  cpu->Init(set.pdb, set.config.out, set.config.sys.step.equil,
            totalSteps);
    int thread_num = omp_get_thread_num();
    cout << "A Thread num " << thread_num << endl;
  // GJS
  if (usingRE){
      re = new barebones_Replica(thread_num, num_replicas, set.config.sys.T.inKelvin, staticValues->mol.count, replica_temps);
      re->init(set.config.sys.step.exchange, set.config.in.prng.seed); 
  }
  // GJS

  //Dump combined PSF
  PSFOutput psfOut(staticValues->mol, *system, set.mol.kindMap,
                   set.pdb.atoms.resKindNames);
  psfOut.PrintPSF(set.config.out.state.files.psf.name);
  std::cout << "Printed combined psf to file "
            << set.config.out.state.files.psf.name << '\n';

}

Simulation::~Simulation()
{
  delete cpu;
  delete system;
  delete staticValues;
}

void Simulation::Init(Setup &set)
{
   //NOTE:
   //IMPORTANT! Keep this order...
   //as system depends on staticValues, and cpu sometimes depends on both.
   //Setup set;
   set.Init(set.config);
   totalSteps = set.config.sys.step.total;
   staticValues = new StaticVals(set);
   system = new System(*staticValues);
   staticValues->Init(set, *system); 
   system->Init(set);
   //recal Init for static value for initializing ewald since ewald is
   //initialized in system
   staticValues->InitOver(set, *system);
   cpu = new CPUSide(*system, *staticValues);
   cpu->Init(set.pdb, set.config.out, set.config.sys.step.equil,
             totalSteps);

   //Dump combined PSF
   PSFOutput psfOut(staticValues->mol, *system, set.mol.kindMap, 
		    set.pdb.atoms.resKindNames);
   psfOut.PrintPSF(set.config.out.state.files.psf.name);
   std::cout << "Printed combined psf to file " 
	     << set.config.out.state.files.psf.name << '\n';

}


Simulation::Simulation()
{
  cpu = nullptr;
  system = nullptr;
  staticValues = nullptr;
}


void Simulation::RunSimulation(void)
{
// GJS

  std::cout << "GJS" << omp_get_thread_num() << std::endl;
// GJS
  double startEnergy = system->potential.totalEnergy.total;
  for (ulong step = 0; step < totalSteps; step++) {
    system->moveSettings.AdjustMoves(step);
// GJS
    if (usingRE) {
        // GJS
        // This epot is either : 1) old , 2) expensive to calculate
        // I need to align the exchange frequency with the potential calc freq
        // GJS
  //      std::cout << step <<"\n" << std::endl;
        //system->moveSettings.ExchangeMoves(step, re, system->potential.totalEnergy.total);
        if(system->moveSettings.ExchangeMoves(step))
            re->replica_exchange((system->calcEnergy.SystemTotal()).totalEnergy.total, step);
    }
// GJS
    system->ChooseAndRunMove(step);
    cpu->Output(step);

    if((step + 1) == cpu->equilSteps) {
      double currEnergy = system->potential.totalEnergy.total;
      if(abs(currEnergy - startEnergy) > 1.0e+10) {
        printf("Info: Performing total energy calculation to preserve the"
               " energy information.\n\n");
        system->calcEwald->Init();
        system->potential = system->calcEnergy.SystemTotal();
      }
    }

#ifndef NDEBUG
    if((step + 1) % 1000 == 0)
      RunningCheck(step);
#endif
  }
  system->PrintTime();
}

void Simulation::RunSimulation(ulong step)
{
  system->ChooseAndRunMove(step);
  cpu->Output(step);
#ifndef NDEBUG
  if ((step + 1) % 1000 == 0)
    RunningCheck(step);
#endif
}


#ifndef NDEBUG
void Simulation::RunningCheck(const uint step)
{
  system->calcEwald->UpdateVectorsAndRecipTerms();
  SystemPotential pot = system->calcEnergy.SystemTotal();

  std::cout
      << "================================================================="
      << std::endl << "-------------------------" << std::endl
      << " STEP: " << step + 1
      << std::endl << "-------------------------" << std::endl
      << "Energy       INTRA B |     INTRA NB |         INTER |           TC |         REAL |         SELF |   CORRECTION |        RECIP"
      << std::endl
      << "System: "
      << std::setw(12) << system->potential.totalEnergy.intraBond << " | "
      << std::setw(12) << system->potential.totalEnergy.intraNonbond << " | "
      << std::setw(12) << system->potential.totalEnergy.inter << " | "
      << std::setw(12) << system->potential.totalEnergy.tc << " | "
      << std::setw(12) << system->potential.totalEnergy.real << " | "
      << std::setw(12) << system->potential.totalEnergy.self << " | "
      << std::setw(12) << system->potential.totalEnergy.correction << " | "
      << std::setw(12) << system->potential.totalEnergy.recip << std::endl
      << "Recalc: "
      << std::setw(12) << pot.totalEnergy.intraBond << " | "
      << std::setw(12) << pot.totalEnergy.intraNonbond << " | "
      << std::setw(12) << pot.totalEnergy.inter << " | "
      << std::setw(12) << pot.totalEnergy.tc << " | "
      << std::setw(12) << pot.totalEnergy.real << " | "
      << std::setw(12) << pot.totalEnergy.self << " | "
      << std::setw(12) << pot.totalEnergy.correction << " | "
      << std::setw(12) << pot.totalEnergy.recip << std::endl
      << "================================================================"
      << std::endl << std::endl;

}
#endif
