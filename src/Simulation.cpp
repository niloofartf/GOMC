/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "Simulation.h"
#include "Setup.h"          //For setup object

#include "EnergyTypes.h"
#include "PSFOutput.h"
#include <iostream>
#include <iomanip>
#include <omp.h>

// GJS
#include "Repl_Ex.cpp"
#include "Replica_State.h"
//

Simulation::Simulation(char const*const configFileName)
{
  //NOTE:
  //IMPORTANT! Keep this order...
  //as system depends on staticValues, and cpu sometimes depends on both.
  Setup set;
  set.Init(configFileName);
// GJS
  int numberOfAtoms;
  usingRE = set.config.sys.usingRE;
  replica_temps = set.config.sys.T.replica_temps;
// GJS
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

  numberOfAtoms = set.pdb.atoms.x.size() ;
  psfOut.PrintPSF(set.config.out.state.files.psf.name);
  std::cout << "Printed combined psf to file "
            << set.config.out.state.files.psf.name << '\n';

  std::cout << "Finished initializing Sim object "
            << set.config.out.state.files.psf.name << '\n';
  
}

Simulation::Simulation(char const*const configFileName, int initiatingLoopIteration, ReplicaExchangeParameters* replExParams)
{
  //NOTE:
  //IMPORTANT! Keep this order...
  //as system depends on staticValues, and cpu sometimes depends on both.
  Setup set;
  set.Init(configFileName, initiatingLoopIteration, replExParams);
  replica_log = set.config.out.statistics.settings.uniqueStr.val;
  replica_log += ".replica_log"; 
// GJS
  usingRE = set.config.sys.usingRE;
  replica_temps = set.config.sys.T.replica_temps;

  #pragma omp single
  {
    replExParams->replica_temps = set.config.sys.T.replica_temps;
  }
// GJS
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
    
Simulation::Simulation(void)
{
  cpu = NULL;
  system = NULL;
  staticValues = NULL;
}


Simulation::~Simulation()
{
  delete cpu;
  delete system;
  delete staticValues;
}

void Simulation::RunSimulation(void)
{
  double startEnergy = system->potential.totalEnergy.total;
  for (ulong step = 0; step < totalSteps; step++) {
    system->moveSettings.AdjustMoves(step);
    system->ChooseAndRunMove(step);
    cpu->Output(step);

    if((step + 1) == cpu->equilSteps) {
      double currEnergy = system->potential.totalEnergy.total;
      printf("energy : %f\n", currEnergy);
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


void Simulation::RunSimulation(ReplicaExchangeParameters* replExParams, Simulation** sim_exchangers)
{

  std::cout << "GJS" << omp_get_thread_num() << std::endl;
  gmx_repl_ex_t     repl_ex = nullptr;

  Replica_State * state_local = new Replica_State();
  Replica_State * state_global = state_local;

  int bFirstStep, bInitStep, bLastStep = false;

  int bDoReplEx, bExchanged;

  FILE *fplog = fopen(replica_log.c_str(), "a");  

  const bool useReplicaExchange = (replExParams->exchangeInterval > 0);
  
  if (useReplicaExchange){
    #pragma omp barrier          
      printf(" calling init w temp of %f\n", staticValues->forcefield.T_in_K);
      repl_ex = init_replica_exchange(fplog, staticValues->forcefield.T_in_K, replExParams);
    }


// GJS
  double startEnergy = system->potential.totalEnergy.total;
    
 
  bFirstStep       = true;

  // Find out where the flag for a restart is located and set this accordingly GJS
  bool startingFromCheckpoint = false;

  bInitStep        = !startingFromCheckpoint;

  for (ulong step = 0; step < totalSteps; step++) {

    system->moveSettings.AdjustMoves(step);
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

  
    bDoReplEx = (useReplicaExchange && (step > 0) && !bLastStep && (step % replExParams->exchangeInterval == 0) && step > cpu->equilSteps);

    if (bDoReplEx) {
    #pragma omp barrier          
            printf("GJS Step : %lu, repl : %d, b4 exch : system->epot = %f\n", step, repl_ex->repl, system->potential.totalEnergy.total);
            bExchanged = replica_exchange(fplog, repl_ex,
                                           state_local, system->potential.totalEnergy.total,
                                           step, replExParams);
   
            if (bExchanged != repl_ex->repl){
                SetTemp(system, sim_exchangers[bExchanged]);        
            }
 
            #pragma omp barrier          
            
            if (bExchanged != repl_ex->repl){
                GetTemp(system, sim_exchangers[repl_ex->repl]); 
            }
            
            #pragma omp barrier          
            
            printf("GJS Step : %lu, repl : %d, after exch : system->epot = %f\n", step, repl_ex->repl, system->potential.totalEnergy.total);
        }
#ifndef NDEBUG
    if((step + 1) % 1000 == 0)
      RunningCheck(step);
#endif
  }
  system->PrintTime();
    if (useReplicaExchange)
     {
         print_replica_exchange_statistics(fplog, repl_ex);
     }
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

void Simulation::GetTemp(System* system, Simulation* sim){

    system->potential           =   sim->system->potential;
    system->com                 =   sim->system->com;
    system->coordinates         =   sim->system->coordinates;
    system->cellList            =   sim->system->cellList;
//    *(system->calcEwald)        =   *(sim->system->calcEwald);

}

void Simulation::SetTemp(System* system_set, Simulation* sim){
    
    sim->system->potential      =   system_set->potential;
    sim->system->com            =   system_set->com;
    sim->system->coordinates    =   system_set->coordinates;
    sim->system->cellList       =   system_set->cellList;
//    *(sim->system->calcEwald)      =   *(system_set->calcEwald);

}
