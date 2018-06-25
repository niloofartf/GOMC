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
  printf("FAST I returned from set.Init\n");
  replica_log = set.config.out.statistics.settings.uniqueStr.val;
  replica_log += ".replica_log"; 
// GJS
  usingRE = set.config.sys.usingRE;

  replica_temps = set.config.sys.T.replica_temps;
  #pragma omp barrier
  #pragma omp single 
  {
    for (int i = 0; i < set.config.sys.T.replica_temps.size(); i++){
        replExParams->replica_temps.push_back(set.config.sys.T.replica_temps[i]);
    }
  }
  #pragma omp barrier
 
#if ENSEMBLE == NPT
  #pragma omp barrier
  #pragma omp single 
  {
    for (int i = 0; i < set.config.sys.gemc.replica_pressures.size(); i++){
        replExParams->replica_pressures.push_back(set.config.sys.gemc.replica_pressures[i]);
    }    
    replExParams->bNPT = 1;
  }
  #pragma omp barrier
#endif

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


void Simulation::RunSimulation(ReplicaExchangeParameters* replExParams)
{

  std::cout << "GJS" << omp_get_thread_num() << std::endl;
  gmx_repl_ex_t     repl_ex = nullptr;

  Replica_State * state_global = new Replica_State();

  int bFirstStep, bInitStep, bLastStep = false;

  int bDoReplEx, bExchanged;

  FILE *fplog = fopen(replica_log.c_str(), "a");  

  const bool useReplicaExchange = (replExParams->exchangeInterval > 0);
 
  if (useReplicaExchange){
    #pragma omp barrier          
      printf(" calling init w temp of %f\n", staticValues->forcefield.T_in_K);
      repl_ex = init_replica_exchange(fplog, staticValues->forcefield.T_in_K, replExParams);
    }

#if ENSEMBLE == NVT
  volume = 0.0;
#endif

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

  #pragma omp barrier          
    bDoReplEx = (useReplicaExchange && (step > 0) && !bLastStep && (step % replExParams->exchangeInterval == 0) && ( step > cpu->equilSteps));

    if (bDoReplEx) {
  #pragma omp barrier          
       
            GetSystem(state_global, system);
           
#if ENSEMBLE == NPT
            // These two prints should be equal.  Its a nice checker that the GetSys works and were keeping track of everything
            //printf("OTTO repl id : %d, vol : %f\n", omp_get_thread_num(), system->boxDimensions->GetTotVolume());
            //printf("OTTO repl id : %d, vol : %f\n", omp_get_thread_num(), (*(state_global->boxDimensions))->GetTotVolume());
            bExchanged = replica_exchange(fplog, repl_ex,
                                           state_global, system->potential.totalEnergy.total,
                                           (*(state_global->boxDimensions))->GetTotVolume(),
                                           step, replExParams);
#endif 
#if ENSEMBLE == NVT
            // These two prints should be equal.  Its a nice checker that the GetSys works and were keeping track of everything
            //printf("OTTO repl id : %d, vol : %f\n", omp_get_thread_num(), system->boxDimensions->GetTotVolume());
//            printf("OTTO repl id : %d, vol : %f\n", omp_get_thread_num(), volume);
            bExchanged = replica_exchange(fplog, repl_ex,
                                           state_global, system->potential.totalEnergy.total,
                                           volume,
                                           step, replExParams);
#endif 
            if (bExchanged){
                #pragma omp barrier
                printf("GJS Step : %lu, before repl : %d ; system epot : %f\n", step, repl_ex->repl, system->potential.totalEnergy.total);
                printf("GJS Step : %lu, before repl : %d ; global array epot : %f\n", step, repl_ex->repl, replExParams->replica_energies[repl_ex->repl]);
                #pragma omp barrier
                state_global = replExParams->replica_states[repl_ex->repl];
                //SetSystem(state_global, system);
                #pragma omp barrier
                printf("GJS Step : %lu, after repl : %d ; system epot : %f\n", step, repl_ex->repl, system->potential.totalEnergy.total);
                printf("GJS Step : %lu, after repl : %d ; global array epot : %f\n", step, repl_ex->repl, replExParams->replica_energies[repl_ex->repl]);
                #pragma omp barrier
            }
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

void Simulation::GetSystem(Replica_State* state_get, System* system_get){

    state_get->potential = &(system_get->potential);
    state_get->com = &(system_get->com);
    state_get->coordinates = &(system_get->coordinates);
    state_get->cellList = &(system_get->cellList);
    state_get->calcEwald = &(system_get->calcEwald);
#if ENSEMBLE == NPT
    state_get->boxDimensions = &(system_get->boxDimensions);
#endif

}

void Simulation::SetSystem(Replica_State* state_set, System* system_set){

    system_set->potential = *(state_set->potential);
    system_set->com = *(state_set->com);
    system_set->coordinates = *(state_set->coordinates);
    system_set->cellList = *(state_set->cellList);
    system_set->calcEwald = *(state_set->calcEwald);
#if ENSEMBLE == NPT
    system_set->boxDimensions = *(state_set->boxDimensions);
#endif

}
