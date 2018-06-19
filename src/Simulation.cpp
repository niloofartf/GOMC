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
  replica_log = set.config.out.statistics.settings.uniqueStr.val;
  replica_log += ".replica_log"; 
// GJS
  usingRE = set.config.sys.usingRE;
  replica_temps = set.config.sys.T.replica_temps;
  replExParams->replica_temps = set.config.sys.T.replica_temps;
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
  // GJS

  std::cout << "GJS" << omp_get_thread_num() << std::endl;
  gmx_repl_ex_t     repl_ex = nullptr;

  Replica_State * state_global = new Replica_State();

// Local state only becomes valid now.
  Replica_State *                state;

  //gmx::ThreeFry2x64<64> rng(replExParams.randomSeed, gmx::RandomDomain::ReplicaExchange);
  //gmx::UniformRealDistribution<real>   uniformRealDist;

  int bFirstStep, bInitStep, bLastStep = false;

  int bDoReplEx, bExchanged;

  printf("beth %s\n", replica_log.c_str());


  FILE *fplog = fopen(replica_log.c_str(), "a");  


 const bool useReplicaExchange = (replExParams->exchangeInterval > 0);
  //if (useReplicaExchange && MASTER(cr))
    // pragma omp master
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

    bDoReplEx = (useReplicaExchange && (step > 0) && !bLastStep && (step % replExParams->exchangeInterval == 0));// &&
                      //( step > cpu->equilSteps));

    if (bDoReplEx) {
  #pragma omp barrier          
       
            for (int i = 0; i < repl_ex->nrepl; i++){
  //              printf("re id : %d , state[%d]->epot : %f\n", repl_ex->repl, i, replExParams->replica_states[i]->potential->totalEnergy.total);    
            }

            GetSystem(state_global, system);
            state = state_global;
            printf("GJS About to call replica_exchange\n");
            printf("GJS PRE re id %d , step: %lu, epot %f\n", repl_ex->repl, step, state->potential->totalEnergy.total);
            bExchanged = replica_exchange(fplog, repl_ex,
                                           state_global, system->potential.totalEnergy.total,
                                           state, step, replExParams);
            printf("GJS Returned from call to replica_exchange\n");
            if (bExchanged){
                state = replExParams->replica_states[repl_ex->repl];
                printf("GJS POST re id %d , step: %lu, epot %f\n", repl_ex->repl, step, state->potential->totalEnergy.total);
                SetSystem(state, system);
            }
            for (int i = 0; i < repl_ex->nrepl; i++){
                //printf("re id : %d , state[%d]->epot : %f\n", repl_ex->repl, i, replExParams->replica_states[i]->potential->totalEnergy.total);    
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

}

void Simulation::SetSystem(Replica_State* state_set, System* system_set){

    system_set->potential = *(state_set->potential);
    system_set->com = *(state_set->com);
    system_set->coordinates = *(state_set->coordinates);
    system_set->cellList = *(state_set->cellList);
    system_set->calcEwald = *(state_set->calcEwald);

}
