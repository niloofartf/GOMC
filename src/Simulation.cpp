/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "Simulation.h"
#include "Setup.h"          //For setup object
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/mdlib/stat.cpp"
#include "gromacs/mdlib/md_support.cpp"

#include "EnergyTypes.h"
#include "PSFOutput.h"
#include <iostream>
#include <iomanip>
#include <omp.h>

// GJS
//#include "repl_ex.h"
#include "gromacs/mdtypes/state.h"
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
  replica_log = set.config.out.statistics.settings.uniqueStr.val;

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
      system->potential = system->calcEnergy.SystemTotal();
      double currEnergy = system->potential.totalEnergy.total;
    system->moveSettings.AdjustMoves(step);
    system->ChooseAndRunMove(step);
    cpu->Output(step);
      printf("energy : %f\n", currEnergy);

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


void Simulation::RunSimulation(const ReplicaExchangeParameters &replExParams)
{
// GJS
  PaddedRVecVector  f {};

  std::cout << "GJS" << omp_get_thread_num() << std::endl;
  gmx_repl_ex_t     repl_ex = nullptr;

  t_state * state_global;
  state_global = new t_state();

  set_state_entries(state_global, inputrec);

// Local state only becomes valid now.
  t_state *                state;


  state_change_natoms(state_global, state_global->natoms);
  /* We need to allocate one element extra, since we might use
   * (unaligned) 4-wide SIMD loads to access rvec entries.
   */
  f.resize(gmx::paddedRVecVectorSize(state_global->natoms));
  /* Copy the pointer to the global state */
  state = state_global;


  gmx::ThreeFry2x64<64> rng(replExParams.randomSeed, gmx::RandomDomain::ReplicaExchange);
  gmx::UniformRealDistribution<real>   uniformRealDist;

  gmx_bool bFirstStep, bInitStep, bLastStep = false;

  gmx_bool bDoReplEx, bExchanged;

  FILE *fplog = fopen (replica_log.c_str(), "a");  


 const bool useReplicaExchange = (replExParams.exchangeInterval > 0);
  //if (useReplicaExchange && MASTER(cr))
    // pragma omp master
  if (useReplicaExchange){

//        printf(" calling init, passing natoms : %d\n", top_global->natoms);
//        repl_ex = init_replica_exchange(fplog, cr->ms, top_global->natoms, ir, replExParams);
        printf(" calling init\n");
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

    bDoReplEx = (useReplicaExchange && (step > 0) && !bLastStep &&
                      do_per_step(step, replExParams.exchangeInterval) && ( step > cpu->equilSteps));

    if (bDoReplEx) {
 
             bExchanged = replica_exchange(fplog, repl_ex,
                                           state_global, system->potential.totalEnergy.total,
                                           state, step);
    }

    //if (useReplicaExchange && MASTER(cr))
    // pragma omp master
    if (useReplicaExchange)
     {
         print_replica_exchange_statistics(fplog, repl_ex);
     }


#ifndef NDEBUG
    if((step + 1) % 1000 == 0)
      RunningCheck(step);
#endif
  }
  system->PrintTime();
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
