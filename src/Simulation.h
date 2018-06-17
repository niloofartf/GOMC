/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef SIMULATION_H
#define SIMULATION_H

//Member vars
#include "CPUSide.h"
#include "System.h"
#include "StaticVals.h"
#include "BasicTypes.h"

#include "Repl_Ex.h"


class Simulation
{
public:
// GJS  
  bool usingRE;
  int numberOfAtoms;
  std::string replica_log;
  std::vector<double> replica_temps;
// GJS
  explicit Simulation(char const*const configFileName);
  explicit Simulation(char const*const configFileName, int initiatingLoopIteration, ReplicaExchangeParameters* replExParams);
  Simulation(void);
  ~Simulation();

  void RunSimulation(void);
  void RunSimulation(ReplicaExchangeParameters* replExParams);

#ifndef NDEBUG
  void RunningCheck(const uint step);
#endif

private:
  StaticVals * staticValues;
  System * system;
  CPUSide * cpu;
  ulong totalSteps;
};

#endif /*SIMULATION_H*/
