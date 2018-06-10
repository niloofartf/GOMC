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


class Simulation
{
public:
// GJS  
  bool usingRE;
  std::vector<double> replica_temps;
  int num_replicas;
// GJS
  explicit Simulation(char const*const configFileName);
  ~Simulation();

  void RunSimulation(void);


#ifndef NDEBUG
  void RunningCheck(const uint step);
#endif

private:
  StaticVals * staticValues;
  System * system;
  CPUSide * cpu;
  ulong totalSteps;
  barebones_Replica* re;
};

#endif /*SIMULATION_H*/
