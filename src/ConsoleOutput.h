/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef CONSOLE_OUTPUT_H
#define CONSOLE_OUTPUT_H

#include "BasicTypes.h" //For uint
#include "OutputAbstracts.h"
#include "Molecules.h"
#include "MoleculeKind.h"
#include "StaticVals.h"
#include "PDBSetup.h"
#include "MoveConst.h"
#include "OutputVars.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>

class System;
namespace config_setup
{
struct Output;
}
class SystemPotential;
class Energy;
class Virial;
class MoveSettings;
class MoleculeLookup;

struct ConsoleOutput : OutputableBase {
public:
  ConsoleOutput(OutputVars & v)
  {
    this->var = &v;
  }

  //Console Output does not need to sample, so does nothing.
  virtual void Sample(const ulong step) {}


  virtual void Init(pdb_setup::Atoms const& atoms,
                    config_setup::Output const& output)
  {
    enableOut = output.console.enable;
    stepsPerOut = output.console.frequency;
    enableEnergy = output.statistics.vars.energy.fluct;
    enablePressure = output.statistics.vars.pressure.fluct;
    enableSurfTension = output.statistics.vars.surfaceTension.fluct;
#ifdef VARIABLE_VOLUME
    enableVolume = output.statistics.vars.volume.fluct;
#else
    enableVolume = false;
#endif

#ifdef VARIABLE_PARTICLE_NUMBER
    enableMol = output.statistics.vars.molNum.fluct;
#else
    enableMol = false;
#endif
    enableDens = output.statistics.vars.density.fluct;
    if (enableVolume || enablePressure || enableMol || enableDens ||
        enableSurfTension) {
      enableStat = true;
    }

    usingRE = output.usingRE;
    writingReplica = output.writingReplica;

    directory_stream << "temp_" << output.temp;
    directory_name = directory_stream.str();
    replica_directory = opendir(directory_name.c_str());

    if(replica_directory){
    
        printf("Directory already exists : %s\n", directory_name.c_str());
        /* Do whatever gromacs does here w the backups #1, #2, ect.
            path_stream << "./" << out.statistics.settings.uniqueStr.val << ".console_out";
            path_string = path_stream.str();
        */
        path_stream << "./" << directory_name << "/" << output.statistics.settings.uniqueStr.val << ".console_out";
        path_string = path_stream.str();
    } else {
        printf("Creating directory : %s\n", directory_name.c_str());
        const int dir_err = mkdir(directory_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (-1 == dir_err){
            printf("Error creating directory! Writing in pwd\n");
            path_stream << "./" << output.statistics.settings.uniqueStr.val << ".console_out";
            path_string = path_stream.str();
        } else {
            path_stream << "./" << directory_name << "/" << output.statistics.settings.uniqueStr.val << ".console_out";
            path_string = path_stream.str();
        }
    }

    const char *path = path_string.c_str();

    this->rep_out = new ofstream(path, std::ofstream::out);

    DoOutput(0);
  }
  virtual void DoOutput(const ulong step);

private:
  const static int elementWidth = 16;
  bool enableEnergy, enablePressure, enableDens, enableVolume, enableMol;
  bool enableSurfTension, enableStat;
  void PrintMove(const uint box, const ulong step) const;
  void PrintMoveStat(const uint box, const ulong step) const;
  void PrintStatistic(const uint box, const ulong step) const;
  void PrintPressureTensor(const uint box, const ulong step) const;
  void PrintEnergy(const uint box, Energy const& en, Virial const& vir,
                   const ulong step) const;
  void PrintEnergyTitle();
  void PrintStatisticTitle();
  void PrintMoveTitle();
  void printElement (const double t, const int width, uint percision = 4) const;
  void printElement (const uint t, const int width) const;
  void printElement (const std::string t, const int width) const;

template <typename T> void printElementStep ( const T t, const ulong step,
                                                const int width) const;

  
  DIR *replica_directory = NULL;
  std::stringstream directory_stream;
  std::string directory_name;
  std::stringstream path_stream;
  std::string path_string;

  ofstream* rep_out;

  bool usingRE = 0;
  bool writingReplica = 0;
};

#endif /*CONSOLE_OUTPUT_H*/
