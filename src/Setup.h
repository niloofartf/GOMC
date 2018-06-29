/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef SETUP_H
#define SETUP_H

#include <string> //for filename
#include <cstdlib>

#include "ConfigSetup.h"
#include "FFSetup.h"
#include "PDBSetup.h"
#include "PRNGSetup.h"
#include "MolSetup.h"
//#include "Repl_Ex.h"

class Setup
{
public:
  //Read order follows each item
  ConfigSetup config;  //1
  PDBSetup pdb;        //2
  FFSetup ff;          //3
  PRNGSetup prng;      //4
  MolSetup mol;        //5
  ReplicaExchangeParameters replExParams;

  void Init(char const*const configFileName)
  {
    //Read in all config data
    config.Init(configFileName);
    //Read in FF data.
    ff.Init(config.in.files.param.name, config.in.ffKind.isCHARMM);
    //Read PDB dat
    pdb.Init(config.in.restart, config.in.files.pdb.name);
    //Read molecule data from psf
    prng.Init(config.in.restart, config.in.prng, config.in.files.seed.name);

    if(mol.Init(config.in.restart, config.in.files.psf.name) != 0) {
      exit(EXIT_FAILURE);
    }
    mol.AssignKinds(pdb.atoms, ff);

  }

  void Init(char const*const configFileName, int initiatingLoopIteration, ReplicaExchangeParameters* replExParams)
  {
    printf("FAST iteration %d entered set.Init\n", initiatingLoopIteration);
    //Read in all config data
    config.Init(configFileName, initiatingLoopIteration, replExParams);
    printf("FAST iteration %d returned from config.Init\n", initiatingLoopIteration);
    //Read in FF data.
    ff.Init(config.in.files.param.name, config.in.ffKind.isCHARMM);
    printf("FAST iteration %d returned from ff.Init\n", initiatingLoopIteration);

    // Add some type of functionailty to try and figure out what the psf and pdf's are named

    std::string inputFileString1 = config.in.files.pdb.name[0];
    std::string inputFileString2 = config.in.files.psf.name[0];
    std::string suffix1 = "_BOX_0_restart.pdb";
    std::string suffix2 = "_merged.psf";

    printf("FAST config.in.files.pdb.name[0] : %s\n", (config.in.files.pdb.name[0]).c_str());
    printf("FAST config.in.files.psf.name[0] : %s\n", (config.in.files.psf.name[0]).c_str());
  
    fstream inputFileReader1;
    //OPEN FILE
    inputFileReader1.open(inputFileString1.c_str(), ios::in | ios::out);

    //CHECK IF FILE IS OPENED...IF NOT OPENED EXCEPTION REASON FIRED
    if (!inputFileReader1.is_open()) {
        printf("FAST I couldn't open %s\n", (config.in.files.pdb.name[0]).c_str());
        std::size_t found = inputFileString1.find(suffix1);
        printf("FAST I tried searching %s for %s\n", (config.in.files.pdb.name[0]).c_str(), suffix1.c_str());
        if (found!=std::string::npos){
            printf("FAST I found %s at position %lu in %s\n", suffix1.c_str(), found, (config.in.files.pdb.name[0]).c_str()); 
            string temp_in = std::to_string((int)(config.sys.T.replica_temps[initiatingLoopIteration]));
            temp_in += 'K';
            printf("FAST temp_in : %s\n", temp_in.c_str());
            inputFileString1.insert(found, temp_in);
            printf("FAST post insertion : %s\n", inputFileString1.c_str());
            inputFileReader1.open(inputFileString1.c_str(), ios::in | ios::out);
            
            printf("FAST replica %d opened : %s\n", initiatingLoopIteration, inputFileString1.c_str());
        }
    }
    
    //CHECK IF FILE IS OPENED...IF NOT OPENED EXCEPTION REASON FIRED
    if (!inputFileReader1.is_open()) {
      std::cout << "Error: Cannot open/find " << inputFileString1 <<
                " in the directory provided!\n";
      exit(EXIT_FAILURE);
    }
    
    //CLOSE FILE TO NOW PASS TO SIMULATION
    inputFileReader1.close();

    printf("FAST time to play with inputFileReader2\n");
    cout.flush();
    #pragma omp barrier

    fstream inputFileReader2;
    //OPEN FILE
    inputFileReader2.open(inputFileString2.c_str(), ios::in | ios::out);

    //CHECK IF FILE IS OPENED...IF NOT OPENED EXCEPTION REASON FIRED
    if (!inputFileReader2.is_open()) {
        std::size_t found = inputFileString2.find(suffix2);
        if (found!=std::string::npos){
            string temp_in = std::to_string((int)(config.sys.T.replica_temps[initiatingLoopIteration]));
            temp_in += 'K';
            inputFileString2.insert(found, temp_in);
            inputFileReader2.open(inputFileString2.c_str(), ios::in | ios::out);
        }
    }

    //CHECK IF FILE IS OPENED...IF NOT OPENED EXCEPTION REASON FIRED
    if (!inputFileReader2.is_open()) {
      std::cout << "Error: Cannot open/find " << inputFileString2 <<
                " in the directory provided!\n";
      exit(EXIT_FAILURE);
    }
    
    //CLOSE FILE TO NOW PASS TO SIMULATION
    inputFileReader2.close();
    
    printf("FAST time to assign the strings back to their config...\n");
    cout.flush();
    #pragma omp barrier

    config.in.files.pdb.name[0] = inputFileString1;

    config.in.files.psf.name[0] = inputFileString2;
    
    //Read PDB dat
    pdb.Init(config.in.restart, config.in.files.pdb.name);
    //Read molecule data from psf
    prng.Init(config.in.restart, config.in.prng, config.in.files.seed.name);

    if(mol.Init(config.in.restart, config.in.files.psf.name) != 0) {
      exit(EXIT_FAILURE);
    }
    mol.AssignKinds(pdb.atoms, ff);

  }

/*


  void Init(char const*const configFileName, int initiatingLoopIteration, ReplicaExchangeParameters* replExParams)
  {
    //Read in all config data
    config.Init(configFileName, initiatingLoopIteration, replExParams);
    //Read in FF data.
    ff.Init(config.in.files.param.name, config.in.ffKind.isCHARMM);
    //Read PDB dat
    pdb.Init(config.in.restart, config.in.files.pdb.name);
    //Read molecule data from psf
    prng.Init(config.in.restart, config.in.prng, config.in.files.seed.name);

    if(mol.Init(config.in.restart, config.in.files.psf.name) != 0) {
      exit(EXIT_FAILURE);
    }
    mol.AssignKinds(pdb.atoms, ff);

  }
*/
};

#endif
