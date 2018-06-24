/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.20
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
//TODO: print weightingFunction to file instead of console/Format console print (see manual, pg 38)
//TODO: find a way to set molKind to the correct kind of molecule (currently only good for single comp liquid/gas interface, no adsorbance)
//TODO: Check surface tension units

//Note: When adding new moves to GOMC, if they might run in a GCMC simulation, the transition matrix 
//must be updated with the acceptance probability. If the move will not end in the number of particles in the 
//main simulation box (BOX0) changing, then the acceptance probability can be captured with the following
//code chunk:
//	if ENSEMBLE == GCMC
//		transitionMatrixRef.addAcceptanceProbToMatrix(0.0, 0);
//	endif
//Files currently affected: Movebase.h, MoleculeTransfer.h, IntraSwap.h, Regrowth.h

#pragma once

#include "StaticVals.h"             //For init
#include "Forcefield.h"
#include "MoleculeLookup.h"
#include "BoxDimensions.h"
#include "BoxDimensionsNonOrth.h"
#include "ConfigSetup.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <string>

#if ENSEMBLE == GCMC

class TransitionMatrix
{
public:
	TransitionMatrix(StaticVals const& stat, MoleculeLookup const& molLookupRef,
                   BoxDimensions const& boxDimRef) : forcefield(stat.forcefield),
                   molLookRef(molLookupRef), currentAxes(boxDimRef)
	{
    biasingOn = false;
  }

	void Init(config_setup::TMMC const& tmmc);
	void AddAcceptanceProbToMatrix(double acceptanceProbability, int move);
	void IncrementAcceptanceProbability(int molKind);
	double CalculateBias(int move);
	void UpdateWeightingFunction(ulong step);
	void PrintTMProbabilityDistribution();
	int nmax, nmin;						//Specify minimum and maximum molecule counts in simulation
										//Important for TMMC because biasing allows extremely unlikely phase space explorations that
										//do not impact on the final result, causing wasted calculations and CPU time.
										//Nmin currently not implemented, but will allow for parallelization via umbrella

	bool biasingOn;						//Config flag, turns biasing on or off.
	ofstream TMfile;

	uint molKind;						//Track what kind of molecule we're interested in; currently defaults to 0 
										//(Current implementation works with single component systems only).
										//If user can select a component to track, adsorption/absorption sims become possible.

private:
  const MoleculeLookup& molLookRef;	//Used to reference number of molecules of interest in the main box
  const BoxDimensions& currentAxes;     //Used for volume
  const Forcefield& forcefield;			//Used for temperature
  ulong biasStep;                       //Biasing steps
  double boxVolume, temperature;// , vaporPressure, surfaceTension, vaporDensity, liquidDensity;
  int INITIAL_WEIGHTINGFUNCTION_VALUE;
  
  

  //std::vector<double> PostProcessTransitionMatrix();
  int GetTMDelIndex(int numMolec);
  int GetTMEtcIndex(int numMolec);
  int GetTMInsIndex(int numMolec);

  std::vector<double> transitionMatrix;			  //Holds flat 2D Transition Matrix array after initialization
  std::vector<double> weightingFunction;          //Holds calculated biasing function
};


inline void TransitionMatrix::Init(config_setup::TMMC const& tmmc) {
	biasingOn = tmmc.enable;
	biasStep = tmmc.step;
	molKind = 0;								//TODO: make input from MoleculeLookup 
	
	boxVolume = currentAxes.GetBoxVolume(mv::BOX0);
	temperature = forcefield.T_in_K;

	INITIAL_WEIGHTINGFUNCTION_VALUE = 1.0;
	int totMolec = molLookRef.NumKindInBox(molKind, mv::BOX0) + molLookRef.NumKindInBox(molKind, mv::BOX1);

	if (tmmc.Nmax != UINT_MAX && nmax < totMolec) {
		nmax = tmmc.Nmax;	//TODO: Print warnings when these changes are made if necessary
	}
	else {
		nmax = totMolec;
	}


	TMfile.open(tmmc.outName + "_TMFile.dat");		//TODO: make this an actual name
	if (TMfile.is_open()) {
		TMfile << setw(16) << left << "Temperature" << setw(16) << left << temperature << endl;
		TMfile << setw(16) << left << "Box Volume" << setw(16) << left << boxVolume << endl;
	}
	//nmin = tmmc.Nmin;
	//if (nmin < 0) {
	//	nmin = 0;
	//}

	//del/etc/ins, del/etc/ins, ... 
	transitionMatrix.resize(3 * (nmax+1));
	for (int i = 0; i < transitionMatrix.size(); i++) { transitionMatrix[i] = 0.0; }

	weightingFunction.resize(totMolec);
	for (int i = 0; i < weightingFunction.size(); i++) { weightingFunction[i] = INITIAL_WEIGHTINGFUNCTION_VALUE; }
}

inline int TransitionMatrix::GetTMDelIndex(int numMolec) { return numMolec * 3; }
inline int TransitionMatrix::GetTMEtcIndex(int numMolec) { return numMolec * 3 + 1; }
inline int TransitionMatrix::GetTMInsIndex(int numMolec) { return numMolec * 3 + 2; }

//Adds acceptance probability from GCMC particle insertion and deletion moves to the transitionMatrix vector
inline void TransitionMatrix::AddAcceptanceProbToMatrix(double acceptanceProbability, int move)
{
	if (!biasingOn)
		return;
	uint numMolecules = molLookRef.NumKindInBox(molKind, mv::BOX0);
	if (acceptanceProbability > 1.0)
		acceptanceProbability = 1.0;

	//move == 0: deletion move; move == 1: insertion move
	if (move == 0)
	{
		transitionMatrix[GetTMDelIndex(numMolecules)] += acceptanceProbability;
		transitionMatrix[GetTMEtcIndex(numMolecules)] += 1.0 - acceptanceProbability;
	}
	else if (move == 1)
	{
		transitionMatrix[GetTMInsIndex(numMolecules)] += acceptanceProbability;
		transitionMatrix[GetTMEtcIndex(numMolecules)] += 1.0 - acceptanceProbability;
	}
}

inline void TransitionMatrix::IncrementAcceptanceProbability(int molKind)
{
	transitionMatrix[GetTMEtcIndex(molLookRef.NumKindInBox(molKind, mv::BOX0))] += 1.0;
}

//Calculates the bias to be applied to GCMC insertion/deletion moves from the transition transitionMatrix.
//Note: Must calculate bias /after/ adding acceptance probabilities
inline double TransitionMatrix::CalculateBias(int move)
{
	//If not running a transition matrix run, don't apply biasing
	if (!biasingOn)
		return 1.0;

	uint numMolecules = molLookRef.NumKindInBox(molKind, mv::BOX0);
	//move == 0: deletion move; move == 1: insertion move
	if (move == 0 && numMolecules != 0)
	{
		return exp(weightingFunction[numMolecules] - weightingFunction[numMolecules - 1]);
	}
	else if (move == 1 && numMolecules != nmax)
	{
		return exp(weightingFunction[numMolecules] - weightingFunction[numMolecules + 1]);
	}
	else 
	{
		return 1.0;
	}
}

//Updates the bias vector using transition matrix data
inline void TransitionMatrix::UpdateWeightingFunction(ulong step)
{
	if (!biasingOn)
		return;

  if((step + 1) % biasStep == 0) {
	  weightingFunction[0] = INITIAL_WEIGHTINGFUNCTION_VALUE;

	  double sumEntry, sumEntryPlusOne, probInsert, probDelete;
	  for (int i = 1; i<weightingFunction.size(); i++) {
		  sumEntry = transitionMatrix[GetTMDelIndex(i - 1)] + transitionMatrix[GetTMEtcIndex(i - 1)] + transitionMatrix[GetTMInsIndex(i - 1)];
		  sumEntryPlusOne = transitionMatrix[GetTMDelIndex(i)] + transitionMatrix[GetTMEtcIndex(i)] + transitionMatrix[GetTMInsIndex(i)];
		  if ((sumEntry > 0.0) && (sumEntryPlusOne > 0.0) && (transitionMatrix[GetTMInsIndex(i - 1)] > 0.0) && (transitionMatrix[GetTMDelIndex(i)] > 0.0)) {
			  probInsert = transitionMatrix[GetTMInsIndex(i - 1)] / sumEntry;
			  probDelete = transitionMatrix[GetTMDelIndex(i)] / sumEntryPlusOne;
		  } else {
			  probInsert = 1.0;
			  probDelete = 1.0;
		  }
		  weightingFunction[i] = weightingFunction[i - 1] + log(probInsert / probDelete);
	  }
	  if(TMfile.is_open()){
		  for (int i = 0; i < nmax; i++) {
			  TMfile << weightingFunction[i] << endl;
		  }
		  TMfile << endl;
	  }
  }
}

inline void TransitionMatrix::PrintTMProbabilityDistribution()
{
	if (!biasingOn)  
		return; 

	UpdateWeightingFunction(biasStep - 1);

	std::cout << "\nTM Particle Number Probability Distribution:\n";
	for (int i = 0; i < nmax; i++) {
		std::cout << weightingFunction[i] << ",";
	}
	std::cout << weightingFunction[nmax] << endl;

	//ofstream TMfile;
	//TMfile.open("TMFile.dat");		//TODO: make this an actual name
	if (TMfile.is_open()) {
		//TMfile << setw(16) << left << "Temperature" << setw(16) << left << temperature << endl;
		//TMfile << setw(16) << left << "Box Volume" << setw(16) << left << boxVolume << endl;
		for (int i = 0; i < nmax; i++) {
			TMfile << weightingFunction[i] << endl;
		}
	}
	
	TMfile.close();
}

#endif