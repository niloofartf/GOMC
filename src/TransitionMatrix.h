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
#include <algorithm>

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

private:
	const MoleculeLookup& molLookRef;	//Used to reference number of molecules of interest in the main box
  const BoxDimensions& currentAxes;     //Used for volume
  const Forcefield& forcefield;			//Used for temperature
	bool biasingOn;						//Config flag, turns biasing on or off
  ulong biasStep;                       //Biasing steps
  double boxVolume, temperature, vaporPressure, surfaceTension, vaporDensity, liquidDensity;
  int INITIAL_WEIGHTINGFUNCTION_VALUE;
  int vaporPeak, midpoint, liquidPeak;
  
  uint molKind;							//Track what kind of molecule we're interested in; currently defaults to 0 
										//(Current implementation works with single component systems only)
										//If user can select a component to track, adsorption/absorption sims become possible

  std::vector<double> PostProcessTransitionMatrix();
  int GetTMDelIndex(int numMolec);
  int GetTMEtcIndex(int numMolec);
  int GetTMInsIndex(int numMolec);

  std::vector<double> transitionMatrix;			  //Holds flat 2D Transition Matrix array after initialization

  //std::vector<double> transitionMatrixDel;        //Tracks sum of acceptance probabilities of deletion moves
  //std::vector<double> transitionMatrixEtc;        //Tracks sum of 1-acceptance probabilities of insert/delete, all attempts of all other moves
  //std::vector<double> transitionMatrixIns;        //Tracks sum of acceptance probabilities of insertion moves
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
	//del/etc/ins, del/etc/ins, ... 
	transitionMatrix.resize(3 * totMolec);
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
	else if (move == 1 && numMolecules != weightingFunction.size()-1)
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
  }
}

inline void TransitionMatrix::PrintTMProbabilityDistribution()
{
	if (!biasingOn)  
		return; 

	UpdateWeightingFunction(biasStep - 1);

	std::cout << "\nTM Particle Number Probability Distribution:\n";
	for (int i = 0; i < weightingFunction.size()-1; i++) {
		std::cout << weightingFunction[i] << ",";
	}
	std::cout << weightingFunction[weightingFunction.size() - 1];
	

	weightingFunction = PostProcessTransitionMatrix();
	
	
	std::cout << "\nEquilibrium Particle Number Probability Distribution:" << "\n";
	for (int i = 0; i < weightingFunction.size()-1; i++) {
		std::cout << weightingFunction[i] << ",";
	}
	std::cout << weightingFunction[weightingFunction.size() - 1] << "\n";
	
	std::cout << "\nVapor peak: " << vaporPeak << "; midpoint: " << midpoint << "; liquid peak: " << liquidPeak;
	std::cout << "\nVapor density: " << vaporDensity << " (mol/m3); Liquid peak: " << liquidDensity << " (mol/m3)";
	std::cout << "\nVapor pressure: " << vaporPressure << " (Pa)";
	std::cout << "\nBox Length (assumes cubic box): " << pow(boxVolume,1/3.0) << " (ang); Surface Tension: " << surfaceTension << " (N/m)\n";
}

inline std::vector<double> TransitionMatrix::PostProcessTransitionMatrix()
{
	std::vector<double> newWeightingFunction(weightingFunction.size());

	int maximumMoleculesSampled = weightingFunction.size() - 1;
	int i = weightingFunction.size() - 2;
	while (weightingFunction[i] == weightingFunction[maximumMoleculesSampled]) {
		maximumMoleculesSampled = i;
		i--;
	}

	int maxMolecules = maximumMoleculesSampled - 1;
	midpoint = maxMolecules / 2;

	double leftArea = 0.0;
	for (int i = 0; i < midpoint; i++) {
		leftArea += weightingFunction[i];
	}
	double rightArea = 0.0;
	for (int i = midpoint; i < maxMolecules; i++) {
		rightArea += weightingFunction[i];
	}
	double areaDif = abs(rightArea - leftArea);

	double oldAreaDif;
	int oldVaporPeak, oldLiquidPeak, oldMidpoint;

	double dChemPot = 0.5;
	double change = 0.1;
	vaporPeak = 0;
	liquidPeak = 0;
	int infLoopPrevention = 0;
	
	do {
		oldVaporPeak = vaporPeak;
		oldMidpoint = midpoint;
		oldLiquidPeak = liquidPeak;

		//Weight macrostate probability function so integral of vapor region == integral of liquid region
		while (areaDif > 0.0001) {
			oldAreaDif = areaDif;
			for (int i = 0; i < maximumMoleculesSampled; i++) {
				newWeightingFunction[i] = weightingFunction[i] + dChemPot*i;
			}

			leftArea = 0.0;
			for (int i = 0; i < midpoint; i++) {
				leftArea += newWeightingFunction[i];
			}
			rightArea = 0.0;
			for (int i = midpoint; i < maxMolecules; i++) {
				rightArea += newWeightingFunction[i];
			}
			areaDif = abs(rightArea - leftArea);

			if (areaDif >= oldAreaDif) {
				change *= -0.5;
			}
			dChemPot += change;
		}

		//Determine new peaks of vapor/liquid regions, midpoint (lowest point between peaks)
		std::cout << "";
		vaporPeak = 0;
		for (int i = 0; i < maxMolecules/2; i++) {
			if (newWeightingFunction[i] > newWeightingFunction[vaporPeak]) {
				vaporPeak = i;
			}
		}
		std::cout << "";
		liquidPeak = midpoint;
		for (int i = midpoint; i < maxMolecules; i++) {
			if (newWeightingFunction[i] > newWeightingFunction[liquidPeak]) {
				liquidPeak = i;
			}
		}
		std::cout << "";
		midpoint = vaporPeak;
		for (int i = vaporPeak; i < liquidPeak; i++) {
			if (newWeightingFunction[i] < newWeightingFunction[midpoint]) {
				midpoint = i;
			}
		}
		std::cout << "";

		//Determine new cutoff for max molecules (see: Errington 2003)
		maxMolecules = liquidPeak;
		for (int i = liquidPeak; i < maximumMoleculesSampled; i++){
			if (newWeightingFunction[liquidPeak] - newWeightingFunction[maxMolecules] <= 10.0)
				maxMolecules++;
		}
		std::cout << "";
		infLoopPrevention += 1;
		//Repeatedly solve until peaks and midpoint stabilize; infLoopPrevention prevents oscillation around a point
	} while ((vaporPeak != oldVaporPeak || midpoint != oldMidpoint || liquidPeak != oldLiquidPeak) && infLoopPrevention<1000);

	//Zero out junk data
	for (int i = maxMolecules + 1; i < newWeightingFunction.size(); i++) { newWeightingFunction[i] = 0.0; }

	//Vapor/liquid-phase densities calculations
	double InvAvogadrosNumTimesTenToTheThirtieth = 1660539.04;		//Inverse Avogadro's Number times 10^30 (box volume from cubic angstroms to cubic meters)
	vaporDensity = InvAvogadrosNumTimesTenToTheThirtieth * vaporPeak / boxVolume;
	liquidDensity = InvAvogadrosNumTimesTenToTheThirtieth * liquidPeak / boxVolume;

	//Vapor Pressure calculation
	double sumFunction = 0.0;
	for (int i = 0; i < maxMolecules; i++) {
		sumFunction += newWeightingFunction[i];
	}
	double KBTimesTenToTheThirtieth = 13806485.2;			//Boltzmann constant times 10^30 (box volume from cubic angstroms to cubic meters)
	vaporPressure = (log(sumFunction) - log(newWeightingFunction[0]) - log(2)) * KBTimesTenToTheThirtieth * temperature / boxVolume;

	//Surface Tension calculation
	double KBTimesTenToTheTwentieth = 0.00138064852;
	surfaceTension = (0.5 * (newWeightingFunction[liquidPeak] + newWeightingFunction[vaporPeak]) - newWeightingFunction[midpoint]) * temperature * KBTimesTenToTheTwentieth / (2 * pow(boxVolume, 2.0 / 3.0));

	return newWeightingFunction;
}
#endif