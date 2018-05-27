/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.20
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
//TODO: Add TM flag to ConfigSetup.cpp to set biasingOn to true/false, set reweighting by config (hardcoded in simulation.cpp)
//TODO: print weightingFunction to file instead of console/Format console print (see manual, pg 38)
//TODO: find a way to set molKind to the correct kind of molecule (currently only good for single comp liquid/gas interface, no adsorbance)
//TODO: Output information (densities, vapor pressure) in correct units

//Note: When adding new moves to GOMC, if they might run in a GCMC simulation, the transition matrix 
//must be updated with the acceptance probability. If the move will not end in the number of particles in the 
//main simulation box (BOX0) changing, then the acceptance probability can be captured with the following
//code chunk:
//	if ENSEMBLE == GCMC
//		transitionMatrixRef.addAcceptanceProbToMatrix(0.0, 0);
//	endif
//Files currently affected: Movebase.h, MoleculeTransfer.h, IntraSwap.h

#ifndef TRANSITIONMATRIX_H
#define TRANSITIONMATRIX_H
#include "EnsemblePreprocessor.h"
#include "System.h"                 //For init
#include "StaticVals.h"             //For init
#include "Forcefield.h"
#include "MoleculeLookup.h"
#include "BoxDimensions.h"
#include "BoxDimensionsNonOrth.h"
#include <vector>
#include <cmath>
#include <iostream>

class System;
class StaticVals;
class Forcefield;
class MoleculeLookup;
class MoleculeKind;
class BoxDimensions;

class TransitionMatrix
{
public:
	TransitionMatrix(StaticVals const& stat, System & sys) : forcefield(stat.forcefield),
    molLookRef(sys.molLookupRef), currentAxes(sys.boxDimRef) {};
    
	void Init(config_setup::Output out);
	void AddAcceptanceProbToMatrix(double acceptanceProbability, int move);
	double CalculateBias(bool isDelMove);
	void UpdateWeightingFunction();
	void PrintTMProbabilityDistribution();

private:
	const MoleculeLookup& molLookRef;			//Used to reference number of molecules of interest in the main box
    const BoxDimensions& currentAxes;       //Used for volume
    const Forcefield& forcefield;
	bool biasingOn;							//Config flag, turns biasing on or off
    ulong biasStep;                          //Biasing steps
    double boxVolume, temperature;
    int vaporPeak, midpoint, liquidPeak;
    uint molKind;                            //Track what kind of molecule we're interested in; currently defaults to 0 (works with single component systems only)
    std::vector<double> PostProcessTransitionMatrix();
    std::vector<double> transitionMatrixDel;        //Tracks sum of acceptance probabilities of deletion moves
    std::vector<double> transitionMatrixEtc;        //Tracks sum of 1-acceptance probabilities of insert/delete, all attempts of all other moves
    std::vector<double> transitionMatrixIns;        //Tracks sum of acceptance probabilities of insertion moves
    std::vector<double> weightingFunction;            //Holds calculated biasing function

};


inline void TransitionMatrix::Init(config_setup::Output out) {
	biasingOn = out.state.files.tmmc.enable;
    biasStep = out.state.files.tmmc.step;
	molKind = 0;
    boxVolume = currentAxes.GetBoxVolume(mv::BOX0);
	temperature = forcefield.T_in_K;
	transitionMatrixDel.push_back((double)0.0);
	transitionMatrixEtc.push_back((double)0.0);
	transitionMatrixIns.push_back((double)0.0);
	weightingFunction.push_back((double)1.0);
}

//Adds acceptance probability from GCMC particle insertion and deletion moves to the transition transitionMatrix
//vector, growing the vector if required.
inline void TransitionMatrix::AddAcceptanceProbToMatrix(double acceptanceProbability, int move)
{
	if (!biasingOn)
		return;
	uint numMolecules = molLookRef.NumKindInBox(molKind, 0);
	while (numMolecules >= transitionMatrixDel.size()) {
		transitionMatrixDel.push_back((double)0.0);
		transitionMatrixEtc.push_back((double)0.0);
		transitionMatrixIns.push_back((double)0.0);
		weightingFunction.push_back(weightingFunction[weightingFunction.size() - 1]);
	}

	if (acceptanceProbability > 1.0)
		acceptanceProbability = 1.0;

	//move == 1: deletion move; move == 2: insertion move; move == anything else: move will keep numMolecules constant
	if (move == 1)
	{
		transitionMatrixDel[numMolecules] = transitionMatrixDel[numMolecules] + acceptanceProbability;
		transitionMatrixEtc[numMolecules] = transitionMatrixEtc[numMolecules] + 1.0 - acceptanceProbability;
	}
	else if (move == 2)
	{
		transitionMatrixIns[numMolecules] = transitionMatrixIns[numMolecules] + acceptanceProbability;
		transitionMatrixEtc[numMolecules] = transitionMatrixEtc[numMolecules] + 1.0 - acceptanceProbability;
	}
	else {
		transitionMatrixEtc[numMolecules] = transitionMatrixEtc[numMolecules] + 1.0;
	}
}

//Calculates the bias to be applied to GCMC insertion/deletion moves from the transition transitionMatrix.
//Note: Must calculate bias /after/ adding acceptance probabilities
inline double TransitionMatrix::CalculateBias(bool isDelMove)
{
	//If not running a transition matrix run, don't apply biasing to move acceptance chance
	if (!biasingOn)
		return 1.0;

	uint numMolecules = molLookRef.NumKindInBox(molKind, 0);
	if (isDelMove && numMolecules != 0)
	{
		return weightingFunction[numMolecules - 1] / weightingFunction[numMolecules];
	}
	else if (numMolecules != weightingFunction.size()-1)
	{
		return weightingFunction[numMolecules] / weightingFunction[numMolecules + 1];
	}

	else 
	{
		return 1.0;
	}
}

//Updates the bias vector using transition matrix data. Should be done around every 10^6 steps.
inline void TransitionMatrix::UpdateWeightingFunction()
{
	if (!biasingOn)
		return;
	weightingFunction[0] = 1.0;

	double sumEntry, sumEntryPlusOne, probInsert, probDelete;
	for (int i = 1; i<weightingFunction.size(); i++)
	{
		std::cout << i;
		sumEntry = transitionMatrixDel[i - 1] + transitionMatrixEtc[i - 1] + transitionMatrixIns[i - 1];
		sumEntryPlusOne = transitionMatrixDel[i] + transitionMatrixEtc[i] + transitionMatrixIns[i];
		if ((sumEntry > 0.0) && (sumEntryPlusOne > 0.0) && (transitionMatrixIns[i - 1] > 0.0) && (transitionMatrixDel[i] > 0.0))
		{
			probInsert = transitionMatrixIns[i - 1] / sumEntry;
			probDelete = transitionMatrixDel[i] / sumEntryPlusOne;
		}
		else
		{
			probInsert = 1.0;
			probDelete = 1.0;
		}
		weightingFunction[i] = weightingFunction[i - 1] + log(probInsert / probDelete);
	}
}

inline void TransitionMatrix::PrintTMProbabilityDistribution()
{
	if (!biasingOn)  
		return; 
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
	
	//Calc/print vapor/liquid densities here (peaks/volume?)
	std::cout << "Box volume: " << boxVolume << "Vapor peak: " << vaporPeak << " Liquid peak: " << liquidPeak << "\n";

	double sumFunction = 0.0;
	for (int i = 0; i < weightingFunction.size(); i++) {
		sumFunction += weightingFunction[i];
	}

	double KBTimesTenToTheThirtieth = 13806485.2;			//Boltzmann constant times 10^30 (box volume from cubic angstroms to cubic meters)
	double vaporPressure = (log(sumFunction) - log(weightingFunction[0]) - log(2)) * KBTimesTenToTheThirtieth * temperature / boxVolume;
	std::cout << "Vapor pressure: " << vaporPressure << "\n";
}

inline std::vector<double> TransitionMatrix::PostProcessTransitionMatrix()
{
	std::vector<double> newWeightingFunction;
	for (int i = 0; i < weightingFunction.size(); i++) {
		newWeightingFunction.push_back(0.0);
	}
	int maxMolecules = newWeightingFunction.size() - 1;
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
			for (int i = 0; i < weightingFunction.size(); i++) {
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
		for (int i = 0; i < midpoint; i++) {
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
		for (int i = liquidPeak; i < newWeightingFunction.size(); i++){
			if (newWeightingFunction[liquidPeak] - newWeightingFunction[maxMolecules] <= 10.0)
				maxMolecules++;
		}
		std::cout << "";
		infLoopPrevention += 1;
		//Repeatedly solve until peaks and midpoint stabilize; infLoopPrevention prevents oscillation around a point
	} while ((vaporPeak != oldVaporPeak || midpoint != oldMidpoint || liquidPeak != oldLiquidPeak) && infLoopPrevention<100);

	return newWeightingFunction;
}
#endif
