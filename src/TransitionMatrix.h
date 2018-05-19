/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.20
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
//TODO: what type does Printer need to be?
//TODO: implement printer
//TODO: find where flags are tracked and add TM flag
//ConfigSetup: config file has 'flags': is that how it works?
//maybe make flag # of steps between weightingFunction updates, w/ -1 'off'?
//Simulation.h has a setup object; maybe initialize TM there?
//TODO: apply AddAccept/CalcBias to MoleculeTransfer.h
//Note: MolTransfer has sourceBox = 0 as delete move, find out where that's set to pass to calcs
//TODO: add bias updating to Simulation.cpp
//TODO: find end of simulation and add postprocess/printer to it

#ifndef TRANSITIONMATRIX_H
#define TRANSITIONMATRIX_H
#include "MoleculeLookup.h"
#include <vector>
#include <cmath>
class TransitionMatrix
{
public:
	TransitionMatrix(MoleculeLookup & molLook) : transitionMatrixDel({ (double)0.0 }), transitionMatrixEtc({ (double)0.0 }),
		transitionMatrixIns({ (double)0.0 }), weightingFunction({ (double)1.0 }), molLookRef(molLook) {}
	void Init();
	void AddAcceptanceProbToMatrix(double acceptanceProbability, int move);
	double CalculateBias(bool isDelMove);
	void UpdateWeightingFunction();
	std::vector<double> PostProcessTransitionMatrix(double beta);
	void TransitionMatrixPrinter();

private:
	std::vector<double> transitionMatrixDel;
	std::vector<double> transitionMatrixEtc;
	std::vector<double> transitionMatrixIns;
	std::vector<double> weightingFunction;
	int vaporPeak, midpoint, liquidPeak;
	uint molKind;							//Track what kind of molecule we're interested in
	MoleculeLookup & molLookRef;			//Used to reference number of molecules of interest in the main box
	bool biasingOn;							//Config flag, turns biasing on or off
};


inline void TransitionMatrix::Init() {}

//Adds acceptance probability from GCMC particle insertion and deletion moves to the transition transitionMatrix
//vector, growing the vector if required.
inline void TransitionMatrix::AddAcceptanceProbToMatrix(double acceptanceProbability, int move)
{
	if (!biasingOn)
		return;
	uint numMolecules = molLookRef.NumKindInBox(molKind, 0);
	if (numMolecules == transitionMatrixDel.size()) {
		transitionMatrixDel.push_back((double)0.0);
		transitionMatrixEtc.push_back((double)0.0);
		transitionMatrixIns.push_back((double)0.0);
		weightingFunction.push_back(weightingFunction[numMolecules - 1]);
	}

	//move == 1: deletion move; move == 2: insertion move; move == anything else: move will keep numMolecules constant
	if (move == 1)
	{
		transitionMatrixDel[numMolecules] += acceptanceProbability;
		transitionMatrixEtc[numMolecules] += 1.0 - acceptanceProbability;
	}
	else if (move == 2)
	{
		transitionMatrixIns[numMolecules] += acceptanceProbability;
		transitionMatrixEtc[numMolecules] += 1.0 - acceptanceProbability;
	}
	else {
		transitionMatrixEtc[numMolecules] += 1.0;
	}
}

//Calculates the bias to be applied to GCMC insertion/deletion moves from the transition transitionMatrix.
//Note: Must calculate bias /after/ adding acceptance probabilities
//TODO: does GOMC increase numMolecules before actually accepting move? If so, change numMolec refs to match prototype
inline double TransitionMatrix::CalculateBias(bool isDelMove)
{
	//If not running a transition matrix run, don't apply biasing to move acceptance chance
	if (!biasingOn)
		return 1.0;

	uint numMolecules = molLookRef.NumKindInBox(molKind, 0);
	if (isDelMove)
	{
		return weightingFunction[numMolecules] / weightingFunction[numMolecules + 1];
	}
	else if (numMolecules != 0)
	{
		return weightingFunction[numMolecules - 1] / weightingFunction[numMolecules];
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
		sumEntry = transitionMatrixDel[i - 1] + transitionMatrixEtc[i - 1] + transitionMatrixIns[i - 1];
		sumEntryPlusOne = transitionMatrixDel[i] + transitionMatrixEtc[i] + transitionMatrixIns[i];

		if (sumEntry > 0 && sumEntryPlusOne > 0 && transitionMatrixIns[i - 1] > 0 && transitionMatrixDel[i] > 0)
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

inline std::vector<double> TransitionMatrix::PostProcessTransitionMatrix(double beta)
{
	if (!biasingOn)
		return;

	int maxMolecules = weightingFunction.size();
	midpoint = maxMolecules / 2;
	double dChemPot = 0.1;
	double change = 0.1;

	std::vector<double> newWeightingFunction;
	newWeightingFunction.reserve(weightingFunction.size());
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
				newWeightingFunction.push_back(weightingFunction[i] + dChemPot * beta*i);
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
		vaporPeak = 0;
		for (int i = 0; i < midpoint; i++) {
			if (newWeightingFunction[i] > newWeightingFunction[vaporPeak]) {
				vaporPeak = i;
			}
		}
		liquidPeak = midpoint;
		for (int i = midpoint; i < maxMolecules + 1; i++) {
			if (newWeightingFunction[i] > newWeightingFunction[liquidPeak]) {
				liquidPeak = i;
			}
		}
		midpoint = vaporPeak;
		for (int i = vaporPeak; i < liquidPeak; i++) {
			if (newWeightingFunction[i] < newWeightingFunction[midpoint]) {
				midpoint = i;
			}
		}

		//Determine new cutoff for max molecules (see: Errington 2003)
		maxMolecules = liquidPeak;
		while (newWeightingFunction[liquidPeak] - newWeightingFunction[maxMolecules] <= 10.0 && maxMolecules <= newWeightingFunction.size()) {
			maxMolecules++;
		}

		infLoopPrevention += 1;
		//Repeatedly solve until peaks and midpoint stabilize; infLoopPrevention prevents oscillation around a point
	} while ((vaporPeak != oldVaporPeak || midpoint != oldMidpoint || liquidPeak != oldLiquidPeak) && infLoopPrevention<1000);

	return newWeightingFunction;
}

inline void TransitionMatrix::TransitionMatrixPrinter()
{
	//write weightingFunction to file
	//Run postprocessing, write to file
}

#endif