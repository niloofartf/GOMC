/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef MOVE_SETTINGS_H
#define MOVE_SETTINGS_H

#include "EnsemblePreprocessor.h" //For BOX_TOTAL
#include "BasicTypes.h"           //for uint
#include "OutputVars.h"
#include "PDBSetup.h" //Primary source of volume.
#include "MoveConst.h"           //For sizes of arrays.


class StaticVals;                 //For various initialization constants.
class BoxDimensions;              //For axis sizes
//class barebones_Replica;              //For re

class barebones_Replica;              //For replica exchange

class MoveSettings
{
public:
  friend class OutputVars;
  MoveSettings(BoxDimensions & dim) : boxDimRef(dim) {}

  MoveSettings& operator=(MoveSettings const& rhs)
  {
    return *this;
  }

  void Init(StaticVals const& statV, pdb_setup::Remarks const& remarks);

  void Update(const bool isAccepted, const uint moveIndex, const uint step);

  void AdjustMoves(const uint step);

// GJS
  //void ExchangeMoves(const uint step, barebones_Replica* re, const float energy);
  bool ExchangeMoves(const uint step);
  //void Exchange(barebones_Replica* re, const float energy, uint step);
// GJS

  void Adjust(const uint majMoveKind, const uint moveIndex, const uint b);

  double Scale(const uint move) const
  {
    return scale[move];
  }
  uint perExchange;

private:
  double scale[mv::SCALEABLE];
  double acceptPercent[mv::COUNT];
  uint accepted[mv::COUNT];
  uint tries[mv::COUNT];
  uint perAdjust;
// GJS
// GJS
  uint tempAccepted[mv::SCALEABLE], tempTries[mv::SCALEABLE];

#if ENSEMBLE == GEMC
  uint GEMC_KIND;
#endif

  BoxDimensions & boxDimRef;

  static const double TARGET_ACCEPT_FRACT;
  static const double TINY_AMOUNT;
};


#endif /*MOVE_SETTINGS_H*/
