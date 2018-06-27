class SystemPotential;
class Coordinates;
class COM;
class CalcEwald;
class CellList;
class Ewald;
#include "EnergyTypes.h"


class Replica_State
{
    public:
        //! Constructor
        Replica_State(){
          potential = nullptr;
        coordinates = nullptr; //ex
                com = nullptr; //ex
          calcEwald = nullptr; //ex
          cellList = nullptr; //ex
        }  
          SystemPotential* potential; //ex
          Coordinates* coordinates; //ex
          COM* com; //ex
          Ewald* *calcEwald; //ex
          CellList* cellList; //ex
};

