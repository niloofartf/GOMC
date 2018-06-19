class SystemPotential;
class Coordinates;
class COM;
class CalcEwald;
class CellList;
class Ewald;

class Replica_State
{
    public:
        //! Constructor
        Replica_State(){
            potential = NULL;
          coordinates = NULL; //ex
          com = NULL; //ex
          calcEwald = NULL; //ex
          cellList = NULL; //ex
        }  
          SystemPotential* potential; //ex
          Coordinates* coordinates; //ex
          COM* com; //ex
          Ewald* *calcEwald; //ex
          CellList* cellList; //ex

};

