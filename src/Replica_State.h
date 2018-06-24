class SystemPotential;
class Coordinates;
class COM;
class CalcEwald;
class CellList;
class Ewald;
class BoxDimensions;

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
#if ENSEMBLE == NPT
          boxDimensions = NULL;
#endif
        }  
          SystemPotential* potential; //ex
          Coordinates* coordinates; //ex
          COM* com; //ex
          Ewald* *calcEwald; //ex
          CellList* cellList; //ex
#if ENSEMBLE == NPT
          BoxDimensions* *boxDimensions;
#endif

};

