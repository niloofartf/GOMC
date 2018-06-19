class SystemPotential;
class Coordinates;
class COM;
class CalcEwald;
class CellList;

class t_state
{
    public:
        //! Constructor
        t_state();
  
          SystemPotential* potential; //ex
          Coordinates* coordinates; //ex
          COM* com; //ex
          Ewald *calcEwald; //ex
          CellList* cellList; //ex

};

