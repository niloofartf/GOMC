

CellList& operator=(CellList const& rhs){

// Use the library vector copy
std::vector<int> list;
std::vector<std::vector<int> > neighbors[BOX_TOTAL];
std::vector<int> head[BOX_TOTAL];

list = rhs.list;

for (int i = 0; i < rhs.BOX_TOTAL; i++){
    neighbors[i] = rhs.neighbors[i];
    head[i]      = rhs.head[i];
    cellSize[i]  = rhs.cellSize[i];
    for (int j = 0; j < 3; j++){
        edgeCells[i][j] = rhs.edgeCells[i][j];
    }
}

// Molecules is kinda complex and doesnt have a cc, rats!
// Here is where molecules cc goes

// I think this will work
dimensions = rhs.dimensions;

cutoff = rhs.cutoff;
isBuilt = rhs.isBuilt;

// Idk re initialize this or maybe its static and can be left alone
const Molecules* mols;


}
