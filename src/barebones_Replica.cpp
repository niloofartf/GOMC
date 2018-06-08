class barebones_Replica {
public:
    int       repl;        /* replica ID */
    int       nrepl;       /* total number of replica */
    float      temp;        /* temperature */
    int      numAtomsInSystem;        /* natoms */
    
    int       type;        /* replica exchange type from ere enum */
    float    **q;           /* quantity, e.g. temperature or lambda; first index is ere, second index is replica ID */
    bool    bNPT;        /* use constant pressure and temperature */
    bool    bNVT;        /* use constant pressure and temperature */
    float     *pres;        /* replica pressures */
    int      *ind;         /* replica indices */
    int      *allswaps;    /* used for keeping track of all the replica swaps */
    int       nst;         /* replica exchange interval (number of steps) */
    int       nex;         /* number of exchanges per interval */
    int       seed;        /* random seed */
    int       nattempt[2]; /* number of even and odd replica change attempts */
    float     *prob_sum;    /* sum of probabilities */
    int     **nmoves;      /* number of moves between replicas i and j */
    int     *nexchange;    /* i-th element of the array is the number of exchanges between replica i-1 and i */
    
    int      *destinations;
    int     **cyclic;
    int     **order;
    int      *tmpswap;
    bool *incycle;
    bool *bEx;
    /* helper arrays to hold the quantities that are exchanged */
    float  *prob;
    float  *Epot;
    float  *beta;
    float  *Vol;
    float **de;
    
    barebones_Replica(  int repl,
                        int nrepl,
                        float temp,
                        int natoms) {
        this->repl = repl;
        this->nrepl = nrepl;
        this->temp = temp;
        this->numAtomsInSystem = natoms;
    }

    void init(int   exchangeInterval
              ){
    
    /* Make an index for increasing replica order */
    /* only makes sense if one or the other is varying, not both!
       if both are varying, we trust the order the person gave. */
        this->ind = (int*)malloc((this->nrepl)*sizeof(int));
        for (int i = 0; i < this->nrepl; i++)
        {
            this->ind[i] = i;
        }
    
/* keep track of all the swaps, starting with the initial placement. */
        this->allswaps = (int*)malloc((this->nrepl)*sizeof(int));
        
        this->prob_sum = (float*)malloc((this->nrepl)*sizeof(float));
        this->nexchange = (int*)malloc((this->nrepl)*sizeof(int));
        this->nmoves = (int**)malloc((this->nrepl)*sizeof(int*));
        

        for (int i = 0; i < this->nrepl; i++){
            this->nmoves[i] = (int*)malloc((this->nrepl)*sizeof(int));
        }
       
        /* generate space for the helper functions so we don't have to snew each time */
 
        this->destinations = (int*)malloc((this->nrepl)*sizeof(int));
        this->incycle      = (bool*)malloc((this->nrepl)*sizeof(bool));
        this->tmpswap      = (int*)malloc((this->nrepl)*sizeof(int));
        this->cyclic = (int**)malloc((this->nrepl)*sizeof(int*));
        this->order = (int**)malloc((this->nrepl)*sizeof(int*));

        for (int i = 0; i < this->nrepl; i++){
            this->cyclic[i] = (int*)malloc((this->nrepl)*sizeof(int));
            this->order[i] = (int*)malloc((this->nrepl)*sizeof(int));
        } 

    
    /* allocate space for the functions storing the data for the replicas */
    /* not all of these arrays needed in all cases, but they don't take
       up much space, since the max size is nrepl**2 */

        this->prob = (float*)malloc((this->nrepl)*sizeof(float));
        this->bEx = (bool*)malloc((this->nrepl)*sizeof(int));
        this->beta = (float*)malloc((this->nrepl)*sizeof(float));
        this->Vol = (float*)malloc((this->nrepl)*sizeof(float));
        this->Epot = (float*)malloc((this->nrepl)*sizeof(float));
        this->de = (float**)malloc((this->nrepl)*sizeof(float*));
        
        for (int i = 0; i < this->nrepl; i++){
            this->de[i] = (float*)malloc((this->nrepl)*sizeof(float));
        }
        this->nst =     exchangeInterval;
    }
};
