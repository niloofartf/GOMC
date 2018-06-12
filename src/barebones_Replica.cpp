#include "System.h"

#define KILO (1e3)
#define AVOGADRO         (6.02214129e23)
#define BOLTZMANN (1.3806488e-23)  /* (J/K, NIST 2010 CODATA */
#define RGAS        (BOLTZMANN*AVOGADRO)   /* (J/(mol K))  */
#define BOLTZ   (RGAS/KILO)    /* (kJ/(mol K)) */
#define PROBABILITYCUTOFF 100
#define BAR_MDUNITS      (1e5*NANO*PICO*PICO/AMU)
#define PRESFAC          (1.0/BAR_MDUNITS)
#define NANO             (1e-9)                            /* A Number  */
#define PICO             (1e-12)                           /* A Number  */
#define AMU              (1.660538921e-27)                 /* kg, NIST 2010 */

class barebones_Replica {
public:
    
enum {
    ereTEMP, erePRESSURE
};
//    System* theSystem;
    int       repl;        /* replica ID */
    int       nrepl;       /* total number of replica */
    float      temp;        /* temperature */
    int      numAtomsInSystem;        /* natoms */

    int i, j, k;
    
    int       type;        /* replica exchange type from ere enum */
    std::vector<double>  replica_temps;           /* temperatures of replicas */
    
    float    **q;           /* quantity, e.g. temperature or lambda; first index is ere, second index is replica ID */
    bool    bNPT;        /* use constant pressure and temperature */
    bool    bVol;        /* use constant pressure and temperature */
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
                        int natoms,
                        std::vector<double> replica_temps){
        
//        this->theSystem = sys;

        // For now, since I'm assuming NVT
        this->bNPT = false;

        this->replica_temps = replica_temps;

        this->repl = repl;
        this->nrepl = nrepl;
        this->temp = temp;
        this->numAtomsInSystem = natoms;
    }

    void init(int   exchangeInterval,
              int   randomSeed)
    {
    
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

        if (this->seed == UINT_MAX){
            // NOTE : Ask Md and Younes how they autoseed            
        } else {
            this->seed = randomSeed;
        }

        printf("\nReplic exchange interval: %d\n", this->nst);
        printf("\nReplica random seed: %d\n", this->seed);

        this->nattempt[0] = 0;
        this->nattempt[1] = 0;

        printf("Replica exchange information below: ex and x = exchange, pr = probability\n");


    }

    //bool replica_exchange(const float energy, PRNG prng, float vol){
    bool replica_exchange(const float energy, uint step){

        int j;
        int replica_id = 0;
        int exchange_partner;
        int maxswap = 0;

        /* Number of rounds of exchanges needed to deal with any multiple
        * exchanges. */
        /* Where each replica ends up after the exchange attempt(s). */
        /* The order in which multiple exchanges will occur. */
    bool bThisReplicaExchanged = false;

    
    // if master thread
        replica_id = this->repl;
        test_for_replica_exchange(energy, step);
        //prepare_to_do_exchange();

        if (bThisReplicaExchanged){
    
            /* There will be only one swap cycle with standard replica
             * exchange, but there may be multiple swap cycles if we
             * allow multiple swaps. */

            for (j = 0; j < maxswap; j++)
            {
                exchange_partner = this->order[replica_id][j];

                if (exchange_partner != replica_id)
                {
                    bool debug = true;
                    /* Exchange the global states between the master nodes */
                    if (debug)
                    {
                        printf("Exchanging %d with %d\n", replica_id, exchange_partner);
                    }
                    exchange_state(exchange_partner);
                }
            }
            /* For temperature-type replica exchange, we need to scale
             * the velocities. */
                scale_velocities();

 
        }
    }

    static void exchange_state(int exchange_partner){}
    void scale_velocities(){}
    void test_for_replica_exchange(const float energy, uint step){

        std::cout << "BATES I am replica number " << this->repl << "and my epot on step " << step << " is " << energy << std::endl;
        return;

        int         m, i, j, a, b, ap, bp, i0, i1, tmp;
        float       delta   = 0.0;
        bool        bPrint, bMultiEx;
        bool        *bEx    = this->bEx;
        float       *prob   = this->prob;
        int         *pind   = this->destinations;   /* Permuted Index */
        bool        bEpot   = false;
        bool        bEDLambda   = false;
        bool        bEVol   = false;
        
        bMultiEx = (this->nex > 1);
        
        if (this-> bNPT){
            for (i = 0; i < this->nrepl ; i++){
                this->Vol[i] = 0.0;
            }
            this->bVol = true;
// NVT only for now
//            this->Vol[this->repl] = vol;
        }

        // ereTEMP
        for (i = 0; i < this->nrepl; i++){
            this->Epot[i] = 0.0;
        }

        bEpot   = true;

        this->Epot[this->repl] = energy; //  POTENTIAL ENERGY
        
        for (i = 0; i < this->nrepl; i++) {
            this->beta[i] = 1.0/(replica_temps[i]*BOLTZ);
        }

        // endERETEMP

    
        if (bVol) {
            // broadcast volumes
        }

        if (bEpot) {
            // broadcast pot eng
        }
    
        // Do some thread com to populate global arrays

        // SPACE HOLDER
        // SPACE HOLDER
        // SPACE HOLDER

        // Standard nearest neighbor replica exchange
/**/
        m = (step / this->nst) % 2;
        for (i = 1; i < this->nrepl; i++)
        {
            a = this->ind[i-1];
            b = this->ind[i];
            bPrint = (this->repl == a || this->repl == b);
            if (i % 2 == m)
            {
                delta = calc_delta(a, b, a, b);
                
                if (delta <= 0)
                {
                    // accepted 
                    prob[i] = 1;
                    bEx[i]  = true;
                }
                else
                {
                    if (delta > PROBABILITYCUTOFF)
                    {
                        prob[i] = 0;
                    }
                    else
                    {
                        prob[i] = exp(-delta);
                    }
                    // roll a number to determine if accepted. For now it is superfluous to
                    // reset, but just in case we ever add more calls in different branches
                    // it is safer to always reset the distribution.
        //            uniformRealDist.reset();
          //          real randVar = uniformRealDist(rng);
                    // generate a random number
                    double randVar = ((double) rand() / (RAND_MAX)) + 1;
                    bEx[i] = randVar < prob[i];
                }
                this->prob_sum[i] += prob[i];

                if (bEx[i])
                {
                    // swap these two 
                    tmp       = pind[i-1];
                    pind[i-1] = pind[i];
                    pind[i]   = tmp;
                    this->nexchange[i]++;  // statistics for back compatibility 
                }
            }
            else
            {
                prob[i] = -1;
                bEx[i]  = false;
            }
        }
        // print some statistics 
//        print_ind(fplog, "ex", this->nrepl, this->ind, bEx);
  //      print_prob(fplog, "pr", this->nrepl, prob);
    //    fprintf(fplog, "\n");
        this->nattempt[m]++;

    // record which moves were made and accepted 
    for (i = 0; i < this->nrepl; i++)
    {
        this->nmoves[this->ind[i]][pind[i]] += 1;
        this->nmoves[pind[i]][this->ind[i]] += 1;
    }
//    fflush(fplog); // make sure we can see what the last exchange was 
    }
    void prepare_to_do_exchange(    const int           replica_id,
                                    int                *maxswap,
                                    bool           *bThisReplicaExchanged){

        int i, j;
        /* Hold the cyclic decomposition of the (multiple) replica
         * exchange. */
        bool bAnyReplicaExchanged = false;
        *bThisReplicaExchanged = false;

        for (i = 0; i < this->nrepl; i++)
        {
            if (this->destinations[i] != this->ind[i])
            {
                /* only mark as exchanged if the index has been shuffled */
                bAnyReplicaExchanged = true;
                break;
            }
        }
        if (bAnyReplicaExchanged)
        {
            /* reinitialize the placeholder arrays */
            for (i = 0; i < this->nrepl; i++)
            {
                for (j = 0; j < this->nrepl; j++)
                {
                    this->cyclic[i][j] = -1;
                    this->order[i][j]  = -1;
                }
            }

            /* Identify the cyclic decomposition of the permutation (very
             * fast if neighbor replica exchange). */
            cyclic_decomposition(this->destinations, this->cyclic, this->incycle, this->nrepl, maxswap);

            /* Now translate the decomposition into a replica exchange
             * order at each step. */
            compute_exchange_order(this->cyclic, this->order, this->nrepl, *maxswap);

            /* Did this replica do any exchange at any point? */
            for (j = 0; j < *maxswap; j++)
            {
                if (replica_id != this->order[replica_id][j])
                {
                    *bThisReplicaExchanged = true;
                    break;
                }
            }
        }

    }

    double calc_delta(int a, int b, int ap, int bp){

        float   ediff, dpV, delta = 0;
        float  *Epot = this->Epot;
        float  *Vol  = this->Vol;
        float **de   = this->de;
        float  *beta = this->beta;

        /* Two cases; we are permuted and not.  In all cases, setting ap = a and bp = b will reduce
           to the non permuted case */

        /*
         * Okabe et. al. Chem. Phys. Lett. 335 (2001) 435-439
         */
       
        ediff = Epot[b] - Epot[a];
        delta = -(beta[bp] - beta[ap])*ediff;
       
        // remove this boolean later - GJS 
        bool bPrint = true;
        if (bPrint)
        {
        //    fprintf(fplog, "Repl %d <-> %d  dE_term = %10.3e (kT)\n", a, b, delta);
        }
        if (this->bNPT)
        {
            /* revist the calculation for 5.0.  Might be some improvements. */
            dpV = (beta[ap]*this->pres[ap]-beta[bp]*this->pres[bp])*(Vol[b]-Vol[a])/PRESFAC;
            if (bPrint)
            {
          //      fprintf(fplog, "  dpV = %10.3e  d = %10.3e\n", dpV, delta + dpV);
            }
            delta += dpV;
        }
        return delta;
    }

    static void
    cyclic_decomposition(const int *destinations,
                         int      **cyclic,
                         bool  *incycle,
                         const int  nrepl,
                         int       *nswap)
    {

        int i, j, c, p;
        int maxlen = 1;
        for (i = 0; i < nrepl; i++)
        {
            incycle[i] = false;
        }
        for (i = 0; i < nrepl; i++)  /* one cycle for each replica */
        {
            if (incycle[i])
            {
                cyclic[i][0] = -1;
                continue;
            }
            cyclic[i][0] = i;
            incycle[i]   = true;
            c            = 1;
            p            = i;
            for (j = 0; j < nrepl; j++) /* potentially all cycles are part, but we will break first */
            {
                p = destinations[p];    /* start permuting */
                if (p == i)
                {
                    cyclic[i][c] = -1;
                    if (c > maxlen)
                    {
                        maxlen = c;
                    }
                    break; /* we've reached the original element, the cycle is complete, and we marked the end. */
                }
                else
                {
                    cyclic[i][c] = p;  /* each permutation gives a new member of the cycle */
                    incycle[p]   = true;
                    c++;
                }
            }
        }
        *nswap = maxlen - 1;
        // remove this boolean later - GJS
        bool debug = true;
        if (debug)
        {
                for (i = 0; i < nrepl; i++)
                {
            //        fprintf(debug, "Cycle %d:", i);
                    for (j = 0; j < nrepl; j++)
                    {
                        if (cyclic[i][j] < 0)
                        {
                            break;
                        }
                        //fprintf(debug, "%2d", cyclic[i][j]);
                    }
                    //fprintf(debug, "\n");
                }
                //fflush(debug);
        }
    }

    static void
    compute_exchange_order(int     **cyclic,
                           int     **order,
                           const int nrepl,
                           const int maxswap)
    {
        int i, j;

        for (j = 0; j < maxswap; j++)
        {
            for (i = 0; i < nrepl; i++)
            {
                if (cyclic[i][j+1] >= 0)
                {
                    order[cyclic[i][j+1]][j] = cyclic[i][j];
                    order[cyclic[i][j]][j]   = cyclic[i][j+1];
                }
            }
            for (i = 0; i < nrepl; i++)
            {
                if (order[i][j] < 0)
                {
                    order[i][j] = i; /* if it's not exchanging, it should stay this round*/
                }
            }
        }
/*
           if (debug)
            {
                fprintf(debug, "Replica Exchange Order\n");
                for (i = 0; i < nrepl; i++)
                {
                    fprintf(debug, "Replica %d:", i);
                    for (j = 0; j < maxswap; j++)
                    {
                        if (order[i][j] < 0)
                        {
                            break;
                        }
                        fprintf(debug, "%2d", order[i][j]);
                    }
                    fprintf(debug, "\n");
                }
                fflush(debug);
            } */
    }


};
