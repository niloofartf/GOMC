/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2011,2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

//#include "gmxpre.h"

#include "Repl_Ex.h"
#include <omp.h>

//#include "config.h"

#include <cmath>

#include <random>
#include <string>
//#include "gromacs/domdec/domdec.h"
//#include "gromacs/gmxlib/network.h"
//#include "gromacs/math/units.h"
//#include "gromacs/math/vec.h"
//#include "gromacs/mdlib/main.h"
//#include "gromacs/mdtypes/commrec.h"
//#include "gromacs/mdtypes/inputrec.h"
//#include "gromacs/mdtypes/md_enums.h"
//#include "gromacs/mdtypes/state.h"
//#include "gromacs/random/threefry.h"
//#include "gromacs/random/uniformintdistribution.h"
//#include "gromacs/random/uniformfloatdistribution.h"
//#include "gromacs/utility/fatalerror.h"
//#include "gromacs/utility/pleasecite.h"
//#include "gromacs/utility/smalloc.h"

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

#define PROBABILITYCUTOFF 100
/* we don't bother evaluating if events are more rare than exp(-100) = 3.7x10^-44 */

// Rank in the multisimulation
#define MSRANK(ms, nodeid)  (nodeid)

enum {
    ereTEMP, ereLAMBDA, ereENDSINGLE, ereTL, ereNR
};
const char *erename[ereNR] = { "temperature", "lambda", "end_single_marker", "temperature and lambda"};
/* end_single_marker merely notes the end of single variable replica exchange. All types higher than
   it are multiple replica exchange methods */
/* Eventually, should add 'pressure', 'temperature and pressure', 'lambda_and_pressure', 'temperature_lambda_pressure'?;
   Let's wait until we feel better about the pressure control methods giving exact ensembles.  Right now, we assume constant pressure  */

typedef struct gmx_repl_ex
{
    int       repl;        /* replica ID */
    int       nrepl;       /* total number of replica */
    float     temp;       /* temperature */
    int       type;        /* replica exchange type from ere enum */
    float    **q;          /* quantity, e.g. temperature or lambda; first index is ere, second index is replica ID */
    int       bNPT;             /* use constant pressure and temperature */
    double     *pres;       /* replica pressures */
    int      *ind;         /* replica indices */
    int      *allswaps;    /* used for keeping track of all the replica swaps */
    int       nst;         /* replica exchange interval (number of steps) */
    int       nex;         /* number of exchanges per interval */
    int       seed;        /* random seed */
    int       nattempt[2]; /* number of even and odd replica change attempts */
    float     *prob_sum;   /* sum of probabilities */
    int     **nmoves;      /* number of moves between replicas i and j */
    int     *nexchange;    /* i-th element of the array is the number of exchanges between replica i-1 and i */
    int      natoms;       /* GJS - REDS - 2/14/18 */

    /* these are helper arrays for replica exchange; allocated here so they
       don't have to be allocated each time */
    int      *destinations;
    int     **cyclic;
    int     **order;
    int      *tmpswap;
    int *incycle;
    int *bEx;

    /* helper arrays to hold the quantities that are exchanged */
    float  *prob;
    float  *Epot;
    float  *beta;
    float  *Vol;
    float **de;

} t_gmx_repl_ex;


static int repl_quantity(struct gmx_repl_ex *re, int ere, float q, ReplicaExchangeParameters* replExParams)
{
    float    *qall;
    int bDiff;
    int      s;

    qall = (float*)malloc(omp_get_num_threads()*sizeof(float));
    qall[re->repl] = q;

//  Fill qall with threads info
//  Broadcast
    for (int i = 0; i < re->nrepl; i++){
        if (i == re->repl){
            if ( q != replExParams->replica_temps[i] ) { 
               std::cout << "Error: There is an issue with assigning replica temperatures\n";
               exit(EXIT_FAILURE);
            } else {
                continue;
            }
        }
        qall[i] = replExParams->replica_temps[i];
    }


    bDiff = false;
    for (s = 1; s < omp_get_num_threads(); s++)
    {
        if (qall[s] != qall[0])
        {
            bDiff = true;
        }
    }

    if (bDiff)
    {
        /* Set the replica exchange type and quantities */
        re->type = ere;
        
        re->q[ere] = (float*)malloc(sizeof(float*)*re->nrepl);
        for (s = 0; s < omp_get_num_threads(); s++)
        {
            re->q[ere][s] = qall[s];
        }
    }
    free(qall);
    return bDiff;
}

gmx_repl_ex_t
init_replica_exchange(FILE                            *fplog,
                      float                             temp,
                      ReplicaExchangeParameters* replExParams)
{
    double                pres;
    int                 i, j, k;
    struct gmx_repl_ex *re;
    int            bTemp;
    int            bLambda = false;
    

    fprintf(fplog, "\nInitializing Replica Exchange\n");
    printf("\nInitializing Replica Exchange\n");

    if (omp_get_num_threads() == 1)
    {
           std::cout << "Error: Nothing to exchange with only one replica, maybe you forgot to set the +pN option of GOMC_XPU_NVT?\n";
           exit(EXIT_FAILURE);
    }

    re = (gmx_repl_ex*)malloc(sizeof(gmx_repl_ex));

    re->bNPT     = 0;
    re->repl     = omp_get_thread_num();
    re->nrepl    = omp_get_num_threads();
    re->q = (float**)malloc(sizeof(float*)*ereENDSINGLE);
    fprintf(fplog, "Repl  There are %d replicas:\n", re->nrepl);

    /* We only check that the number of atoms in the systms match.
     * This, of course, do not guarantee that the systems are the same,
     * but it does guarantee that we can perform replica exchange.
     */
    const int nst = replExParams->exchangeInterval;
    re->bNPT = replExParams->bNPT;

    re->temp = temp;

    // Does GOMC use temperature coupling groups???

    re->type = -1;
    bTemp    = repl_quantity(re, ereTEMP, re->temp, replExParams);
    
    if (re->type == -1)  /* nothing was assigned */
    {
        printf("Error: The properties of the %d systems are all the same, there is nothing to exchange\n", re->nrepl);;
             exit(EXIT_FAILURE);
    }


    if (re->bNPT){
        re->pres = (double*)malloc((re->nrepl)*sizeof(double));
        if (replExParams->replica_temps.size() != replExParams->replica_pressures.size()){
            printf("Error : Number of temperatures %lu does not equal number of pressures %lu", replExParams->replica_temps.size(), replExParams->replica_pressures.size());
            exit(EXIT_FAILURE);
        }
      
        #pragma omp critical
        { 
            for ( int i = 0; i < re->nrepl; i++ ){
                re->pres[i] = replExParams->replica_pressures[i];
            }
        }
        #pragma omp barrier
    }

    

    /* Make an index for increasing replica order */
    /* only makes sense if one or the other is varying, not both!
       if both are varying, we trust the order the person gave. */
    
    re->ind = (int*)malloc((re->nrepl)*sizeof(int));
    for (int i = 0; i < re->nrepl; i++)
    {
        re->ind[i] = i;
    }

    if (re->type < ereENDSINGLE)
    {

        for (i = 0; i < re->nrepl; i++)
        {
            for (j = i+1; j < re->nrepl; j++)
            {
                if (re->q[re->type][re->ind[j]] < re->q[re->type][re->ind[i]])
                {
                    /* Unordered replicas are supposed to work, but there
                     * is still an issues somewhere.
                     * Note that at this point still re->ind[i]=i.
                     */
                    printf("Error : Replicas with indices %d < %d have %ss %g > %g, please order your replicas on increasing %s",
                              i, j,
                              erename[re->type],
                              re->q[re->type][i], re->q[re->type][j],
                              erename[re->type]);

                    k          = re->ind[i];
                    re->ind[i] = re->ind[j];
                    re->ind[j] = k;
                    exit(EXIT_FAILURE);
                    
                }
                else if (re->q[re->type][re->ind[j]] == re->q[re->type][re->ind[i]])
                {
                    printf("Error: Two replicas have identical %ss", erename[re->type]);
                    exit(EXIT_FAILURE);
                }
            }
        }
    }

    /* keep track of all the swaps, starting with the initial placement. */
    re->allswaps = (int*)malloc((re->nrepl)*sizeof(int));
    for (i = 0; i < re->nrepl; i++)
    {
        re->allswaps[i] = re->ind[i];
    }

    switch (re->type)
    {
        case ereTEMP:
            printf("\nReplica exchange in temperature\n");
            fprintf(fplog, "\nReplica exchange in temperature\n");
            for (i = 0; i < re->nrepl; i++)
            {
                fprintf(fplog, " %5.1f", re->q[re->type][re->ind[i]]);
            }
            fprintf(fplog, "\n");
            break;
        default:
            printf("Unknown replica exchange quantity");
            exit(EXIT_FAILURE);
            
    }
    
    if (re->bNPT)
    {
        fprintf(fplog, "\nRepl  p");
        for (i = 0; i < re->nrepl; i++)
        {
            fprintf(fplog, " %5.2f", re->pres[re->ind[i]]);
        }

        for (i = 0; i < re->nrepl; i++)
        {
            if ((i > 0) && (re->pres[re->ind[i]] < re->pres[re->ind[i-1]]))
            {
                fprintf(fplog, "\nWARNING: The reference pressures decrease with increasing temperatures\n\n");
                fprintf(stderr, "\nWARNING: The reference pressures decrease with increasing temperatures\n\n");
            }
        }
    }


    re->nst = replExParams->exchangeInterval;

    re->seed = replExParams->randomSeed;


    fprintf(fplog, "\nReplica exchange interval: %d\n", re->nst);
    fprintf(fplog, "\nReplica random seed: %d\n", re->seed);

    re->nattempt[0] = 0;
    re->nattempt[1] = 0;

    re->prob_sum = (float*)malloc((re->nrepl)*sizeof(float));
    re->nexchange = (int*)malloc((re->nrepl)*sizeof(int));
    re->nmoves = (int**)malloc((re->nrepl)*sizeof(int*));
    
    for (i = 0; i < re->nrepl; i++)
    {
        re->nmoves[i] = (int*)malloc((re->nrepl)*sizeof(int));
    }
    fprintf(fplog, "Replica exchange information below: ex and x = exchange, pr = probability\n");

    /* generate space for the helper functions so we don't have to snew each time */


    re->destinations = (int*)malloc((re->nrepl)*sizeof(int));
    re->incycle      = (int*)malloc((re->nrepl)*sizeof(bool));
    re->tmpswap      = (int*)malloc((re->nrepl)*sizeof(int));
    re->cyclic = (int**)malloc((re->nrepl)*sizeof(int*));
    re->order = (int**)malloc((re->nrepl)*sizeof(int*));

    for (int i = 0; i < re->nrepl; i++){
        re->cyclic[i] = (int*)malloc((re->nrepl)*sizeof(int));
        re->order[i] = (int*)malloc((re->nrepl)*sizeof(int));
    } 

    /* allocate space for the functions storing the data for the replicas */
    /* not all of these arrays needed in all cases, but they don't take
       up much space, since the max size is nrepl**2 */

    re->prob = (float*)malloc((re->nrepl)*sizeof(float));
    re->bEx = (int*)malloc((re->nrepl)*sizeof(int));
    re->beta = (float*)malloc((re->nrepl)*sizeof(float));
    re->Vol = (float*)malloc((re->nrepl)*sizeof(float));
    re->Epot = (float*)malloc((re->nrepl)*sizeof(float));

    #pragma omp barrier

    #pragma omp single 
    { 
        replExParams->replica_states = (Replica_State**)malloc(sizeof(Replica_State*)*re->nrepl);
        for (int i = 0; i < re->nrepl; i++){
            replExParams->replica_energies.push_back(0.0);
        }
    }

    re->de = (float**)malloc((re->nrepl)*sizeof(float*));
    
    for (int i = 0; i < re->nrepl; i++){
        re->de[i] = (float*)malloc((re->nrepl)*sizeof(float));

    }
    
    for (int i = 0; i < re->nrepl; i++){
        re->nexchange[i] = 0;
    }
    
    for (int i = 0; i < re->nrepl; i++){
        for (int j = 0; j < re->nrepl; j++){
            re->nmoves[i][j] = 0;
        }
    }

 
    re->nex = replExParams->numExchanges;
    printf("I finished intializing\n");
    return re;
}

/*
static void exchange_floats(const gmx_multisim_t gmx_unused *ms, int gmx_unused b, float *v, int n)
{
    float *buf;
    int   i;

    if (v)
    {
        snew(buf, n);
        for (i = 0; i < n; i++)
        {
            v[i] = buf[i];
        }
        sfree(buf);
    }
}

static void exchange_zetas(const gmx_multisim_t gmx_unused *ms, int gmx_unused b, float *v, int n)
{
    float *buf;
    int   i;

    if (v)
    {
        snew(buf, n);
        for (i = 0; i < n; i++)
        {
            v[i] = buf[i];
        }
        sfree(buf);
    }
}

static void exchange_doubles(const gmx_multisim_t gmx_unused *ms, int gmx_unused b, double *v, int n)
{
    double *buf;
    int     i;

    if (v)
    {
        snew(buf, n);
        for (i = 0; i < n; i++)
        {
            v[i] = buf[i];
        }
        sfree(buf);
    }
}
static void exchange_rvecs(const gmx_multisim_t gmx_unused *ms, int gmx_unused b, rvec *v, int n)
{
    rvec *buf;
    int   i;

    if (v)
    {
        snew(buf, n);
        for (i = 0; i < n; i++)
        {
            copy_rvec(buf[i], v[i]);
        }
        sfree(buf);
    }
}
*/

static void copy_state_serial(const Replica_State *src, Replica_State *dest)
{
    if (dest != src)
    {
        // Currently the local state is always a pointer to the global
        //  in serial, so we should never end up here.
        // TODO: Implement a (trivial) Replica_State copy once converted to C++.
        //
    }
}

/*
static void scale_velocities(Replica_State *state, float fac)
{
    int i;

    if (as_rvec_array(state->v.data()))
    {
        for (i = 0; i < state->natoms; i++)
        {
            svmul(fac, state->v[i], state->v[i]);
        }
    }
}
*/
static void print_transition_matrix(FILE *fplog, int n, int **nmoves, int *nattempt)
{
    int   i, j, ntot;
    float Tprint;

    ntot = nattempt[0] + nattempt[1];
    printf("YOGA %d\n", ntot);
    fprintf(fplog, "\n");
    fprintf(fplog, "Repl");
    for (i = 0; i < n; i++)
    {
        fprintf(fplog, "    ");  /* put the title closer to the center */
    }
    fprintf(fplog, "Empirical Transition Matrix\n");

    fprintf(fplog, "Repl");
    for (i = 0; i < n; i++)
    {
        fprintf(fplog, "%8d", (i+1));
    }
    fprintf(fplog, "\n");

    for (i = 0; i < n; i++)
    {
        fprintf(fplog, "Repl");
        for (j = 0; j < n; j++)
        {
            Tprint = 0.0;
            if (nmoves[i][j] > 0)
            {
                Tprint = nmoves[i][j]/(2.0*ntot);
            }
            fprintf(fplog, "%8.4f", Tprint);
        }
        fprintf(fplog, "%3d\n", i);
    }
}

static void print_ind(FILE *fplog, const char *leg, int n, int *ind, int *bEx)
{
    int i;

    fprintf(fplog, "Repl %2s %2d", leg, ind[0]);
    for (i = 1; i < n; i++)
    {
        fprintf(fplog, " %c %2d", (bEx != nullptr && bEx[i]) ? 'x' : ' ', ind[i]);
    }
    fprintf(fplog, "\n");
}

static void print_allswitchind(FILE *fplog, int n, int *pind, int *allswaps, int *tmpswap)
{
    int i;

    for (i = 0; i < n; i++)
    {
        tmpswap[i] = allswaps[i];
    }
    for (i = 0; i < n; i++)
    {
        allswaps[i] = tmpswap[pind[i]];
    }

    fprintf(fplog, "\nAccepted Exchanges:   ");
    for (i = 0; i < n; i++)
    {
        fprintf(fplog, "%d ", pind[i]);
    }
    fprintf(fplog, "\n");

    /* the "Order After Exchange" is the state label corresponding to the configuration that
       started in state listed in order, i.e.

       3 0 1 2

       means that the:
       configuration starting in simulation 3 is now in simulation 0,
       configuration starting in simulation 0 is now in simulation 1,
       configuration starting in simulation 1 is now in simulation 2,
       configuration starting in simulation 2 is now in simulation 3
     */
    fprintf(fplog, "Order After Exchange: ");
    for (i = 0; i < n; i++)
    {
        fprintf(fplog, "%d ", allswaps[i]);
    }
    fprintf(fplog, "\n\n");
}

static void print_prob(FILE *fplog, const char *leg, int n, float *prob)
{
    int  i;
    char buf[8];

    fprintf(fplog, "Repl %2s ", leg);
    for (i = 1; i < n; i++)
    {
        if (prob[i] >= 0)
        {
            sprintf(buf, "%4.2f", prob[i]);
            fprintf(fplog, "  %3s", buf[0] == '1' ? "1.0" : buf+1);
        }
        else
        {
            fprintf(fplog, "     ");
        }
    }
    fprintf(fplog, "\n");
}

static void print_count(FILE *fplog, const char *leg, int n, int *count)
{
    int i;

    fprintf(fplog, "Repl %2s ", leg);
    for (i = 1; i < n; i++)
    {
        fprintf(fplog, " %4d", count[i]);
    }
    fprintf(fplog, "\n");
}

static float calc_delta(FILE *fplog, int bPrint, struct gmx_repl_ex *re, int a, int b, int ap, int bp)
{
    float   ediff, dpV, delta = 0;
    float  *Epot = re->Epot;
    float  *Vol  = re->Vol;
    float **de   = re->de;
    float  *beta = re->beta;
    
    /* Two cases; we are permuted and not.  In all cases, setting ap = a and bp = b will reduce
       to the non permuted case */

    switch (re->type)
    {
        case ereTEMP:
            /*
             * Okabe et. al. Chem. Phys. Lett. 335 (2001) 435-439
             */

            ediff = Epot[b] - Epot[a];
    
            delta = -(beta[bp] - beta[ap])*ediff;

            break;
        case ereLAMBDA:
            /* two cases:  when we are permuted, and not.  */
            /* non-permuted:
               ediff =  E_new - E_old
                     =  [H_b(x_a) + H_a(x_b)] - [H_b(x_b) + H_a(x_a)]
                     =  [H_b(x_a) - H_a(x_a)] + [H_a(x_b) - H_b(x_b)]
                     =  de[b][a] + de[a][b] */

            /* permuted:
               ediff =  E_new - E_old
                     =  [H_bp(x_a) + H_ap(x_b)] - [H_bp(x_b) + H_ap(x_a)]
                     =  [H_bp(x_a) - H_ap(x_a)] + [H_ap(x_b) - H_bp(x_b)]
                     =  [H_bp(x_a) - H_a(x_a) + H_a(x_a) - H_ap(x_a)] + [H_ap(x_b) - H_b(x_b) + H_b(x_b) - H_bp(x_b)]
                     =  [H_bp(x_a) - H_a(x_a)] - [H_ap(x_a) - H_a(x_a)] + [H_ap(x_b) - H_b(x_b)] - H_bp(x_b) - H_b(x_b)]
                     =  (de[bp][a] - de[ap][a]) + (de[ap][b] - de[bp][b])    */
            /* but, in the current code implementation, we flip configurations, not indices . . .
               So let's examine that.
                     =  [H_b(x_ap) - H_a(x_a)] - [H_a(x_ap) - H_a(x_a)] + [H_a(x_bp) - H_b(x_b)] - H_b(x_bp) - H_b(x_b)]
                     =  [H_b(x_ap) - H_a(x_ap)]  + [H_a(x_bp) - H_b(x_pb)]
                     = (de[b][ap] - de[a][ap]) + (de[a][bp] - de[b][bp]
                     So, if we exchange b<=> bp and a<=> ap, we return to the same result.
                     So the simple solution is to flip the
                     position of perturbed and original indices in the tests.
             */

            ediff = (de[bp][a] - de[ap][a]) + (de[ap][b] - de[bp][b]);
            delta = ediff*beta[a]; /* assume all same temperature in this case */
            break;
        case ereTL:
            /* not permuted:  */
            /* delta =  reduced E_new - reduced E_old
                     =  [beta_b H_b(x_a) + beta_a H_a(x_b)] - [beta_b H_b(x_b) + beta_a H_a(x_a)]
                     =  [beta_b H_b(x_a) - beta_a H_a(x_a)] + [beta_a H_a(x_b) - beta_b H_b(x_b)]
                     =  [beta_b dH_b(x_a) + beta_b H_a(x_a) - beta_a H_a(x_a)] +
                        [beta_a dH_a(x_b) + beta_a H_b(x_b) - beta_b H_b(x_b)]
                     =  [beta_b dH_b(x_a) + [beta_a dH_a(x_b) +
                        beta_b (H_a(x_a) - H_b(x_b)]) - beta_a (H_a(x_a) - H_b(x_b))
                     =  beta_b dH_b(x_a) + beta_a dH_a(x_b) - (beta_b - beta_a)(H_b(x_b) - H_a(x_a) */
            /* delta = beta[b]*de[b][a] + beta[a]*de[a][b] - (beta[b] - beta[a])*(Epot[b] - Epot[a]; */
            /* permuted (big breath!) */
            /*   delta =  reduced E_new - reduced E_old
                     =  [beta_bp H_bp(x_a) + beta_ap H_ap(x_b)] - [beta_bp H_bp(x_b) + beta_ap H_ap(x_a)]
                     =  [beta_bp H_bp(x_a) - beta_ap H_ap(x_a)] + [beta_ap H_ap(x_b) - beta_bp H_bp(x_b)]
                     =  [beta_bp H_bp(x_a) - beta_ap H_ap(x_a)] + [beta_ap H_ap(x_b) - beta_bp H_bp(x_b)]
                        - beta_pb H_a(x_a) + beta_ap H_a(x_a) + beta_pb H_a(x_a) - beta_ap H_a(x_a)
                        - beta_ap H_b(x_b) + beta_bp H_b(x_b) + beta_ap H_b(x_b) - beta_bp H_b(x_b)
                     =  [(beta_bp H_bp(x_a) - beta_bp H_a(x_a)) - (beta_ap H_ap(x_a) - beta_ap H_a(x_a))] +
                        [(beta_ap H_ap(x_b)  - beta_ap H_b(x_b)) - (beta_bp H_bp(x_b) - beta_bp H_b(x_b))]
             + beta_pb H_a(x_a) - beta_ap H_a(x_a) + beta_ap H_b(x_b) - beta_bp H_b(x_b)
                     =  [beta_bp (H_bp(x_a) - H_a(x_a)) - beta_ap (H_ap(x_a) - H_a(x_a))] +
                        [beta_ap (H_ap(x_b) - H_b(x_b)) - beta_bp (H_bp(x_b) - H_b(x_b))]
             + beta_pb (H_a(x_a) - H_b(x_b))  - beta_ap (H_a(x_a) - H_b(x_b))
                     =  ([beta_bp de[bp][a] - beta_ap de[ap][a]) + beta_ap de[ap][b]  - beta_bp de[bp][b])
             + (beta_pb-beta_ap)(H_a(x_a) - H_b(x_b))  */
            delta = beta[bp]*(de[bp][a] - de[bp][b]) + beta[ap]*(de[ap][b] - de[ap][a]) - (beta[bp]-beta[ap])*(Epot[b]-Epot[a]);
            break;
        default:
            printf("Unknown replica exchange quantity");
            exit(EXIT_FAILURE);
    }
    
    if (bPrint)
    {

        fprintf(fplog, "Repl %d <-> %d  dE_term = %10.3e (kT)\n", a, b, delta);
    
    }

    if (re->bNPT)
    {
        /* revist the calculation for 5.0.  Might be some improvements. */
        dpV = (beta[ap]*re->pres[ap]-beta[bp]*re->pres[bp])*(Vol[b]-Vol[a])/PRESFAC;
        if (bPrint)
        {
            fprintf(fplog, "  dpV = %10.3e  d = %10.3e\n", dpV, delta + dpV);
        }
        delta += dpV;
    }
    
    return delta;
}


void
test_for_replica_exchange(FILE                          *fplog,
                          struct gmx_repl_ex            *re,
                          float                         enerd,
                          double                         vol,
                          int                           step,
                          ReplicaExchangeParameters*    replExParams,
                          Replica_State*                      state_global)
{
    printf("GJS inside test_for_replica_exchange\n");
    int                                  m, i, j, a, b, ap, bp, i0, i1, tmp;
    float                                 delta = 0;
    int                             bPrint, bMultiEx;
    int                            *bEx      = re->bEx;
    float                                *prob     = re->prob;
    int                                 *pind     = re->destinations; /* permuted index */
    int                             bEpot    = false;
    int                             bDLambda = false;
    int                             bVol     = false;
    //gmx::ThreeFry2x64<64>                rng(re->seed, gmx::RandomDomain::ReplicaExchange);
    //gmx::UniformRealDistribution<float>   uniformRealDist;
    //gmx::UniformIntDistribution<int>     uniformNreplDist(0, re->nrepl-1);

    bMultiEx = (re->nex > 1);  /* multiple exchanges at each state */
    fprintf(fplog, "Replica exchange at step %d\n", step);

    if (re->bNPT)
    {
        for (i = 0; i < re->nrepl; i++)
        {
            re->Vol[i] = 0;
        }
        bVol               = true;
        re->Vol[re->repl]  = vol;
    }
    if ((re->type == ereTEMP || re->type == ereTL))
    {
        for (i = 0; i < re->nrepl; i++)
        {
            re->Epot[i] = 0;
        }
        bEpot              = true;
       
        /* temperatures of different states*/
        for (i = 0; i < re->nrepl; i++)
        {
            re->beta[i] = 1.0/(re->q[ereTEMP][i]*BOLTZ);
        }
    }
    else{
    
        for (i = 0; i < re->nrepl; i++)
        {
            re->beta[i] = 1.0/(re->temp*BOLTZ);  /* we have a single temperature */
        }
    }
    
    /* now actually do the communication */
    if (bVol)
    {
        #pragma omp barrier
        #pragma omp critical 
        {        
            replExParams->replica_volumes[re->repl] = vol;
        }
        #pragma omp barrier
        
        #pragma omp critical
        {
            for (int i = 0; i < re->nrepl; i++) 
                re->Vol[i] = replExParams->replica_volumes[i];        
        }
        #pragma omp barrier
    }
    
    if (bEpot)
    {
        #pragma omp barrier
        #pragma omp critical 
        {        
            replExParams->replica_energies[re->repl] = enerd;

            replExParams->replica_states[re->repl] = state_global;
        }
        #pragma omp barrier
        
        #pragma omp critical
        {
            for (int i = 0; i < re->nrepl; i++) 
                re->Epot[i] = replExParams->replica_energies[i];        
        }
        #pragma omp barrier
    }
    /* make a duplicate set of indices for shuffling */
    for (i = 0; i < re->nrepl; i++)
    {
        pind[i] = re->ind[i];
    }
    
    /* make a duplicate set of indices for shuffling */
    for (i = 0; i < re->nrepl; i++)
    {
        re->prob[i] = 0.0;
        re->prob_sum[i] = 0.0;
        re->bEx[i] = 0;
    }
    
    /* standard nearest neighbor replica exchange */

    m = (step / re->nst) % 2;
    
    for (i = 1; i < re->nrepl; i++)
    {
        a = re->ind[i-1];
        b = re->ind[i];

        bPrint = (re->repl == a || re->repl == b);
        
        if (i % 2 == m)
        {
        
        delta = calc_delta(fplog, bPrint, re, a, b, a, b);
            
            if (delta <= 0)
            {
                /* accepted */
                prob[i] = 1;
                bEx[i]  = 1;
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
                //uniformRealDist.reset();
                float randVar = ((float) rand() / (RAND_MAX)) + 1;
                bEx[i] = (int)(randVar < prob[i]);
            }
            
            re->prob_sum[i] += prob[i];

            if (bEx[i])
            {
                /* swap these two */
                tmp       = pind[i-1];
                pind[i-1] = pind[i];
                pind[i]   = tmp;
                re->nexchange[i]++;  /* statistics for back compatibility */
                printf("YOGA repl : %d, i : %d, numexch : %d\n", re->repl, i, re->nexchange[i]);  /* statistics for back compatibility */
            }
        }
        else
        {
            prob[i] = -1;
            bEx[i]  = 0;
        }
    }
    
    /* print some statistics */
    print_ind(fplog, "ex", re->nrepl, re->ind, bEx);
    print_prob(fplog, "pr", re->nrepl, prob);
    fprintf(fplog, "\n");
    re->nattempt[m]++;

    /* record which moves were made and accepted */
    for (i = 0; i < re->nrepl; i++)
    {
        re->nmoves[re->ind[i]][pind[i]] += 1;
        re->nmoves[pind[i]][re->ind[i]] += 1;
    }
    fflush(fplog); /* make sure we can see what the last exchange was */
}

static void
cyclic_decomposition(const int *destinations,
                     int      **cyclic,
                     int  *incycle,
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
/*
    if (debug)
    {
            for (i = 0; i < nrepl; i++)
            {
                fprintf(debug, "Cycle %d:", i);
                for (j = 0; j < nrepl; j++)
                {
                    if (cyclic[i][j] < 0)
                    {
                        break;
                    }
                    fprintf(debug, "%2d", cyclic[i][j]);
                }
                fprintf(debug, "\n");
            }
            fflush(debug);
    }
*/
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
        }
*/
}

static void
prepare_to_do_exchange(struct gmx_repl_ex *re,
                       const int           replica_id,
                       int                *maxswap,
                       int           *bThisReplicaExchanged)
{
    int i, j;
    /* Hold the cyclic decomposition of the (multiple) replica
     * exchange. */
    int bAnyReplicaExchanged = false;
    *bThisReplicaExchanged = false;

    for (i = 0; i < re->nrepl; i++)
    {
        if (re->destinations[i] != re->ind[i])
        {
            /* only mark as exchanged if the index has been shuffled */
            bAnyReplicaExchanged = true;
            break;
        }
    }
    if (bAnyReplicaExchanged)
    {
        /* reinitialize the placeholder arrays */
        for (i = 0; i < re->nrepl; i++)
        {
            for (j = 0; j < re->nrepl; j++)
            {
                re->cyclic[i][j] = -1;
                re->order[i][j]  = -1;
            }
        }

        /* Identify the cyclic decomposition of the permutation (very
         * fast if neighbor replica exchange). */
        cyclic_decomposition(re->destinations, re->cyclic, re->incycle, re->nrepl, maxswap);

        /* Now translate the decomposition into a replica exchange
         * order at each step. */
        compute_exchange_order(re->cyclic, re->order, re->nrepl, *maxswap);

        /* Did this replica do any exchange at any point? */
        for (j = 0; j < *maxswap; j++)
        {
            if (replica_id != re->order[replica_id][j])
            {
                *bThisReplicaExchanged = true;
                break;
            }
        }
    }
}

int replica_exchange(FILE *fplog, struct gmx_repl_ex *re,
                          Replica_State *state_global, float enerd, double vol_par,
                          int step, ReplicaExchangeParameters* replExParams)
{
    printf("GJS inside replica_exchange\n");
    int j;
    int replica_id = 0;
    int exchange_partner;
    int maxswap = 0;
    /* Number of rounds of exchanges needed to deal with any multiple
     * exchanges. */
    /* Where each replica ends up after the exchange attempt(s). */
    /* The order in which multiple exchanges will occur. */
    int bThisReplicaExchanged = false;

//    if (MASTER(cr))
        replica_id  = re->repl;

    // GJS figure out where vol is stored
    double vol = vol_par;

    test_for_replica_exchange(fplog, re, enerd, vol, step, replExParams, state_global);
    
    prepare_to_do_exchange(re, replica_id, &maxswap, &bThisReplicaExchanged);

    if (bThisReplicaExchanged)
    {
          //copy_state_serial(state_local, state);

//        if (MASTER(cr))
            /* There will be only one swap cycle with standard replica
             * exchange, but there may be multiple swap cycles if we
             * allow multiple swaps. */
            for (j = 0; j < maxswap; j++)
            {
                exchange_partner = re->order[replica_id][j];

                if (exchange_partner != replica_id && (replica_id % 2 == 0))
                {
                    /* Exchange the global states between the master nodes */
                    printf("GJS Exchanging %d with %d\n", replica_id, exchange_partner);
                    //exchange_state(state, exchange_partner, replExParams);
                    state_global = replExParams->replica_states[exchange_partner];
                    std::swap(replExParams->replica_states[replica_id], replExParams->replica_states[exchange_partner]);

                }
            }
            /* For temperature-type replica exchange, we need to scale
             * the velocities. */
            if (re->type == ereTEMP || re->type == ereTL)
            {
                //scale_velocities(state, sqrt(re->q[ereTEMP][replica_id]/re->q[ereTEMP][re->destinations[replica_id]]));
            }

            /* Copy the global state to the local state data structure */
            //copy_state_serial(state, state_local);
    }

    #pragma omp barrier

    return bThisReplicaExchanged;
}

void print_replica_exchange_statistics(FILE *fplog, struct gmx_repl_ex *re)
{
    int  i;

    fprintf(fplog, "\nReplica exchange statistics\n");

    if (re->nex == 0)
    {
        fprintf(fplog, "Repl  %d attempts, %d odd, %d even\n",
                re->nattempt[0]+re->nattempt[1], re->nattempt[1], re->nattempt[0]);

        fprintf(fplog, "Repl  average probabilities:\n");
        for (i = 1; i < re->nrepl; i++)
        {
            if (re->nattempt[i%2] == 0)
            {
                re->prob[i] = 0;
            }
            else
            {
                re->prob[i] =  re->prob_sum[i]/re->nattempt[i%2];
            }
        }
        print_ind(fplog, "", re->nrepl, re->ind, nullptr);
        print_prob(fplog, "", re->nrepl, re->prob);

        fprintf(fplog, "Repl  number of exchanges:\n");
        print_ind(fplog, "", re->nrepl, re->ind, nullptr);
        print_count(fplog, "", re->nrepl, re->nexchange);

        fprintf(fplog, "Repl  average number of exchanges:\n");
        for (i = 1; i < re->nrepl; i++)
        {
            if (re->nattempt[i%2] == 0)
            {
                re->prob[i] = 0;
            }
            else
            {
                re->prob[i] =  ((float)re->nexchange[i])/re->nattempt[i%2];
            }
        }
        print_ind(fplog, "", re->nrepl, re->ind, nullptr);
        print_prob(fplog, "", re->nrepl, re->prob);

        fprintf(fplog, "\n");
    }
    /* print the transition matrix */
    print_transition_matrix(fplog, re->nrepl, re->nmoves, re->nattempt);
}

