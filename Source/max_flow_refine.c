#include "core.h"

#define ENQUEUE(i, q, qend) ( (q)[++(qend)] = (i) )
#define DEQUEUE(q, qstart) ((q)[(qstart)++])
#define QUEUE_NONEMPTY(qstart, qend) ( (qstart) <= (qend) )
#define OTHER_SIDE(i) (1 - (i))

// Gives 0, 1, or 2 if i is in side 0 (X side), side 1 (Y side), or separator
#define SIDE(i, sides) (2 - 2*(sides[0])[i] -  (sides[1])[i])


#define SS_INSERT(i, Adense, Asparse, numA)          \
    do { (Asparse)[i] = (numA) ;                     \
         (Adense)[(numA)++] = (i) ; } while (0)



                                            

static bool In_Sep(VSPINT i, VSPINT *sides[2])
{
    return sides[0][i] == 0 && sides[1][i] == 0 ;
}

static VSPLONG Degree(VSPINT i, VSPLONG *Ap)
{
    return Ap[i+1] - Ap[i] ;
}

// Return the number of vertices in A. A return value of 0 means that
// either the separator is empty, or its degree exceeds deg_budget
static VSPINT Get_A
(
    VSPINT  *Ai,
    VSPLONG *Ap,
    VSPINT   n,
    VSPINT   u,
    VSPINT  *c,
    VSPINT  *sides[2],
    VSPINT   wgts[3],     // assumed to be accurate
    VSPLONG  deg_budget,
    VSPINT  *Adense,      // can be garbage
    VSPINT  *Asparse      // assumed to be all -1's
)
{
    // Initialize variables
    VSPINT qstart = 0 ;
    VSPINT numA = 0 ;
    VSPINT side_budgets[2] = {u-wgts[2]-wgts[1],
                              u-wgts[2]-wgts[0]  } ;
    int side_enabled[3] = {TRUE, TRUE, TRUE} ;
    int deg_budget_remains = TRUE ;

    // Uninitialized
    int bi, bj ;
    VSPINT i, j, cost, deg ;
    VSPLONG p, end ;

    // Add S to A. Return early if separator exhausts all degree budget
    for(i = 0 ; i < n ; i++)
    {
        if(In_Sep(i, sides))
        {
            deg_budget -= Degree(i, Ap) ;
            if(deg_budget < 0)
            {
                // Clear Asparse and return
                for(j = 0 ; j < numA ; j++) Asparse[ Adense[j] ] = EMPTY ;
                return 0 ;
            }

            // Insert i into A and the queue
            Asparse[i] = numA ;
            Adense[numA++] = i ;
        }
    }

    // Do the BFS
    while( qstart < numA
           && (side_enabled[0] || side_enabled[1])
           && deg_budget_remains )
    {
        i = Adense[qstart++] ;
        bi = SIDE(i, sides) ;
        end = Ap[i+1] ;
        
        // If i is in the 0 or 1 side and that side is disabled, none of its
        // neighbors (which are in its own side or in the separator) can be
        // added to A. So skip i. i can however be in the separator, so
        // The side_enabled[2] is always TRUE, and this will make sure
        // we don't skip it, since sep can have neighbors on both sides
        if(side_enabled[bj] == FALSE) continue ;

        // Process neighbors of i
        for(p = Ap[i] ; p < end ; p++)
        {
            j = Ai[p] ;
            bj = SIDE(j, sides) ; // i may be in the sep and j on a side

            // Barrier 1: skip nodes of A (all separator nodes are in A already)
            if(Asparse[j] >= 0) continue ;



            // Permanently disable nodes entering from j's side if
            // the cost of j surpasses the remaining budget for that side
            cost = c[j] ;
            if(side_budgets[bj] - cost < 0) side_enabled[bj] = FALSE ;

            // Barrier 2: skip j, and all remaining neighbors of i, if j's
            // side became disabled
            if(side_enabled[bj] == FALSE)
            {
                if(bj < 2) break ;
                else      continue ;
            }

            // Barrier 3: stop building A if about to surpass the degree budget
            deg = Degree(j, Ap) ;
            if(deg_budget - deg < 0)
            {
                deg_budget_remains = FALSE ;
                break ;
            }

            // Insertion into A: finally, j passed all the barriers so it's
            // time to insert it in the queue, and decrease the budgets
            Asparse[j] = numA ;
            Adense[numA++] = j ;
            side_budgets[bj] -= cost ;
            deg_budget -= deg;
        }
    }
    return numA ;
}

static void Generate_Max_Flow_Graph
(
    VSPINT  *Ai,
    VSPLONG *Ap,
    VSPINT   n,
    VSPINT  *Fi,
    VSPINT  *Fp,
    VSPINT  *numF,   // return: number of vertices in F
    VSPINT  *Adense,
    VSPINT  *Asparse,
    VSPINT   numA
)
{

}