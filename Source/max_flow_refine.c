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
    VSPINT   wgts[3], // assumed to be accurate
    VSPLONG  deg_budget,
    VSPINT  *Adense,
    VSPINT  *Asparse // assumed to be all -1's
)
{
    // Initialize variables
    VSPINT qstart = 0 ;
    VSPINT numA = 0 ;
    VSPINT side_budgets[2] = {u-wgts[2]-wgts[1],
                              u-wgts[2]-wgts[0]  } ;
    int side_enabled[2] = {TRUE, TRUE} ;
    int deg_budget_remains = TRUE ;

    // Uninitialized
    int b ;
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
        i = DEQUEUE(q, qstart) ;
        b = SIDE(i, sides) ;
        end = Ap[i+1] ;
        
        // If i's side is disabled, none of its neighbors (which are in its own
        // side or in the separator) can be added to A. So skip i
        if(side_enabled[b] == FALSE) continue ;

        // Process neighbors of i
        for(p = Ap[i] ; p < end ; p++)
        {
            j = Ai[p] ;

            // Barrier 1: skip nodes of A (all separator nodes are in A already)
            if(Asparse[j] >= 0) continue ;

            // Barrier 2: stop building A if about to surpass the degree budget
            deg = Degree(j, Ap) ;
            if(deg_budget - deg < 0)
            {
                deg_budget_remains = FALSE ;
                break ;
            }

            // Permanently disable nodes entering from j's side if
            // the cost of j surpasses the remaining budget for that side
            cost = c[j] ;
            if(side_budgets[b] - cost < 0) side_enabled[b] = FALSE ;

            // Barrier 3: skip j, and all remaining neighbors of i, if j's
            // side became disabled
            if(side_enabled[b] == FALSE) break ;

            // Insertion into A: finally, j passed all the barriers so it's
            // time to insert it in the queue, and decrease the budgets
            Asparse[j] = numA ;
            Adense[numA++] = j ;
            q[++qend] = j ;
            side_budgets[b] -= cost ;
            deg_budget -= deg;
        }
    }
    return numA ;
}