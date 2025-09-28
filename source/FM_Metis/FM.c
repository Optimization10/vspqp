#include "FM_internal.h"
#include "ds_api.h"
#include "shared_api.h"

/*----------------------------------------------------------------------------*/
/*      Static Functions Prototypes                                           */
/*----------------------------------------------------------------------------*/

static void CheckAllFMDataNoLock
(
    VSPINT              *where,
    VSPINT            pwgts[3],
    QueueFM         *queues[2],
    ExtNeighWgts      *extWgts,
    VSPINT          badmaxpwgt,

    VSPINT             *bndind,
    VSPINT             *bndptr,
    VSPINT                nbnd,

    VSP_problem       *problem
) ;

static void Where_To_States
(
    VSPINT        *where,
    VSP_state         *Z,
    VSP_state         *W,
    VSPINT      pwgts[3],
    VSP_problem *problem
) ;

static void Check_Ext_Neigh_Wgts
(
    ExtNeighWgts     *extWgts,
    VSPINT                *Ai,
    VSPINT                *Ap,
    VSPINT                 *c,
    VSPINT                  n,
    VSPINT             *where 
) ;

static void BND_Insert
(
    VSPINT       *nbnd,
    VSPINT     *bndind,
    VSPINT     *bndptr,
    VSPINT    toInsert
) ;

static void BND_Delete
(
    VSPINT       *nbnd,
    VSPINT     *bndind,
    VSPINT     *bndptr,
    VSPINT   nodeToDel
) ;

static void Check_FM_Inputs
(
    VSP_state         *Z,
    VSP_state         *W,
    VSP_problem *problem
) ;

static void Fill_Int_Array_With_Value
(
    VSPINT     val,
    VSPINT  *array,
    VSPINT       n
) ;

static void Init_Where_And_Sep
(
    VSPINT           *sep,
    VSPINT   *sep_locator,
    VSPINT          *nsep,
    VSPINT         *where,
    VSPINT             *z,
    VSPINT             *w,
    VSPINT              n
) ;

static void Set_Ext_Neigh_Wgts
(
    ExtNeighWgts     *extWgts,
    VSPINT                *Ai,
    VSPINT                *Ap,
    VSPINT                 *c,
    VSPINT               *sep,
    VSPINT               nsep,
    VSPINT             *where
) ;

// Do FM moves where the gain is >= 0, and perhaps some neutral moves as well
void MovesImprovingCost
(
    VSPINT              *where,
    VSPINT            pwgts[3],
    QueueFM         *queues[2],
    ExtNeighWgts      *extWgts,
    VSPINT          badmaxpwgt,
    VSPINT             *bndind,
    VSPINT             *bndptr,
    VSPINT               *nbnd,
    VSPINT                pass,
    VSPINT allowedNeutralMoves,
    VSP_problem       *problem
)
{
    SET_FUNCTION_COUNTER ;

    VSPINT topVals[2], topGains[2], to, other, higain, j, k, oldgain,
    numNeutralMoves = 0,
    n = problem->Graph->n,
    *c = problem->Graph->c,
    *Ap = problem->Graph->Ap,
    *Ai = problem->Graph->Ai,
    nflips = 0 ; // nflips accurate at beginning of loop


    while( 1 )
    {
        #ifdef VSP_DEBUG
            CheckAllFMDataNoLock
            (
                where, pwgts, queues, extWgts, badmaxpwgt, bndind,
                bndptr, *nbnd, problem
            ) ;
        #endif

        /*--------------------------------------------------------------------*/
        // Get the best "to"
        /*--------------------------------------------------------------------*/
        topVals[0] = QueueFM_See_Top_Value(queues[0]);
        topVals[1] = QueueFM_See_Top_Value(queues[1]);

        if (topVals[0] != EMPTY && topVals[1] != EMPTY)
        {
            topGains[0] = c[topVals[0]] - extWgts[topVals[0]].wgts[1];
            topGains[1] = c[topVals[1]] - extWgts[topVals[1]].wgts[0];

            // to = (topGains[0] > topGains[1] ? 0 : (topGains[0] < topGains[1] ? 1 : pass % 2));

            if(topGains[0] > topGains[1]) to = 0 ;
            else if(topGains[0] < topGains[1]) to = 1 ;
            else
            {
                if(pwgts[0] > pwgts[1]) to = 1 ;
                else if(pwgts[0] < pwgts[1]) to = 0 ;
                else to = rand() % 2 ;
            } 
            other = 1 - to ;

            VSP_ASSERT( IN_RANGE(topVals[to], 0, n ) ) ;
            if (pwgts[to] + c[topVals[to]] > badmaxpwgt)
            {
                to = (to + 1) % 2 ;
                other = 1 - to ;
                if (pwgts[to] + c[topVals[to]] > badmaxpwgt)
                {
                    // printf("\nBlocked infeasible FM move. Problem: %s  Level: %i", problem->name, problem->level) ;
                    break ;
                }
            }

            // to is the final selection. Only let it go through if the move is
            // beneficial, or if it's neutral and we are still under the limit
            // for these moves
            // Guard agains bad moves
            if(topGains[to] == 0)
            {
                if(numNeutralMoves < allowedNeutralMoves)
                {
                    numNeutralMoves++ ;
                }
                else break ;
            }
            else if(topGains[to] < 0) break ;
        }
        else if (topVals[0] == EMPTY && topVals[1] == EMPTY)
        {
            break;
        }
        else if (topVals[0] != EMPTY && pwgts[0] + c[topVals[0]] <= badmaxpwgt)
        {
            to = 0 ;
            other = 1 - to ;

            // Guard agains bad moves
            if(topGains[to] == 0)
            {
                if(numNeutralMoves < allowedNeutralMoves)
                {
                    numNeutralMoves++ ;
                }
                else break ;
            }
            else if(topGains[to] < 0) break ;
        }
        else if (topVals[1] != EMPTY && pwgts[1] + c[topVals[1]] <= badmaxpwgt)
        {
            to = 1 ;
            other = 1 - to ;

            // Guard agains bad moves
            if(topGains[to] == 0)
            {
                if(numNeutralMoves < allowedNeutralMoves)
                {
                    numNeutralMoves++ ;
                }
                else break ;
            }
            else if(topGains[to] < 0) break ;
        }
        else break;

        other = (to + 1) % 2 ;

        /*--------------------------------------------------------------------*/
        // Apply updates due to primary move
        /*--------------------------------------------------------------------*/
        higain = QueueFM_Get_Top_Value(queues[to]) ;
        QueueFM_Delete(queues[other], higain) ;

        VSP_ASSERT( IN_RANGE(higain, 0, n) ) ;
        VSP_ASSERT(bndptr[higain] != -1);

        pwgts[2] -= (c[higain] - extWgts[higain].wgts[other]) ;

        /* Move node higain from the separator into to */
        BND_Delete(nbnd, bndind, bndptr, higain) ;
        pwgts[to]     += c[higain] ;
        where[higain]  = to ;

        #ifdef VSP_DEBUG
            BND_Check(*nbnd, bndind, bndptr, where, n) ;
        #endif

        /*--------------------------------------------------------------------*/
        // Apply updates due to secondary moves
        /*--------------------------------------------------------------------*/
        for (j = Ap[higain]; j < Ap[higain + 1]; j++)
        {
            #ifdef VSP_DEBUG
                BND_Check(*nbnd, bndind, bndptr, where, n) ;
            #endif
            k = Ai[j];

            if (where[k] == 2) /* If k is in the separator */
            { 
                /* For the in-separator vertices modify their edegree[to] */
                oldgain = c[k] - extWgts[k].wgts[to];
                extWgts[k].wgts[to] += c[higain];

                // Update gain of moving k into other
                QueueFM_Update(queues[other], k, oldgain - c[higain]) ;
            }
            else if (where[k] == other) /* If k is in other */
            { 
                /* Pull k into the separator */
                where[k] = 2 ;
                BND_Insert(nbnd, bndind, bndptr, k) ;
                pwgts[other] -= c[k] ;

                // Why not pwgts[2] += c[k]? Because pwgts[2] was fully updated
                // with pwgts[2] -= (c[higain] - extWgts[higain].wgts[other]) ;

                /* Now that k is in the separator, compute its external wgts */
                /* TODO: I think we only need to move this if k is not locked */
                (extWgts[k].wgts)[0] = 0 ;
                (extWgts[k].wgts)[1] = 0 ;
                for (VSPINT jj = Ap[k]; jj < Ap[k + 1]; jj++)
                {
                    VSPINT kk = Ai[jj];
                    if (where[kk] != 2)
                    {
                        (extWgts[k].wgts)[where[kk]] += c[kk];
                    }
                    else // kk is in the separator
                    {
                        /* kk is in the separator and is a neighbor of k. Since
                        k just moved from other into the sep, there is a
                        change in the external weight of kk */
                        oldgain = c[kk] - (extWgts[kk].wgts)[other] ;
                        (extWgts[kk].wgts)[other] -= c[k] ;

                        /* If kk were to move to to, it no longer has to
                        drag its neighbor k from other, because its
                        neighbor k is in the separator now. So the gain
                        for pushing kk into to has improved by c[k] */
                        QueueFM_Update(queues[to], kk, oldgain + c[k]);
                        
                    }
                }

                // Insert k into queues
                QueueFM_Insert(queues[to], k, c[k] - (extWgts[k].wgts)[other]) ;
                QueueFM_Insert(queues[other], k, c[k] - (extWgts[k].wgts)[to]) ;

            }
        }
        nflips++ ;
    }


    #ifdef VSP_DEBUG
        CheckAllFMDataNoLock
        (
            where, pwgts, queues, extWgts, badmaxpwgt, bndind,
            bndptr, *nbnd, problem
        ) ;
    #endif
    return ;
}



/*----------------------------------------------------------------------------*/
/*      Driver Function                                                       */
/*----------------------------------------------------------------------------*/
/* Assume we start with a bona fide separator */
void FM_Metis
(
    VSP_state             *Z,
    VSP_state             *W,
    VSP_problem     *problem,
    UNUSED VSP_options     *options,
    VSPINT          n_passes
)
{

    TIMER_START ;

    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count++ ;

        Check_FM_Inputs(Z, W, problem) ;
        
        VSPFLOAT objin = VSP_compute_obj(Z, W, problem, problem->Graph->max_vtx_cost) ;
    #endif
    /*--- Parameters (temporary, to be put in options struct -----------------*/
    VSPINT *moved = malloc( problem->Graph->n * sizeof(*moved) ) ;
    /*------------------------------------------------------------------------*/

    /* Declare local variables */
    VSPINT i, j, to, other, topGains[2], topVals[2], nswaps ;
    VSPINT k, oldgain ;

    /* Unpack structs */
    VSPINT *z  = Z->point->int_part ;
    VSPINT *w  = W->point->int_part ;
    VSPINT  n  = problem->Graph->n ;
    VSPINT *c  = problem->Graph->c ;
    VSPINT *Ai = problem->Graph->Ai ;
    VSPINT *Ap = problem->Graph->Ap ;

    /* Initialize swaps array */
    VSPINT *swaps = malloc( n * sizeof(*swaps) ) ;

    /* Initialize partition weights */
    VSPINT pwgts[3] ;
    pwgts[0] = Z->weight ;
    pwgts[1] = W->weight ;
    pwgts[2] = problem->Graph->total_vtx_cost - pwgts[0] - pwgts[1] ;

    /* Initalize the mind structure */
    VSPINT *mind    = malloc(2 * n * sizeof(*mind)) ;
    VSPINT *mptr    = malloc( (n + 1) * sizeof(*mptr) ) ;
    VSPINT  nmind   = 0 ;
    mptr[0]         = 0 ;
    

    /* Initialize separator structure and where array */
    VSPINT nbnd         = 0;
    VSPINT *bndind      = malloc(n * sizeof(*bndind)) ;
    VSPINT *bndptr      = malloc(n * sizeof(*bndptr)) ;
    VSPINT *where       = malloc(n * sizeof(*where)) ;
    Init_Where_And_Sep(bndind, bndptr, &nbnd, where, z, w, n) ;


    /* Initialize the external neighbor weights */
    ExtNeighWgts *extWgts = malloc(n * sizeof(*extWgts) ) ;
    Set_Ext_Neigh_Wgts(extWgts, Ai, Ap, c, bndind, nbnd, where) ;

    /* Initialize the queues */
    QueueFM *queues[2] ;
    queues[0] = QueueFM_Construct(n) ;
    queues[1] = QueueFM_Construct(n) ;

    #ifdef VSP_DEBUG
        Check_Ext_Neigh_Wgts(extWgts, Ai, Ap, c, n, where) ;
        BND_Check(nbnd, bndind, bndptr, where, n) ;
        QueueFM_Check_Heap_Light(queues[0]) ;
        QueueFM_Check_Heap_Light(queues[1]) ;
    #endif

    /*------------------------------------------------------------------------*/
    /* FM Loop                                                                */
    /*------------------------------------------------------------------------*/


    /* higain is the vertex with the highest gain */
    VSPINT higain ;

    /* mindiff is the diff in c'x and x'y at the best x, y so far */
    VSPINT mindiff    = abs( pwgts[0] - pwgts[1] ) ; 
    to                = (pwgts[0] < pwgts[1] ? 0 : 1) ;

    /* This is the current difference between c'x and c'y at the current
       nswap iteration */
    VSPINT newdiff ;

    /* badmaxpwgt is the upper balance constraint */
    VSPINT badmaxpwgt = problem->u ;

    /* limit is used to quit FM early if it's taking too long */
    VSPINT limit      = options->FM_Limit ; // In Metis, limit = (ctrl->compress ? gk_min(5*nbnd, 400) : gk_min(2*nbnd, 300));





    /* TODO: pass is the number of FM passes. It is used in the code to throw
       in some randomness */
    VSPINT pass = 0 ;

    // n_passes = 1 ;
    // n_passes = problem->level + 4 ;
    for(pass = 0 ; pass < n_passes ; pass++)
    {

        /* Initialize moved array */
        Fill_Int_Array_With_Value(EMPTY, moved, n) ;

        /* Initialize queues */
        Fill_Int_Array_With_Value(EMPTY, queues[0]->locator, n) ;
        Fill_Int_Array_With_Value(EMPTY, queues[1]->locator, n) ;
        queues[0]->numNodes = 0 ;
        queues[1]->numNodes = 0 ;
        for(i = 0 ; i < nbnd ; i++)
        {
            j = bndind[i] ;
            QueueFM_Insert(queues[0], j, c[j] - extWgts[j].wgts[1] ) ;
            QueueFM_Insert(queues[1], j, c[j] - extWgts[j].wgts[0] ) ;
        }


        MovesImprovingCost(where, pwgts, queues, extWgts, badmaxpwgt, bndind, bndptr, &nbnd, pass, 0, problem) ;

        mindiff = abs( pwgts[0] - pwgts[1] ) ;

        /* mincut is the weight of the best separator so far */
        VSPINT mincut = pwgts[2] ;


        /* mincutorder is the nswap at which the best separator is achieved */
        VSPINT mincutorder = EMPTY ;

        // TODO: Careful with abs here if VSPINT is not int
        mptr[0]           = 0 ;
        nmind             = 0 ;

        // Start an FM pass. Moves here may not improve cost
        for (nswaps = 0; nswaps < n; nswaps++)
        {
            topVals[0] = QueueFM_See_Top_Value(queues[0]) ;
            topVals[1] = QueueFM_See_Top_Value(queues[1]) ;

            /* This if structure is all about figuring out the best "to" */
            if (topVals[0] != EMPTY && topVals[1] != EMPTY)
            {
                topGains[0] = c[topVals[0]] - extWgts[topVals[0]].wgts[1];
                topGains[1] = c[topVals[1]] - extWgts[topVals[1]].wgts[0];

                to = (topGains[0] > topGains[1] ? 0 : (topGains[0] < topGains[1] ? 1 : pass % 2));
                other = 1 - to ;

                VSP_ASSERT( IN_RANGE(topVals[to], 0, n ) ) ;

                // Check if the chosen move "to" is feasible. If not, chose the other move
                if (pwgts[to] + c[topVals[to]] > badmaxpwgt)
                {
                    to = (to + 1) % 2 ;
                    other = 1 - to ;
                    if (pwgts[to] + c[topVals[to]] > badmaxpwgt)
                    {
                        break ;
                    }
                }

                VSP_ASSERT( IN_RANGE(topVals[to], 0, n ) ) ;
                if (pwgts[to] + c[topVals[to]] > badmaxpwgt)
                {
                    to = (to + 1) % 2 ;
                    other = 1 - to ;
                    if (pwgts[to] + c[topVals[to]] > badmaxpwgt)
                    {
                        break ;
                    }
                } 
            }
            else if (topVals[0] == EMPTY && topVals[1] == EMPTY)
            {
                break;
            }
            else if (topVals[0] != EMPTY && pwgts[0] + c[topVals[0]] <= badmaxpwgt)
            {
                to = 0 ;
                other = 1 - to ;
            }
            else if (topVals[1] != EMPTY && pwgts[1] + c[topVals[1]] <= badmaxpwgt)
            {
                to = 1 ;
                other = 1 - to ;
            }
            else break;

            VSP_ASSERT(to + other == 1) ;

            /* Ovi: Other is the opposite of to. If to is 0, other is 1, and v.v. */
            other = (to + 1) % 2 ;

            /* Retrieve top value and delete it from the queue */
            // higain is node we want to move
            higain = QueueFM_Get_Top_Value(queues[to]) ;
            VSP_ASSERT( IN_RANGE(higain, 0, n) ) ;
            VSP_ASSERT(bndptr[higain] != -1);

            /* Delete if it was in the separator originally */ // Ovi: why delte this conditionally, and not always?
            if (moved[higain] == EMPTY) 
            {
                QueueFM_Delete(queues[other], higain) ;
            }

            /* The following check is to ensure we break out if there is a possibility
            of over-running the mind array. Ovi: Recall mind has size 2n */
            if (nmind + Ap[higain + 1] - Ap[higain] >= 2 * n - 1)
            {
                break;
            }

            pwgts[2] -= (c[higain] - extWgts[higain].wgts[other]);

            /* Calculate new weight difference between z and w */
            newdiff = abs(pwgts[to] + c[higain] - (pwgts[other] - extWgts[higain].wgts[other])) ;

            /* Determine if the current paritition is best so far */
            if (pwgts[2] < mincut || (pwgts[2] == mincut && newdiff < mindiff))
            {
                mincut      = pwgts[2];
                mincutorder = nswaps;
                mindiff     = newdiff;
            }
            else
            {
                /* Determine if we should quit early */
                bool ran_too_long = (nswaps - mincutorder > 2 * limit) ;
                bool too_long_since_last_improvement = (nswaps - mincutorder > limit && pwgts[2] > 1.10 * mincut) ;

                if (ran_too_long || too_long_since_last_improvement)
                {
                    pwgts[2] += (c[higain] - extWgts[higain].wgts[other]);
                    break; /* No further improvement, break out */
                }
            }

            /* Move node higain from the separator into to */
            BND_Delete(&nbnd, bndind, bndptr, higain);
            pwgts[to]     += c[higain];
            where[higain]  = to;
            moved[higain]  = nswaps;
            swaps[nswaps]  = higain;

            #ifdef VSP_DEBUG
                BND_Check(nbnd, bndind, bndptr, where, n) ;
            #endif

            /**********************************************************
             * Update the degrees of the affected nodes
             ***********************************************************/
            for (j = Ap[higain]; j < Ap[higain + 1]; j++)
            {
                #ifdef VSP_DEBUG
                    BND_Check(nbnd, bndind, bndptr, where, n) ;
                #endif
                k = Ai[j];

                if (where[k] == 2) /* If k is in the separator */
                { 
                    /* For the in-separator vertices modify their edegree[to] */
                    oldgain = c[k] - extWgts[k].wgts[to];
                    extWgts[k].wgts[to] += c[higain];

                    /* TODO: What is this? I think I know the answer. k SHOULD only
                    be in the queue if it started in the sep and never moved,
                    or if it moved at some point in the past from other into sep.

                    if moved[k] == EMPTY, we should update the other queue (there is no update to make to the to queue)
                    if moved[k] signifies that k was dragged in from to, it can no longer go back to to but it can go to other, and the queue[other] should be updated
                    if moved[k] signifies k was dragged in from other, it cannot go back to other, it can only go to to. But the to queue was no affected, so there is
                                no need to make any updated
                    if moved[k] signifies that k moved as a vertex from the sep of its own volition,
                                then it's cannot move independently anymore so shouldn't be in any queue


                    If neither condition is met, then k should not be in the to
                    queue to begin with. */
                    if (moved[k] == EMPTY || moved[k] == -(2 + other))
                    {
                        QueueFM_Update(queues[other], k, oldgain - c[higain]);
                    }
                }
                else if (where[k] == other) /* If k is in other */
                { 
                    /* Pull k into the separator */
                    where[k] = 2 ;
                    BND_Insert(&nbnd, bndind, bndptr, k) ;
                    pwgts[other] -= c[k] ;

                    /* Update mind array with this move */
                    mind[nmind++] = k ;

                    /* Now that k is in the separator, compute its external wgts */
                    /* TODO: I think we only need to move this if k is not locked */
                    (extWgts[k].wgts)[0] = 0 ;
                    (extWgts[k].wgts)[1] = 0 ;
                    for (VSPINT jj = Ap[k]; jj < Ap[k + 1]; jj++)
                    {
                        VSPINT kk = Ai[jj];
                        if (where[kk] != 2)
                        {
                            (extWgts[k].wgts)[where[kk]] += c[kk];
                        }
                        else // kk is in the separator
                        {
                            /* kk is in the separator and is a neighbor of k. Since
                            k just moved from other into the sep, there is a
                            change in the external weight of kk */
                            oldgain = c[kk] - (extWgts[kk].wgts)[other] ;
                            (extWgts[kk].wgts)[other] -= c[k] ;

                            /* kk is in the queue only if it hasn't moved, or if
                            if it was dragged into the separator */
                            if (moved[kk] == -1 || moved[kk] == -(2 + to))
                            {
                                /* If kk were to move to to, it no longer has to
                                drag its neighbor k from other, because its
                                neighbor k is in the separator now. So the gain
                                for pushing kk into to has improved by c[k] */
                                QueueFM_Update(queues[to], kk, oldgain + c[k]);
                            }
                        }
                    }

                    /* Insert the new vertex into the priority queue. Only one side! */
                    if (moved[k] == -1)
                    {
                        QueueFM_Insert(queues[to], k, c[k] - (extWgts[k].wgts)[other]);
                        moved[k] = -(2 + to);
                    }
                }
            }
            mptr[nswaps + 1] = nmind;

            #ifdef VSP_DEBUG
                Check_Ext_Neigh_Wgts(extWgts, Ai, Ap, c, n, where) ;
                BND_Check(nbnd, bndind, bndptr, where, n) ;
                QueueFM_Check_Heap_Light(queues[0]) ;
                QueueFM_Check_Heap_Light(queues[1]) ;
                if(pwgts[0] > problem->u)
                {
                    VSP_error("Part 0 violated upper bound", DBG_INFO) ;
                }
                if(pwgts[1] > problem->u)
                {
                    VSP_error("Part 1 violated upper bound", DBG_INFO) ;
                }

            #endif
        }

        // printf("\n\n Fraction of swaps: %f\n\n", (nswaps / (double) n) ) ;

        /****************************************************************
         * Roll back computation
         *****************************************************************/
        for (nswaps--; nswaps > mincutorder; nswaps--)
        {
            higain = swaps[nswaps] ;

            to = where[higain];
            other = (to + 1) % 2;

            /* Revert the weights of the partitions */
            // TODO: If we just remember the weights, it should be unneccesary to
            // revert them step by step like this
            INC_DEC(pwgts[2], pwgts[to], c[higain]);

            /* Move higain back into the separator */
            // TODO: on the last FM pass, we don't have to keep the bnd structure
            // updated anymore
            where[higain] = 2;
            BND_Insert(&nbnd, bndind, bndptr, higain);

            /* Compute edegrees of higain */
            extWgts[higain].wgts[0] = 0 ;
            extWgts[higain].wgts[1] = 0 ;
            for (j = Ap[higain] ; j < Ap[higain + 1] ; j++)
            {
                k = Ai[j];
                if (where[k] == 2) // Update ext wgt of neighbors of higain in sep
                {
                    extWgts[k].wgts[to] -= c[higain] ;
                }
                else // Update ext wgt of higain with the cost of k
                {
                    extWgts[higain].wgts[where[k]] += c[k] ;
                }
            }

            /* Push nodes out of the separator */
            for (j = mptr[nswaps]; j < mptr[nswaps + 1]; j++)
            {
                k = mind[j];
                VSP_ASSERT(where[k] == 2) ;
                #ifdef VSP_DEBUG
                    if(where[k] != 2)
                    {
                        VSP_error("Node not in sep", DBG_INFO) ;
                    }
                #endif

                /* Move k back to other */
                where[k] = other;
                INC_DEC(pwgts[other], pwgts[2], c[k]);
                BND_Delete(&nbnd, bndind, bndptr, k);
                for(VSPINT jj = Ap[k]; jj < Ap[k + 1]; jj++)
                {
                    VSPINT kk = Ai[jj];
                    if (where[kk] == 2)
                    {
                        extWgts[kk].wgts[other] += c[k] ;
                    }
                }
            }
        }

        #ifdef VSP_DEBUG
            if(mincut != pwgts[2])
            {
                VSP_error("Mincut doesn't match size of separator", DBG_INFO) ;
            }
        #endif


        // Flips between FM passes
        // Where_To_States(where, Z, W, pwgts, problem) ;
        // Flips(Z, W, problem, 5) ;

        // // Prepare data structures for next pass of FM
        // pwgts[0] = Z->weight ;
        // pwgts[1] = W->weight ;
        // pwgts[2] = problem->Graph->total_vtx_cost - pwgts[0] - pwgts[1] ;
        
        // nbnd = 0;
        // Init_Where_And_Sep(bndind, bndptr, &nbnd, where, z, w, n) ;
        
        // Set_Ext_Neigh_Wgts(extWgts, Ai, Ap, c, bndind, nbnd, where) ;
    }

    

    /* Convert FM info back to states */
    Where_To_States(where, Z, W, pwgts, problem) ;


    /* Free memory */
    free(swaps) ;
    free(bndind) ;
    free(bndptr) ;
    free(where) ;
    free(mind) ;
    free(mptr) ;
    free(moved) ;
    free(extWgts) ;
    QueueFM_Destruct(queues[0]) ;
    QueueFM_Destruct(queues[1]) ;

    #ifdef VSP_DEBUG
        VSP_check_states(Z, W, problem) ;
        VSPFLOAT objout = VSP_compute_obj(Z, W, problem, problem->Graph->max_vtx_cost) ;
        if(objout < objin)
        {
            VSP_error("Objective decreased during FM", DBG_INFO) ;
        }
    #endif

    TIMER_END ;
    PERFORMANCE_UPDATE_TIME(FM) ;

    return ;
}


/*----------------------------------------------------------------------------*/
/*      Helper Functions                                                      */
/*----------------------------------------------------------------------------*/

void BND_Check
(
    VSPINT        nbnd,
    VSPINT     *bndind,
    VSPINT     *bndptr,
    VSPINT      *where,
    VSPINT           n  // size of where array
)
{
    static long long unsigned count = 0 ;
    count++ ;

    VSPINT i, j ;
    
    /* Check that everything stored in bndind is in the separator */
    for(i = 0 ; i < nbnd ; i++)
    {
        if(bndind[i] < 0 || bndind[i] > n - 1)
        {
            VSP_error("Boundary vertex out of bounds", DBG_INFO) ;
        }

        if( where[bndind[i]] != 2 )
        {
            VSP_error("Boundary vertex not in separator", DBG_INFO) ;
            break ;
        }
    }

    /* Check that everything in the separator is in the bndind array.
       Also, check bndptr */
    for(i = 0 ; i < n ; i++)
    {
        if(where[i] < 0 || where[i] > 2)
        {
            VSP_error("Partition number out of bounds", DBG_INFO) ;
        }

        if(where[i] == 2)
        {
            /* Check bndptr */
            j = bndptr[i] ;
            if( !IN_RANGE(j, 0, n) )
            {
                VSP_error("Out of bounds", DBG_INFO) ;
                break ;
            }

            if(bndind[j] != i || j >= nbnd)
            {
                VSP_error("Bndptr incorrect", DBG_INFO) ;
            }


            /* Do linear search through bndind for i */
            bool found = false ;
            for(j = 0 ; j < nbnd ; j++)
            {
                if( bndind[j] == i)
                {
                    found = true ;
                    break ;
                }
            }
            if(!found)
            {
                VSP_error("Vertex in separator but not in bnd structure", DBG_INFO) ;
            }
        }
        else // if where[i] != 2
        {
            if(bndptr[i] != EMPTY)
            {
                VSP_error("Bndptr claims vertex is in sep when it's not", DBG_INFO) ;
                break ;
            }

            /* Linear search through bndind to make sure i is not there */
            bool found = false ;
            for(j = 0 ; j < nbnd ; j++)
            {
                if( bndind[j] == i)
                {
                    found = true ;
                    break ;
                }
            }
            if(found)
            {
                VSP_error("Vertex in separator but not in bnd structure", DBG_INFO) ;
            }
        }
    }
    return ;
}

/* ASSUMED: the node is NOT in the list */
static void BND_Insert
(
    VSPINT       *nbnd,
    VSPINT     *bndind,
    VSPINT     *bndptr,
    VSPINT    toInsert
)
{
    bndind[*nbnd] = toInsert ;
    bndptr[toInsert] = *nbnd ;
    (*nbnd)++ ;

    return ;
}

/* ASSUMED: The node IS in the list */
static void BND_Delete
(
    VSPINT       *nbnd,
    VSPINT     *bndind,
    VSPINT     *bndptr,
    VSPINT   nodeToDel
)
{
    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count++ ;
    #endif

    #ifdef VSP_DEBUG
        // Validate that there are nodes to delete
        if (*nbnd <= 0)
        {
            char msg[ERR_MSG_LEN];
            snprintf(msg, sizeof(msg), "Error: No nodes to delete. *nbnd = %d", *nbnd);
            VSP_error(msg, DBG_INFO);
        }

        // Validate nodeToDel
        if (nodeToDel < 0)
        {
            char msg[ERR_MSG_LEN];
            snprintf(msg, sizeof(msg), "Deleting invalid node") ;
            VSP_error(msg, DBG_INFO);
        }

        // Check if the node is present in list
        if (bndptr[nodeToDel] == EMPTY)
        {
            char msg[ERR_MSG_LEN];
            snprintf(msg, sizeof(msg), "Node not present in list");
            VSP_error(msg, DBG_INFO);
        }
    #endif

    /* Note: programming trick here to avoid branching. Normally, you would
       think you have to check if the node to be deleted is the last node,
       and in that case, you do not need to overwrite anything. In fact, if
       you do things in the wrong order, you can write the last node back after
       you deleted it. Doing it the way I have it below does do a bit of extra
       work when the node to be deleted is the last node (rare) but it avoids
       branching in all cases. */

    /* Overwrite the node to be deleted with the last node */
    VSPINT locDelete = bndptr[nodeToDel] ;
    VSPINT lastNode = bndind[--*nbnd] ;
    bndind[locDelete] = lastNode ;
    bndptr[lastNode] = locDelete ;
 
    /* Delete the node from pointers */
    bndptr[nodeToDel] = EMPTY ;


    #ifdef VSP_DEBUG
        // Validate locDelete
        if (locDelete < 0)
        {
            char msg[ERR_MSG_LEN];
            snprintf(msg, sizeof(msg), "Invalid location") ;
            VSP_error(msg, DBG_INFO);
        }

        if (lastNode < 0)
        {
            char msg[ERR_MSG_LEN];
            snprintf(msg, sizeof(msg), "Invalid node") ;
            VSP_error(msg, DBG_INFO);
        }
    #endif
    return ;
}

/* Check that the points are feasible, there are no fractional entries, and 
   we have a bona fide separator (i.e. the cross term is 0)                   */
static void Check_FM_Inputs
(
    VSP_state         *Z,
    VSP_state         *W,
    VSP_problem *problem
)
{
    static long long unsigned count = 0 ;
    count++ ;


    VSPINT *z     = Z->point->int_part ;
    VSPINT *w     = W->point->int_part ;

    VSPINT *Az    = Z->Az->int_part ;
    VSPINT *Aw    = W->Az->int_part ;

    VSPINT  n     = problem->Graph->n ;

    /* Check feasibility, among other things */
    VSP_check_states(Z, W, problem) ;

    /* Check for no fractional component */
    if(Z->point->ind_frac != EMPTY || W->point->ind_frac != EMPTY)
    {
        VSP_error("Fractional component in FM", DBG_INFO) ;
    }

    /* Check that we have a bona fide separator by checking the cross term's 0*/
    VSPINT zAw = 0, wAz = 0 ;
    for(VSPINT i = 0 ; i < n ; i++)
    {
        zAw += z[i] * Aw[i] ;
        wAz += w[i] * Az[i] ;
    }
    if(zAw != 0 || wAz != 0)
    {
        VSP_error("Cross term isn't zero in FM", DBG_INFO) ;
    }
    return ;
}

static void Fill_Int_Array_With_Value
(
    VSPINT     val,
    VSPINT  *array,
    VSPINT       n
)
{
    for(VSPINT i = 0 ; i < n ; i++) array[i] = val ;
}

// TODO: Consider if initializing the external neighbor weights here as well
//       would be beneficial */
static void Init_Where_And_Sep
(
    VSPINT           *sep,
    VSPINT   *sep_locator,
    VSPINT          *nsep,
    VSPINT         *where,
    VSPINT             *z,
    VSPINT             *w,
    VSPINT              n
)
{
    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count++ ;
    #endif

    VSPINT i, nsep_local = 0 ;
    for(i = 0 ; i < n ; i++)
    {
        if( z[i] )
        {
            where[i]  = 0 ;
            sep_locator[i] = EMPTY ;
        }
        else if( w[i] )
        {
            where[i]  = 1 ;
            sep_locator[i] = EMPTY ;
        }
        else
        {
            where[i] = 2 ;

            sep[nsep_local] = i ;
            sep_locator[i] = nsep_local ;
            nsep_local++ ;
        }
    }
    *nsep = nsep_local ;

    #ifdef VSP_DEBUG
        BND_Check(*nsep, sep, sep_locator, where, n) ;
    #endif

    return ;
}

static void Set_Ext_Neigh_Wgts
(
    ExtNeighWgts     *extWgts,
    VSPINT                *Ai,
    VSPINT                *Ap,
    VSPINT                 *c,
    VSPINT               *sep,
    VSPINT               nsep,
    VSPINT             *where
)
{
    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count++ ;
    #endif

    VSPINT i, j, k, p, end, which_set ;
    for(k = 0 ; k < nsep ; k++)
    {
        i = sep[k] ;
        extWgts[i].wgts[0] = 0 ;
        extWgts[i].wgts[1] = 0 ;
        end = Ap[i + 1] ;
        for(p = Ap[i] ; p < end ; p++)
        {
            j = Ai[p] ;
            which_set = where[j] ;
            VSP_ASSERT(0 <= which_set && which_set <= 2) ;
            if(which_set != 2)
            {
                extWgts[i].wgts[which_set] += c[j] ;
            }
        }
    }
}


static void Check_Ext_Neigh_Wgts
(
    ExtNeighWgts     *extWgts,
    VSPINT                *Ai,
    VSPINT                *Ap,
    VSPINT                 *c,
    VSPINT                  n,
    VSPINT             *where 
)
{
    static long long unsigned count = 0 ;
    count++ ;

    for(VSPINT i = 0 ; i < n ; i++)
    {
        if(where[i] == 2)
        {
            VSPINT ext0 = 0 ;
            VSPINT ext1 = 0 ;
            
            /* Go through neighbors and compute external weights */
            for(VSPINT p = Ap[i] ; p < Ap[i+1] ; p++)
            {
                VSPINT k = Ai[p] ;
                if(where[k] == 0)      ext0 += c[k] ;
                else if(where[k] == 1) ext1 += c[k] ;
            }

            /* Check the result */
            if(ext0 != extWgts[i].wgts[0])
            {
                VSP_error("External weight into partition 0 incorrect", DBG_INFO) ;
                break ;
            }

            if(ext1 != extWgts[i].wgts[1])
            {
                VSP_error("External weight into partition 1 incorrect", DBG_INFO) ;
                break ;
            }
        }
    }
    return ;
}

/* Convert where array into VSP_states */
static void Where_To_States
(
    VSPINT        *where,
    VSP_state         *Z,
    VSP_state         *W,
    VSPINT      pwgts[3],
    VSP_problem *problem
)
{
    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count++ ;

        VSP_check_states(Z, W, problem) ;
        VSPFLOAT objin = VSP_compute_obj(Z, W, problem, problem->Graph->max_vtx_cost) ;
    #endif

    VSPINT n = problem->Graph->n ;

    VSPINT *points[2] ;
    points[0] = Z->point->int_part ;
    points[1] = W->point->int_part ;

    VSP_state *states[2] ;
    states[0] = Z ;
    states[1] = W ;

    VSPINT *Azs[2] ;
    Azs[0] = Z->Az->int_part ;
    Azs[1] = W->Az->int_part ;

    /* Update weights */
    states[0]->weight = pwgts[0] ;
    states[1]->weight = pwgts[1] ;

    /* Update points and Azs */
    for(VSPINT i = 0 ; i < n ; i++)
    {
        VSP_ASSERT( IN_RANGE(where[i], 0, 3) ) ;
        
        /* If i is in the separator, both points should be 0 */
        if(where[i] == 2)
        {
            if(points[0][i] != 0)
            {
                points[0][i] = 0 ;
                VSP_Az_pm_col_of_A(Azs[0], i, MINUS, problem) ;
            }

            if(points[1][i] != 0)
            {
                points[1][i] = 0 ;
                VSP_Az_pm_col_of_A(Azs[1], i, MINUS, problem) ;
            }
        }
        else // if i is in Z or W
        {
            VSPINT to    = where[i] ;
            VSPINT other = (to + 1) % 2 ;

            /* to point should be 1 */
            if(points[to][i] == 0)
            {
                points[to][i] = 1 ;
                VSP_Az_pm_col_of_A(Azs[to], i, PLUS, problem) ;
            }

            /* other point should be 0 */
            if(points[other][i] != 0)
            {
                points[other][i] = 0 ;
                VSP_Az_pm_col_of_A(Azs[other], i, MINUS, problem) ;
            }  
        }      
    }

    #ifdef VSP_DEBUG
        VSP_check_states(Z, W, problem) ;

        VSPFLOAT objout = VSP_compute_obj(Z, W, problem, problem->Graph->max_vtx_cost) ;
        if(objout < objin)
        {
            VSP_error("Objective decreased after FM", DBG_INFO) ;
        }

        VSPINT objFM = pwgts[0] + pwgts[1] ;
        if(objout != objFM)
        {
            VSP_error("Objective doesn't match", DBG_INFO) ;
        }

        /* Cross term should be 0 */
        VSPINT zAw = VSP_penalty_term(Z->point->int_part, W->Az->int_part, n) ;
        VSPINT wAz = VSP_penalty_term(W->point->int_part, Z->Az->int_part, n) ;

        if(zAw != 0 || wAz != 0)
        {
            VSP_error("Penalty term nonzero", DBG_INFO) ;
        }
    #endif

    return ;
}


// This is meant for the phase where we do no locking whatsoever
// Check the following:
// parition weights are correct and partition is feasible. All verties are in part 0, 1, or 2
// If i is in the separator, then its extWgts info is correct
// bndind equals the separator
// Each queue equals the separator, and heap proeprty is satisfied
static void CheckAllFMDataNoLock
(
    VSPINT              *where,
    VSPINT            pwgts[3],
    QueueFM         *queues[2],
    ExtNeighWgts      *extWgts,
    VSPINT          badmaxpwgt,

    VSPINT             *bndind,
    VSPINT             *bndptr,
    VSPINT                nbnd,

    VSP_problem       *problem
)
{
    static long long unsigned count = 0 ; 
    count++ ;
    WARN_IF_NOT_DEBUG_MODE ;

    VSPINT i, p, j, k, key,
    n = problem->Graph->n,
    *c = problem->Graph->c,
    *Ai = problem->Graph->Ai,
    *Ap = problem->Graph->Ap ;
    /*------------------------------------------------------------------------*/
    // - parition weights are correct
    // - partition is feasible
    // - part[i] is 0, 1, or 2 for all vertices i
    /*------------------------------------------------------------------------*/
    VSPINT localpwgts[3] = {0} ;
    for(i = 0 ; i < n ; i++)
    {
        VSP_ASSERT(IN_RANGE(where[i], 0, 3)) ;
        localpwgts[ where[i] ] += c[i] ;
    }

    for(p = 0 ; p < 3 ; p++)
    {
        if(pwgts[p] != localpwgts[p]) ERR("Wrong partition weights") ;
    }

    if( pwgts[0] > badmaxpwgt || pwgts[1] > badmaxpwgt)
    {
        ERR("Infeasible partition") ;
    }

    /*------------------------------------------------------------------------*/
    // If i is in the separator, then its extWgts info is correct
    /*------------------------------------------------------------------------*/
    for(i = 0 ; i < n ; i++)
    {
        if(where[i] != 2) continue ;

        ExtNeighWgts localExtWgt = { .wgts = {0} } ; 
        for(p = Ap[i] ; p < Ap[i+1] ; p++)
        {
            j = Ai[p] ;
            if(where[j] != 2) localExtWgt.wgts[ where[j] ] += c[j] ;
        }

        if(localExtWgt.wgts[0] != extWgts[i].wgts[0])
        {
            ERR("External weights incorrect") ;
        }

        if(localExtWgt.wgts[1] != extWgts[i].wgts[1])
        {
            ERR("External weights incorrect") ;
        }
    }

    /*------------------------------------------------------------------------*/
    // Check bnd structure. What I would like it to do is the following,
    // altough I can't guarantee that the checker does all of this:
    // - bndind equals the separator
    // - bndptr[i] != EMPTY iff i is in the sep, and if i is in the sep, then
    //   bndptr[i] = its location in bndind
    /*------------------------------------------------------------------------*/
    BND_Check(nbnd, bndind, bndptr, where, n ) ;


    /*------------------------------------------------------------------------*/
    // - Each queue equals the separator
    // - queue keys are correct
    // - Call the queue checker
    /*------------------------------------------------------------------------*/
    // Check that all nodes in each quere are in the separator, and their
    // keys are correct
    for(k = 0 ; k < 2 ; k++)
    {
        QueueFM *q = queues[k] ;
        Node Node ;
        for(p = 1 ; p <= q->numNodes ; p++)
        {
            i   = q->heap[p].value ;
            key = q->heap[p].key ;

            if(where[i] != 2)
            {
                ERR("i is in queue but not in separator") ;
            }

            VSPINT localKey = c[i] - extWgts[i].wgts[ 1 - k ] ;
            if (key != localKey)
            {
                ERR("Heap key incorrect") ;
            }
        }
    }

    // Check that all nodes in the separator are in the queues
    for(i = 0 ; i < n ; i++)
    {
        if(where[i] != 2) continue ;
        if(queues[0]->locator[i] == EMPTY)
        {
            ERR("Node in separator but not in queue") ;
        }

        if(queues[1]->locator[i] == EMPTY)
        {
            ERR("Node in separator but not in queue") ;
        }
    }

    QueueFM_Check_Heap_Light(queues[0] ) ;
    QueueFM_Check_Heap_Light(queues[1] ) ;
}