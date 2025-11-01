#include "FM_Flips_Internal.h"
#include "limits.h"

#define GAIN(dir, ci, gamma, Awi) ( (dir) * ( (ci) - (gamma) * (Awi) ) )

void FM_Flips_Struct_Check
(
    FM_Flips_Struct *data,
    bool         check_Az
)
{
    static long long unsigned count = 0 ;
    count++ ;

    WARN_IF_NOT_DEBUG_MODE ;

    // Unpack struct
    VSPINT *z      = data->z ;
    VSPINT *Aw     = data->Aw ;
    VSPINT *Az     = data->Az ;
    VSPINT   n     = data->n ;
    VSPINT *locked = data->locked ;
    VSPINT *c      = data->c ;

    VSPINT i, j, entry ;

    // Check that the weight is correct
    VSPINT wgt = 0 ;
    for(i = 0 ; i < n ; i++)
    {
        wgt += z[i] * c[i] ;
    }

    if(wgt != data->weight)
    {
        VSP_error("Wrong weight", DBG_INFO) ;
    }

    // Check the heap property of A_up
    Max_Queue_v1 *A_up = data->A_up ;
    for(i = 1 ; i <= data->A_up->n_queue ; i++)
    {
        VSPINT val_par = A_up->queue[i] ;
        VSPINT key_par = data->c[val_par] ;

        // Check left child
        if(2 * i <= A_up->n_queue)
        {
            VSPINT val_left = A_up->queue[2 * i] ;
            VSPINT key_left = data->c[val_left] ;
            if(key_par < key_left)
            {
                VSP_error("Max heap property violated", DBG_INFO) ;
                break ;
            }
        }

        if(2 * i + 1 <= A_up->n_queue)
        {
            VSPINT val_right = A_up->queue[2 * i + 1] ;
            VSPINT key_right = data->c[val_right] ;
            if(key_par < key_right)
            {
                VSP_error("Max heap property violated", DBG_INFO) ;
                break ;
            }
        }
    }


    // Check the heap property of A_dn
    Min_Queue_v1 *A_dn = data->A_dn ;
    for(i = 1 ; i <= data->A_dn->n_queue ; i++)
    {
        VSPINT val_par = A_dn->queue[i] ;
        VSPINT key_par = data->c[val_par] ;

        // Check left child
        if(2 * i <= A_dn->n_queue)
        {
            VSPINT val_left = A_dn->queue[2 * i] ;
            VSPINT key_left = data->c[val_left] ;
            if(key_par > key_left)
            {
                VSP_error("Min heap property violated", DBG_INFO) ;
                break ;
            }
        }

        if(2 * i + 1 <= A_dn->n_queue)
        {
            VSPINT val_right = A_dn->queue[2 * i + 1] ;
            VSPINT key_right = data->c[val_right] ;
            if(key_par > key_right)
            {
                VSP_error("Min heap property violated", DBG_INFO) ;
                break ;
            }
        }
    }

    // Check that Az is correct
    if(check_Az)
    {
        VSPINT *Ai = data->Ai ;
        VSPINT *Ap = data->Ap ;
        VSPINT *Az_true = malloc(n * sizeof(VSPINT)) ;
        for(i = 0; i < n; i++) Az_true[i] = 0 ;

        for (VSPINT col = 0; col < n; col++)
        {
            if (z[col])
            {
                VSPINT start = Ap[col] ;
                VSPINT end   = Ap[col + 1] ;
                for (i = start ; i < end ; i++) Az_true[ Ai[i] ] ++ ;
                /* Take identity into account */
                Az_true[col] ++ ;
            }
        }
        for(i = 0 ; i < n ; i++)
        {
            VSPINT Azi = Az[i] ;
            VSPINT Azi_true = Az_true[i] ;
            if(Azi != Azi_true)
            {
                VSP_error("Az is incorrect", DBG_INFO) ;
            }
        }
        free(Az_true) ;
    }


    // Check that zi = 0 --> i is not locked
    for(i = 0 ; i < n ; i++)
    {
        if(z[i] == 0 && locked[i])
        {
            VSP_error("i is locked while zi is 0", DBG_INFO) ;
            break ;
        }
    }

    // Check p_lo and p_hi
    VSPINT p_lo = n + 1 ;
    VSPINT p_hi = -1 ;
    for(i = 0 ; i < n ; i++)
    {
        if( z[i] == 1 && !locked[i] )
        {
                if(Aw[i] > p_hi) p_hi = Aw[i] ;
        }
        else if (z[i] == 0)
        {
            if(Aw[i] < p_lo) p_lo = Aw[i] ;
        }
    }
    if(p_lo != data->p_lo)
    {

        VSP_error("Wrong p_lo", DBG_INFO) ;
    }
    if(p_hi != data->p_hi)
    {
        VSP_error("Wrong p_hi", DBG_INFO) ;
    }

    // Check c_lo and c_hi
    VSPINT c_lo = data->gamma + 1 ;
    VSPINT c_hi =  0 ;
    UNUSED VSPINT c_lo_witness = EMPTY ;
    UNUSED VSPINT c_hi_witness = EMPTY ;

    for(i = 0 ; i < n ; i++)
    {
        if(z[i] == 1 && Aw[i] == p_hi && !locked[i] )
        {
            if(c[i] < c_lo)
            {
                c_lo = c[i] ;
                c_lo_witness = i ;
            }
        }
        else if(z[i] == 0 && Aw[i] == p_lo)
        {
             if(c[i] > c_hi)
             {
                c_hi = c[i] ;
                c_hi_witness = i ;
             }
        }
    }
    if(c_lo != data->c_lo)
    {
        VSP_error("Wrong c_lo", DBG_INFO) ;
    }
    if(c_hi != data->c_hi)
    {
        VSP_error("Wrong c_hi", DBG_INFO) ;
    }

    // Check that each member of O_up has p_lo and c_hi
    Direct_Access_List *O_up = data->O_up ;
    for(i = 0 ; i < O_up->size ; i++)
    {
        j = O_up->list[i] ;
        if(Aw[j] != p_lo)
        {
            VSP_error("Member of O_up does not have p_lo", DBG_INFO) ;
        }

        if(c[j] != c_hi)
        {
            VSP_error("Member of O_up does not have c_hi", DBG_INFO) ;
        }
    }

    // Check that each member of A_up has p_lo and a c less than c_hi
    for(i = 1 ; i <= A_up->n_queue ; i++)
    {
        j = A_up->queue[i] ;
        if(Aw[j] != p_lo)
        {
            VSP_error("Member of A_up does not have p_lo", DBG_INFO) ;
        }

        if(c[j] >= c_hi)
        {
            VSP_error("Member of A_up has too high c", DBG_INFO) ;
        }
    }

    // Check that each member of O_dn has p_hi and c_lo
    Direct_Access_List *O_dn = data->O_dn ;
    for(i = 0 ; i < O_dn->size ; i++)
    {
        j = O_dn->list[i] ;
        if(Aw[j] != p_hi)
        {
            VSP_error("Member of O_dn does not have p_hi", DBG_INFO) ;
        }

        if(c[j] != c_lo)
        {
            VSP_error("Member of O_dn does not have c_lo", DBG_INFO) ;
        }
    }

    // Check that each member of A_dn has p_hi and a c greater than c_lo
    for(i = 1 ; i <= A_dn->n_queue ; i++)
    {
        j = A_dn->queue[i] ;
        if(Aw[j] != p_hi)
        {
            VSP_error("Member of A_dn does not have p_hi", DBG_INFO) ;
        }

        if(c[j] <= c_lo)
        {
            VSP_error("Member of A_dn has too low c", DBG_INFO) ;
        }
    }

    // Go through all entries and make sure they're in the right structure
    for(i = 0 ; i < n ; i++)
    {
        if(z[i] == 1 && Aw[i] == p_hi && !locked[i])
        {
            if(c[i] > c_lo)
            {
                if(A_dn->lookup[i] == EMPTY)
                {
                    VSP_error("Bug", DBG_INFO) ;
                }
            }
            else if(c[i] == c_lo)
            {
                if(O_dn->lookup[i] == EMPTY)
                {
                    VSP_error("Bug", DBG_INFO) ;
                }
            }
            else
            {
                VSP_error("Bug", DBG_INFO) ;
            }
        }
        else if(z[i] == 0 && Aw[i] == p_lo)
        {
            if(c[i] < c_hi)
            {
                if(A_up->lookup[i] == EMPTY)
                {
                    VSP_error("i should be in A_up but it's not", DBG_INFO) ;
                }
            }
            else if(c[i] == c_hi)
            {
                if(O_up->lookup[i] == EMPTY)
                {
                    VSP_error("Bug", DBG_INFO) ;
                }
            }
            else
            {
                VSP_error("Bug", DBG_INFO) ;
            }
        }
    }

    // Go through all the linked lists and make sure there are no duplicates,
    // and that 0 ... n-1 are partitioned into locked nodes and nodes in the
    // linked lists
    VSPINT *L0spread = Doubly_LLs_Unsized_Spread(data->L0) ;
    VSPINT *L1spread = Doubly_LLs_Unsized_Spread(data->L1) ;

    for(i = 0 ; i < n ; i++)
    {
        VSPINT in_L0 = (L0spread[i] == 2) ;
        VSPINT in_L1 = (L1spread[i] == 2) ;
        VSPINT in_lk = (locked[i]   == 1) ;

        if(in_L0 + in_L1 + in_lk != 1)
        {
            VSP_error("Error in FM flips struct", DBG_INFO) ;
            break ;
        }
    }
    free(L0spread) ;
    free(L1spread) ;

    // Go through all the lists and verify that the zi and Awi values are 
    // correct
    VSPINT *fwd = data->L0->fwd ;
    for(VSPINT whichll = 0 ; whichll <= n ; whichll++)
    {
        entry = data->L0->entries[whichll] ;
        if(entry == EMPTY) continue ;

        i = entry ;
        do
        {
            if(z[i] != 0 || Aw[i] != whichll)
            {
                VSP_error("Problem in linked lists", DBG_INFO) ;
            }
        } while( (i = fwd[i]) != entry ) ;
    }

    fwd = data->L1->fwd ;
    for(VSPINT whichll = 0 ; whichll <= n ; whichll++)
    {
        entry = data->L1->entries[whichll] ;
        if(entry == EMPTY) continue ;

        i = entry ;
        do
        {
            if(z[i] != 1 || Aw[i] != whichll)
            {
                VSP_error("Problem in linked lists", DBG_INFO) ;
            }
        } while( (i = fwd[i]) != entry ) ;
    }
    return ;
}



// NOTE: This function looks for best flip, but recall that if we run up against
// constraint we don't actually do the best flip.
// Check if the given flip is indeed feasible and optimal. Take locked array
// into account. Lower balance is IGNORED
// Inputs left as VSPINT*'s so that I have flexibility in which structs I use
// Varying the structs I use can help catch bugs
// ASSUME: the flip has not been locked yet
void FM_Flips_Check_Best_Flip
(
    VSPINT              *z,
    VSPINT              *w,
    VSPINT             *Az,
    VSPINT             *Aw,
    VSPINT       *locked_z,
    VSPINT       *locked_w,
    VSP_problem   *problem,
    VSPINT        ind_flip, // EMPTY if none
    char          var_flip
)
{
    static long long unsigned count = 0 ;
    count++ ;
    WARN_IF_NOT_DEBUG_MODE ;

    // Local vars
    VSPINT i ;
    VSPINT n = problem->Graph->n ;
    VSPINT u = problem->u ;
    VSPINT *c = problem->Graph->c ;
    VSPINT gamma = problem->Graph->max_vtx_cost ;

    // Compute weights
    VSPINT wgt_z = 0, wgt_w = 0 ;
    for(i = 0 ; i < n ; i++)
    {
        wgt_z += z[i] * c[i] ;
        wgt_w += w[i] * c[i] ;
    }

    // Check given flip satisfies upper balance constraint. Get obj change
    VSPINT dobj_given, dir_flip ;
    switch(var_flip)
    {
        case 'z':
            dir_flip    = 1 - 2 * z[ind_flip] ;
            dobj_given  = c[ind_flip] - gamma * Aw[ind_flip] ;
            dobj_given *= dir_flip ;

            if(locked_z[ind_flip])
            {
                VSP_error("Trying to flip locked entry", DBG_INFO) ;
            }

            if(dir_flip > 0 && wgt_z + c[ind_flip] > u)
            {
                VSP_error("Flip violates upper balance", DBG_INFO) ;
            }
            break ;
        case 'w':
            dir_flip    = 1 - 2 * w[ind_flip] ;
            dobj_given  = c[ind_flip] - gamma * Az[ind_flip] ;
            dobj_given *= dir_flip ;

            if(locked_w[ind_flip])
            {
                VSP_error("Trying to flip locked entry", DBG_INFO) ;
            }

            if(dir_flip > 0 && wgt_w + c[ind_flip] > u)
            {
                VSP_error("Flip violates upper balance", DBG_INFO) ;
            }
            break ;
        default:
            VSP_error("Invalid variable given", DBG_INFO) ;
            return ;
    }

    // Search for the best flip
    VSPINT ind_best_flip = EMPTY ;
    VSPINT gain_best = - gamma * n - 1 ;
    VSPINT gain_current ;
    UNUSED char var_best_flip ;
    VSPINT dir ;
    for(i = 0 ; i < n ; i++)
    {
        // Consider the z flip
        if(!(locked_z[i]) )
        {
            dir = 1 - 2 * z[i] ;
            gain_current = dir * (c[i] - gamma * Aw[i]) ;
            if(dir < 0 || wgt_z + c[i] <= u)
            {
                if(gain_current > gain_best)
                {
                    gain_best = gain_current ;
                    ind_best_flip = i ;
                    var_best_flip = 'z' ;
                }
            }
        }

        // Consider the w flip
        if( !(locked_w[i]) )
        {
            dir = 1 - 2 * w[i] ;
            gain_current = dir * (c[i] - gamma * Az[i]) ;
            if(dir < 0 || wgt_w + c[i] <= u)
            {
                if(gain_current > gain_best)
                {
                    gain_best = gain_current ;
                    ind_best_flip = i ;
                    var_best_flip = 'w' ;
                }
            }
        }
    }

    if(ind_flip == EMPTY)
    {
        if(ind_best_flip != EMPTY)
        {
            VSP_error("Best flip not found", DBG_INFO) ;
        }
    }
    else
    {
        if(dobj_given != gain_best)
        {
            VSP_error("Given flip is not optimal", DBG_INFO) ;
        }
    }
    return ;
}