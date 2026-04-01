#include "shared_api.h"
#include "ds_api.h" // IWYU pragma: keep

static void VSP_check_penalties_one_state
(
    VSP_state *Z,
    VSP_state *W,
    VSP_penalties *Z_pen,
    VSPINT n
) ;

static void VSP_check_a_linked_list(
    VSPINT     entry,
    VSPINT        *z,
    VSPINT *fwd_list,
    VSPINT *rev_list,
    VSPINT     z_val,
    VSPINT    Aw_val,
    VSPINT       *Aw,
    VSPINT  *witness
) ;

/* Check that the linked lists are all correct */
void VSP_check_linked_lists
(
    VSPINT* z,
    VSPINT *Aw,
    VSP_penalties *Z_pen
)
{

    static long long unsigned count = 0 ;
    count++ ;

    WARN_IF_NOT_DEBUG_MODE ;

    VSPINT entry ;
    VSPINT n = Z_pen->n ;

    VSPINT *z0_entry    = Z_pen->z0_entry ;
    VSPINT *z0_fwd_list = Z_pen->z0_fwd_list ;
    VSPINT *z0_rev_list = Z_pen->z0_rev_list ;

    VSPINT *z1_entry    = Z_pen->z1_entry ;
    VSPINT *z1_fwd_list = Z_pen->z1_fwd_list ;
    VSPINT *z1_rev_list = Z_pen->z1_rev_list ;


    /* Allocate memory */
    VSPINT *witness = calloc(n, sizeof(*witness) ) ;

    /* Check the lists */
    for(VSPINT Aw_val = 0 ; Aw_val < n + 1 ; Aw_val++) {

        /* Transverse the z0 linked list corresponding to Aw_val */
        entry = z0_entry[Aw_val] ;
        if(entry != EMPTY) {
            VSP_check_a_linked_list(entry , z , z0_fwd_list, z0_rev_list, 0,
                                    Aw_val, Aw, witness                      ) ;
        }

        /* Transverse the z1 linked list corresponding to Aw_val */
        entry = z1_entry[Aw_val] ;
        if(entry != EMPTY) {
            VSP_check_a_linked_list(entry , z , z1_fwd_list, z1_rev_list, 1,
                                    Aw_val, Aw, witness                      ) ;
        }
    }

    /* Check that each node is present in precisely one linked list */
    for(VSPINT i = 0 ; i < n ; i++) {
        if(witness[i] != 1) {
            VSP_error("Problem in linked lists", DBG_INFO) ;
            break ;
        }
    }

    /* Free memory */
    free(witness) ;

    return ;
}

static void VSP_check_a_linked_list(
    VSPINT     entry,
    VSPINT        *z,
    VSPINT *fwd_list,
    VSPINT *rev_list,
    VSPINT     z_val,
    VSPINT    Aw_val,
    VSPINT       *Aw,
    VSPINT  *witness
) {
    static long long unsigned count = 0 ;
    count++ ;

    WARN_IF_NOT_DEBUG_MODE ;

    VSPINT p = entry ;
    do{
        witness[p]++ ;
        if( rev_list[ fwd_list[p] ] != p)
        {
            VSP_error("Linked list not correct", DBG_INFO) ;
            return ;
        }
        if( z[p] != z_val )
        {
            VSP_error("Linked list not correct", DBG_INFO) ;
            return ;
        }
        if( Aw[p] != Aw_val ) {
            VSP_error("Linked list not correct", DBG_INFO) ;
            return ;
        }

        p = fwd_list[p] ;
    } while(p != entry) ;
}

void VSP_check_penalties
(
    VSP_state *X,
    VSP_state *Y,
    VSP_penalties *X_pen,
    VSP_penalties *Y_pen,
    VSPINT n
)
{
    WARN_IF_NOT_DEBUG_MODE ;
    VSP_check_penalties_one_state(X, Y, X_pen, n) ;
    VSP_check_penalties_one_state(Y, X, Y_pen, n) ;

    return ;
}


static void VSP_check_penalties_one_state
(
    VSP_state *Z,
    VSP_state *W,
    VSP_penalties *Z_pen,
    VSPINT n
)
{

    static long long unsigned count = 0 ;
    count++ ;

    WARN_IF_NOT_DEBUG_MODE ;

    /* 1. Check the linked lists */
    VSP_check_linked_lists
    (
        Z->point->int_part,
        W->Az->int_part,
        Z_pen
    ) ;


    /* 2. Check the min and max penalties */
    VSPINT i ;
    VSPINT local_min_pen = n + 1 ;
    VSPINT local_max_pen = -1 ;
    VSPINT *z = Z->point->int_part ;
    VSPINT *Aw = W->Az->int_part ;

    for(i = 0 ; i < n ; i++)
    {
        if(z[i])
        {
            if(Aw[i] > local_max_pen)
            {
                local_max_pen = Aw[i] ;
            }
        }
        else
        {
            if(Aw[i] < local_min_pen)
            {
                local_min_pen = Aw[i] ;
            }
        }
    }

    if(local_min_pen != Z_pen->min_pen)
    {
        VSP_error("Error in min penalty.", DBG_INFO) ;
    }

    if(local_max_pen != Z_pen->max_pen)
    {
        VSP_error("Error in max penalty.", DBG_INFO) ;
    }

    return ;
}
