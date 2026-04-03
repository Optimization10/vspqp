#include "core.h"
#include "shared_api.h"
#include "ds_api.h"


Penalties *construct_penalties(VSPINT n)
{
    Penalties *pens = malloc(sizeof(*pens)) ;

    pens->z0_entry = malloc( (n + 1) * sizeof(VSPINT) ) ;
    pens->z1_entry = malloc( (n + 1) * sizeof(VSPINT) ) ;

    pens->z0_fwd_list = malloc( n * sizeof(VSPINT) ) ;
    pens->z1_fwd_list = malloc( n * sizeof(VSPINT) ) ;

    pens->z0_rev_list = malloc( n * sizeof(VSPINT) ) ;
    pens->z1_rev_list = malloc( n * sizeof(VSPINT) ) ;

    pens->n = n ;

    return pens ;
}

void destruct_penalties(Penalties *pens) {
    free(pens->z0_entry) ;
    free(pens->z1_entry) ;

    free(pens->z0_fwd_list) ;
    free(pens->z1_fwd_list) ;

    free(pens->z0_rev_list) ;
    free(pens->z1_rev_list) ;

    free(pens) ;

    return ;
}


/* Add an element to a linked list in a Penalties structure */
void Add_To_List
(
    VSPINT      *entry, /* Heads of all the linked lists                      */
    VSPINT   *fwd_list, /* Contains all the forward linked lists              */
    VSPINT   *rev_list, /* Contains all the reverse linked lists              */
    VSPINT         new, /* Entry to be added to the linked list               */
    VSPINT      Aw_val  /* Remove from the ll at heads[Aw_val]                */
)
{
    /* If the list is empty, initialize it */
    VSPINT i = entry[Aw_val] ;
    if(i == EMPTY) {
        entry[Aw_val] = new ;
        fwd_list[new] = new ;
        rev_list[new] = new ;
        return ;
    }

    /* Change the linked list from "i <--> next" to "i <--> new <--> next" */
    VSPINT next = fwd_list[i] ;

    /* Link i to new */
    fwd_list[i]   = new ;
    rev_list[new] = i ;

    /* Link new to next */
    fwd_list[new] = next ;
    rev_list[next] = new ;

    return ;
}

void build_penalties
(
    VSPINT                 *z,
    Penalties     *Z_pens,
    VSPINT                *Aw,
    Problem      *problem
)
{
    #ifdef DEBUG
        static long long unsigned count = 0 ;
        count ++ ;
    #endif

    /* Local vars */
    VSPINT  i ;
    VSPINT  n           = problem->Graph->n ;

    VSPINT *z0_entry    = Z_pens->z0_entry ;
    VSPINT *z0_fwd_list = Z_pens->z0_fwd_list ;
    VSPINT *z0_rev_list = Z_pens->z0_rev_list ;

    VSPINT *z1_entry    = Z_pens->z1_entry ;
    VSPINT *z1_fwd_list = Z_pens->z1_fwd_list ;
    VSPINT *z1_rev_list = Z_pens->z1_rev_list ;


    /* Initialize the heads of the linked lists */
    for (i = 0 ; i <=n  ; i++)
    {
        z0_entry [i] = EMPTY ;
        z1_entry [i] = EMPTY ;
    }

    /* Build the linked lists */
    VSPINT max_Aw_z1 = -1 ;
    VSPINT min_Aw_z0 = n + 1 ;
    for (i = 0; i < n; i++) {

        VSPINT Awi = Aw[i] ;
        if (z[i])
        {
            if ( Awi > max_Aw_z1 ) max_Aw_z1 = Awi ;

            Add_To_List(z1_entry, z1_fwd_list, z1_rev_list, i, Awi) ;
        }
        else
        {
            if ( Awi < min_Aw_z0 ) min_Aw_z0 = Awi ;

            Add_To_List(z0_entry, z0_fwd_list, z0_rev_list, i, Awi) ;
        }
    }

    Z_pens->max_pen = max_Aw_z1 ;
    Z_pens->min_pen = min_Aw_z0 ;

    #ifdef DEBUG
        Check_Linked_Lists(z, Aw, Z_pens) ;
    #endif

    return ;
}

static void check_a_linked_list
(
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
            Error("Linked list not correct", DBG_INFO) ;
            return ;
        }
        if( z[p] != z_val )
        {
            Error("Linked list not correct", DBG_INFO) ;
            return ;
        }
        if( Aw[p] != Aw_val ) {
            Error("Linked list not correct", DBG_INFO) ;
            return ;
        }

        p = fwd_list[p] ;
    } while(p != entry) ;
}

/* Check that the linked lists are all correct */
void Check_Linked_Lists
(
    VSPINT            *z,
    VSPINT           *Aw,
    Penalties *Z_pen
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
            check_a_linked_list(entry , z , z0_fwd_list, z0_rev_list, 0,
                                    Aw_val, Aw, witness                      ) ;
        }

        /* Transverse the z1 linked list corresponding to Aw_val */
        entry = z1_entry[Aw_val] ;
        if(entry != EMPTY) {
            check_a_linked_list(entry , z , z1_fwd_list, z1_rev_list, 1,
                                    Aw_val, Aw, witness                      ) ;
        }
    }

    /* Check that each node is present in precisely one linked list */
    for(VSPINT i = 0 ; i < n ; i++) {
        if(witness[i] != 1) {
            Error("Problem in linked lists", DBG_INFO) ;
            break ;
        }
    }

    /* Free memory */
    free(witness) ;

    return ;
}






static void check_penalties_one_state
(
    State *Z,
    State *W,
    Penalties *Z_pen,
    VSPINT n
)
{
    static long long unsigned count = 0 ;
    count++ ;

    WARN_IF_NOT_DEBUG_MODE ;

    /* 1. Check the linked lists */
    Check_Linked_Lists
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
        Error("Error in min penalty.", DBG_INFO) ;
    }

    if(local_max_pen != Z_pen->max_pen)
    {
        Error("Error in max penalty.", DBG_INFO) ;
    }

    return ;
}

void check_penalties
(
    State *X,
    State *Y,
    Penalties *X_pen,
    Penalties *Y_pen,
    VSPINT n
)
{
    WARN_IF_NOT_DEBUG_MODE ;
    check_penalties_one_state(X, Y, X_pen, n) ;
    check_penalties_one_state(Y, X, Y_pen, n) ;

    return ;
}


void copy_penalties
(
    Penalties *from,
    Penalties *to
)
{
    VSPINT n = from->n ;
    
    to->n = n ;
    to->max_pen = from->max_pen ;
    to->min_pen = from->min_pen ;
    for(VSPINT i = 0 ; i < n ; i++)
    {
        to->z0_entry[i]    = from->z0_entry[i] ;
        to->z0_fwd_list[i] = from->z0_fwd_list[i] ;
        to->z0_rev_list[i] = from->z0_rev_list[i] ;

        to->z1_entry[i]    = from->z1_entry[i] ;
        to->z1_fwd_list[i] = from->z1_fwd_list[i] ;
        to->z1_rev_list[i] = from->z1_rev_list[i] ;
    }
    to->z0_entry[n]    = from->z0_entry[n] ;
    to->z1_entry[n]    = from->z1_entry[n] ;
    return ;
}