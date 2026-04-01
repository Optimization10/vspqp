#include "shared_api.h"
// #include "flips_internal_v1.h"
#include "ds_api.h"

/***** Declarations ***********************************************************/
void VSP_add_to_list_v1
(
    VSPINT      *entry, /* Heads of all the linked lists                      */
    VSPINT   *fwd_list, /* Contains all the forward linked lists              */
    VSPINT   *rev_list, /* Contains all the reverse linked lists              */
    VSPINT         new, /* Entry to be added to the linked list               */
    VSPINT      Aw_val  /* Remove from the ll at heads[Aw_val]                */
) ;

void VSP_check_linked_lists
(
    VSPINT* z,
    VSPINT *Aw,
    VSP_penalties *Z_pen
) ;

/******** Definitions *********************************************************/
/* Fill out the fields in Z_pen */


/**************** VSP_penalties ********************************************/
VSP_penalties *VSP_construct_penalties(VSPINT n)
{
    VSP_penalties *pens = malloc(sizeof(*pens)) ;

    pens->z0_entry = malloc( (n + 1) * sizeof(VSPINT) ) ;
    pens->z1_entry = malloc( (n + 1) * sizeof(VSPINT) ) ;

    pens->z0_fwd_list = malloc( n * sizeof(VSPINT) ) ;
    pens->z1_fwd_list = malloc( n * sizeof(VSPINT) ) ;

    pens->z0_rev_list = malloc( n * sizeof(VSPINT) ) ;
    pens->z1_rev_list = malloc( n * sizeof(VSPINT) ) ;

    pens->n = n ;

    return pens ;
}

void VSP_destruct_penalties(VSP_penalties *pens) {
    free(pens->z0_entry) ;
    free(pens->z1_entry) ;

    free(pens->z0_fwd_list) ;
    free(pens->z1_fwd_list) ;

    free(pens->z0_rev_list) ;
    free(pens->z1_rev_list) ;

    free(pens) ;

    return ;
}

void VSP_build_penalties
(
    VSPINT                 *z,
    VSP_penalties     *Z_pens,
    VSPINT                *Aw,
    VSP_problem      *problem
)
{
    #ifdef VSP_DEBUG
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

            VSP_add_to_list_v1(z1_entry, z1_fwd_list, z1_rev_list, i, Awi) ;
        }
        else
        {
            if ( Awi < min_Aw_z0 ) min_Aw_z0 = Awi ;

            VSP_add_to_list_v1(z0_entry, z0_fwd_list, z0_rev_list, i, Awi) ;
        }
    }

    Z_pens->max_pen = max_Aw_z1 ;
    Z_pens->min_pen = min_Aw_z0 ;

    #ifdef VSP_DEBUG
        VSP_check_linked_lists(z, Aw, Z_pens) ;
    #endif

    return ;
}

/* Add an element to a linked list in a VSP_penalties structure */
void VSP_add_to_list_v1
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