#include "FM_Flips_Internal.h"

FM_Flips_Struct *Construct_FM_Flips_Struct
(
    VSPINT            *z,
    VSPINT           *Az,
    VSPINT           *Aw,
    VSPINT            *c,
    VSPINT             n,
    VSPINT        weight,
    VSPINT         gamma,
    VSPINT           *Ai,
    VSPLONG *Ap,
    VSPINT             u
)
{
    FM_Flips_Struct *ret = malloc( sizeof(*ret) ) ;
    ret->z  = z  ;
    ret->Az = Az ;
    ret->Aw = Aw ;
    ret->c  = c  ;
    ret->n  = n  ;
    ret->weight = weight ;
    ret->gamma  = gamma ;

    ret->Ai = Ai ;
    ret->Ap = Ap ;

    ret->L0 = Doubly_LLs_Unsized_Init(n) ;
    ret->L1 = Doubly_LLs_Unsized_Init(n) ;

    ret->O_up = Construct_Direct_Access_List(n) ;
    ret->O_dn = Construct_Direct_Access_List(n)  ;

    ret->A_up = VSP_construct_max_queue_v1(n, c) ;
    ret->A_dn = VSP_construct_min_queue_v1(n, c) ;

    ret->u = u ;

    ret->locked = calloc(n, sizeof(VSPINT)) ;
    ret->var = '\0' ;

    return ret ;
}

// Free everything except z, Aw, c
void Destruct_FM_Flips_Struct
(
    FM_Flips_Struct *data
)
{
    Doubly_LLs_Unsized_Destroy(data->L0) ;
    Doubly_LLs_Unsized_Destroy(data->L1) ;

    Destruct_Direct_Access_List(data->O_up) ;
    Destruct_Direct_Access_List(data->O_dn) ;

    VSP_destruct_max_queue_v1(data->A_up) ;
    VSP_destruct_min_queue_v1(data->A_dn) ;

    free(data->locked) ;
    free(data) ;
}

// If I want to use this to reset the struct, I just need to reset the locked
// arrays manually. This function does not touch the locked array
void FM_Flips_Struct_Build
(
    FM_Flips_Struct *data,
    VSPINT          gamma
)
{
    SET_FUNCTION_COUNTER ;

    /* Local vars */
    VSPINT *z  = data->z ;
    VSPINT *Aw = data->Aw ;
    VSPINT *c  = data->c ;
    VSPINT  n  = data->n ;
    VSPINT i ;

    /*------------------------------------------------------------------------*/
    /*      Initialize L0, L1, p_lo, p_hi                                     */
    /*------------------------------------------------------------------------*/
    Doubly_LLs_Unsized *L0 = data->L0 ;
    Doubly_LLs_Unsized *L1 = data->L1 ;
    VSPINT  p_hi            = -1 ;
    VSPINT  p_lo            = n + 1 ;

    /* Build the linked lists */
    for (i = 0; i < n; i++) {

        VSPINT Awi = Aw[i] ;
        if (z[i])
        {
            Doubly_LLs_Unsized_Insert(L1, Awi, i) ;
            if ( Awi > p_hi ) p_hi = Awi ;
        }
        else
        {
            Doubly_LLs_Unsized_Insert(L0, Awi, i) ;
            if ( Awi < p_lo ) p_lo = Awi ;
        }
    }
    data->p_lo = p_lo ;
    data->p_hi = p_hi ;

    /*------------------------------------------------------------------------*/
    /*      Initialize O_up, A_up, c_hi                                       */
    /*------------------------------------------------------------------------*/
    Direct_Access_List *O_up = data->O_up ;
    Max_Queue_v1       *A_up = data->A_up ;
    VSPINT c_hi   = 0 ;
    VSPINT entry0 = L0->entries[p_lo] ;
    VSPINT n_queue ;

    // Get c_hi
    i = entry0 ;
    do
    {
        if(c[i] > c_hi) c_hi = c[i] ;
    } while( (i = L0->fwd[i]) != entry0 ) ;
    data->c_hi = c_hi ;

    // Build O_up and A_up
    n_queue = 0 ;
    i = entry0 ; // TODO: This line not needed, but I think it may give compiler
                 // a hint for optimization
    do
    {
        if(c[i] == c_hi)
        {
            Direct_Access_List_Insert(O_up, i) ;
        }
        else
        {
            A_up->queue[++n_queue] = i ;
            A_up->lookup[i] = n_queue ;
        }
    } while( (i = L0->fwd[i]) != entry0 ) ;
    A_up->n_queue = n_queue ;
    VSP_max_queue_build_v1(A_up) ;

    /*------------------------------------------------------------------------*/
    /*      Initialize O_dn, A_dn, c_lo                                       */
    /*------------------------------------------------------------------------*/
    Direct_Access_List *O_dn = data->O_dn ;
    Min_Queue_v1       *A_dn = data->A_dn ;
    VSPINT c_lo   = gamma + 1 ;
    VSPINT entry1 = L1->entries[p_hi] ;

    if(p_hi == EMPTY)
    {
        data->c_lo = c_lo ;
        data->p_hi = p_hi ;
        A_dn->n_queue = 0 ;
        O_dn->size    = 0 ;
        return ;
    }

    // Get c_lo
    i = entry1 ;
    do
    {
        if(c[i] < c_lo) c_lo = c[i] ;
    } while( (i = L1->fwd[i]) != entry1 ) ;
    data->c_lo = c_lo ;

    // Build O_dn and A_dn
    n_queue = 0 ;
    i = entry1 ; // TODO: This line not needed, but I think it may give compiler
                 // a hint for optimization
    do
    {
        if(c[i] == c_lo)
        {
            Direct_Access_List_Insert(O_dn, i) ;
        }
        else
        {
            A_dn->queue[++n_queue] = i ;
            A_dn->lookup[i] = n_queue ;
        }
    } while( (i = L1->fwd[i]) != entry1 ) ;
    A_dn->n_queue = n_queue ;
    VSP_min_queue_build_v1(A_dn) ;

    // Locked array should'be been set to all 0's by calloc

    return ;
}


// Reset linked lists, heaps, direct access lists, and locked. Everything
// else is left as is. After this reset, the struct is ready to go into build
void FM_Flips_Struct_Reset
(
    FM_Flips_Struct *data
)
{
    Doubly_LLs_Unsized_Reset(data->L0) ;
    Doubly_LLs_Unsized_Reset(data->L1) ;

    Max_Queue_Reset_v1(data->A_up) ;
    Min_Queue_Reset_v1(data->A_dn) ;

    Direct_Access_List_Reset(data->O_dn) ;
    Direct_Access_List_Reset(data->O_up) ;

    VSPINT n = data->n ;
    VSPINT *locked = data->locked ;
    for(VSPINT i = 0 ; i < n ; i++)
    {
        locked[i] = 0 ;
    }

}