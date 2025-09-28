#include "FM_Flips_Internal.h"


// Take the linked list corresponding to which_ll and builds O_up, A_up, and c_hi
// The function resets O_up and A_up so it is not neccesary to reset them before
void Build_O_up_And_A_up_And_c_hi_From_Linked_List
(
    FM_Flips_Struct    *data,
    VSPINT          which_ll
)
{
    SET_FUNCTION_COUNTER ;

    VSPINT *c                = data->c ;
    Direct_Access_List *O_up = data->O_up ;
    Max_Queue_v1       *A_up = data->A_up ;
    Doubly_LLs_Unsized *L0  = data->L0 ;
    VSPINT c_hi              = 0 ;
    VSPINT entry0            = data->L0->entries[which_ll] ;
    VSPINT n_queue ;

    // Get c_hi
    VSPINT i = entry0 ;
    do
    {
        if(c[i] > c_hi) c_hi = c[i] ;
    } while( (i = L0->fwd[i]) != entry0 ) ;
    data->c_hi = c_hi ;

    // Build O_up and A_up
    Direct_Access_List_Reset(O_up) ;
    Max_Queue_Reset_v1(A_up) ;
    n_queue = 0 ;
    i       = entry0 ; // TODO: Not needed but maybe it gives compiler a hint 
    do
    {
        if(c[i] == c_hi)
        {
            Direct_Access_List_Insert(O_up, i) ;
        }
        else
        {
            A_up->queue[++n_queue] = i ;
            A_up->lookup[i]        = n_queue ;
        }
    } while( (i = L0->fwd[i]) != entry0 ) ;
    A_up->n_queue = n_queue ;
    VSP_max_queue_build_v1(A_up) ;
}



// Take the linked list corresponding to which_ll and builds O_up, A_up, and c_hi
// The function resets O_up and A_up so it is not neccesary to reset them before
void Build_O_dn_And_A_dn_And_c_lo_From_Linked_List
(
    FM_Flips_Struct    *data,
    VSPINT          which_ll
)
{
    SET_FUNCTION_COUNTER ;
    VSPINT *c                = data->c ;
    Direct_Access_List *O_dn = data->O_dn ;
    Min_Queue_v1       *A_dn = data->A_dn ;
    Doubly_LLs_Unsized *L1   = data->L1 ;
    VSPINT c_lo              = data->gamma + 1 ;
    VSPINT entry1            = data->L1->entries[which_ll] ;
    VSPINT n_queue ;

    // Get c_hi
    VSPINT i = entry1 ;
    VSP_ASSERT(i != EMPTY) ;
    do
    {
        #ifdef VSP_DEBUG
            if(i < 0 || i >= data->n)
            {
                VSP_error("Out of bounds", DBG_INFO) ;
            }
        #endif
        if(c[i] < c_lo) c_lo = c[i] ;
    } while( (i = L1->fwd[i]) != entry1 ) ;
    data->c_lo = c_lo ;

    // Build O_up and A_up
    Direct_Access_List_Reset(O_dn) ;
    Min_Queue_Reset_v1(A_dn) ;
    i       = entry1 ; // TODO: Not needed but maybe it gives compiler a hint 
    n_queue = 0 ;
    do
    {
        if(c[i] == c_lo)
        {
            Direct_Access_List_Insert(O_dn, i) ;
        }
        else
        {
            A_dn->queue[++n_queue] = i ;
            A_dn->lookup[i]        = n_queue ;
        }
    } while( (i = L1->fwd[i]) != entry1 ) ;
    A_dn->n_queue = n_queue ;
    VSP_min_queue_build_v1(A_dn) ;
}