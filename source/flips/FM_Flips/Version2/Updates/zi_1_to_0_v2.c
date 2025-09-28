#include "FM_Flips_Internal.h"

void FM_Flips_Struct_Update_zi_1_To_0_Where_i_Is_In_O_dn_v2
(
    FM_Flips_Struct *data,
    VSPINT              i,
    VSPINT          idx_i  // index of i in O_dn
)
{
    SET_FUNCTION_COUNTER ;
    VSP_ASSERT( !data->locked[i] ) ;
    VSP_ASSERT( data->z[i] == 1 ) ;

    #ifdef VSP_DEBUG
        FM_Flips_Struct_Check_v2(data, true) ;
    #endif
    data->z[i] = 0 ;
    data->locked[i] = 1 ;
    data->weight -= data->c[i] ;

    Min_Queue_v1 *A_dn = data->A_dn ;
    Direct_Access_List *O_dn = data->O_dn ;
    VSPINT Awi = data->Aw[i] ; 
    VSPINT  *c = data->c ;
    VSPINT j, k ;

    // Delete i from the linked lists
    Doubly_LLs_Unsized_Delete(data->L1, Awi, i) ;


    // Update O_dn, A_dn, p_hi, c_lo if neccesary
    Direct_Access_List_Delete_By_Idx(O_dn, idx_i) ;
    if(O_dn->size == 0)
    {
        if(A_dn->n_queue)
        {
            Transfer_Nodes_From_Nonempty_A_dn_To_The_Empty_O_dn(A_dn, O_dn) ;
            data->c_lo = data->c[ data->O_dn->list[0] ] ;
        }
        else // No more nodes with p_hi
        {
            VSP_ASSERT(data->L1->entries[ data->p_hi] == EMPTY) ;
             
             // Search for new p_hi
            VSPINT p_hi = data->p_hi ;

            // REMARK: This code assumed that the flip corresponded to p_hi
            while(--p_hi >= 0 && data->L1->entries[p_hi] == EMPTY) ;
            data->p_hi = p_hi ;

            if(p_hi != -1)
            {
                Build_O_dn_And_A_dn_And_c_lo_From_Linked_List(data, p_hi) ;
            }
            else
            {
                data->c_lo = data->gamma + 1 ;
            }
        }
    }

    // NOTE: Not needed since downflips can get locked now
    // Place i in the upflip structures if neccesary
    // Max_Queue_v1 *A_up = data->A_up ;
    // Direct_Access_List *O_up = data->O_up ;
    // if(Awi == data->p_lo)
    // {
    //     if(c[i] > data->c_hi)
    //     {
    //         // New O_up consists of just i. Every old member is transferred
    //         // to A_up
    //         for(j = 0 ; j < data->O_up->size ; j++)
    //         {
    //             k = O_up->list[j] ;
    //             Max_Queue_Insert_v1(A_up, k) ;
    //         }
    //         Direct_Access_List_Reset(data->O_up) ;
    //         Direct_Access_List_Insert(data->O_up, i) ;
    //         data->c_hi = c[i] ;

    //     }
    //     else if (c[i] == data->c_hi)
    //     {
    //         // Insert i into O_up
    //         Direct_Access_List_Insert(data->O_up, i) ;
    //     }
    //     else // ci < c_hi
    //     {
    //         // Insert i into A_up
    //         Max_Queue_Insert_v1(A_up, i) ;
    //     }
    // }
    // else if(Awi < data->p_lo)
    // {
    //     Direct_Access_List_Reset(O_up) ;
    //     Max_Queue_Reset_v1(A_up) ;
    //     Direct_Access_List_Insert(O_up, i) ;
    //     data->p_lo = Awi ;
    //     data->c_hi = c[i] ;
    // }
    #ifdef VSP_DEBUG
        FM_Flips_Struct_Check_v2(data, false) ;
    #endif
}


// Take all nodes corresponding whose c-value matches the top of the heap, and
// place those nodes into O_dn
void Transfer_Nodes_From_Nonempty_A_dn_To_The_Empty_O_dn
(
    Min_Queue_v1        *A_dn,
    Direct_Access_List  *O_dn
)
{
    VSPINT *c     = A_dn->priorities ;
    VSPINT  c_dn  = c[ A_dn->queue[1] ] ;
    VSPINT *queue = A_dn->queue ;

    // Keep taking the root of the heap and placing it in O_dn
    while(c[queue[1]] == c_dn && A_dn->n_queue > 0)
    {
        VSPINT i = A_dn->queue[1] ;

        /* Delete the root from the heap and place it in O_dn */
        Min_Queue_Delete_v1(A_dn, i) ; // TODO: We know we're deleting the root
        Direct_Access_List_Insert(O_dn, i) ;
    }
}