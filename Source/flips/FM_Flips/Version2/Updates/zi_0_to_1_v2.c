#include "FM_Flips_Internal.h"




void FM_Flips_Struct_Update_zi_0_To_1_Where_i_Is_In_O_up_v2
(
    FM_Flips_Struct *data,
    VSPINT              i,
    VSPINT          idx_i  // index of i in O_up
)
{
    SET_FUNCTION_COUNTER ;
    VSP_ASSERT( data->z[i] == 0 ) ;

    #ifdef VSP_DEBUG
        FM_Flips_Struct_Check_v2(data, true) ;
    #endif

    data->z[i] = 1 ;
    data->locked[i] = 1 ;
    data->weight += data->c[i] ;


    VSPINT *Aw = data->Aw ; 
    data->locked[i] = 1 ;
    Max_Queue_v1 *A_up = data->A_up ;
    Direct_Access_List *O_up = data->O_up ;


    // Delete i from the lls. Since we lock upflips, no need to add to another
    Doubly_LLs_Unsized_Delete(data->L0, Aw[i], i)   ;

    // Update O_up, A_up, p_lo, c_hi if neccesary
    Direct_Access_List_Delete_By_Idx(data->O_up, idx_i) ;
    if(O_up->size == 0)
    {
        if(A_up->n_queue)
        {
            Transfer_Nodes_From_Nonempty_A_up_To_The_Empty_O_Up(A_up, O_up) ;
            data->c_hi = data->c[ data->O_up->list[0] ] ;
        }
        else // No more unlocked nodes corresponding to p_lo.
        {
            VSP_ASSERT(data->L0->entries[ data->p_lo] == EMPTY) ;
            // Set_New_p_lo_And_c_hi_And_O_up_And_A_up_When_zi_0_To_1(data, data->p_lo) ;
            // Search for new p_lo
            VSPINT p_lo = data->p_lo ;
            VSPINT n    = data->n ;
            while( ++p_lo <= n && data->L0->entries[p_lo] == EMPTY ) ;
            data->p_lo = p_lo ;

            if(p_lo != n + 1) // If found unlocked upflip
            {
                Build_O_up_And_A_up_And_c_hi_From_Linked_List(data, p_lo) ;
            }
            else
            {
                data->c_hi = 0 ;
            }
        }
    }

    #ifdef VSP_DEBUG
        FM_Flips_Struct_Check_v2(data, false) ;
    #endif
    return ;
}


void Transfer_Nodes_From_Nonempty_A_up_To_The_Empty_O_Up
(
    Max_Queue_v1        *A_up,
    Direct_Access_List  *O_up
)
{
    VSPINT *c     = A_up->priorities ;
    VSPINT  c_up  = c[ A_up->queue[1] ] ;
    VSPINT *queue = A_up->queue ;

    // Keep taking the root of the heap and placing it in O_up
    while(c[queue[1]] == c_up && A_up->n_queue > 0)
    {
        VSPINT i = A_up->queue[1] ;

        /* Delete the root from the heap and place it in O_up */
        Max_Queue_Delete_v1(A_up, i) ; // TODO: We know we're deleting the root
        Direct_Access_List_Insert(O_up, i) ;
    }
}