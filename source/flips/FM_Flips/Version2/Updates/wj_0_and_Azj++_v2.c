#include "FM_Flips_Internal.h"

void FM_Flips_Struct_Update_wj_0_And_Azj_Increases_By_1_v2
(
    FM_Flips_Struct *data,
    VSPINT              j,
    VSPINT        Azj_old
)
{
    SET_FUNCTION_COUNTER ;
    #ifdef VSP_DEBUG
        FM_Flips_Struct_Check_v2(data, true) ;
    #endif
    data->Aw[j]++ ;

    // If j is locked, it has no effect on the data structures
    if(data->locked[j]) return ;

    VSPINT Azj_new = Azj_old + 1 ;
    Max_Queue_v1 *A_up = data->A_up ;
    Direct_Access_List *O_up = data->O_up ;

    // Move j through the linked lists
    Doubly_LLs_Unsized_Delete(data->L0, Azj_old, j) ;
    Doubly_LLs_Unsized_Insert(data->L0, Azj_new, j) ;


    if(A_up->lookup[j] != EMPTY) // j is in A_up
    {
        // TODO: these two steps can probably be combined into one
        Max_Queue_Delete_v1(A_up, j) ;
    }
    else if(O_up->lookup[j] != EMPTY) // j is in O_up
    {
        // Remove j from O_up
        // TODO: These two steps can probably be combined into one
        Direct_Access_List_Delete_By_Value(O_up, j) ;
        if(O_up->size == 0)
        {
            if(A_up->n_queue) // If A_up is not empty
            {
                // TODO: This can be done more efficiently because A_up
                // already exists
                Build_O_up_And_A_up_And_c_hi_From_Linked_List(data, data->p_lo) ;
            }
            else // If A_up is empty also
            {
                data->p_lo++ ;
                Build_O_up_And_A_up_And_c_hi_From_Linked_List(data, data->p_lo) ;
            }
        }
    }
    #ifdef VSP_DEBUG
        FM_Flips_Struct_Check_v2(data, false) ;
    #endif
}