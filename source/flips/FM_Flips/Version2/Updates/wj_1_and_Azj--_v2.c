#include "FM_Flips_Internal.h"

/*
For the variable w:
    p_lo = min{ Azj : wj = 0 }
    p_hi = max{ Azj : wj = 1, j unlocked }
    c_lo = min{ cj s.t. j Azj = p_hi, wj = 1, j unlocked }
    c_hi = max{ cj s.t. j Azj = p_lo, wj = 0 }

If Azj_old = p_hi, then p_hi needs to be updated, and with it, O_dn, A_dn and c_lo
It is not possible for Azj_old to be > p_hi

*/


// Should be identical to v1
void FM_Flips_Struct_Update_wj_1_And_Azj_Decreases_By_1_v2
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
    data->Aw[j]-- ;

    // If j is locked, it has no effect on the data structures
    if(data->locked[j]) return ;

    Min_Queue_v1 *A_dn = data->A_dn ;
    Direct_Access_List *O_dn = data->O_dn ;

    VSPINT Azj_new = Azj_old - 1 ;

    // Move j through the linked lists
    Doubly_LLs_Unsized_Delete(data->L1, Azj_old, j) ;
    Doubly_LLs_Unsized_Insert(data->L1, Azj_new, j) ;

    if(Azj_old == data->p_hi)
    {
        if(A_dn->lookup[j] != EMPTY)
        {
            Min_Queue_Delete_v1(A_dn, j) ;
        }
        else // In this case, j should be in O_dn
        {
            VSP_ASSERT(O_dn->lookup[j] != EMPTY) ;
            Direct_Access_List_Delete_By_Value(O_dn, j) ;
            if(O_dn->size == 0)
            {
                if(A_dn->n_queue) // If A_dn is not empty
                {
                    // TODO: This can be done more efficiently
                    Build_O_dn_And_A_dn_And_c_lo_From_Linked_List(data, data->p_hi) ;
                }
                else // if A_dn is empty
                {
                    // TODO: This can be done more efficiently
                    data->p_hi-- ;
                    Build_O_dn_And_A_dn_And_c_lo_From_Linked_List(data, data->p_hi) ;
                }
            }
        }
    }
    #ifdef VSP_DEBUG
        FM_Flips_Struct_Check_v2(data, false) ;
    #endif
}