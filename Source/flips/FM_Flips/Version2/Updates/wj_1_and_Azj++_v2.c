#include "FM_Flips_Internal.h"

// Should be identical to v1
void FM_Flips_Struct_Update_wj_1_And_Azj_Increases_By_1_v2
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

    if(data->locked[j]) return ;

    VSPINT Azj_new = Azj_old + 1 ;
    Min_Queue_v1 *A_dn = data->A_dn ;
    Direct_Access_List *O_dn = data->O_dn ;
    VSPINT *c = data->c ;

    // Move j through the linked lists
    Doubly_LLs_Unsized_Delete(data->L1, Azj_old, j) ;
    Doubly_LLs_Unsized_Insert(data->L1, Azj_new, j) ;

    if(Azj_old == data->p_hi)
    {
        // TODO: Could be done more efficiently here
        data->p_hi++ ;
        Build_O_dn_And_A_dn_And_c_lo_From_Linked_List(data, Azj_new) ;
    }
    else if(Azj_new == data->p_hi)
    {
        if(c[j] < data->c_lo)
        {
            // Set new c_lo
            data->c_lo = c[j] ;

            // TODO: Resetting O_dn and inserting can be done manually in one
            // loop

            // Transfer everything from O_dn to A_dn
            for(VSPINT i = 0 ; i < O_dn->size ; i++)
            {
                Min_Queue_Insert_v1( A_dn, O_dn->list[i] ) ;
            }

            // O_dn contains only j now
            Direct_Access_List_Reset(O_dn) ;
            Direct_Access_List_Insert(O_dn, j) ;
            
        }
        else if(c[j] == data->c_lo)
        {
            // Insert j into O_dn
            Direct_Access_List_Insert(O_dn, j) ;
        }
        else // cj > c_lo
        {
            // Insert j into A_dn
            Min_Queue_Insert_v1(A_dn, j) ;
        }
    }
    #ifdef VSP_DEBUG
        FM_Flips_Struct_Check_v2(data, false) ;
    #endif
}