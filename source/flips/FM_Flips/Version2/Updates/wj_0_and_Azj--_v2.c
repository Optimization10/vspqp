#include "FM_Flips_Internal.h"


// Only difference from v1 is that nothing is done if wj is locked
void FM_Flips_Struct_Update_wj_0_And_Azj_Decreases_By_1_v2
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

    VSPINT Azj_new = Azj_old - 1 ;
    Max_Queue_v1 *A_up = data->A_up ;
    Direct_Access_List *O_up = data->O_up ;
    VSPINT *c = data->c ;

    // Move j through the linked lists
    Doubly_LLs_Unsized_Delete(data->L0, Azj_old, j) ;
    Doubly_LLs_Unsized_Insert(data->L0, Azj_new, j) ;

    if( Azj_old == data->p_lo )
    {
        data->p_lo-- ;
        data->c_hi = c[j] ;
        Direct_Access_List_Reset(O_up) ;
        Max_Queue_Reset_v1(A_up) ;
        Direct_Access_List_Insert(O_up, j) ;
    }
    else if( Azj_new == data->p_lo )
    {
        if(c[j] < data->c_hi)
        {
            Max_Queue_Insert_v1(A_up, j) ;
        }
        else if(c[j] == data->c_hi)
        {
            Direct_Access_List_Insert(O_up, j) ;
        }
        else // if c[j] > data->c_hi
        {
            // Transfer everything from O_up to A_up
            data->c_hi = c[j] ;
            for(VSPINT i = 0 ; i < O_up->size ; i++)
            {
                VSPINT node = O_up->list[i] ;
                Max_Queue_Insert_v1(A_up, node) ;
            }
            Direct_Access_List_Reset(O_up) ;
            Direct_Access_List_Insert(O_up, j) ;
        }
    }
    #ifdef VSP_DEBUG
        FM_Flips_Struct_Check_v2(data, false) ;
    #endif
}