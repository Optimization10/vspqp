#include "shared_api.h"
#include "mc_internal.h"


/* ============================================================================
   ================ VSP_compare_value_index_floats ============================
   ============================================================================

 PURPOSE: 
 Function to compare the float values in two value_index struct. This is used
 when calling qsort to sort these structs.
*/
int VSP_compare_value_index_floats(const void *a, const void *b) 
{
    VSP_value_index* A = (VSP_value_index*) a;
    VSP_value_index* B = (VSP_value_index*) b;

    if (A->float_value > B->float_value)     return -1; /* a < b  */
    else if (A->float_value < B->float_value) return 1; /* a > b  */
    else return 0;                                      /* a == b */
}


/* ============================================================================
   ================ VSP_compare_value_index_floats ============================
   ============================================================================

 PURPOSE: 
 Function to compare the float values in two value_index struct. This is used
 when calling qsort to sort these structs.
*/
int VSP_compare_val_ind_structs_inc(const void *a, const void *b) 
{
    VSP_value_index* A = (VSP_value_index*) a;
    VSP_value_index* B = (VSP_value_index*) b;

    if (A->float_value > B->float_value)     return   1; /* a > b  */
    else if (A->float_value < B->float_value) return -1; /* a < b  */
    else return 0;                                       /* a == b */
}

