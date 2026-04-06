#include "shared_api.h"
#include "mc_internal.h"


/* ============================================================================
   ================ compare_value_index_floats ============================
   ============================================================================

 PURPOSE: 
 Function to compare the float values in two value_index struct. This is used
 when calling qsort to sort these structs.
*/
int compare_value_index_floats(const void *a, const void *b) 
{
    value_index* A = (value_index*) a;
    value_index* B = (value_index*) b;

    if (A->float_value > B->float_value)     return -1; /* a < b  */
    else if (A->float_value < B->float_value) return 1; /* a > b  */
    else return 0;                                      /* a == b */
}