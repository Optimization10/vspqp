#include "core.h"
#include "shared_api.h"
#include "ds_api.h"

Doubly_LLs_Sized *Doubly_LLs_Sized_Init(VSPINT n)
{
    SET_FUNCTION_COUNTER ;
    TMP_ASSERT(n > 0) ;

    Doubly_LLs_Sized *ret = malloc( sizeof(*ret) ) ;
    ret->entries = malloc((n+1) * sizeof(VSPINT)) ;
    ret->sizes   = calloc( n+1  , sizeof(VSPINT)) ;
    ret->fwd     = malloc( n    * sizeof(VSPINT)) ;
    ret->rev     = malloc( n    * sizeof(VSPINT)) ;
    ret->n       = n ;

    VSPINT *entries = ret->entries ;
    for(VSPINT i = 0 ; i <= n ; i++) entries[i] = EMPTY ;
    return ret ;
}

void Doubly_LLs_Sized_Destroy(Doubly_LLs_Sized *lls)
{
    free(lls->entries) ;
    free(lls->sizes) ;
    free(lls->fwd) ;
    free(lls->rev) ;
    free(lls) ;
    return ;
}

// ASSUMEDL vtx is not in the list
void Doubly_LLs_Sized_Insert
(
    Doubly_LLs_Sized *lls,
    VSPINT        whichll,
    VSPINT            vtx
)
{
    Doubly_Linked_Lists_Insert
    (
        lls->entries, lls->fwd, lls->rev, whichll, vtx, lls->n
    ) ;
    lls->sizes[whichll]++ ;
}

// ASSUMED: vtx is in the list
void Doubly_LLs_Sized_Delete
(
    Doubly_LLs_Sized *lls,
    VSPINT        whichll,
    VSPINT            vtx
)
{
    Doubly_Linked_Lists_Delete
    (
        lls->entries, lls->fwd, lls->rev, whichll, vtx, lls->n
    ) ;
    lls->sizes[whichll]-- ;
}

// Check the following:
// The unlocked vertices are contained in one and only one ll
// All vertices in each ll have the correct zi and Awi value
// The linked list sizes are correct
// The n's in lls 0 and lls1 match
void Doubly_LLs_Sized_Check
(
    Doubly_LLs_Sized     *lls0,
    Doubly_LLs_Sized     *lls1,
    VSPINT                  *z,
    VSPINT                 *Aw,
    VSPINT             *locked
)
{
    SET_FUNCTION_COUNTER ;

    Doubly_Linked_Lists_Check
    (
        lls0->entries,
        lls0->fwd,
        lls0->rev,
        lls1->entries,
        lls1->fwd,
        lls1->rev,
        z, Aw, locked, lls0->n
    ) ;

    // Check that n's agree
    TMP_ASSERT(lls0->n == lls1->n) ;

    /*------------------------------------------------------------------------*/
    // Check that the sizes are correct
    /*------------------------------------------------------------------------*/
    VSPINT n = lls0->n, entry, i ;
    VSPINT *sizes0Dbg = calloc(n+1, sizeof(VSPINT)) ;
    for(VSPINT whichll = 0 ; whichll <= n ; whichll++)
    {
        entry = lls0->entries[whichll] ;
        if(entry == EMPTY) continue ;

        i = entry ;
        do
        {
            if(locked == NULL || !locked[i])
            {
                sizes0Dbg[whichll]++ ;
            }
        } while( (i = lls0->fwd[i]) != entry ) ;
    }

    // Get sizes of 1 lists
    VSPINT *sizes1Dbg = calloc(n+1, sizeof(VSPINT)) ;
    for(VSPINT whichll = 0 ; whichll <= n ; whichll++)
    {
        entry = lls1->entries[whichll] ;
        if(entry == EMPTY) continue ;

        i = entry ;
        do
        {
            if(locked == NULL || !locked[i])
            {
                sizes1Dbg[whichll]++ ;
            }
        } while( (i = lls1->fwd[i]) != entry ) ;
    }

    VSP_check_arrays_are_equal(sizes0Dbg, lls0->sizes, n) ;
    VSP_check_arrays_are_equal(sizes1Dbg, lls1->sizes, n) ;

    free(sizes0Dbg) ;
    free(sizes1Dbg) ;
}

// Set all entries  to EMPTY and all sizes to 0
void Doubly_LLs_Sized_Reset(Doubly_LLs_Sized *lls)
{
    VSPINT
    n = lls->n,
    *entries = lls->entries,
    *sizes   = lls->sizes ;
    for(VSPINT i = 0 ; i <= n ; i++)
    {
        entries[i] = EMPTY ;
        sizes  [i] = 0 ;
    }
    return ;
}

void Doubly_LLs_Sized_See_Entry
(
    Doubly_LLs_Sized *lls, VSPINT whichll, VSPINT *ret
)
{
    SET_FUNCTION_COUNTER ;
    TMP_ASSERT( IN_RANGE(whichll, 0, lls->n+1) ) ;
    Doubly_Linked_Lists_See_Entry(lls->entries, whichll, ret) ;
    return ;
}

VSPINT Doubly_LLs_Sized_Get_Size(Doubly_LLs_Sized *lls, VSPINT whichll)
{
    return lls->sizes[whichll] ;
}