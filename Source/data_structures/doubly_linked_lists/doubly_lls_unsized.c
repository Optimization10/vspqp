#include "core.h"
#include "shared_api.h"
#include "ds_api.h"



Doubly_LLs_Unsized *Doubly_LLs_Unsized_Init(VSPINT n)
{
    SET_FUNCTION_COUNTER;
    VSP_ASSERT(n > 0);

    Doubly_LLs_Unsized *ret = malloc(sizeof(*ret));
    ret->entries = malloc((n+1) * sizeof(VSPINT));
    ret->fwd     = malloc(n * sizeof(VSPINT));
    ret->rev     = malloc(n * sizeof(VSPINT));
    ret->n       = n;

    VSPINT *entries = ret->entries;
    for(VSPINT i = 0; i <= n; i++) entries[i] = EMPTY;
    return ret;
}

void Doubly_LLs_Unsized_Destroy(Doubly_LLs_Unsized *lls)
{
    free(lls->entries);
    free(lls->fwd);
    free(lls->rev);
    free(lls);
}

void Doubly_LLs_Unsized_Insert
(
    Doubly_LLs_Unsized *lls,
    VSPINT        whichll,
    VSPINT            vtx
)
{
    Doubly_Linked_Lists_Insert
    (
        lls->entries, lls->fwd, lls->rev, whichll, vtx, lls->n
    );
}

void Doubly_LLs_Unsized_Delete
(
    Doubly_LLs_Unsized *lls,
    VSPINT        whichll,
    VSPINT            vtx
)
{
    Doubly_Linked_Lists_Delete
    (
        lls->entries, lls->fwd, lls->rev, whichll, vtx, lls->n
    );
}

void Doubly_LLs_Unsized_Check
(
    Doubly_LLs_Unsized *lls0,
    Doubly_LLs_Unsized *lls1,
    VSPINT            *z,
    VSPINT           *Aw,
    VSPINT       *locked
)
{
    SET_FUNCTION_COUNTER;

    Doubly_Linked_Lists_Check
    (
        lls0->entries,
        lls0->fwd,
        lls0->rev,
        lls1->entries,
        lls1->fwd,
        lls1->rev,
        z, Aw, locked, lls0->n
    );

    // Check that n's agree
    VSP_ASSERT(lls0->n == lls1->n);
}

void Doubly_LLs_Unsized_Reset(Doubly_LLs_Unsized *lls)
{
    VSPINT n = lls->n;
    VSPINT *entries = lls->entries;

    for(VSPINT i = 0; i <= n; i++)
    {
        entries[i] = EMPTY;
    }
}

void Doubly_LLs_Unsized_See_Entry
(
    Doubly_LLs_Unsized *lls,
    VSPINT whichll,
    VSPINT *ret
)
{
    SET_FUNCTION_COUNTER;
    VSP_ASSERT(IN_RANGE(whichll, 0, lls->n+1));
    Doubly_Linked_Lists_See_Entry(lls->entries, whichll, ret);
}



VSPINT *Doubly_LLs_Unsized_Spread
(
    Doubly_LLs_Unsized *lls
)
{
    VSPINT *ret = Doubly_Linked_Lists_Spread
                    (
                        lls->entries,
                        lls->fwd,
                        lls->rev,
                        lls->n
                    ) ;
    return ret ;
}