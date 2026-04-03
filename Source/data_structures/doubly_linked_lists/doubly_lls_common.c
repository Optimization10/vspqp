#include "core.h"
#include "ds_api.h"
// Functions that take primitive data types and can be re-used in many places

void Doubly_Linked_Lists_Insert
(
    VSPINT   *entries,
    VSPINT       *fwd,
    VSPINT       *rev,
    VSPINT   which_ll,
    VSPINT  to_insert,
    VSPINT          n
)
{
    SET_FUNCTION_COUNTER ;
    #ifdef VSP_DEBUG
        if(Doubly_Linked_Lists_Is_Member_Any_List(entries, fwd, rev, to_insert, n))
        {
            VSP_error("Inserting member which is already there", DBG_INFO) ;
        }
    #endif

    /* If the list is empty, initialize it */
    VSPINT i = entries[which_ll] ;
    if(i == EMPTY)
    {
        entries[which_ll] = to_insert ;
        fwd    [to_insert] = to_insert ;
        rev    [to_insert] = to_insert ;
        return ;
    }

    // Change the linked list from "i <--> next" to 
    // "i <--> to_insert <--> next"
    VSPINT next = fwd[i] ;

    /* Link i to to_insert */
    fwd[i]         = to_insert ;
    rev[to_insert] = i ;

    /* Link to_insert to next */
    fwd[to_insert] = next ;
    rev[next]      = to_insert ;
}

void Doubly_Linked_Lists_Delete
(
    VSPINT   *entries,
    VSPINT       *fwd,
    VSPINT       *rev,
    VSPINT   which_ll,
    VSPINT  to_delete,
    VSPINT          n
)
{
    SET_FUNCTION_COUNTER ;
    
    #ifdef VSP_DEBUG
        if(!Doubly_Linked_Lists_Is_Member_Specific_List(entries, fwd, rev, which_ll, to_delete, n))
        {
            VSP_error("Attempting to delete something that's not in the ll", DBG_INFO) ;
        }
    #endif

    VSPINT next = fwd[to_delete] ;
    if(next != to_delete) // size of list is > 1
    {
        entries[which_ll] = next ;
        VSPINT prev = rev[to_delete] ;
        fwd[prev]   = next ;   
        rev[next]   = prev ;
    }
    else /* ind_flip is the only entry in its ll */
    {
        entries[which_ll] = EMPTY ;
    }
}

bool Doubly_Linked_Lists_Is_Member_Specific_List
(
    VSPINT *entries,
    VSPINT     *fwd,
    VSPINT     *rev,
    VSPINT which_ll,
    VSPINT   member,
    VSPINT        n
)
{
    SET_FUNCTION_COUNTER ;
    VSP_ASSERT(0 <= which_ll && which_ll <= n) ;
    VSP_ASSERT(0 <= member   && member   <  n) ;

    VSPINT entry = entries[which_ll] ;
    if(entry == EMPTY) return false ;

    bool ret = false ;
    VSPINT i = entry ;
    do
    {
        if(i == member)
        {
            ret = true ;
            break ;
        }
        VSP_ASSERT(0 <= i && i < n) ;
    } while( (i = fwd[i]) != entry) ;
    
    return ret ;
}

bool Doubly_Linked_Lists_Is_Member_Any_List
(
    VSPINT *entries,
    VSPINT     *fwd,
    VSPINT     *rev,
    VSPINT   member,
    VSPINT        n
)
{
    for(VSPINT i = 0 ; i <= n ; i++)
    {
        if( Doubly_Linked_Lists_Is_Member_Specific_List(entries, fwd, rev, i, member, n) )
        {
            return true ;
        }
    }
    return false ;
}

// Increment "spread" twice in all entries i for i in whichll. The first
// increment is done by following the fwd list, and the second increment is done
// by following the rev list
void Doubly_Linked_Lists_Spread_One_List
(
    VSPINT *entries,
    VSPINT     *fwd,
    VSPINT     *rev,
    VSPINT  whichll,
    VSPINT  *spread,
    VSPINT        n
)
{
    SET_FUNCTION_COUNTER ;
    VSPINT entry = entries[whichll] ;
    if(entry == EMPTY) return ;
    
    VSPINT i = entry ;
    do
    {
        spread[i]++ ;
        VSP_ASSERT(0 <= i && i < n) ;
    } while( (i = fwd[i]) != entry) ;

    i = entry ;
    do
    {
        spread[i]++ ;
        VSP_ASSERT(0 <= i && i < n) ;
    } while( (i = rev[i]) != entry) ;
    return ;
}

// Return an array called "spread". spread[i] = 0 if i is not in any of the
// linked lists. spread[i] = 2 if it is in a linked list
VSPINT *Doubly_Linked_Lists_Spread
(
    VSPINT *entries,
    VSPINT     *fwd,
    VSPINT     *rev,
    VSPINT        n
)
{
    SET_FUNCTION_COUNTER ;
    VSPINT *spread = calloc(n, sizeof(*spread)) ;

    for(VSPINT whichll = 0 ; whichll < n ; whichll++)
    {
        Doubly_Linked_Lists_Spread_One_List(entries, fwd, rev, whichll, spread, n) ;
    }

    #ifdef VSP_DEBUG
        for(VSPINT i = 0 ; i < n ; i++)
        {
            if(spread[i] != 0 && spread[i] != 2)
            {
                VSP_error("Problem in linked lists", DBG_INFO) ;
            }
        }
    #endif
    return spread ;
}

// Check the following:
// The unlocked vertices are contained in one and only one ll
// All vertices in each ll have the correct zi and Awi value
void Doubly_Linked_Lists_Check
(
    VSPINT *entries0,
    VSPINT     *fwd0,
    VSPINT     *rev0,
    VSPINT *entries1,
    VSPINT     *fwd1,
    VSPINT     *rev1,
    VSPINT        *z,
    VSPINT       *Aw,
    VSPINT   *locked,  // If NULL, ignore
    VSPINT         n
)
{
    static long long unsigned count = 0 ;
    count ++ ;

    WARN_IF_NOT_DEBUG_MODE ;
    VSPINT i, entry ;

    /*------------------------------------------------------------------------*/
    // Check that each vertex is contained in one and only one list
    /*------------------------------------------------------------------------*/
    VSPINT *L0spread = Doubly_Linked_Lists_Spread(entries0, fwd0, rev0, n) ;
    VSPINT *L1spread = Doubly_Linked_Lists_Spread(entries1, fwd1, rev1, n) ;

    for(i = 0 ; i < n ; i++)
    {
        VSPINT in_L0  = (L0spread[i] == 2) ;
        VSPINT in_L1  = (L1spread[i] == 2) ;

        if(locked == NULL)
        {
            if(in_L0 + in_L1 != 1)
            {
                VSP_error("Error in linked lists", DBG_INFO) ;
                break ;
            }
        }
        else
        {
            VSPINT is_lkd = (locked[i]   == 1) ;
            if(in_L0 + in_L1 + is_lkd != 1)
            {
                VSP_error("Error in linked lists", DBG_INFO) ;
                break ;
            }
        }
    }
    free(L0spread) ;
    free(L1spread) ;

    /*------------------------------------------------------------------------*/
    // Verify that in each list, the zi and Awi vals of the vertices are correct
    /*------------------------------------------------------------------------*/
    // Check the 0 lists
    for(VSPINT whichll = 0 ; whichll <= n ; whichll++)
    {
        entry = entries0[whichll] ;
        if(entry == EMPTY) continue ;

        i = entry ;
        do
        {
            if(z[i] != 0 || Aw[i] != whichll)
            {
                VSP_error("Problem in linked lists", DBG_INFO) ;
            }

            if(locked != NULL && locked[i])
            {
                ERR("Locked vertex in ll") ;
            }
        } while( (i = fwd0[i]) != entry ) ;
    }

    // Check the 1 lists
    for(VSPINT whichll = 0 ; whichll <= n ; whichll++)
    {
        entry = entries1[whichll] ;
        if(entry == EMPTY) continue ;

        i = entry ;
        do
        {
            if(z[i] != 1 || Aw[i] != whichll)
            {
                VSP_error("Problem in linked lists", DBG_INFO) ;
            }

            if(locked != NULL && locked[i])
            {
                ERR("Locked vertex in ll") ;
            }
        } while( (i = fwd1[i]) != entry ) ;
    }
}

void Doubly_Linked_Lists_See_Entry
(
    VSPINT *entries,
    VSPINT  whichll,
    VSPINT     *ret
)
{
    SET_FUNCTION_COUNTER ;
    VSP_ASSERT(whichll >= 0) ;
    *ret = entries[whichll] ;
    return ;
}