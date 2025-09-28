#include "ds_api.h"

Direct_Access_List *Construct_Direct_Access_List(VSPINT max_size)
{
    Direct_Access_List *ret = malloc( sizeof(*ret) ) ;
    ret->list     = malloc(max_size * sizeof(VSPINT)) ;
    ret->lookup   = malloc(max_size * sizeof(VSPINT)) ;
    for(VSPINT i = 0 ; i < max_size ; i++) ret->lookup[i] = EMPTY ; // TODO: See if it's better to use a local variable for lookup
    ret->size     = 0 ;
    ret->max_size = max_size ;
    return ret ;
}


void Destruct_Direct_Access_List(Direct_Access_List *list)
{
    free(list->list) ;
    free(list->lookup) ;
    free(list) ;
}

// This function does not check that there is still room in the list. It is the
// responsibility of the user to not put too many things in there
void Direct_Access_List_Insert
(
    Direct_Access_List      *list,
    VSPINT              new_entry
)
{
    SET_FUNCTION_COUNTER ;
    VSP_ASSERT(list->size < list->max_size) ;

    list->list[ list->size ]  = new_entry ;
    list->lookup[ new_entry ] = list->size ;
    list->size++ ; 
}


// Delete entry in a given index
void Direct_Access_List_Delete_By_Idx
(
    Direct_Access_List       *list,
    VSPINT              idx_to_del
)
{
    SET_FUNCTION_COUNTER ;
    VSP_ASSERT(0 <= idx_to_del && idx_to_del < list->size) ;

    // Overwrite the index to be deleted with the last entry
    VSPINT node_to_del = list->list[idx_to_del] ;
    VSPINT last = list->list[ list->size - 1 ] ;
    list->list[ idx_to_del ] = last ;
    list->lookup[ last ]     = idx_to_del ;

    // TODO: This whole situation can be done better
    list->lookup[ node_to_del ] = EMPTY ;

    list->size-- ;
}

// Delete the value in a list, at unspecified location.
// ASSUMED: Value is in the list
void Direct_Access_List_Delete_By_Value
(
    Direct_Access_List       *list,
    VSPINT              val_to_del
)
{
    SET_FUNCTION_COUNTER ;
    VSP_ASSERT(0 <= val_to_del && val_to_del < list->max_size) ;
    VSP_ASSERT(0 <= list->lookup[val_to_del]
                 && list->lookup[val_to_del] < list->max_size) ;

    // TODO: Probably it's better to do this manually instead of calling
    // this function
    Direct_Access_List_Delete_By_Idx(list, list->lookup[val_to_del]) ;
}


// Reset the list so it contains nothing
void Direct_Access_List_Reset
(
    Direct_Access_List *list
)
{
    SET_FUNCTION_COUNTER ;
    VSPINT size = list->size ;
    // VSPINT *lookup->l
    for(VSPINT i = 0 ; i < size ; i++)
    {
        (list->lookup)[ list->list[i] ] = EMPTY ;
    }
    list->size = 0 ;

    #ifdef VSP_DEBUG
        for(VSPINT i = 0 ; i < list->max_size ; i++)
        {
            if(list->lookup[i] != EMPTY)
            {
                VSP_error("Error in resetting list", DBG_INFO) ;
                break ;
            }
        }
    #endif

    return ;
}

void Direct_Access_List_Print
(
    Direct_Access_List *list,
    VSPINT                *c, // Assumed that list contains indices into c
    char                *msg  // Ignored if NULL
)
{
    FILE *file = fopen("output/log_debug.txt", "a") ;

    // Print the title
    fprintf(file, "\n    Printing Direct Access List:  ") ;
    if(msg != NULL)
    {
        fprintf(file, "%s", msg) ;
    }
    fprintf(file, "\n") ;

    // Print the list
    if(list->size == 0)
    {
        fprintf(file, "List is Empty") ;
    }
    else
    {
        fprintf(file, "| ") ;
        for(VSPINT i = 0 ; i < list->size ; i++)
        {
            fprintf(file, "i = %i,  ci = %i  | ", list->list[i], c[list->list[i]] ) ;
        }
    }
    fclose(file) ;
}

void Direct_Access_List_Check
(
    Direct_Access_List *List
)
{
    SET_FUNCTION_COUNTER ;
    WARN_IF_NOT_DEBUG_MODE ;

    #ifdef VSP_DEBUG
        // Check that for each entry in the list, the pointer array agrees
        VSPINT *list   = List->list ;
        VSPINT *lookup = List->lookup ;
        VSPINT size    = List->size ;
        VSPINT maxSize = List->max_size ;
        VSPINT i, p, k ;


        if(size > maxSize)
            {
                ERR("Direct acces list too full") ;
            }

        for(p = 0 ; p < size ; p++)
        {
            k = list[p] ;
            if(k < 0)
            {
                ERR("Negative entry in direct access list") ;
            }

            if(k >= maxSize)
            {
                ERR("Too large of an entry in direct access list") ;
            }

            if(lookup[k] != p)
            {
                ERR("Lookup array incorrect in direct access list") ;
            }
        }
    #endif
}