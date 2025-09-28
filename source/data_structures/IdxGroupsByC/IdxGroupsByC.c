#include "ds_api.h"

IdxGroupsByC *IdxGroupsByC_Init
(
    VSPINT  gamma,
    VSPINT *maxGroupSizes,
    VSPINT *c,
    VSPINT nvtx, // Number of vertices at current level
    VSPINT targetzi,
    VSPINT targetAwi,
    VSPINT *z,
    VSPINT *Aw,
    char   var,
    MinMax_Heap_Flag flag
)
{
    SET_FUNCTION_COUNTER ;
    IdxGroupsByC *Groups = malloc( sizeof(*Groups) ) ;
    Groups->var       = var ;
    Groups->targetzi  = targetzi ;
    Groups->targetAwi = targetAwi ;
    Groups->z         = z ;
    Groups->Aw        = Aw ;

    Groups->flag = flag ;
    if(flag == MIN_HEAP_FLAG)
    {
        Groups->minHeap = Min_Heap_Init_v4(gamma, c, nvtx) ;
        Groups->maxHeap = NULL ;
    }
    else if(flag == MAX_HEAP_FLAG)
    {
        Groups->minHeap = NULL ;
        Groups->maxHeap = Max_Heap_Init_v4(gamma, c, nvtx) ;
    }
    else
    {
        #ifdef VSP_DEBUG
            ERR("Invalid heap choice") ;
        #endif
    }
    

    Groups->groups = malloc( (gamma + 1) * sizeof( VSPINT* ) ) ;
    Groups->sizes  = calloc( (gamma + 1), sizeof(VSPINT) ) ;

    for(VSPINT cval = 0 ; cval <= gamma ; cval++)
    {
        if(maxGroupSizes[cval])
        {
            Groups->groups[cval] = malloc( maxGroupSizes[cval] * sizeof(VSPINT) ) ;
        }
        else
        {
            Groups->groups[cval] = NULL ;
        }
    }

    Groups->ptrs = malloc( nvtx * sizeof(VSPINT) ) ;
    for(VSPINT i = 0 ; i < nvtx ; i++)
    {
        Groups->ptrs[i] = EMPTY ;
    }
    Groups->c         = c ;
    Groups->gamma     = gamma ;

    return Groups ;
}

void IdxGroupsByC_Destroy(IdxGroupsByC *Groups)
{
    SET_FUNCTION_COUNTER;

    // Free the heap based on the flag
    if (Groups->flag == MIN_HEAP_FLAG)
    {
        Min_Heap_Destroy_v4(Groups->minHeap);
    }
    else if (Groups->flag == MAX_HEAP_FLAG)
    {
        Max_Heap_Destroy_v4(Groups->maxHeap);
    }

    // Free each allocated group
    if (Groups->groups)
    {
        for (VSPINT cval = 0; cval <= Groups->gamma; cval++)
        {
            if(Groups->groups[cval] != NULL)
            {
                free(Groups->groups[cval]) ;
            }
        }
        free(Groups->groups);
    }

    // Free other allocated arrays
    free(Groups->sizes);
    free(Groups->ptrs);

    // Finally, free the main structure
    free(Groups);
}



void IdxGroupsByC_Insert
(
    IdxGroupsByC      *Groups,
    VSPINT        idxToInsert
)
{
    VSPINT cval = Groups->c[idxToInsert] ;

    // Insert cval into the heap (the insert functions check for membership)
    if(Groups->flag == MIN_HEAP_FLAG)
    {
        Min_Heap_Insert_v4(Groups->minHeap, cval) ;
    }
    else
    {
        Max_Heap_Insert_v4(Groups->maxHeap, cval) ;
    }

    // Insert idx into its group
    VSPINT sizeGroup = Groups->sizes[cval] ;
    Groups->groups[cval][sizeGroup] = idxToInsert ;
    Groups->ptrs[idxToInsert] = sizeGroup ;
    Groups->sizes[cval]++ ;

    return ;
}


void IdxGroupsByC_Delete
(
    IdxGroupsByC    *Groups,
    VSPINT         idxToDel
)
{
    SET_FUNCTION_COUNTER ;
    VSPINT  cval  = Groups->c[idxToDel] ;
    VSPINT *group = Groups->groups[cval] ;

    VSPINT *ptrs  = Groups->ptrs ;
    VSPINT  loc   = ptrs[idxToDel] ;

    ptrs[idxToDel] = EMPTY ;
    Groups->sizes[cval]-- ; 

    if(Groups->sizes[cval])
    {
        // In the group, overwrite idxToDel with the last thing in the group
        VSPINT k   = group[ Groups->sizes[cval] ] ;
        group[loc] = k ;
        ptrs[k]    = loc ;
    }
    else
    {
        if(Groups->flag == MIN_HEAP_FLAG)
        {
            Min_Heap_Delete_v4(Groups->minHeap, cval) ;
        }
        else
        {
            Max_Heap_Delete_v4(Groups->maxHeap, cval) ;
        }
    }
    return ;
}

// Get a random entry from the group associated with the cvalue at the top of
// the heap. If the struct is empty, or if this flip is not feasible, then we
// idxReturn is set to EMPTY
void IdxGroupsByC_See_Top_Cval_If_Feasible
(
    IdxGroupsByC    *Groups,
    VSPINT              wgt, // Current weight
    VSPINT                u,
    VSPINT         *cvalRet  // EMPTY if not feasible
)
{
    VSPINT cval ;
    if(Groups->flag == MIN_HEAP_FLAG)
    {
        if(Groups->minHeap->sizeHeap)
        {
            cval = Groups->minHeap->heap[1] ;
        }
        else
        {
            *cvalRet = EMPTY ;
            return ;
        }

        if(wgt != cval)
        {
            *cvalRet = cval ;
            return ;
        }
        else
        {
            *cvalRet = EMPTY ;
        }
    }
    else
    {
        if(Groups->maxHeap->sizeHeap)
        {
            cval = Groups->maxHeap->heap[1] ;
        }
        else
        {
            *cvalRet = EMPTY ;
            return ;
        }

        if(wgt + cval <= u)
        {
            *cvalRet = cval ;
        }
        else
        {
            *cvalRet = EMPTY ;
        }
    }
    return ;
}


// Get a random entry from the group associated with cval at the top of the heap
// ASSUMED: The group is non-empty
void IdxGroupsByC_Get_Random_Idx_Top_C
(
    IdxGroupsByC    *Groups,
    VSPINT          randInt, // Used for randInt % sizeGroup
    VSPINT       *idxReturn  // Return
)
{
    VSPINT cval ;
    if(Groups->flag == MIN_HEAP_FLAG)
    {
        cval = Groups->minHeap->heap[1] ;
    }
    else
    {
        cval = Groups->maxHeap->heap[1] ;
    }

    VSPINT *group   = Groups->groups[cval] ;
    VSPINT  size    = Groups->sizes[cval] ;
    VSPINT  randPtr = randInt % size ;
    *idxReturn      = group[randPtr] ;
    IdxGroupsByC_Delete(Groups, *idxReturn) ;

    return ;
}


void IdxGroupsByC_See_Random_Idx_Top_C
(
    IdxGroupsByC    *Groups,
    VSPINT          randInt, // Used for randInt % sizeGroup
    VSPINT       *idxReturn  // Return
)
{
    VSPINT cval ;
    if(Groups->flag == MIN_HEAP_FLAG)
    {
        cval = Groups->minHeap->heap[1] ;
    }
    else
    {
        cval = Groups->maxHeap->heap[1] ;
    }

    VSPINT *group   = Groups->groups[cval] ;
    VSPINT  size    = Groups->sizes[cval] ;
    VSPINT  randPtr = randInt % size ;
    *idxReturn      = group[randPtr] ;

    return ;
}


void IdxGroupsByC_Check
(
    IdxGroupsByC *Groups,
    VSPINT             n  // Size of z and Aw arrays
)
{
    static long long unsigned count = 0 ;
    count++ ;

    WARN_IF_NOT_DEBUG_MODE ;
    VSPINT i, p, cval, size, *group ;

    VSPINT gamma     = Groups->gamma ;
    VSPINT **groups  = Groups->groups ;
    VSPINT *sizes    = Groups->sizes ;
    const VSPINT *c        = Groups->c ;
    VSPINT *maxSizes = Groups->maxSizes ;

    const VSPINT *z        = Groups->z ;
    const VSPINT *Aw       = Groups->Aw ;
    VSPINT targetzi  = Groups->targetzi ;
    VSPINT targetAwi = Groups->targetAwi ;

    /*------------------------------------------------------------------------*/
    // Check that everything in each group has the correct zi, Awi, and ci value
    // Count the number of times an index occurs in any group (witnessGroups)
    // Record all the c's which occur in the group (csInGroup)
    /*------------------------------------------------------------------------*/
    // Counts the number of times each index is seen in the group
    VSPINT *witnessGroups = calloc(n, sizeof(*witnessGroups)) ;

    // Put all the c values seen. They will be added even if they're not
    // distinct
    VSPINT *csInGroup     = calloc(n, sizeof(*csInGroup)) ;
    VSPINT  numCsInGroup  = 0 ;

    for(cval = 0 ; cval <= gamma ; cval++)
    {
        group = groups[cval] ;
        size  = sizes[cval] ;
        if(group == NULL)
        {
            if(sizes[cval] != 0)
            {
                ERR("Group size not zero in empty group") ;
            }
            continue ;
        }

        for(p = 0 ; p < size ; p++)
        {
            i = group[p] ;
            witnessGroups[i]++ ;
            csInGroup[numCsInGroup++] = c[i] ;
            if(z[i] != targetzi)
            {
                ERR("zi doesn't match target zi") ;
            }
            if(Aw[i] != targetAwi)
            {
                ERR("Awi doesn't match target Awi") ;
            }

            if(c[i] != cval)
            {
                ERR("ci doesn't match group of i") ;
            }
        }
    }

    /*------------------------------------------------------------------------*/
    // Heap Checks:
    // Check that every c value in csInGroups occurs in the heap, and v.v
    // Call the heap checker, which checks for duplicates as well
    /*------------------------------------------------------------------------*/
    // Check that every c value in csInGroups is in the heap
    VSPINT *heap ;
    VSPINT  sizeHeap ;
    if(Groups->flag == MIN_HEAP_FLAG)
    {
        heap     = Groups->minHeap->heap ;
        sizeHeap = Groups->minHeap->sizeHeap ;
    }
    else
    {
        heap     = Groups->maxHeap->heap ;
        sizeHeap = Groups->maxHeap->sizeHeap ;
    }


    for(p = 0 ; p < numCsInGroup ; p++)
    {
        cval = csInGroup[p] ;
        if( !Is_Array_Member(cval, heap, 1, sizeHeap + 1) )
        {
            ERR("c value in groups but not in heap") ;
        }
    }


    // Check that every c value in the heap is in the groups
    for(p = 1 ; p <= sizeHeap ; p++)
    {
        cval = heap[p] ;
        if( !Is_Array_Member(cval, csInGroup, 0, numCsInGroup) )
        {
            ERR("c value in groups but not in heap") ;
        }
    }

    // Do extensive checks with the heap checker (includes verifying heap
    // property and uniqueness of entries)
    if(Groups->flag == MIN_HEAP_FLAG)
    {
        Min_Heap_Check_v4(Groups->minHeap) ;
    }
    else
    {
        Max_Heap_Check_v4(Groups->maxHeap) ;
    }

    /*------------------------------------------------------------------------*/
    // Ensure that all i with the target zi and Awi are in the groups, and there
    // are no duplicates in the groups
    /*------------------------------------------------------------------------*/
    VSPINT *witnessAccurate = calloc(n, sizeof(*witnessAccurate)) ;
    for(i = 0 ; i < n ; i++)
    {
        if(z[i] == targetzi && Aw[i] == targetAwi)
        {
            witnessAccurate[i]++ ;
        }
    }

    for(i = 0 ; i < n ; i++)
    {
        if(witnessAccurate[i] != witnessGroups[i])
        {
            ERR("Mismatch between groups and i which"
                "have target zi and Awi values") ;
        }

        if(witnessGroups[i] != 0 && witnessGroups[i] != 1)
        {
            ERR("An i is in a group more than once") ;
        }

        if(witnessAccurate[i] != 0 && witnessAccurate[i] != 1)
        {
            ERR("Very strange error") ;
        }
    }

    free(witnessGroups) ;
    free(witnessAccurate) ;
    free(csInGroup) ;
}

// ASSUMED: Groups is not empty
void IdxGroupsByC_Num_Idx_With_Top_C
(
    IdxGroupsByC    *Groups,
    VSPINT       *retNumIdx
)
{
    VSPINT topCval ;
    if(Groups->flag == MIN_HEAP_FLAG)
    {
        topCval = Groups->minHeap->heap[1] ;
    }
    else
    {
        topCval = Groups->maxHeap->heap[1] ;
    }
    *retNumIdx = Groups->sizes[topCval] ;


    return ;
}



void IdxGroupsByC_See_Random_Idx_Top_C_v5
(
    IdxGroupsByC    *Groups,
    VSPINT          randInt, // Used for randInt % sizeGroup
    VSPINT       *idxReturn  // Return
)
{
    VSPINT cval ;
    if(Groups->flag == MIN_HEAP_FLAG)
    {
        cval = Groups->minHeap->heap[1] ;
    }
    else
    {
        cval = Groups->maxHeap->heap[1] ;
    }

    VSPINT *group   = Groups->groups[cval] ;
    VSPINT  size    = Groups->sizes[cval] ;
    VSPINT  randPtr = randInt % size ;
    *idxReturn      = group[randPtr] ;

    return ;
}


void IdxGroupsByC_See_Random_Idx_Given_C
(
    IdxGroupsByC    *Groups,
    VSPINT          randInt, // Used for randInt % sizeGroup
    VSPINT        givenCval,
    VSPINT       *idxReturn  // Return
)
{
    VSPINT *group   = Groups->groups[givenCval] ;
    VSPINT  size    = Groups->sizes[givenCval] ;
    VSPINT  randPtr = randInt % size ;
    *idxReturn      = group[randPtr] ;

    return ;
}



// ASSUMED: Groups is not empty
void IdxGroupsByC_Num_Idx_With_Given_C
(
    IdxGroupsByC    *Groups,
    VSPINT        givenCval,
    VSPINT       *retNumIdx
)
{
    *retNumIdx = Groups->sizes[givenCval] ;
    return ;
}