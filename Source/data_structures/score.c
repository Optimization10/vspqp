#include "core.h"
#include "score.h"

/*----------------------------------------------------------------------------*/
// module score

// Purpose: score a point based on objective value and imbalance
// Scores are totally ordered lexicographically, first by obj and then by imb
/*----------------------------------------------------------------------------*/

// Return 1 if score 1 is highest, return 2 if score 2 is higher, and return 0
// if equal
VSPINT PointScore_Compare(PointScore *Score1, PointScore *Score2)
{
    VSPINT ret ;

    if( equal(Score1->obj, Score2->obj) )
    {
        if(Score1->imbal < Score2->imbal) ret = 1 ;
        else if(Score1->imbal > Score2->imbal) ret = 2 ;
        else ret = 0 ;
    }
    else if(Score1->obj > Score2->obj) ret = 1 ;
    else ret = 2 ;
    return ret ;
}

void PointScore_Update(PointScore *Score, State *X, State *Y, Problem *problem)
{
    VSPFLOAT obj = compute_obj(X, Y, problem, problem->Graph->max_vtx_cost) ;
    VSPINT   imb = IMBALANCE(X->weight, Y->weight) ;
    PointScore_SetObj(Score, obj) ;
    PointScore_SetImb(Score, imb) ;
}

void PointScore_SetObj(PointScore *Score, VSPFLOAT obj)
{
    Score->obj = obj ;
}

void PointScore_SetImb(PointScore *Score, VSPINT imb)
{
    Score->imbal = imb ;
}

bool PointScore_IsStrictlyBetter(PointScore *New, PointScore *Base)
{
    return (PointScore_Compare(New, Base) == 1) ;
}