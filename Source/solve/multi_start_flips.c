#include "shared_api.h"
#include "solve_internal.h"
#include "greedy_flips_api.h"


/*----------------------------------------------------------------------------*/
// Start from the given point and apply greedy flips numFlipsRuns number of
// times. Return the best point.
/*----------------------------------------------------------------------------*/
void Multi_Start_Flips
(
    VSPINT                 **xBestHome,
    VSPINT                **AxBestHome,
    VSPINT               *wgtxBestHome,
    VSPINT                 **yBestHome,
    VSPINT                **AyBestHome,
    VSPINT               *wgtyBestHome,
    Problem                   *problem,
    Options                   *options,
    VSPINT                 numFlipRuns,
    VSPINT            numPassesInFlips
)
{
    SET_FUNCTION_COUNTER ;

    /*------------------------------------------------------------------------*/
    // Check inputs
    /*------------------------------------------------------------------------*/
    #ifdef DEBUG
        Check_Part
        (
            *xBestHome   , *AxBestHome,
            *yBestHome   , *AyBestHome,
            *wgtxBestHome, *wgtyBestHome, // Dereferencing to pass by value
            problem
        ) ;
    #endif

    /*------------------------------------------------------------------------*/
    // Compute current best objective and imbalance
    /*------------------------------------------------------------------------*/
    VSPINT objBest ;
    Compute_ObjValSimple(*xBestHome, *AxBestHome, *wgtxBestHome, *yBestHome, *AyBestHome, *wgtyBestHome, problem, &objBest) ;
    VSPINT imbBest = IMBALANCE(*wgtxBestHome, *wgtyBestHome) ;

    /*------------------------------------------------------------------------*/
    // Variables and memory
    /*------------------------------------------------------------------------*/
    VSPINT i ;
    VSPINT n = problem->graph->n ;

    VSPINT objStar, imbStar ;
    VSPINT wgtxstar, wgtystar ;
    VSPINT *xstar  = malloc(n * sizeof(VSPINT)) ;
    VSPINT *Axstar = malloc(n * sizeof(VSPINT)) ;
    VSPINT *ystar  = malloc(n * sizeof(VSPINT)) ;
    VSPINT *Aystar = malloc(n * sizeof(VSPINT)) ;

    VSPINT wgtxBest, wgtyBest ;
    VSPINT *xBest  = malloc(n * sizeof(VSPINT)) ;
    VSPINT *AxBest = malloc(n * sizeof(VSPINT)) ;
    VSPINT *yBest  = malloc(n * sizeof(VSPINT)) ;
    VSPINT *AyBest = malloc(n * sizeof(VSPINT)) ;



    // At the beginning of the iteration, (xstar, ystar) is the input point
    // will and it is improved with flips. If the objective is then better than
    // objBest, 
    bool foundBetterThanOriginal = false ;
    for(VSPINT p = 0 ; p < numFlipRuns ; p++)
    {
        // Copy the original into the star point
        Copy_Data
        (
            *xBestHome, *AxBestHome, wgtxBestHome,
            *yBestHome, *AyBestHome, wgtyBestHome,
            xstar, Axstar, &wgtxstar,
            ystar, Aystar, &wgtystar,
            n
        ) ;




// Apply flips to the star point
        State Xstar, Ystar ;
        Point xstarPoint, ystarPoint ;
        Az    xstarAz,    ystarAz ;

        Xstar.point = &xstarPoint ;
        Xstar.Az    = &xstarAz ;
        Ystar.point = &ystarPoint ;
        Ystar.Az    = &ystarAz ;

        Xstar.point->int_part = xstar ;
        Xstar.point->ind_frac = EMPTY ;
        Xstar.point->n = problem->graph->n ;
        Xstar.Az->int_part = Axstar ;
        Xstar.Az->col_frac_z = NULL ;
        Xstar.Az->ind_frac_z = EMPTY ;
        Xstar.Az->n = problem->graph->n ;
        Xstar.weight = wgtxstar ;

        Ystar.point->int_part = ystar ;
        Ystar.point->ind_frac = EMPTY ;
        Ystar.point->n = problem->graph->n ;
        Ystar.Az->int_part = Aystar ;
        Ystar.Az->col_frac_z = NULL ;
        Ystar.Az->ind_frac_z = EMPTY ;
        Ystar.Az->n = problem->graph->n ;
        Ystar.weight = wgtystar ;

        Greedy_Flips(&Xstar, &Ystar, problem, numPassesInFlips, true) ;

        wgtxstar = Xstar.weight ;
        wgtystar = Ystar.weight ;

        #ifdef DEBUG
            Check_Part
            (
                xstar, Axstar, 
                ystar, Aystar, 
                wgtxstar, wgtystar, // Dereferencing to pass by value
                problem
            ) ;
        #endif        

        Compute_ObjValSimple(xstar, Axstar, wgtxstar, ystar, Aystar, wgtystar, problem, &objStar) ;
        imbStar = IMBALANCE(wgtxstar, wgtystar) ;
        if( (objStar > objBest) || (objStar == objBest && imbStar < imbBest) )
        {
            foundBetterThanOriginal = true ;

            imbBest = IMBALANCE(wgtxstar, wgtystar) ;
            objBest = objStar ;

            SWAP(xstar, xBest, VSPINT*) ;
            SWAP(Axstar, AxBest, VSPINT*) ;
            SWAP(wgtxstar, wgtxBest, VSPINT) ;

            SWAP(ystar, yBest, VSPINT*) ;
            SWAP(Aystar, AyBest, VSPINT*) ;
            SWAP(wgtystar, wgtyBest, VSPINT) ;
        }
    }

    // Copy the best solution into the star point
    if(foundBetterThanOriginal)
    {
        SWAP(*xBestHome, xBest, VSPINT*) ;
        SWAP(*AxBestHome, AxBest, VSPINT*) ;
        SWAP(*wgtxBestHome, wgtxBest, VSPINT) ;

        SWAP(*yBestHome, yBest, VSPINT*) ;
        SWAP(*AyBestHome, AyBest, VSPINT*) ;
        SWAP(*wgtyBestHome, wgtyBest, VSPINT) ;
    }

    /*------------------------------------------------------------------------*/
    // Cleanup
    /*------------------------------------------------------------------------*/
    free(xstar);
    free(Axstar);
    free(ystar);
    free(Aystar);
    
    free(xBest);
    free(AxBest);
    free(yBest);
    free(AyBest);

    return ;
}


void Copy_Data
(
    VSPINT *xfrom,
    VSPINT *Axfrom,
    VSPINT *wgtxfrom,
    VSPINT *yfrom,
    VSPINT *Ayfrom,
    VSPINT *wgtyfrom,
    VSPINT *xto,
    VSPINT *Axto,
    VSPINT *wgtxto,
    VSPINT *yto,
    VSPINT *Ayto,
    VSPINT *wgtyto,
    VSPINT n
)
{
    VSPINT i ;
    *wgtxto = *wgtxfrom ;
    *wgtyto = *wgtyfrom ;

    for(i = 0 ; i < n ; i++)
    {
        xto[i]  = xfrom[i];
        Axto[i] = Axfrom[i];
        yto[i]  = yfrom[i];
        Ayto[i] = Ayfrom[i];
    }
}