#include "shared_api.h"
#include "flips_public_v4.h"
#include "flips_public_v5.h"
#include "flips_public_v7.h"
#include "flips_api.h"


/*----------------------------------------------------------------------------*/
// Start from the given point and apply greedy flips numFlipsRuns number of
// times. Return the best point.
/*----------------------------------------------------------------------------*/
void TryFlipsManyStartingPoints
(
    VSPINT                 **xBestHome,
    VSPINT                **AxBestHome,
    VSPINT               *wgtxBestHome,
    VSPINT                 **yBestHome,
    VSPINT                **AyBestHome,
    VSPINT               *wgtyBestHome,
    VSP_problem               *problem,
    VSP_options               *options,
    VSPINT                 numFlipRuns,
    VSPINT            numPassesInFlips
)
{
    SET_FUNCTION_COUNTER ;

    /*------------------------------------------------------------------------*/
    // Check inputs
    /*------------------------------------------------------------------------*/
    #ifdef VSP_DEBUG
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
    VSPINT n = problem->Graph->n ;

    VSPINT objStar, imbStar ;
    VSPINT wgtxstar, wgtystar ;
    VSPINT *xstar  = malloc(n * sizeof(VSPINT)) ;
    VSPINT *Axstar = malloc(n * sizeof(VSPINT)) ;
    VSPINT *ystar  = malloc(n * sizeof(VSPINT)) ;
    VSPINT *Aystar = malloc(n * sizeof(VSPINT)) ;

    VSPINT wgtxBest, wgtyBest ;
    VSPINT *xBest  = malloc(n * sizeof(VSPINT)) ;
    VSPINT *AxBest  = malloc(n * sizeof(VSPINT)) ;
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
        // Flips() 
        Flips_v5
        (
            xstar, Axstar, &wgtxstar,
            ystar, Aystar, &wgtystar,
            problem, numPassesInFlips
        ) ;



        #ifdef VSP_DEBUG
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