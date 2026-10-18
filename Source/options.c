#include "core.h"
#include "options.h"
#include "VSPQP.h"
#include "stdio.h"




Levels* construct_levels
(
    Levels_selection_type selection_type
)
{
    #ifdef DEBUG
        static long long unsigned count = 0 ;
        count++ ;
    #endif

    Levels *levels = malloc(sizeof(*levels)) ;
    levels->top_levels                = false ;
    levels->bot_levels                = false ;
    levels->all_levels                = false ;
    levels->no_levels                 = false ;
    levels->specify_selected_levels   = false ;
    levels->specify_unselected_levels = false ;

    levels->levels_selected       = NULL ;
    levels->levels_not_selected   = NULL ;

    if(selection_type == TOP_LEVELS)
    {
        levels->top_levels     = true ;
        levels->num_top_levels = 0 ;
    }
    else if(selection_type == BOT_LEVELS)
    {
        levels->bot_levels     = true ;
        levels->num_bot_levels = 0 ;
    }
    else if(selection_type == ALL_LEVELS)
    {
        levels->all_levels = true ;
    }
    else if(selection_type == NO_LEVELS)
    {
        levels->no_levels = true ;
    }
    else if(selection_type == SPECIFY_SELECTED_LEVELS)
    {
        levels->specify_selected_levels = true ;
        levels->levels_selected = construct_dyn_stk() ;
    }
    else if(selection_type == SPECIFY_UNSELECTED_LEVELS)
    {
        levels->specify_unselected_levels = true ;
        levels->levels_not_selected = construct_dyn_stk() ;
    }
    else
    {
        #ifdef DEBUG
            Error("Invalid selection type", DBG_INFO) ;
        #endif
    }

    return levels ;
}


void destruct_levels
(
    Levels* levels
)
{
    if(levels->specify_selected_levels)
    {
        destruct_dyn_stk(levels->levels_selected) ;
    }
    else if(levels->specify_unselected_levels)
    {
        destruct_dyn_stk(levels->levels_not_selected) ;
    }

    free(levels) ;
}


bool Levels_is_member
(
    Problem *problem,
    Levels  *levels
)
{
    #ifdef DEBUG
        static long long unsigned count = 0 ;

        /* Check that only one  value is true */
        VSPINT sum =    levels->top_levels +
                        levels->bot_levels +
                        levels->all_levels +
                        levels->no_levels  +
                        levels->specify_selected_levels +
                        levels->specify_unselected_levels ;
        if(sum != 1)
        {
            Error("Invalid selection type", DBG_INFO) ;
        }
    #endif
    VSPINT problem_level = problem->level ;

    if(levels->all_levels)
    {
        return true ;
    }
    else if(levels->no_levels)
    {
        return false ;
    }
    else if(levels->top_levels)
    {
        return (problem_level < levels->num_top_levels) ;
    }
    else if(levels->bot_levels)
    {
        /* Find the bottom level */
        VSPINT bot_level = problem->level ;
        Problem *p = problem ;
        while(p->child_problem != NULL)
        {
            p = p->child_problem ;
            bot_level = p->level ;
        }
        return (problem_level > bot_level - levels->num_bot_levels) ;
    }
    else if(levels->specify_selected_levels)
    {
        bool ret = dyn_stk_is_member
                    (
                        problem_level,
                        levels->levels_selected
                    ) ;
        return ret ;
    }
    else if(levels->specify_unselected_levels)
    {
        bool ret = !dyn_stk_is_member
                        (
                            problem_level,
                            levels->levels_not_selected
                        ) ;
        return ret ;
    }
    else
    {
        #ifdef DEBUG
            Error("Invalid selection type", DBG_INFO) ;
        #endif
    }
    return false ;
}

void Levels_select
(
    VSPINT       level_num,
    Levels  *levels
)
{
    #ifdef DEBUG
        static long long unsigned count = 0 ;
        count++ ;
        if( levels->specify_selected_levels == false)
        {
            Error("Using wrong function to select levels", DBG_INFO) ;
        }
    #endif

    #ifdef DEBUG
        if(dyn_stk_is_member(level_num, levels->levels_selected))
        {
            Error("Level already selected", DBG_INFO) ;
        }
    #endif

    dyn_stk_push(levels->levels_selected, level_num) ;
}


void Levels_unselect
(
    VSPINT      level_num,
    Levels    *levels
)
{
    #ifdef DEBUG
        static long long unsigned count = 0 ;
        count++ ;
        if( levels->specify_unselected_levels == false)
        {
            Error("Using wrong function to select levels", DBG_INFO) ;
        }
    #endif

    #ifdef DEBUG
        if( dyn_stk_is_member
            (
                level_num,
                levels->levels_not_selected
            ) )
        {
            Error("Level already unselected", DBG_INFO) ;
        }
    #endif

    dyn_stk_push(levels->levels_not_selected, level_num) ;

    return ;
}


void Levels_select_top
(
    VSPINT n_top_levels,
    Levels *levels
)
{
    #ifdef DEBUG
        static long long unsigned count = 0 ;
        count++ ;
        if(levels->top_levels == false)
        {
            Error("Using wrong function to select levels", DBG_INFO) ;
        }
    #endif

    levels->num_top_levels = n_top_levels;

    return ;
}

void Levels_select_bot
(
    VSPINT n_bot_levels,
    Levels *levels
)
{
    #ifdef DEBUG
        static long long unsigned count = 0 ;
        count++ ;
        if(levels->bot_levels == false)
        {
            Error("Using wrong function to select levels", DBG_INFO) ;
        }
    #endif

    levels->num_bot_levels = n_bot_levels;

    return ;
}




/*----------------------------------------------------------------------------*/
/* Dynamic Stack Definitions */
/*----------------------------------------------------------------------------*/

/* Construct a dynamic stack. Initial size is 32 elements */
dyn_stk *construct_dyn_stk()
{
    dyn_stk *ds  = malloc(sizeof(*ds)) ;
    ds->max_size = 32 ;
    ds->stack     = malloc( (ds->max_size) * sizeof(VSPINT) ) ;
    ds->n_stack   = 0 ;

    return ds ;
}

void destruct_dyn_stk(dyn_stk* ds)
{
    free(ds->stack) ;
    free(ds) ;
}

void dyn_stk_push(dyn_stk* ds, VSPINT to_push)
{
    if (ds->n_stack == ds->max_size)
    {
        ds->max_size *= 2 ;
        ds->stack = realloc(ds->stack, ds->max_size * sizeof(VSPINT)) ;
    }

    ds->stack[ds->n_stack] = to_push ;
    ds->n_stack++ ;
}

bool dyn_stk_is_member(VSPINT i, dyn_stk* ds)
{
    for (VSPINT j = 0 ; j < ds->n_stack ; j++)
    {
        if (ds->stack[j] == i)
        {
            return true ;
        }
    }

    return false ;
}

/* Create a deep copy */
dyn_stk *dyn_stk_copy(dyn_stk *stack_in)
{
    dyn_stk *stack_out = construct_dyn_stk() ;

    for (VSPINT i = 0 ; i < stack_in->n_stack ; i++)
    {
        dyn_stk_push(stack_out, stack_in->stack[i]) ;
    }

    return stack_out ;
}





/* Just construct and set everything to empty */
restart_decision *construct_restart_decision()
{
    restart_decision *ret = malloc(sizeof(*ret)) ;

    /* Levels */
    ret->use_levels = false ;
    ret->levels     = NULL ;

    /* Rule */
    ret->use_rule             = false ;
    ret->restarts_done        = 0 ;
    ret->restarts_allowed     = 0 ;
    ret->rule_param           = 0 ;
    ret->finest_level_allowed = 0 ;

    return ret ;
}

void destruct_restart_decision(restart_decision *dec)
{
    if (dec->levels != NULL)
    {
        destruct_levels(dec->levels) ;
    }

    free(dec) ;

    return ;
}

/* If using the rule, and if it returns true, then it also increments
   dec->restarts_done. */
bool restart_decide
(
    restart_decision       *dec,
    VSPFLOAT                  objIn, /* If not using rule, this is ignored */
    VSPFLOAT                 objOut, /* If not using rule, this is ignored */
    Problem                *problem  /* If using levels, this is ignored   */
)
{
    #ifdef DEBUG
        static long long unsigned count = 0 ;
        count++ ;

        if(dec->use_levels + dec->use_rule != 1)
        {
            Error("Error in restart decision", DBG_INFO) ;
        }

        if(dec->use_rule)
        {
            if(dec->restarts_done > dec->restarts_allowed)
            {
                Error("Error in restart decision", DBG_INFO) ;
            }

            if(objOut < objIn)
            {
                Error("Probable error. Final objective worse.", DBG_INFO) ;
            }
        }
        else
        {
            if(dec->levels == NULL)
            {
                Error("Error in restart decision", DBG_INFO) ;
            }
        }
    #endif



    if (dec->use_levels)
    {
        return Levels_is_member(problem, dec->levels) ;
    }
    else if (dec->use_rule)
    {
        if(dec->restarts_done >= dec->restarts_allowed
            || problem->level < dec->finest_level_allowed)
        {
            return false ;
        }

        /* Check if we should restart */
        VSPFLOAT rel_change = (objOut - objIn) / objIn ;
        if (rel_change <= dec->rule_param)
        {
            dec->restarts_done++ ;
            return true ;
        }
        else
        {
            return false ;
        }
    }

    return false ;
}

void restart_decide_reset(restart_decision *dec)
{
    dec->restarts_done = 0 ;
    return ;
}



/***************** Options ************************************************/
Options* construct_options()
{
    Options* options = malloc(sizeof(Options));
    options->restart_decision     = NULL ;
    return options;
}


void destruct_options(Options* options)
{
    if(options->restart_decision != NULL)
    {
        destruct_restart_decision(options->restart_decision) ;
    }


    /* Free the struct itself */
    free(options);
}


void set_default_options(Options *options)
{
    /*------------------------------------------------------------------------*/
    /* Restart Level Options                                                  */
    /*------------------------------------------------------------------------*/
    /* All levels in [2, infty) */
    options->restart_decision             = construct_restart_decision() ;
    options->restart_decision->use_levels = true ;
    options->restart_decision->levels     = construct_levels(SPECIFY_UNSELECTED_LEVELS) ;
    Levels_unselect(0, options->restart_decision->levels) ;
    Levels_unselect(1, options->restart_decision->levels) ;

    /*------------------------------------------------------------------------*/
    /* Gamma Reduction Options                                                */
    /*------------------------------------------------------------------------*/
    options->max_num_vert_small_graph = 8000 ;
    options->num_gr_for_small_graphs  = 4 ;
    options->num_gr_for_large_graphs  = 1 ;

    /*------------------------------------------------------------------------*/
    /* Relax Upper Bound Options                                              */
    /*------------------------------------------------------------------------*/
    options->tighten_u        = true ;
    options->u_ratio_finest   = 0.6 ;
    options->u_ratio_coarsest = 0.5 ;
    options->ub_based_on_bisection = false ;

    /*------------------------------------------------------------------------*/
    /* Starting Guess Options                                                 */
    /*------------------------------------------------------------------------*/
    options->initial_guess                  = centroid_initial_guess ;
    options->centroidVal                    = 0.35 ;
    options->random_point_weight_proportion = 0.38 ;

    /*------------------------------------------------------------------------*/
    /* FM Options                                                             */
    /*------------------------------------------------------------------------*/
    options->FM_Limit = 300 ;

    /*------------------------------------------------------------------------*/
    /* Mountain Climbing Options                                              */
    /*------------------------------------------------------------------------*/
    options->tol_LP_opt_checker = 1E-10 ;
    options->start_mc_both_vars = TRUE ;

    /*------------------------------------------------------------------------*/
    /* Simple Coarsening Options                                              */
    /*------------------------------------------------------------------------*/
    options->compression_threshold = 0.51 ;
    options->node_cutoff           = 32 ;
    options->density_cutoff        = 0.5 ;

    options->initial_project       = FALSE ;

    /* ====================================================================== */
    // Vertex Cover Options
    /* ====================================================================== */
    options->cover_type            = MIN_CARD_COVER ;
    options->try_both_covers       = TRUE ;

    /* ====================================================================== */
    // Bisection ngp options
    /* ====================================================================== */
    options->bisection_ngp = 3 ;

    /* ====================================================================== */
    // FM subspace options
    /* ====================================================================== */
    options->use_FM_subspace_grad_proj = FALSE ;
    options->noFM = 0.2 ;


    options->attempt_second_start = TRUE ;
    options->min_imbal_ratio_second_start = 0.07 ;
    options->min_level_second_start = 4 ;

    /* ====================================================================== */
    // Gradient projection applied to vsp
    /* ====================================================================== */
    options->vsp_grad_proj = false ;
    options->vsp_ngp = 3 ;

    /* ====================================================================== */
    // Hub removal
    /* ====================================================================== */
    options-> hub_treshold = 0.5 ;

    return ;
}

static const char *initial_guess_to_str(initial_guess g)
{
    switch (g)
    {
        case centroid_initial_guess:  return "centroid";
        case random_initial_guess:    return "random";
        case manual_initial_guess:    return "manual";
        case bisection_initial_guess: return "bisection";
        default:                      return "unknown";
    }
}

static const char *cover_type_to_str(cover_type c)
{
    switch (c)
    {
        case MIN_COST_COVER:          return "min-cost-cover" ;
        case MIN_CARD_COVER:          return "min-cardinality-cover" ;
        case APPROX_MIN_CARD_COVER:   return "approx-min-cardinality-cover" ;
        default:                      return "unknown";
    }
}

void Print_Options(const Options *opts, const char *path)
{
    if (path == NULL)
    {
        printf("Print_Options: path is NULL\n");
        return;
    }
    FILE *f = fopen(path, "w");
    if (f == NULL)
    {
        printf("Print_Options: could not open %s for writing\n", path);
        return;
    }

    if (opts == NULL)
    {
        fprintf(f, "Options: (null)\n");
        fclose(f);
        return;
    }

    fprintf(f, "Initial guess\n") ;
    fprintf(f, "  initial_guess                  = %s\n",   initial_guess_to_str(opts->initial_guess));
    fprintf(f, "  random_point_weight_proportion = %g\n",   opts->random_point_weight_proportion);
    fprintf(f, "  centroidVal                    = %g\n",   opts->centroidVal);

    /* Gamma reduction */
    fprintf(f, "Gamma reduction\n") ;
    fprintf(f, "  num_gr_for_small_graphs        = %lld\n", PRI(opts->num_gr_for_small_graphs));
    fprintf(f, "  num_gr_for_large_graphs        = %lld\n", PRI(opts->num_gr_for_large_graphs));
    fprintf(f, "  max_num_vert_small_graph       = %lld\n", PRI(opts->max_num_vert_small_graph));

    /* Miscellaneous */
    fprintf(f, "Miscellaneous\n") ;
    fprintf(f, "  tighten_u                      = %s\n",   opts->tighten_u ? "true" : "false");
    fprintf(f, "  u_ratio_finest                 = %f\n",   opts->u_ratio_finest);
    fprintf(f, "  u_ratio_coarsest               = %f\n",   opts->u_ratio_coarsest);

    /* FM Options */
    fprintf(f, "  FM_Limit                       = %lld\n", PRI(opts->FM_Limit));

    /* Mountain Climbing Options */
    fprintf(f, "Mountain climbing\n") ;
    fprintf(f, "  tol_LP_opt_checker             = %g\n",   opts->tol_LP_opt_checker);
    fprintf(f, "  Initial mc var strategy        = %s\n",   opts->start_mc_both_vars ? "try_both_vars" : "use_lower_cost_var");

    /* Simple Coarsening Options */
    fprintf(f, "Coarsening\n") ;
    fprintf(f, "  compression_threshold          = %g\n",   opts->compression_threshold);
    fprintf(f, "  node_cutoff                    = %i\n",   opts->node_cutoff);
    fprintf(f, "  density_cutoff                 = %f\n",   opts->density_cutoff);

    /* Cover type options */
    fprintf(f, "  cover_type                     = %s\n",   cover_type_to_str(opts->cover_type) ) ;
    fprintf(f, "  try_both_covers                = %s\n",   (opts->try_both_covers) ? "yes" : "no") ;

    /* vs options */
    fprintf(f, "VS startup\n") ;
    fprintf(f, "  initial_projection             = %s\n",   (opts->initial_project) ? "yes" : "no") ;
    fprintf(f, "  bisection_ngp                  = %i\n",   opts->bisection_ngp) ;

    /* FM subspace gradient projection */
    fprintf(f, "FM Subspace\n") ;
    fprintf(f, "  Use FM subpsace grad proj      = %s\n",   (opts->use_FM_subspace_grad_proj) ? "yes" : "no") ;
    fprintf(f, "  Max subspace size              = %f\n",   opts->noFM) ;


    /* Second start options */
    fprintf(f, "Second start\n") ;
    fprintf(f, "  attempt_second_start           = %s\n",   opts->attempt_second_start ? "true" : "false");
    fprintf(f, "  min_imbal_ratio_second_start   = %g\n",   opts->min_imbal_ratio_second_start);
    fprintf(f, "  min_level_second_start         = %lld\n", PRI(opts->min_level_second_start));

    /* VSP gradient projection options */
    fprintf(f, "VSP gradient projection\n") ;
    fprintf(f, "  vsp_grad_proj                  = %s\n",   (opts->vsp_grad_proj) ? "yes" : "no") ;
    fprintf(f, "  vsp_ngp                        = %i\n",   opts->vsp_ngp ) ;

    /* hub removal */
    fprintf(f, "Hub removal\n") ;
    fprintf(f, "  hub_treshold                   = %f\n",   opts->hub_treshold ) ;

    fclose(f);
}