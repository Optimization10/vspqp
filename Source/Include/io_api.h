#ifndef IO_API_H
#define IO_API_H

#include "shared_api.h"


void VSP_reset_output_files() ;
void VSP_reset_file(char* name) ;
VSP_problem* VSP_read_problem(char* file_name, VSP_options* options) ;
void VSP_print_stats(VSP_stats *stats) ;

/* Dump graph to a file in coordinate format */
void VSP_dump_graph
(
    VSP_state *X,
    VSP_state *Y,
    VSP_problem *problem
) ;



void VSP_write_to_csv
(
    char                *graph_name,
    char                *col_title_given, /* If empty string, a title is
                                             generated from the options       */
    VSP_solution_stats  *stats,
    VSP_options         *options
) ;

void VSP_write_to_csv_2
(
    char                *graph_name,
    char                *col_title_given, /* If empty string, a title is
                                             generated from the options       */
    VSP_solution_stats  *solStats,
    VSP_problem          *problem,
    VSP_options          *options
) ;

void VSP_print_solution_stats
(
    VSPINT            *x,
    VSPINT            *y,
    VSP_problem *problem,
    char           *name,  /* name of the problem */
    VSPFLOAT        time
) ;

VSP_solution_stats* VSP_compute_solution_stats(
    VSPINT            *x,
    VSPINT            *y,
    VSP_problem *problem,
    VSPFLOAT        time
) ;

#endif