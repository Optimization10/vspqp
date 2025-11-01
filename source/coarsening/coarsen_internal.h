#ifndef COARSEN_INTERNAL_H
#define COARSEN_INTERNAL_H

#include "shared_api.h"
#include "coarsen_api.h" // IWYU pragma: keep

VSP_problem* VSP_coarsen_all_the_way_down
(
    VSP_problem* problem,
    VSP_options* options
) ;

VSPINT* VSP_get_all_matches
(
    VSPINT    vertex,
    VSPINT *matching
) ;

void VSP_match(VSP_problem* problem, VSP_options* options) ;

void VSP_matching_HEMSR
(
    VSP_problem* problem, 
    VSP_options* options
);

VSPINT VSP_get_match
(
    VSPINT vertex,
    VSPINT* matching
) ;

VSPINT VSP_get_edge_weight
(
    VSPINT      vertex_a,
    VSPINT      vertex_b,
    VSPINT           *Ai,
    VSPINT           *Ap,
    VSPINT *Ae
) ;


VSPINT* VSP_get_neighbors(VSPINT v, VSPINT* Ap, VSPINT* Ai) ;

VSPINT* VSP_get_preimage_of_vertex
(
    VSPINT vertex_coarse,         /* A vertex in the coarse problem */
    VSPINT     *matching,
    VSPINT  *invmatchmap
) ;

VSPINT VSP_there_is_an_edge
(
    VSPINT vertex_a,
    VSPINT vertex_b,
    VSPINT      *Ai,
    VSPINT      *Ap
) ;

void VSP_check_matrix_is_symmetric(
    VSPINT           *Ai,
    VSPINT           *Ap,
    VSPINT           *Ae,
    VSPINT            nz,
    VSPINT             n,
    VSP_options* options
) ;

VSPINT VSP_check_for_edge_between_arrays_of_vertices
(
    VSPINT *arr1,
    VSPINT *arr2,
    VSPINT  size,
    VSPINT   *Ai,
    VSPINT   *Ap
) ;

VSPINT VSP_check_matrix
(
    VSPINT* Ai,
    VSPINT* Ap,
    VSPINT* Ae,
    VSPLONG nz,
    VSPINT n,
    VSP_options* options
) ;

void Heavy_Edge_Matching(VSP_problem* problem, VSP_options* options) ;

bool Is_Matched(VSPINT vertex, VSPINT* matching) ;

void VSP_create_match
(
    VSP_problem* problem, 
    VSPINT vertexA, VSPINT vertexB, 
    VSP_match_type matchtype
) ;


bool IsMatchedMatchtype
(
    VSPINT vertex,
    VSPINT* matchingArray,
    VSP_match_type *matchtypeArray,
    VSP_match_type desiredType
) ;

void Match_Check
(
    VSPINT                     n,
    VSPINT                    cn,
    VSPINT             *matching,
    VSPINT             *matchmap,
    VSPINT          *invmatchmap,
    VSP_match_type    *matchtype
) ;

#endif