#ifndef COARSEN_INTERNAL_H
#define COARSEN_INTERNAL_H

#include "shared_api.h"
#include "coarsen_api.h"


VSPINT* get_all_matches
(
    VSPINT    vertex,
    VSPINT *matching
) ;

VSPINT get_match
(
    VSPINT vertex,
    VSPINT* matching
) ;

VSPINT get_edge_weight
(
    VSPINT    vertex_a,
    VSPINT    vertex_b,
    VSPINT         *Ai,
    VSPLONG        *Ap,
    VSPINT         *Ae
) ;


VSPINT* get_neighbors(VSPINT v, VSPLONG* Ap, VSPINT* Ai) ;

VSPINT* get_preimage_of_vertex
(
    VSPINT vertex_coarse,         /* A vertex in the coarse problem */
    VSPINT     *matching,
    VSPINT  *invmatchmap
) ;

VSPINT there_is_an_edge
(
    VSPINT vertex_a,
    VSPINT vertex_b,
    VSPINT      *Ai,
    VSPLONG     *Ap
) ;

void check_matrix_is_symmetric
(
    VSPINT           *Ai,
    VSPLONG          *Ap,
    VSPINT           *Ae,
    VSPLONG           nz,
    VSPINT             n,
    Options* options
) ;

VSPINT check_for_edge_between_arrays_of_vertices
(
    VSPINT *arr1,
    VSPINT *arr2,
    VSPINT  size,
    VSPINT   *Ai,
    VSPLONG  *Ap
) ;

VSPINT check_matrix
(
    VSPINT*           Ai,
    VSPLONG*          Ap,
    VSPINT*           Ae,
    VSPLONG           nz,
    VSPINT             n,
    Options* options
) ;

bool Is_Matched(VSPINT vertex, VSPINT* matching) ;

bool IsMatchedMatchtype
(
    VSPINT vertex,
    VSPINT* matchingArray,
    match_type *matchtypeArray,
    match_type desiredType
) ;

void Match_Check
(
    VSPINT                     n,
    VSPINT                    cn,
    VSPINT             *matching,
    VSPINT             *matchmap,
    VSPINT          *invmatchmap,
    match_type    *matchtype
) ;

#endif