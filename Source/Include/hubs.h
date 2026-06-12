#ifndef HUBS_H
#define HUBS_H
#include "core.h"

void Get_Hubs
(
    VSPINT      n,
    VSPLONG    *Ap,
    VSPFLOAT    t,
    VSPINT    **hubs,    // Return
    VSPINT     *num_hubs // Return
) ;

void Remove_Hubs
(
    VSPINT    n_input,
    VSPINT   *Ai_input,
    VSPLONG  *Ap_input,
    VSPINT   *hubs,
    VSPINT    num_hubs,
    VSPINT   *n,    // output
    VSPINT  **Ai,   // output
    VSPLONG **Ap    // output
) ;

void Restore_Hubs
(
    State    *state,     // modified in place (its point only)
    VSPINT   *hubs,
    VSPINT    num_hubs
) ;

#endif