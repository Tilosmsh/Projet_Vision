/* -------------------- */
/* --- sigmadelta.c --- */
/* -------------------- */

/*
 * Copyright (c) 2004 - 2013, Lionel Lacassagne, All rights reserved
 * University of Paris Sud, Laboratoire de Recherche en Informatique 
 * Creation: 2004-05-18 :
 * Creation: 2021-01-06 : version line pour pipeline
 */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#include "nrtype.h"
#include "nrdef.h"
#include "nrutil.h"

#include "sigmadelta.h"

/* --------------- */
int findPower2(int x)
/* --------------- */
{
    int p = 0;
    if(!x) return 0;
    
    while(!(x & 1)) {
        x = x / 2;
        p = p + 1;
    }
    return p;
}
// -----------------------------------------------------------------------------------------
void SigmaDelta_Step0_line(uint8 *I, uint8 *M, uint8 *O, uint8 *V, uint8 *E, int j0, int j1)
// -----------------------------------------------------------------------------------------
{
}
// ------------------------------------------------------------------------------------------------
void SigmaDelta_1Step_line(uint8 *I, uint8 *M, uint8 *O, uint8 *V, uint8 *E, int k, int j0, int j1)
// ------------------------------------------------------------------------------------------------
{
}
// ---------------------------------------------------------------------------------------------------------
void SigmaDelta_Step0(uint8 **I, uint8 **M, uint8 **O, uint8 **V, uint8 **E, int i0, int i1, int j0, int j1)
// ---------------------------------------------------------------------------------------------------------
{
    for(int i=i0; i<=i1; i++) {
        SigmaDelta_Step0_line(I[i], M[i], O[i], V[i], E[i], j0, j1);
    }
}
// ----------------------------------------------------------------------------------------------------------------
void SigmaDelta_1Step(uint8 **I, uint8 **M, uint8 **O, uint8 **V, uint8 **E, int k, int i0, int i1, int j0, int j1)
// ----------------------------------------------------------------------------------------------------------------
{
    for(int i=i0; i<=i1; i++) {
        SigmaDelta_1Step_line(I[i], M[i], O[i], V[i], E[i], k, j0, j1);
    }
}
