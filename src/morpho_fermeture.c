/* -------------------------- */
/* --- morpho_fermeture.c --- */
/* -------------------------- */

/*
 * Copyright (c) 2004 - 2013, Lionel Lacassagne, All rights reserved
 * University of Paris Sud, Laboratoire de Recherche en Informatique 
 
 */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#include "macro.h"
#include "nrtype.h"
#include "nrdef.h"
#include "nrutil.h"
//#include "sequence.h"

#include "swp.h"

#include "morpho.h"
#include "morpho_erosion.h"
#include "morpho_dilatation.h"

#include "morpho_fermeture.h"

// -------------------------------------------------------------------------------
void line_fermeture3_ui8matrix_fusion_basic(uint8 **X, int i, int j0, int j1, uint8 **Y)
// -------------------------------------------------------------------------------------
{
    uint8 x[5][5];
    uint8 z[3][3];
    uint8 max;

    for(int j = j0 ; j <= j1 ; j++ ){
        x[0][0]=load2(X,i-2,j-2);
        x[1][0]=load2(X,i-1,j-2);
        x[2][0]=load2(X,i  ,j-2);
        x[3][0]=load2(X,i+1,j-2);
        x[4][0]=load2(X,i+2,j-2);

        x[0][1]=load2(X,i-2,j-1);
        x[1][1]=load2(X,i-1,j-1);
        x[2][1]=load2(X,i  ,j-1);
        x[3][1]=load2(X,i+1,j-1);
        x[4][1]=load2(X,i+2,j-1);

        x[0][2]=load2(X,i-2,j  );
        x[1][2]=load2(X,i-1,j  );
        x[2][2]=load2(X,i  ,j  );
        x[3][2]=load2(X,i+1,j  );
        x[4][2]=load2(X,i+2,j  );

        x[0][3]=load2(X,i-2,j+1);
        x[1][3]=load2(X,i-1,j+1);
        x[2][3]=load2(X,i  ,j+1);
        x[3][3]=load2(X,i+1,j+1);
        x[4][3]=load2(X,i+2,j+1);

        x[0][4]=load2(X,i-2,j+2);
        x[1][4]=load2(X,i-1,j+2);
        x[2][4]=load2(X,i  ,j+2);
        x[3][4]=load2(X,i+1,j+2);
        x[4][4]=load2(X,i+2,j+2);

        z[0][0]=or9_mat33(x,0,0);
        z[0][1]=or9_mat33(x,0,1);
        z[0][2]=or9_mat33(x,0,2);

        z[1][0]=or9_mat33(x,1,0);
        z[1][1]=or9_mat33(x,1,1);
        z[1][2]=or9_mat33(x,1,2);

        z[2][0]=or9_mat33(x,2,0);
        z[2][1]=or9_mat33(x,2,1);
        z[2][2]=or9_mat33(x,2,2);

        max=and9_mat33(z,0,0);
        store2(Y,i,j,max);
    }
}
// -----------------------------------------------------------------------------------
void line_fermeture3_ui8matrix_fusion_ilu5_red(uint8 **X, int i, int j0, int j1, uint8 **Y)
// ----------------------------------------------------------------------------------------
{
}
// ---------------------------------------------------------------------------------------------
void line_fermeture3_ui8matrix_fusion_ilu5_elu2_red(uint8 **X, int i, int j0, int j1, uint8 **Y)
// ---------------------------------------------------------------------------------------------
{
}
// ----------------------------------------------------------------------------------------------------
void line_fermeture3_ui8matrix_fusion_ilu5_elu2_red_factor(uint8 **X, int i, int j0, int j1, uint8 **Y)
// ----------------------------------------------------------------------------------------------------
{
}
// -----------------------------------------------------------------------------------------
void line_fermeture3_ui8matrix_fusion_ilu15_red(uint8 **X, int i, int j0, int j1, uint8 **Y)
// -----------------------------------------------------------------------------------------
{
}
// ---------------------------------------------------------------------------------------------
void fermeture3_ui8matrix_basic(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y, uint8 **Z)
// ---------------------------------------------------------------------------------------------
{
}
// -----------------------------------------------------------------------------------
void fermeture3_ui8matrix_fusion(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// -----------------------------------------------------------------------------------
{
}
// --------------------------------------------------------------------------------------------
void fermeture3_ui8matrix_fusion_ilu5_red(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// --------------------------------------------------------------------------------------------
{
}
// -------------------------------------------------------------------------------------------------
void fermeture3_ui8matrix_fusion_ilu5_elu2_red(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// -------------------------------------------------------------------------------------------------
{
}
// --------------------------------------------------------------------------------------------------------
void fermeture3_ui8matrix_fusion_ilu5_elu2_red_factor(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// --------------------------------------------------------------------------------------------------------
{
}
// --------------------------------------------------------------------------------------------
void fermeture3_ui8matrix_fusion_ilu15_red(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// --------------------------------------------------------------------------------------------
{
}
// ------------------------------------------------------------------------------------------------------
void fermeture3_ui8matrix_pipeline_basic(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y, uint8 **Z)
// ------------------------------------------------------------------------------------------------------
{
}
// ---------------------------------------------------------------------------------------------------------
void fermeture3_ui8matrix_pipeline_ilu3_red(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y, uint8 **Z)
// ---------------------------------------------------------------------------------------------------------
{
}
// ----------------------------------------------------------------------------------------------------------------
void fermeture3_ui8matrix_pipeline_elu2_red_factor(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y, uint8 **Z)
// ----------------------------------------------------------------------------------------------------------------
{
}
