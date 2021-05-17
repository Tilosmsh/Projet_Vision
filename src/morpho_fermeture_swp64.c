/* -------------------------------- */
/* --- morpho_fermeture_swp64.c --- */
/* -------------------------------- */

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
//#include "nrutil_ext.h" // printfM8

//#include "sequence.h"

#include "swp.h"  // left right
#include "morpho.h"

#include "morpho_dilatation_swp64.h"
#include "morpho_erosion_swp64.h"
#include "morpho_fermeture_swp64.h"
// -----------------------------------------------------------------------------------------------
void line_fermeture3_ui64matrix_swp64_fusion_basic(uint64 **X, int i, int j0, int j1, uint64 **Y)
// -----------------------------------------------------------------------------------------------
{
    uint64 x[5][5];
    uint64 z[3][3];
    uint64 max;

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

        z[0][0]=i64dilatation_mat33(x,1,1);
        z[0][1]=i64dilatation_mat33(x,1,2);
        z[0][2]=i64dilatation_mat33(x,1,3);

        z[1][0]=i64dilatation_mat33(x,2,1);
        z[1][1]=i64dilatation_mat33(x,2,2);
        z[1][2]=i64dilatation_mat33(x,2,3);

        z[2][0]=i64dilatation_mat33(x,3,1);
        z[2][1]=i64dilatation_mat33(x,3,2);
        z[2][2]=i64dilatation_mat33(x,3,3);

        max=i64erosion_mat33(z,1,1);
        store2(Y,i,j,max);
    }
}
// --------------------------------------------------------------------------------------------
void line_fermeture3_ui64matrix_swp64_fusion_red(uint64 **X, int i, int j0, int j1, uint64 **Y)
// --------------------------------------------------------------------------------------------
{
    uint64  x00,x01,x02,x03,x04,
            x10,x11,x12,x13,x14,
            x20,x21,x22,x23,x24,
            x30,x31,x32,x33,x34,
            x40,x41,x42,x43,x44;
    uint64  z00,z01,z02,
            z10,z11,z12,
            z20,z21,z22;
    uint64 max;

    for(int j = j0 ; j <= j1 ; j++ ){
        x00=load2(X,i-2,j-2);
        x10=load2(X,i-1,j-2);
        x20=load2(X,i  ,j-2);
        x30=load2(X,i+1,j-2);
        x40=load2(X,i+2,j-2);

        x01=load2(X,i-2,j-1);
        x11=load2(X,i-1,j-1);
        x21=load2(X,i  ,j-1);
        x31=load2(X,i+1,j-1);
        x41=load2(X,i+2,j-1);

        x02=load2(X,i-2,j  );
        x12=load2(X,i-1,j  );
        x22=load2(X,i  ,j  );
        x32=load2(X,i+1,j  );
        x42=load2(X,i+2,j  );

        x03=load2(X,i-2,j+1);
        x13=load2(X,i-1,j+1);
        x23=load2(X,i  ,j+1);
        x33=load2(X,i+1,j+1);
        x43=load2(X,i+2,j+1);

        x04=load2(X,i-2,j+2);
        x14=load2(X,i-1,j+2);
        x24=load2(X,i  ,j+2);
        x34=load2(X,i+1,j+2);
        x44=load2(X,i+2,j+2);

        z00 = i64dilatation9_red(  x00,x01,x02,
                            x10,x11,x12,
                            x20,x21,x22);

        z10 = i64dilatation9_red(  x10,x11,x12,
                            x20,x21,x22,
                            x30,x31,x32);

        z20 = i64dilatation9_red(  x20,x21,x22,
                            x30,x31,x32,
                            x40,x41,x42);

        z01 = i64dilatation9_red(  x01,x02,x03,
                            x11,x12,x13,
                            x21,x22,x23);

        z11 = i64dilatation9_red(  x11,x12,x13,
                            x21,x22,x23,
                            x31,x32,x33);

        z21 = i64dilatation9_red(  x21,x22,x23,
                            x31,x32,x33,
                            x41,x42,x43);

        z02 = i64dilatation9_red(  x02,x03,x04,
                            x12,x13,x14,
                            x22,x23,x24);

        z12 = i64dilatation9_red(  x12,x13,x14,
                            x22,x23,x24,
                            x32,x33,x34);

        z22 = i64dilatation9_red(  x22,x23,x24,
                            x32,x33,x34,
                            x42,x43,x44);

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);

        store2(Y,i,j,max);
    }
}
// -------------------------------------------------------------------------------------------------
void line_fermeture3_ui64matrix_swp64_fusion_ilu3_red(uint64 **X, int i, int j0, int j1, uint64 **Y)
// -------------------------------------------------------------------------------------------------
//pas encore fait ilu3 mais ilu5
{
    uint64  n00,n01,n02,n03,n04,
            n10,n11,n12,n13,n14,
            n20,n21,n22,n23,n24,
            n30,n31,n32,n33,n34,
            n40,n41,n42,n43,n44;

    uint64  z00,z01,z02,
            z10,z11,z12,
            z20,z21,z22;
    uint64 max;

    int j, r = (j1 - j0 + 1) % 5;

    //prologue
        n00 = load2(X,i-2,j0-2);
        n10 = load2(X,i-1,j0-2);
        n20 = load2(X,i  ,j0-2);
        n30 = load2(X,i+1,j0-2);
        n40 = load2(X,i+2,j0-2);
        
        n01 = load2(X,i-2,j0-1);
        n11 = load2(X,i-1,j0-1);
        n21 = load2(X,i  ,j0-1);
        n31 = load2(X,i+1,j0-1);
        n41 = load2(X,i+2,j0-1);
        
        n02 = load2(X,i-2,j0  );
        n12 = load2(X,i-1,j0  );
        n22 = load2(X,i  ,j0  );
        n32 = load2(X,i+1,j0  );
        n42 = load2(X,i+2,j0  );

        n03 = load2(X,i-2,j0+1);
        n13 = load2(X,i-1,j0+1);
        n23 = load2(X,i  ,j0+1);
        n33 = load2(X,i+1,j0+1);
        n43 = load2(X,i+2,j0+1);
        
        
    //boucle
    for(j = j0 ; j < j1 - r ; j = j + 5 ){
        n04 = load2(X,i-2,j+2);
        n14 = load2(X,i-1,j+2);
        n24 = load2(X,i  ,j+2);
        n34 = load2(X,i+1,j+2);
        n44 = load2(X,i+2,j+2);
//-----------------0---------------------

        z00 = i64dilatation9_red( n00,n01,n02,
                    n10,n11,n12,
                    n20,n21,n22);

        z10 = i64dilatation9_red( n10,n11,n12,
                    n20,n21,n22,
                    n30,n31,n32);

        z20 = i64dilatation9_red( n20,n21,n22,
                    n30,n31,n32,
                    n40,n41,n42);
//--------------------------------------

        z01 = i64dilatation9_red(  n01,n02,n03,
                            n11,n12,n13,
                            n21,n22,n23);

        z11 = i64dilatation9_red(  n11,n12,n13,
                            n21,n22,n23,
                            n31,n32,n33);

        z21 = i64dilatation9_red(  n21,n22,n23,
                            n31,n32,n33,
                            n41,n42,n43);
//--------------------------------------

        z02 = i64dilatation9_red(  n02,n03,n04,
                            n12,n13,n14,
                            n22,n23,n24);

        z12 = i64dilatation9_red(  n12,n13,n14,
                            n22,n23,n24,
                            n32,n33,n34);

        z22 = i64dilatation9_red(  n22,n23,n24,
                            n32,n33,n34,
                            n42,n43,n44);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);

        store2(Y,i,j,max);
//-----------------------1--------------------
        n00 = load2(X,i-2,j+3);
        n10 = load2(X,i-1,j+3);
        n20 = load2(X,i  ,j+3);
        n30 = load2(X,i+1,j+3);
        n40 = load2(X,i+2,j+3);
        
//--------------------------------------

        z00 = i64dilatation9_red( n01,n02,n03,
                    n11,n12,n13,
                    n21,n22,n23);

        z10 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z20 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);
//--------------------------------------

        z01 = i64dilatation9_red( n02,n03,n04,
                    n12,n13,n14,
                    n22,n23,n24);

        z11 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z21 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);
//--------------------------------------

        z02 = i64dilatation9_red( n03,n04,n00,
                    n13,n14,n10,
                    n23,n24,n20);

        z12 = i64dilatation9_red( n13,n14,n10,
                    n23,n24,n20,
                    n33,n34,n30);

        z22 = i64dilatation9_red( n23,n24,n20,
                    n33,n34,n30,
                    n43,n44,n40);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);


        store2(Y,i,j+1,max);

        n01 = load2(X,i-2,j+4);
        n11 = load2(X,i-1,j+4);
        n21 = load2(X,i  ,j+4);
        n31 = load2(X,i+1,j+4);
        n41 = load2(X,i+2,j+4);

        
//---------------2-----------------------

        z00 = i64dilatation9_red( n02,n03,n04,
                    n12,n13,n14,
                    n22,n23,n24);

        z10 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z20 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);
//--------------------------------------

        z01 = i64dilatation9_red( n03,n04,n00,
                    n13,n14,n10,
                    n23,n24,n20);

        z11 = i64dilatation9_red( n13,n14,n10,
                    n23,n24,n20,
                    n33,n34,n30);

        z21 = i64dilatation9_red( n23,n24,n20,
                    n33,n34,n30,
                    n43,n44,n40);
//--------------------------------------

        z02 = i64dilatation9_red( n04,n00,n01,
                    n14,n10,n11,
                    n24,n20,n21);

        z12 = i64dilatation9_red( n14,n10,n11,
                    n24,n20,n21,
                    n34,n30,n31);

        z22 = i64dilatation9_red( n24,n20,n21,
                    n34,n30,n31,
                    n44,n40,n41);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);
        
        

        store2(Y,i,j+2,max);

        n02 = load2(X,i-2,j+5);
        n12 = load2(X,i-1,j+5);
        n22 = load2(X,i  ,j+5);
        n32 = load2(X,i+1,j+5);
        n42 = load2(X,i+2,j+5);
//------------------3--------------------

        z00 = i64dilatation9_red( n03,n04,n00,
                    n13,n14,n10,
                    n23,n24,n20);

        z10 = i64dilatation9_red( n13,n14,n10,
                    n23,n24,n20,
                    n33,n34,n30);

        z20 = i64dilatation9_red( n23,n24,n20,
                    n33,n34,n30,
                    n43,n44,n40);
//--------------------------------------

        z01 = i64dilatation9_red( n04,n00,n01,
                    n14,n10,n11,
                    n24,n20,n21);

        z11 = i64dilatation9_red( n14,n10,n11,
                    n24,n20,n21,
                    n34,n30,n31);

        z21 = i64dilatation9_red( n24,n20,n21,
                    n34,n30,n31,
                    n44,n40,n41);
//--------------------------------------

        z02 = i64dilatation9_red( n00,n01,n02,
                    n10,n11,n12,
                    n20,n21,n22);

        z12 = i64dilatation9_red( n10,n11,n12,
                    n20,n21,n22,
                    n30,n31,n32);

        z22 = i64dilatation9_red( n20,n21,n22,
                    n30,n31,n32,
                    n40,n41,n42);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);
       

        store2(Y,i,j+3,max);

        n03 = load2(X,i-2,j+6);
        n13 = load2(X,i-1,j+6);
        n23 = load2(X,i  ,j+6);
        n33 = load2(X,i+1,j+6);
        n43 = load2(X,i+2,j+6);
//-----------------4---------------------

        z00 = i64dilatation9_red( n04,n00,n01,
                    n14,n10,n11,
                    n24,n20,n21);

        z10 = i64dilatation9_red( n14,n10,n11,
                    n24,n20,n21,
                    n34,n30,n31);

        z20 = i64dilatation9_red( n24,n20,n21,
                    n34,n30,n31,
                    n44,n40,n41);
//--------------------------------------

        z01 = i64dilatation9_red( n00,n01,n02,
                    n10,n11,n12,
                    n20,n21,n22);

        z11 = i64dilatation9_red( n10,n11,n12,
                    n20,n21,n22,
                    n30,n31,n32);

        z21 = i64dilatation9_red( n20,n21,n22,
                    n30,n31,n32,
                    n40,n41,n42);
//--------------------------------------

        z02 = i64dilatation9_red( n01,n02,n03,
                    n11,n12,n13,
                    n21,n22,n23);

        z12 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z22 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);

        store2(Y,i,j+4,max);
    }

    //epilogue
    switch (r)
    {
        case 1:
        n04 = load2(X,i-2,j+2);
        n14 = load2(X,i-1,j+2);
        n24 = load2(X,i  ,j+2);
        n34 = load2(X,i+1,j+2);
        n44 = load2(X,i+2,j+2);
//-----------------0---------------------

        z00 = i64dilatation9_red( n00,n01,n02,
                    n10,n11,n12,
                    n20,n21,n22);

        z10 = i64dilatation9_red( n10,n11,n12,
                    n20,n21,n22,
                    n30,n31,n32);

        z20 = i64dilatation9_red( n20,n21,n22,
                    n30,n31,n32,
                    n40,n41,n42);
//--------------------------------------

        z01 = i64dilatation9_red( n01,n02,n03,
                    n11,n12,n13,
                    n21,n22,n23);

        z11 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z21 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);
//--------------------------------------

        z02 = i64dilatation9_red( n02,n03,n04,
                    n12,n13,n14,
                    n22,n23,n24);

        z12 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z22 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);
        

        store2(Y,i,j,max);

        break;

        case 2:
        n04 = load2(X,i-2,j+2);
        n14 = load2(X,i-1,j+2);
        n24 = load2(X,i  ,j+2);
        n34 = load2(X,i+1,j+2);
        n44 = load2(X,i+2,j+2);
//-----------------0---------------------

        z00 = i64dilatation9_red( n00,n01,n02,
                    n10,n11,n12,
                    n20,n21,n22);

        z10 = i64dilatation9_red( n10,n11,n12,
                    n20,n21,n22,
                    n30,n31,n32);

        z20 = i64dilatation9_red( n20,n21,n22,
                    n30,n31,n32,
                    n40,n41,n42);
//--------------------------------------

        z01 = i64dilatation9_red( n01,n02,n03,
                    n11,n12,n13,
                    n21,n22,n23);

        z11 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z21 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);
//--------------------------------------

        z02 = i64dilatation9_red( n02,n03,n04,
                    n12,n13,n14,
                    n22,n23,n24);

        z12 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z22 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);

        store2(Y,i,j,max);
//-----------------------1--------------------
        n00 = load2(X,i-2,j+3);
        n10 = load2(X,i-1,j+3);
        n20 = load2(X,i  ,j+3);
        n30 = load2(X,i+1,j+3);
        n40 = load2(X,i+2,j+3);
        
//--------------------------------------

        z00 = i64dilatation9_red( n01,n02,n03,
                    n11,n12,n13,
                    n21,n22,n23);

        z10 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z20 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);
//--------------------------------------

        z01 = i64dilatation9_red( n02,n03,n04,
                    n12,n13,n14,
                    n22,n23,n24);

        z11 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z21 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);
//--------------------------------------

        z02 = i64dilatation9_red( n03,n04,n00,
                    n13,n14,n10,
                    n23,n24,n20);

        z12 = i64dilatation9_red( n13,n14,n10,
                    n23,n24,n20,
                    n33,n34,n30);

        z22 = i64dilatation9_red( n23,n24,n20,
                    n33,n34,n30,
                    n43,n44,n40);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);


        store2(Y,i,j+1,max);

        break;
    
        case 3:
        n04 = load2(X,i-2,j+2);
        n14 = load2(X,i-1,j+2);
        n24 = load2(X,i  ,j+2);
        n34 = load2(X,i+1,j+2);
        n44 = load2(X,i+2,j+2);
//-----------------0---------------------

        z00 = i64dilatation9_red( n00,n01,n02,
                    n10,n11,n12,
                    n20,n21,n22);

        z10 = i64dilatation9_red( n10,n11,n12,
                    n20,n21,n22,
                    n30,n31,n32);

        z20 = i64dilatation9_red( n20,n21,n22,
                    n30,n31,n32,
                    n40,n41,n42);
//--------------------------------------

        z01 = i64dilatation9_red( n01,n02,n03,
                    n11,n12,n13,
                    n21,n22,n23);

        z11 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z21 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);
//--------------------------------------

        z02 = i64dilatation9_red( n02,n03,n04,
                    n12,n13,n14,
                    n22,n23,n24);

        z12 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z22 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);

        store2(Y,i,j,max);
//-----------------------1--------------------
        n00 = load2(X,i-2,j+3);
        n10 = load2(X,i-1,j+3);
        n20 = load2(X,i  ,j+3);
        n30 = load2(X,i+1,j+3);
        n40 = load2(X,i+2,j+3);
        
//--------------------------------------

        z00 = i64dilatation9_red( n01,n02,n03,
                    n11,n12,n13,
                    n21,n22,n23);

        z10 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z20 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);
//--------------------------------------

        z01 = i64dilatation9_red( n02,n03,n04,
                    n12,n13,n14,
                    n22,n23,n24);

        z11 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z21 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);
//--------------------------------------

        z02 = i64dilatation9_red( n03,n04,n00,
                    n13,n14,n10,
                    n23,n24,n20);

        z12 = i64dilatation9_red( n13,n14,n10,
                    n23,n24,n20,
                    n33,n34,n30);

        z22 = i64dilatation9_red( n23,n24,n20,
                    n33,n34,n30,
                    n43,n44,n40);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);


        store2(Y,i,j+1,max);

        n01 = load2(X,i-2,j+4);
        n11 = load2(X,i-1,j+4);
        n21 = load2(X,i  ,j+4);
        n31 = load2(X,i+1,j+4);
        n41 = load2(X,i+2,j+4);

        
//---------------2-----------------------

        z00 = i64dilatation9_red( n02,n03,n04,
                    n12,n13,n14,
                    n22,n23,n24);

        z10 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z20 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);
//--------------------------------------

        z01 = i64dilatation9_red( n03,n04,n00,
                    n13,n14,n10,
                    n23,n24,n20);

        z11 = i64dilatation9_red( n13,n14,n10,
                    n23,n24,n20,
                    n33,n34,n30);

        z21 = i64dilatation9_red( n23,n24,n20,
                    n33,n34,n30,
                    n43,n44,n40);
//--------------------------------------

        z02 = i64dilatation9_red( n04,n00,n01,
                    n14,n10,n11,
                    n24,n20,n21);

        z12 = i64dilatation9_red( n14,n10,n11,
                    n24,n20,n21,
                    n34,n30,n31);

        z22 = i64dilatation9_red( n24,n20,n21,
                    n34,n30,n31,
                    n44,n40,n41);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);
        
        

        store2(Y,i,j+2,max);

        break;

        case 4:
        n04 = load2(X,i-2,j+2);
        n14 = load2(X,i-1,j+2);
        n24 = load2(X,i  ,j+2);
        n34 = load2(X,i+1,j+2);
        n44 = load2(X,i+2,j+2);
//-----------------0---------------------

        z00 = i64dilatation9_red( n00,n01,n02,
                    n10,n11,n12,
                    n20,n21,n22);

        z10 = i64dilatation9_red( n10,n11,n12,
                    n20,n21,n22,
                    n30,n31,n32);

        z20 = i64dilatation9_red( n20,n21,n22,
                    n30,n31,n32,
                    n40,n41,n42);
//--------------------------------------

        z01 = i64dilatation9_red( n01,n02,n03,
                    n11,n12,n13,
                    n21,n22,n23);

        z11 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z21 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);
//--------------------------------------

        z02 = i64dilatation9_red( n02,n03,n04,
                    n12,n13,n14,
                    n22,n23,n24);

        z12 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z22 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);

        store2(Y,i,j,max);
//-----------------------1--------------------
        n00 = load2(X,i-2,j+3);
        n10 = load2(X,i-1,j+3);
        n20 = load2(X,i  ,j+3);
        n30 = load2(X,i+1,j+3);
        n40 = load2(X,i+2,j+3);
        
//--------------------------------------

        z00 = i64dilatation9_red( n01,n02,n03,
                    n11,n12,n13,
                    n21,n22,n23);

        z10 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z20 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);
//--------------------------------------

        z01 = i64dilatation9_red( n02,n03,n04,
                    n12,n13,n14,
                    n22,n23,n24);

        z11 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z21 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);
//--------------------------------------

        z02 = i64dilatation9_red( n03,n04,n00,
                    n13,n14,n10,
                    n23,n24,n20);

        z12 = i64dilatation9_red( n13,n14,n10,
                    n23,n24,n20,
                    n33,n34,n30);

        z22 = i64dilatation9_red( n23,n24,n20,
                    n33,n34,n30,
                    n43,n44,n40);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);


        store2(Y,i,j+1,max);

        n01 = load2(X,i-2,j+4);
        n11 = load2(X,i-1,j+4);
        n21 = load2(X,i  ,j+4);
        n31 = load2(X,i+1,j+4);
        n41 = load2(X,i+2,j+4);

        
//---------------2-----------------------

        z00 = i64dilatation9_red( n02,n03,n04,
                    n12,n13,n14,
                    n22,n23,n24);

        z10 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z20 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);
//--------------------------------------

        z01 = i64dilatation9_red( n03,n04,n00,
                    n13,n14,n10,
                    n23,n24,n20);

        z11 = i64dilatation9_red( n13,n14,n10,
                    n23,n24,n20,
                    n33,n34,n30);

        z21 = i64dilatation9_red( n23,n24,n20,
                    n33,n34,n30,
                    n43,n44,n40);
//--------------------------------------

        z02 = i64dilatation9_red( n04,n00,n01,
                    n14,n10,n11,
                    n24,n20,n21);

        z12 = i64dilatation9_red( n14,n10,n11,
                    n24,n20,n21,
                    n34,n30,n31);

        z22 = i64dilatation9_red( n24,n20,n21,
                    n34,n30,n31,
                    n44,n40,n41);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);
        
        

        store2(Y,i,j+2,max);

        n02 = load2(X,i-2,j+5);
        n12 = load2(X,i-1,j+5);
        n22 = load2(X,i  ,j+5);
        n32 = load2(X,i+1,j+5);
        n42 = load2(X,i+2,j+5);
//------------------3--------------------

        z00 = i64dilatation9_red( n03,n04,n00,
                    n13,n14,n10,
                    n23,n24,n20);

        z10 = i64dilatation9_red( n13,n14,n10,
                    n23,n24,n20,
                    n33,n34,n30);

        z20 = i64dilatation9_red( n23,n24,n20,
                    n33,n34,n30,
                    n43,n44,n40);
//--------------------------------------

        z01 = i64dilatation9_red( n04,n00,n01,
                    n14,n10,n11,
                    n24,n20,n21);

        z11 = i64dilatation9_red( n14,n10,n11,
                    n24,n20,n21,
                    n34,n30,n31);

        z21 = i64dilatation9_red( n24,n20,n21,
                    n34,n30,n31,
                    n44,n40,n41);
//--------------------------------------

        z02 = i64dilatation9_red( n00,n01,n02,
                    n10,n11,n12,
                    n20,n21,n22);

        z12 = i64dilatation9_red( n10,n11,n12,
                    n20,n21,n22,
                    n30,n31,n32);

        z22 = i64dilatation9_red( n20,n21,n22,
                    n30,n31,n32,
                    n40,n41,n42);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);
       

        store2(Y,i,j+3,max);

        break;
    }
}
// -------------------------------------------------------------------------------------------------
void line_fermeture3_ui64matrix_swp64_fusion_elu2_red(uint64 **X, int i, int j0, int j1, uint64 **Y)
// -------------------------------------------------------------------------------------------------
{
    uint64 x[6][5];
    uint64 z[3][3];
    uint64 max;

    for(int j = j0 ; j <= j1 ; j++ ){
        x[0][0]=load2(X,i-2,j-2);
        x[1][0]=load2(X,i-1,j-2);
        x[2][0]=load2(X,i  ,j-2);
        x[3][0]=load2(X,i+1,j-2);
        x[4][0]=load2(X,i+2,j-2);
        x[5][0]=load2(X,i+3,j-2);

        x[0][1]=load2(X,i-2,j-1);
        x[1][1]=load2(X,i-1,j-1);
        x[2][1]=load2(X,i  ,j-1);
        x[3][1]=load2(X,i+1,j-1);
        x[4][1]=load2(X,i+2,j-1);
        x[5][1]=load2(X,i+3,j-1);

        x[0][2]=load2(X,i-2,j  );
        x[1][2]=load2(X,i-1,j  );
        x[2][2]=load2(X,i  ,j  );
        x[3][2]=load2(X,i+1,j  );
        x[4][2]=load2(X,i+2,j  );
        x[5][2]=load2(X,i+3,j  );

        x[0][3]=load2(X,i-2,j+1);
        x[1][3]=load2(X,i-1,j+1);
        x[2][3]=load2(X,i  ,j+1);
        x[3][3]=load2(X,i+1,j+1);
        x[4][3]=load2(X,i+2,j+1);
        x[5][3]=load2(X,i+3,j+1);

        x[0][4]=load2(X,i-2,j+2);
        x[1][4]=load2(X,i-1,j+2);
        x[2][4]=load2(X,i  ,j+2);
        x[3][4]=load2(X,i+1,j+2);
        x[4][4]=load2(X,i+2,j+2);
        x[5][4]=load2(X,i+3,j+2);

        z[0][0]=i64dilatation_mat33(x,1,1);
        z[0][1]=i64dilatation_mat33(x,1,2);
        z[0][2]=i64dilatation_mat33(x,1,3);

        z[1][0]=i64dilatation_mat33(x,2,1);
        z[1][1]=i64dilatation_mat33(x,2,2);
        z[1][2]=i64dilatation_mat33(x,2,3);

        z[2][0]=i64dilatation_mat33(x,3,1);
        z[2][1]=i64dilatation_mat33(x,3,2);
        z[2][2]=i64dilatation_mat33(x,3,3);

        max=i64erosion_mat33(z,1,1);
        store2(Y,i,j,max);

        z[0][0]=i64dilatation_mat33(x,2,1);
        z[0][1]=i64dilatation_mat33(x,2,2);
        z[0][2]=i64dilatation_mat33(x,2,3);

        z[1][0]=i64dilatation_mat33(x,3,1);
        z[1][1]=i64dilatation_mat33(x,3,2);
        z[1][2]=i64dilatation_mat33(x,3,3);

        z[2][0]=i64dilatation_mat33(x,4,1);
        z[2][1]=i64dilatation_mat33(x,4,2);
        z[2][2]=i64dilatation_mat33(x,4,3);

        max=i64erosion_mat33(z,1,1);
        store2(Y,i+1,j,max);
    }
}
// ------------------------------------------------------------------------------------------------------
void line_fermeture3_ui64matrix_swp64_fusion_ilu3_elu2_red(uint64 **X, int i, int j0, int j1, uint64 **Y)
// ------------------------------------------------------------------------------------------------------
{
    //ilu5elu2...
    uint64  n00,n01,n02,n03,n04,
            n10,n11,n12,n13,n14,
            n20,n21,n22,n23,n24,
            n30,n31,n32,n33,n34,
            n40,n41,n42,n43,n44,
            n50,n51,n52,n53,n54;

    uint64   z00,z01,z02,
            z10,z11,z12,
            z20,z21,z22;

    uint64 max;

    int j, r = (j1 - j0 + 1) % 5;

    //prologue
        n00 = load2(X,i-2,j0-2);
        n10 = load2(X,i-1,j0-2);
        n20 = load2(X,i  ,j0-2);
        n30 = load2(X,i+1,j0-2);
        n40 = load2(X,i+2,j0-2);
        n50 = load2(X,i+3,j0-2);
        
        n01 = load2(X,i-2,j0-1);
        n11 = load2(X,i-1,j0-1);
        n21 = load2(X,i  ,j0-1);
        n31 = load2(X,i+1,j0-1);
        n41 = load2(X,i+2,j0-1);
        n51 = load2(X,i+3,j0-1);
        
        n02 = load2(X,i-2,j0  );
        n12 = load2(X,i-1,j0  );
        n22 = load2(X,i  ,j0  );
        n32 = load2(X,i+1,j0  );
        n42 = load2(X,i+2,j0  );
        n52 = load2(X,i+3,j0  );

        n03 = load2(X,i-2,j0+1);
        n13 = load2(X,i-1,j0+1);
        n23 = load2(X,i  ,j0+1);
        n33 = load2(X,i+1,j0+1);
        n43 = load2(X,i+2,j0+1);
        n53 = load2(X,i+3,j0+1);

    //boucle
    //boucle
    for(j = j0 ; j < j1 - r ; j = j + 5 ){
        n04 = load2(X,i-2,j+2);
        n14 = load2(X,i-1,j+2);
        n24 = load2(X,i  ,j+2);
        n34 = load2(X,i+1,j+2);
        n44 = load2(X,i+2,j+2);
        n54 = load2(X,i+3,j+2);
//-----------------0-0-------------------

        z00 = i64dilatation9_red( n00,n01,n02,
                    n10,n11,n12,
                    n20,n21,n22);

        z10 = i64dilatation9_red( n10,n11,n12,
                    n20,n21,n22,
                    n30,n31,n32);

        z20 = i64dilatation9_red( n20,n21,n22,
                    n30,n31,n32,
                    n40,n41,n42);
//--------------------------------------

        z01 = i64dilatation9_red( n01,n02,n03,
                    n11,n12,n13,
                    n21,n22,n23);

        z11 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z21 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);
//--------------------------------------

        z02 = i64dilatation9_red( n02,n03,n04,
                    n12,n13,n14,
                    n22,n23,n24);

        z12 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z22 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);

        store2(Y,i,j,max);

//-------------------0-1-------------------
        z00 = i64dilatation9_red( n10,n11,n12,
                    n20,n21,n22,
                    n30,n31,n32);

        z10 = i64dilatation9_red( n20,n21,n22,
                    n30,n31,n32,
                    n40,n41,n42);

        z20 = i64dilatation9_red( n30,n31,n32,
                    n40,n41,n42,
                    n50,n51,n52);
//--------------------------------------

        z01 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z11 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);

        z21 = i64dilatation9_red( n31,n32,n33,
                    n41,n42,n43,
                    n51,n52,n53);
//--------------------------------------

        z02 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z12 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);

        z22 = i64dilatation9_red( n32,n33,n34,
                    n42,n43,n44,
                    n52,n53,n54);
//--------------------------------------
        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);

        store2(Y,i+1,j,max);

//-----------------------1-0--------------------
        n00 = load2(X,i-2,j+3);
        n10 = load2(X,i-1,j+3);
        n20 = load2(X,i  ,j+3);
        n30 = load2(X,i+1,j+3);
        n40 = load2(X,i+2,j+3);
        n50 = load2(X,i+3,j+3);
        
//--------------------------------------

        z00 = i64dilatation9_red( n01,n02,n03,
                    n11,n12,n13,
                    n21,n22,n23);

        z10 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z20 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);
//--------------------------------------

        z01 = i64dilatation9_red( n02,n03,n04,
                    n12,n13,n14,
                    n22,n23,n24);

        z11 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z21 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);
//--------------------------------------

        z02 = i64dilatation9_red( n03,n04,n00,
                    n13,n14,n10,
                    n23,n24,n20);

        z12 = i64dilatation9_red( n13,n14,n10,
                    n23,n24,n20,
                    n33,n34,n30);

        z22 = i64dilatation9_red( n23,n24,n20,
                    n33,n34,n30,
                    n43,n44,n40);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);


        store2(Y,i,j+1,max);

//-------------------1-1------------------

        z00 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z10 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);

        z20 = i64dilatation9_red( n31,n32,n33,
                    n41,n42,n43,
                    n51,n52,n53);
//--------------------------------------

        z01 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z11 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);

        z21 = i64dilatation9_red( n32,n33,n34,
                    n42,n43,n44,
                    n52,n53,n54);
//--------------------------------------

        z02 = i64dilatation9_red( n13,n14,n10,
                    n23,n24,n20,
                    n33,n34,n30);

        z12 = i64dilatation9_red( n23,n24,n20,
                    n33,n34,n30,
                    n43,n44,n40);

        z22 = i64dilatation9_red( n33,n34,n30,
                    n43,n44,n40,
                    n53,n54,n50);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);


        store2(Y,i+1,j+1,max);


        n01 = load2(X,i-2,j+4);
        n11 = load2(X,i-1,j+4);
        n21 = load2(X,i  ,j+4);
        n31 = load2(X,i+1,j+4);
        n41 = load2(X,i+2,j+4);
        n51 = load2(X,i+3,j+4);
        
//---------------2-0----------------------

        z00 = i64dilatation9_red( n02,n03,n04,
                    n12,n13,n14,
                    n22,n23,n24);

        z10 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z20 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);
//--------------------------------------

        z01 = i64dilatation9_red( n03,n04,n00,
                    n13,n14,n10,
                    n23,n24,n20);

        z11 = i64dilatation9_red( n13,n14,n10,
                    n23,n24,n20,
                    n33,n34,n30);

        z21 = i64dilatation9_red( n23,n24,n20,
                    n33,n34,n30,
                    n43,n44,n40);
//--------------------------------------

        z02 = i64dilatation9_red( n04,n00,n01,
                    n14,n10,n11,
                    n24,n20,n21);

        z12 = i64dilatation9_red( n14,n10,n11,
                    n24,n20,n21,
                    n34,n30,n31);

        z22 = i64dilatation9_red( n24,n20,n21,
                    n34,n30,n31,
                    n44,n40,n41);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);
        
        store2(Y,i,j+2,max);

//---------------2-1----------------------

        z00 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z10 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);

        z20 = i64dilatation9_red( n32,n33,n34,
                    n42,n43,n44,
                    n52,n53,n54);
//--------------------------------------

        z01 = i64dilatation9_red( n13,n14,n10,
                    n23,n24,n20,
                    n33,n34,n30);

        z11 = i64dilatation9_red( n23,n24,n20,
                    n33,n34,n30,
                    n43,n44,n40);

        z21 = i64dilatation9_red( n33,n34,n30,
                    n43,n44,n40,
                    n53,n54,n50);
//--------------------------------------

        z02 = i64dilatation9_red( n14,n10,n11,
                    n24,n20,n21,
                    n34,n30,n31);

        z12 = i64dilatation9_red( n24,n20,n21,
                    n34,n30,n31,
                    n44,n40,n41);

        z22 = i64dilatation9_red( n34,n30,n31,
                    n44,n40,n41,
                    n54,n50,n51);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);
        
        store2(Y,i+1,j+2,max);

        n02 = load2(X,i-2,j+5);
        n12 = load2(X,i-1,j+5);
        n22 = load2(X,i  ,j+5);
        n32 = load2(X,i+1,j+5);
        n42 = load2(X,i+2,j+5);
        n52 = load2(X,i+3,j+5);
//------------------3-0--------------------

        z00 = i64dilatation9_red( n03,n04,n00,
                    n13,n14,n10,
                    n23,n24,n20);

        z10 = i64dilatation9_red( n13,n14,n10,
                    n23,n24,n20,
                    n33,n34,n30);

        z20 = i64dilatation9_red( n23,n24,n20,
                    n33,n34,n30,
                    n43,n44,n40);
//--------------------------------------

        z01 = i64dilatation9_red( n04,n00,n01,
                    n14,n10,n11,
                    n24,n20,n21);

        z11 = i64dilatation9_red( n14,n10,n11,
                    n24,n20,n21,
                    n34,n30,n31);

        z21 = i64dilatation9_red( n24,n20,n21,
                    n34,n30,n31,
                    n44,n40,n41);
//--------------------------------------

        z02 = i64dilatation9_red( n00,n01,n02,
                    n10,n11,n12,
                    n20,n21,n22);

        z12 = i64dilatation9_red( n10,n11,n12,
                    n20,n21,n22,
                    n30,n31,n32);

        z22 = i64dilatation9_red( n20,n21,n22,
                    n30,n31,n32,
                    n40,n41,n42);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);
       
        store2(Y,i,j+3,max);

//------------------3-1--------------------

        z00 = i64dilatation9_red( n13,n14,n10,
                    n23,n24,n20,
                    n33,n34,n30);

        z10 = i64dilatation9_red( n23,n24,n20,
                    n33,n34,n30,
                    n43,n44,n40);

        z20 = i64dilatation9_red( n33,n34,n30,
                    n43,n44,n40,
                    n53,n54,n50);
//--------------------------------------

        z01 = i64dilatation9_red( n14,n10,n11,
                    n24,n20,n21,
                    n34,n30,n31);

        z11 = i64dilatation9_red( n24,n20,n21,
                    n34,n30,n31,
                    n44,n40,n41);

        z21 = i64dilatation9_red( n34,n30,n31,
                    n44,n40,n41,
                    n54,n50,n51);
//--------------------------------------

        z02 = i64dilatation9_red( n10,n11,n12,
                    n20,n21,n22,
                    n30,n31,n32);

        z12 = i64dilatation9_red( n20,n21,n22,
                    n30,n31,n32,
                    n40,n41,n42);

        z22 = i64dilatation9_red( n30,n31,n32,
                    n40,n41,n42,
                    n50,n51,n52);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);
       
        store2(Y,i+1,j+3,max);

        n03 = load2(X,i-2,j+6);
        n13 = load2(X,i-1,j+6);
        n23 = load2(X,i  ,j+6);
        n33 = load2(X,i+1,j+6);
        n43 = load2(X,i+2,j+6);
        n53 = load2(X,i+3,j+6);
//-----------------4-0--------------------

        z00 = i64dilatation9_red( n04,n00,n01,
                    n14,n10,n11,
                    n24,n20,n21);

        z10 = i64dilatation9_red( n14,n10,n11,
                    n24,n20,n21,
                    n34,n30,n31);

        z20 = i64dilatation9_red( n24,n20,n21,
                    n34,n30,n31,
                    n44,n40,n41);
//--------------------------------------

        z01 = i64dilatation9_red( n00,n01,n02,
                    n10,n11,n12,
                    n20,n21,n22);

        z11 = i64dilatation9_red( n10,n11,n12,
                    n20,n21,n22,
                    n30,n31,n32);

        z21 = i64dilatation9_red( n20,n21,n22,
                    n30,n31,n32,
                    n40,n41,n42);
//--------------------------------------

        z02 = i64dilatation9_red( n01,n02,n03,
                    n11,n12,n13,
                    n21,n22,n23);

        z12 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z22 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);

        store2(Y,i,j+4,max);

    //-----------------4-1--------------------

        z00 = i64dilatation9_red( n14,n10,n11,
                    n24,n20,n21,
                    n34,n30,n31);

        z10 = i64dilatation9_red( n24,n20,n21,
                    n34,n30,n31,
                    n44,n40,n41);

        z20 = i64dilatation9_red( n34,n30,n31,
                    n44,n40,n41,
                    n54,n50,n51);
//--------------------------------------

        z01 = i64dilatation9_red( n10,n11,n12,
                    n20,n21,n22,
                    n30,n31,n32);

        z11 = i64dilatation9_red( n20,n21,n22,
                    n30,n31,n32,
                    n40,n41,n42);

        z21 = i64dilatation9_red( n30,n31,n32,
                    n40,n41,n42,
                    n50,n51,n52);
//--------------------------------------

        z02 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z12 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);

        z22 = i64dilatation9_red( n31,n32,n33,
                    n41,n42,n43,
                    n51,n52,n53);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);

        store2(Y,i+1,j+4,max);
    }

    //epilogue
    switch (r)
    {
        case 1:
        n04 = load2(X,i-2,j+2);
        n14 = load2(X,i-1,j+2);
        n24 = load2(X,i  ,j+2);
        n34 = load2(X,i+1,j+2);
        n44 = load2(X,i+2,j+2);
        n54 = load2(X,i+3,j+2);
//-----------------0-0-------------------

        z00 = i64dilatation9_red( n00,n01,n02,
                    n10,n11,n12,
                    n20,n21,n22);

        z10 = i64dilatation9_red( n10,n11,n12,
                    n20,n21,n22,
                    n30,n31,n32);

        z20 = i64dilatation9_red( n20,n21,n22,
                    n30,n31,n32,
                    n40,n41,n42);
//--------------------------------------

        z01 = i64dilatation9_red( n01,n02,n03,
                    n11,n12,n13,
                    n21,n22,n23);

        z11 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z21 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);
//--------------------------------------

        z02 = i64dilatation9_red( n02,n03,n04,
                    n12,n13,n14,
                    n22,n23,n24);

        z12 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z22 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);

        store2(Y,i,j,max);

//-------------------0-1-------------------
        z00 = i64dilatation9_red( n10,n11,n12,
                    n20,n21,n22,
                    n30,n31,n32);

        z10 = i64dilatation9_red( n20,n21,n22,
                    n30,n31,n32,
                    n40,n41,n42);

        z20 = i64dilatation9_red( n30,n31,n32,
                    n40,n41,n42,
                    n50,n51,n52);
//--------------------------------------

        z01 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z11 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);

        z21 = i64dilatation9_red( n31,n32,n33,
                    n41,n42,n43,
                    n51,n52,n53);
//--------------------------------------

        z02 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z12 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);

        z22 = i64dilatation9_red( n32,n33,n34,
                    n42,n43,n44,
                    n52,n53,n54);
//--------------------------------------
        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);

        store2(Y,i+1,j,max);

        break;

        case 2:
        n04 = load2(X,i-2,j+2);
        n14 = load2(X,i-1,j+2);
        n24 = load2(X,i  ,j+2);
        n34 = load2(X,i+1,j+2);
        n44 = load2(X,i+2,j+2);
        n54 = load2(X,i+3,j+2);
//-----------------0-0-------------------

        z00 = i64dilatation9_red( n00,n01,n02,
                    n10,n11,n12,
                    n20,n21,n22);

        z10 = i64dilatation9_red( n10,n11,n12,
                    n20,n21,n22,
                    n30,n31,n32);

        z20 = i64dilatation9_red( n20,n21,n22,
                    n30,n31,n32,
                    n40,n41,n42);
//--------------------------------------

        z01 = i64dilatation9_red( n01,n02,n03,
                    n11,n12,n13,
                    n21,n22,n23);

        z11 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z21 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);
//--------------------------------------

        z02 = i64dilatation9_red( n02,n03,n04,
                    n12,n13,n14,
                    n22,n23,n24);

        z12 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z22 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);

        store2(Y,i,j,max);

//-------------------0-1-------------------
        z00 = i64dilatation9_red( n10,n11,n12,
                    n20,n21,n22,
                    n30,n31,n32);

        z10 = i64dilatation9_red( n20,n21,n22,
                    n30,n31,n32,
                    n40,n41,n42);

        z20 = i64dilatation9_red( n30,n31,n32,
                    n40,n41,n42,
                    n50,n51,n52);
//--------------------------------------

        z01 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z11 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);

        z21 = i64dilatation9_red( n31,n32,n33,
                    n41,n42,n43,
                    n51,n52,n53);
//--------------------------------------

        z02 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z12 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);

        z22 = i64dilatation9_red( n32,n33,n34,
                    n42,n43,n44,
                    n52,n53,n54);
//--------------------------------------
        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);

        store2(Y,i+1,j,max);

//-----------------------1-0--------------------
        n00 = load2(X,i-2,j+3);
        n10 = load2(X,i-1,j+3);
        n20 = load2(X,i  ,j+3);
        n30 = load2(X,i+1,j+3);
        n40 = load2(X,i+2,j+3);
        n50 = load2(X,i+3,j+3);
        
//--------------------------------------

        z00 = i64dilatation9_red( n01,n02,n03,
                    n11,n12,n13,
                    n21,n22,n23);

        z10 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z20 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);
//--------------------------------------

        z01 = i64dilatation9_red( n02,n03,n04,
                    n12,n13,n14,
                    n22,n23,n24);

        z11 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z21 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);
//--------------------------------------

        z02 = i64dilatation9_red( n03,n04,n00,
                    n13,n14,n10,
                    n23,n24,n20);

        z12 = i64dilatation9_red( n13,n14,n10,
                    n23,n24,n20,
                    n33,n34,n30);

        z22 = i64dilatation9_red( n23,n24,n20,
                    n33,n34,n30,
                    n43,n44,n40);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);


        store2(Y,i,j+1,max);

//-------------------1-1------------------

        z00 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z10 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);

        z20 = i64dilatation9_red( n31,n32,n33,
                    n41,n42,n43,
                    n51,n52,n53);
//--------------------------------------

        z01 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z11 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);

        z21 = i64dilatation9_red( n32,n33,n34,
                    n42,n43,n44,
                    n52,n53,n54);
//--------------------------------------

        z02 = i64dilatation9_red( n13,n14,n10,
                    n23,n24,n20,
                    n33,n34,n30);

        z12 = i64dilatation9_red( n23,n24,n20,
                    n33,n34,n30,
                    n43,n44,n40);

        z22 = i64dilatation9_red( n33,n34,n30,
                    n43,n44,n40,
                    n53,n54,n50);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);


        store2(Y,i+1,j+1,max);

        break;
    
        case 3:
        n04 = load2(X,i-2,j+2);
        n14 = load2(X,i-1,j+2);
        n24 = load2(X,i  ,j+2);
        n34 = load2(X,i+1,j+2);
        n44 = load2(X,i+2,j+2);
        n54 = load2(X,i+3,j+2);
//-----------------0-0-------------------

        z00 = i64dilatation9_red( n00,n01,n02,
                    n10,n11,n12,
                    n20,n21,n22);

        z10 = i64dilatation9_red( n10,n11,n12,
                    n20,n21,n22,
                    n30,n31,n32);

        z20 = i64dilatation9_red( n20,n21,n22,
                    n30,n31,n32,
                    n40,n41,n42);
//--------------------------------------

        z01 = i64dilatation9_red( n01,n02,n03,
                    n11,n12,n13,
                    n21,n22,n23);

        z11 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z21 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);
//--------------------------------------

        z02 = i64dilatation9_red( n02,n03,n04,
                    n12,n13,n14,
                    n22,n23,n24);

        z12 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z22 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);

        store2(Y,i,j,max);

//-------------------0-1-------------------
        z00 = i64dilatation9_red( n10,n11,n12,
                    n20,n21,n22,
                    n30,n31,n32);

        z10 = i64dilatation9_red( n20,n21,n22,
                    n30,n31,n32,
                    n40,n41,n42);

        z20 = i64dilatation9_red( n30,n31,n32,
                    n40,n41,n42,
                    n50,n51,n52);
//--------------------------------------

        z01 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z11 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);

        z21 = i64dilatation9_red( n31,n32,n33,
                    n41,n42,n43,
                    n51,n52,n53);
//--------------------------------------

        z02 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z12 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);

        z22 = i64dilatation9_red( n32,n33,n34,
                    n42,n43,n44,
                    n52,n53,n54);
//--------------------------------------
        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);

        store2(Y,i+1,j,max);

//-----------------------1-0--------------------
        n00 = load2(X,i-2,j+3);
        n10 = load2(X,i-1,j+3);
        n20 = load2(X,i  ,j+3);
        n30 = load2(X,i+1,j+3);
        n40 = load2(X,i+2,j+3);
        n50 = load2(X,i+3,j+3);
        
//--------------------------------------

        z00 = i64dilatation9_red( n01,n02,n03,
                    n11,n12,n13,
                    n21,n22,n23);

        z10 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z20 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);
//--------------------------------------

        z01 = i64dilatation9_red( n02,n03,n04,
                    n12,n13,n14,
                    n22,n23,n24);

        z11 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z21 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);
//--------------------------------------

        z02 = i64dilatation9_red( n03,n04,n00,
                    n13,n14,n10,
                    n23,n24,n20);

        z12 = i64dilatation9_red( n13,n14,n10,
                    n23,n24,n20,
                    n33,n34,n30);

        z22 = i64dilatation9_red( n23,n24,n20,
                    n33,n34,n30,
                    n43,n44,n40);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);


        store2(Y,i,j+1,max);

//-------------------1-1------------------

        z00 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z10 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);

        z20 = i64dilatation9_red( n31,n32,n33,
                    n41,n42,n43,
                    n51,n52,n53);
//--------------------------------------

        z01 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z11 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);

        z21 = i64dilatation9_red( n32,n33,n34,
                    n42,n43,n44,
                    n52,n53,n54);
//--------------------------------------

        z02 = i64dilatation9_red( n13,n14,n10,
                    n23,n24,n20,
                    n33,n34,n30);

        z12 = i64dilatation9_red( n23,n24,n20,
                    n33,n34,n30,
                    n43,n44,n40);

        z22 = i64dilatation9_red( n33,n34,n30,
                    n43,n44,n40,
                    n53,n54,n50);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);


        store2(Y,i+1,j+1,max);


        n01 = load2(X,i-2,j+4);
        n11 = load2(X,i-1,j+4);
        n21 = load2(X,i  ,j+4);
        n31 = load2(X,i+1,j+4);
        n41 = load2(X,i+2,j+4);
        n51 = load2(X,i+3,j+4);
        
//---------------2-0----------------------

        z00 = i64dilatation9_red( n02,n03,n04,
                    n12,n13,n14,
                    n22,n23,n24);

        z10 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z20 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);
//--------------------------------------

        z01 = i64dilatation9_red( n03,n04,n00,
                    n13,n14,n10,
                    n23,n24,n20);

        z11 = i64dilatation9_red( n13,n14,n10,
                    n23,n24,n20,
                    n33,n34,n30);

        z21 = i64dilatation9_red( n23,n24,n20,
                    n33,n34,n30,
                    n43,n44,n40);
//--------------------------------------

        z02 = i64dilatation9_red( n04,n00,n01,
                    n14,n10,n11,
                    n24,n20,n21);

        z12 = i64dilatation9_red( n14,n10,n11,
                    n24,n20,n21,
                    n34,n30,n31);

        z22 = i64dilatation9_red( n24,n20,n21,
                    n34,n30,n31,
                    n44,n40,n41);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);
        
        store2(Y,i,j+2,max);

//---------------2-1----------------------

        z00 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z10 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);

        z20 = i64dilatation9_red( n32,n33,n34,
                    n42,n43,n44,
                    n52,n53,n54);
//--------------------------------------

        z01 = i64dilatation9_red( n13,n14,n10,
                    n23,n24,n20,
                    n33,n34,n30);

        z11 = i64dilatation9_red( n23,n24,n20,
                    n33,n34,n30,
                    n43,n44,n40);

        z21 = i64dilatation9_red( n33,n34,n30,
                    n43,n44,n40,
                    n53,n54,n50);
//--------------------------------------

        z02 = i64dilatation9_red( n14,n10,n11,
                    n24,n20,n21,
                    n34,n30,n31);

        z12 = i64dilatation9_red( n24,n20,n21,
                    n34,n30,n31,
                    n44,n40,n41);

        z22 = i64dilatation9_red( n34,n30,n31,
                    n44,n40,n41,
                    n54,n50,n51);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);
        
        store2(Y,i+1,j+2,max);

        break;

        case 4:
        n04 = load2(X,i-2,j+2);
        n14 = load2(X,i-1,j+2);
        n24 = load2(X,i  ,j+2);
        n34 = load2(X,i+1,j+2);
        n44 = load2(X,i+2,j+2);
        n54 = load2(X,i+3,j+2);
//-----------------0-0-------------------

        z00 = i64dilatation9_red( n00,n01,n02,
                    n10,n11,n12,
                    n20,n21,n22);

        z10 = i64dilatation9_red( n10,n11,n12,
                    n20,n21,n22,
                    n30,n31,n32);

        z20 = i64dilatation9_red( n20,n21,n22,
                    n30,n31,n32,
                    n40,n41,n42);
//--------------------------------------

        z01 = i64dilatation9_red( n01,n02,n03,
                    n11,n12,n13,
                    n21,n22,n23);

        z11 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z21 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);
//--------------------------------------

        z02 = i64dilatation9_red( n02,n03,n04,
                    n12,n13,n14,
                    n22,n23,n24);

        z12 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z22 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);

        store2(Y,i,j,max);

//-------------------0-1-------------------
        z00 = i64dilatation9_red( n10,n11,n12,
                    n20,n21,n22,
                    n30,n31,n32);

        z10 = i64dilatation9_red( n20,n21,n22,
                    n30,n31,n32,
                    n40,n41,n42);

        z20 = i64dilatation9_red( n30,n31,n32,
                    n40,n41,n42,
                    n50,n51,n52);
//--------------------------------------

        z01 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z11 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);

        z21 = i64dilatation9_red( n31,n32,n33,
                    n41,n42,n43,
                    n51,n52,n53);
//--------------------------------------

        z02 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z12 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);

        z22 = i64dilatation9_red( n32,n33,n34,
                    n42,n43,n44,
                    n52,n53,n54);
//--------------------------------------
        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);

        store2(Y,i+1,j,max);

//-----------------------1-0--------------------
        n00 = load2(X,i-2,j+3);
        n10 = load2(X,i-1,j+3);
        n20 = load2(X,i  ,j+3);
        n30 = load2(X,i+1,j+3);
        n40 = load2(X,i+2,j+3);
        n50 = load2(X,i+3,j+3);
        
//--------------------------------------

        z00 = i64dilatation9_red( n01,n02,n03,
                    n11,n12,n13,
                    n21,n22,n23);

        z10 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z20 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);
//--------------------------------------

        z01 = i64dilatation9_red( n02,n03,n04,
                    n12,n13,n14,
                    n22,n23,n24);

        z11 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z21 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);
//--------------------------------------

        z02 = i64dilatation9_red( n03,n04,n00,
                    n13,n14,n10,
                    n23,n24,n20);

        z12 = i64dilatation9_red( n13,n14,n10,
                    n23,n24,n20,
                    n33,n34,n30);

        z22 = i64dilatation9_red( n23,n24,n20,
                    n33,n34,n30,
                    n43,n44,n40);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);


        store2(Y,i,j+1,max);

//-------------------1-1------------------

        z00 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z10 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);

        z20 = i64dilatation9_red( n31,n32,n33,
                    n41,n42,n43,
                    n51,n52,n53);
//--------------------------------------

        z01 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z11 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);

        z21 = i64dilatation9_red( n32,n33,n34,
                    n42,n43,n44,
                    n52,n53,n54);
//--------------------------------------

        z02 = i64dilatation9_red( n13,n14,n10,
                    n23,n24,n20,
                    n33,n34,n30);

        z12 = i64dilatation9_red( n23,n24,n20,
                    n33,n34,n30,
                    n43,n44,n40);

        z22 = i64dilatation9_red( n33,n34,n30,
                    n43,n44,n40,
                    n53,n54,n50);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);


        store2(Y,i+1,j+1,max);


        n01 = load2(X,i-2,j+4);
        n11 = load2(X,i-1,j+4);
        n21 = load2(X,i  ,j+4);
        n31 = load2(X,i+1,j+4);
        n41 = load2(X,i+2,j+4);
        n51 = load2(X,i+3,j+4);
        
//---------------2-0----------------------

        z00 = i64dilatation9_red( n02,n03,n04,
                    n12,n13,n14,
                    n22,n23,n24);

        z10 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z20 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);
//--------------------------------------

        z01 = i64dilatation9_red( n03,n04,n00,
                    n13,n14,n10,
                    n23,n24,n20);

        z11 = i64dilatation9_red( n13,n14,n10,
                    n23,n24,n20,
                    n33,n34,n30);

        z21 = i64dilatation9_red( n23,n24,n20,
                    n33,n34,n30,
                    n43,n44,n40);
//--------------------------------------

        z02 = i64dilatation9_red( n04,n00,n01,
                    n14,n10,n11,
                    n24,n20,n21);

        z12 = i64dilatation9_red( n14,n10,n11,
                    n24,n20,n21,
                    n34,n30,n31);

        z22 = i64dilatation9_red( n24,n20,n21,
                    n34,n30,n31,
                    n44,n40,n41);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);
        
        store2(Y,i,j+2,max);

//---------------2-1----------------------

        z00 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z10 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);

        z20 = i64dilatation9_red( n32,n33,n34,
                    n42,n43,n44,
                    n52,n53,n54);
//--------------------------------------

        z01 = i64dilatation9_red( n13,n14,n10,
                    n23,n24,n20,
                    n33,n34,n30);

        z11 = i64dilatation9_red( n23,n24,n20,
                    n33,n34,n30,
                    n43,n44,n40);

        z21 = i64dilatation9_red( n33,n34,n30,
                    n43,n44,n40,
                    n53,n54,n50);
//--------------------------------------

        z02 = i64dilatation9_red( n14,n10,n11,
                    n24,n20,n21,
                    n34,n30,n31);

        z12 = i64dilatation9_red( n24,n20,n21,
                    n34,n30,n31,
                    n44,n40,n41);

        z22 = i64dilatation9_red( n34,n30,n31,
                    n44,n40,n41,
                    n54,n50,n51);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);
        
        store2(Y,i+1,j+2,max);

        n02 = load2(X,i-2,j+5);
        n12 = load2(X,i-1,j+5);
        n22 = load2(X,i  ,j+5);
        n32 = load2(X,i+1,j+5);
        n42 = load2(X,i+2,j+5);
        n52 = load2(X,i+3,j+5);
//------------------3-0--------------------

        z00 = i64dilatation9_red( n03,n04,n00,
                    n13,n14,n10,
                    n23,n24,n20);

        z10 = i64dilatation9_red( n13,n14,n10,
                    n23,n24,n20,
                    n33,n34,n30);

        z20 = i64dilatation9_red( n23,n24,n20,
                    n33,n34,n30,
                    n43,n44,n40);
//--------------------------------------

        z01 = i64dilatation9_red( n04,n00,n01,
                    n14,n10,n11,
                    n24,n20,n21);

        z11 = i64dilatation9_red( n14,n10,n11,
                    n24,n20,n21,
                    n34,n30,n31);

        z21 = i64dilatation9_red( n24,n20,n21,
                    n34,n30,n31,
                    n44,n40,n41);
//--------------------------------------

        z02 = i64dilatation9_red( n00,n01,n02,
                    n10,n11,n12,
                    n20,n21,n22);

        z12 = i64dilatation9_red( n10,n11,n12,
                    n20,n21,n22,
                    n30,n31,n32);

        z22 = i64dilatation9_red( n20,n21,n22,
                    n30,n31,n32,
                    n40,n41,n42);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);
       
        store2(Y,i,j+3,max);

//------------------3-1--------------------

        z00 = i64dilatation9_red( n13,n14,n10,
                    n23,n24,n20,
                    n33,n34,n30);

        z10 = i64dilatation9_red( n23,n24,n20,
                    n33,n34,n30,
                    n43,n44,n40);

        z20 = i64dilatation9_red( n33,n34,n30,
                    n43,n44,n40,
                    n53,n54,n50);
//--------------------------------------

        z01 = i64dilatation9_red( n14,n10,n11,
                    n24,n20,n21,
                    n34,n30,n31);

        z11 = i64dilatation9_red( n24,n20,n21,
                    n34,n30,n31,
                    n44,n40,n41);

        z21 = i64dilatation9_red( n34,n30,n31,
                    n44,n40,n41,
                    n54,n50,n51);
//--------------------------------------

        z02 = i64dilatation9_red( n10,n11,n12,
                    n20,n21,n22,
                    n30,n31,n32);

        z12 = i64dilatation9_red( n20,n21,n22,
                    n30,n31,n32,
                    n40,n41,n42);

        z22 = i64dilatation9_red( n30,n31,n32,
                    n40,n41,n42,
                    n50,n51,n52);
//--------------------------------------

        max = i64erosion9_red(z00,z01,z02,z10,z11,z12,z20,z21,z22);
       
        store2(Y,i+1,j+3,max);

        break;
    }
}
// -------------------------------------------------------------------------------------------------------------
void line_fermeture3_ui64matrix_swp64_fusion_ilu3_elu2_red_factor(uint64 **X, int i, int j0, int j1, uint64 **Y)
// -------------------------------------------------------------------------------------------------------------
{
    uint64 n00,n01,n02,n03,n04,
           n10,n11,n12,n13,n14,
           n20,n21,n22,n23,n24,
           n30,n31,n32,n33,n34, 
           n40,n41,n42,n43,n44,
           n50,n51,n52,n53,n54;

    uint64  z00,z01,z02,
            z10,z11,z12,
            z20,z21,z22;

    uint64 factor;

    uint64 max;

    int j, r = (j1 - j0 + 1) % 5;

    //prologue
        n00 = load2(X,i-2,j0-2);
        n10 = load2(X,i-1,j0-2);
        n20 = load2(X,i  ,j0-2);
        n30 = load2(X,i+1,j0-2);
        n40 = load2(X,i+2,j0-2);
        n50 = load2(X,i+3,j0-2);
        
        n01 = load2(X,i-2,j0-1);
        n11 = load2(X,i-1,j0-1);
        n21 = load2(X,i  ,j0-1);
        n31 = load2(X,i+1,j0-1);
        n41 = load2(X,i+2,j0-1);
        n51 = load2(X,i+3,j0-1);
        
        n02 = load2(X,i-2,j0  );
        n12 = load2(X,i-1,j0  );
        n22 = load2(X,i  ,j0  );
        n32 = load2(X,i+1,j0  );
        n42 = load2(X,i+2,j0  );
        n52 = load2(X,i+3,j0  );

        n03 = load2(X,i-2,j0+1);
        n13 = load2(X,i-1,j0+1);
        n23 = load2(X,i  ,j0+1);
        n33 = load2(X,i+1,j0+1);
        n43 = load2(X,i+2,j0+1);
        n53 = load2(X,i+3,j0+1);

    //boucle
    for(j = j0 ; j < j1 - r ; j = j + 5 ){
        n04 = load2(X,i-2,j+2);
        n14 = load2(X,i-1,j+2);
        n24 = load2(X,i  ,j+2);
        n34 = load2(X,i+1,j+2);
        n44 = load2(X,i+2,j+2);
        n54 = load2(X,i+3,j+2);

//-----------------calc factor-------------------
        z10 = i64dilatation9_red( n10,n11,n12,
                    n20,n21,n22,
                    n30,n31,n32);

        z20 = i64dilatation9_red( n20,n21,n22,
                    n30,n31,n32,
                    n40,n41,n42);

        z11 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z21 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);

        z12 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z22 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);

        factor = i64erosion6_red(z10,z11,z12,z20,z21,z22);

//-----------------0-0-------------------

        z00 = i64dilatation9_red( n00,n01,n02,
                    n10,n11,n12,
                    n20,n21,n22);

        z01 = i64dilatation9_red( n01,n02,n03,
                    n11,n12,n13,
                    n21,n22,n23);

        z02 = i64dilatation9_red( n02,n03,n04,
                    n12,n13,n14,
                    n22,n23,n24);

        max = i64erosion4(z00,z01,z02,factor);

        store2(Y,i,j,max);

//-------------------0-1-------------------

        z20 = i64dilatation9_red( n30,n31,n32,
                    n40,n41,n42,
                    n50,n51,n52);

        z21 = i64dilatation9_red( n31,n32,n33,
                    n41,n42,n43,
                    n51,n52,n53);

        z22 = i64dilatation9_red( n32,n33,n34,
                    n42,n43,n44,
                    n52,n53,n54);

        max = i64erosion4(z20,z21,z22,factor);

        store2(Y,i+1,j,max);

//-----------------------1-0--------------------
        n00 = load2(X,i-2,j+3);
        n10 = load2(X,i-1,j+3);
        n20 = load2(X,i  ,j+3);
        n30 = load2(X,i+1,j+3);
        n40 = load2(X,i+2,j+3);
        n50 = load2(X,i+3,j+3);
        
//-----------------calc factor-------------------
        z10 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z20 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);
                    
        z11 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z21 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);

        z12 = i64dilatation9_red( n13,n14,n10,
                    n23,n24,n20,
                    n33,n34,n30);

        z22 = i64dilatation9_red( n23,n24,n20,
                    n33,n34,n30,
                    n43,n44,n40);

        factor = i64erosion6_red(z10,z11,z12,z20,z21,z22);
//----------------1-0---------------------

        z00 = i64dilatation9_red( n01,n02,n03,
                    n11,n12,n13,
                    n21,n22,n23);

        z01 = i64dilatation9_red( n02,n03,n04,
                    n12,n13,n14,
                    n22,n23,n24);

        z02 = i64dilatation9_red( n03,n04,n00,
                    n13,n14,n10,
                    n23,n24,n20);

        max = i64erosion4(z00,z01,z02,factor);

        store2(Y,i,j+1,max);

//-------------------1-1------------------

        z20 = i64dilatation9_red( n31,n32,n33,
                    n41,n42,n43,
                    n51,n52,n53);


        z21 = i64dilatation9_red( n32,n33,n34,
                    n42,n43,n44,
                    n52,n53,n54);


        z22 = i64dilatation9_red( n33,n34,n30,
                    n43,n44,n40,
                    n53,n54,n50);


        max = i64erosion4(z20,z21,z22,factor);


        store2(Y,i+1,j+1,max);


        n01 = load2(X,i-2,j+4);
        n11 = load2(X,i-1,j+4);
        n21 = load2(X,i  ,j+4);
        n31 = load2(X,i+1,j+4);
        n41 = load2(X,i+2,j+4);
        n51 = load2(X,i+3,j+4);
        
//-----------------calc factor-------------------
        z10 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z20 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);
        
        z11 = i64dilatation9_red( n13,n14,n10,
                    n23,n24,n20,
                    n33,n34,n30);

        z21 = i64dilatation9_red( n23,n24,n20,
                    n33,n34,n30,
                    n43,n44,n40);

        z12 = i64dilatation9_red( n14,n10,n11,
                    n24,n20,n21,
                    n34,n30,n31);

        z22 = i64dilatation9_red( n24,n20,n21,
                    n34,n30,n31,
                    n44,n40,n41);
        
        factor = i64erosion6_red(z10,z11,z12,z20,z21,z22);
//---------------2-0----------------------

        z00 = i64dilatation9_red( n02,n03,n04,
                    n12,n13,n14,
                    n22,n23,n24);

        z01 = i64dilatation9_red( n03,n04,n00,
                    n13,n14,n10,
                    n23,n24,n20);

        z02 = i64dilatation9_red( n04,n00,n01,
                    n14,n10,n11,
                    n24,n20,n21);

        max = i64erosion4(z00,z01,z02,factor);
        
        store2(Y,i,j+2,max);

//---------------2-1----------------------

        z20 = i64dilatation9_red( n32,n33,n34,
                    n42,n43,n44,
                    n52,n53,n54);

        z21 = i64dilatation9_red( n33,n34,n30,
                    n43,n44,n40,
                    n53,n54,n50);

        z22 = i64dilatation9_red( n34,n30,n31,
                    n44,n40,n41,
                    n54,n50,n51);
//--------------------------------------

        max = i64erosion4(z20,z21,z22,factor);
        
        store2(Y,i+1,j+2,max);

        n02 = load2(X,i-2,j+5);
        n12 = load2(X,i-1,j+5);
        n22 = load2(X,i  ,j+5);
        n32 = load2(X,i+1,j+5);
        n42 = load2(X,i+2,j+5);
        n52 = load2(X,i+3,j+5);

//-----------------calc factor-------------------
        z10 = i64dilatation9_red( n13,n14,n10,
                    n23,n24,n20,
                    n33,n34,n30);

        z20 = i64dilatation9_red( n23,n24,n20,
                    n33,n34,n30,
                    n43,n44,n40);

        z11 = i64dilatation9_red( n14,n10,n11,
                    n24,n20,n21,
                    n34,n30,n31);

        z21 = i64dilatation9_red( n24,n20,n21,
                    n34,n30,n31,
                    n44,n40,n41);

        z12 = i64dilatation9_red( n10,n11,n12,
                    n20,n21,n22,
                    n30,n31,n32);

        z22 = i64dilatation9_red( n20,n21,n22,
                    n30,n31,n32,
                    n40,n41,n42);

        factor = i64erosion6_red(z10,z11,z12,z20,z21,z22);
//------------------3-0--------------------

        z00 = i64dilatation9_red( n03,n04,n00,
                    n13,n14,n10,
                    n23,n24,n20);

        z01 = i64dilatation9_red( n04,n00,n01,
                    n14,n10,n11,
                    n24,n20,n21);

        z02 = i64dilatation9_red( n00,n01,n02,
                    n10,n11,n12,
                    n20,n21,n22);

        max = i64erosion4(z00,z01,z02,factor);
       
        store2(Y,i,j+3,max);

//------------------3-1--------------------


        z20 = i64dilatation9_red( n33,n34,n30,
                    n43,n44,n40,
                    n53,n54,n50);

        z21 = i64dilatation9_red( n34,n30,n31,
                    n44,n40,n41,
                    n54,n50,n51);

        z22 = i64dilatation9_red( n30,n31,n32,
                    n40,n41,n42,
                    n50,n51,n52);
//--------------------------------------

        max = i64erosion4(z20,z21,z22,factor);
       
        store2(Y,i+1,j+3,max);

        n03 = load2(X,i-2,j+6);
        n13 = load2(X,i-1,j+6);
        n23 = load2(X,i  ,j+6);
        n33 = load2(X,i+1,j+6);
        n43 = load2(X,i+2,j+6);
        n53 = load2(X,i+3,j+6);

//-----------------calc factor-------------------
        z10 = i64dilatation9_red( n14,n10,n11,
                    n24,n20,n21,
                    n34,n30,n31);

        z20 = i64dilatation9_red( n24,n20,n21,
                    n34,n30,n31,
                    n44,n40,n41);

        z11 = i64dilatation9_red( n10,n11,n12,
                    n20,n21,n22,
                    n30,n31,n32);

        z21 = i64dilatation9_red( n20,n21,n22,
                    n30,n31,n32,
                    n40,n41,n42);

        z12 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z22 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);

        factor = i64erosion6_red(z10,z11,z12,z20,z21,z22);

//-----------------4-0--------------------

        z00 = i64dilatation9_red( n04,n00,n01,
                    n14,n10,n11,
                    n24,n20,n21);

//--------------------------------------

        z01 = i64dilatation9_red( n00,n01,n02,
                    n10,n11,n12,
                    n20,n21,n22);

//--------------------------------------

        z02 = i64dilatation9_red( n01,n02,n03,
                    n11,n12,n13,
                    n21,n22,n23);

//--------------------------------------

        max = i64erosion4(z00,z01,z02,factor);

        store2(Y,i,j+4,max);

    //-----------------4-1--------------------

        z20 = i64dilatation9_red( n34,n30,n31,
                    n44,n40,n41,
                    n54,n50,n51);

        z21 = i64dilatation9_red( n30,n31,n32,
                    n40,n41,n42,
                    n50,n51,n52);

        z22 = i64dilatation9_red( n31,n32,n33,
                    n41,n42,n43,
                    n51,n52,n53);
//--------------------------------------

        max = i64erosion4(z20,z21,z22,factor);

        store2(Y,i+1,j+4,max);
    }

    //epilogue
    switch (r)
    {
        case 1:
        n04 = load2(X,i-2,j+2);
        n14 = load2(X,i-1,j+2);
        n24 = load2(X,i  ,j+2);
        n34 = load2(X,i+1,j+2);
        n44 = load2(X,i+2,j+2);
        n54 = load2(X,i+3,j+2);

//-----------------calc factor-------------------
        z10 = i64dilatation9_red( n10,n11,n12,
                    n20,n21,n22,
                    n30,n31,n32);

        z20 = i64dilatation9_red( n20,n21,n22,
                    n30,n31,n32,
                    n40,n41,n42);

        z11 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z21 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);

        z12 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z22 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);

        factor = i64erosion6_red(z10,z11,z12,z20,z21,z22);

//-----------------0-0-------------------

        z00 = i64dilatation9_red( n00,n01,n02,
                    n10,n11,n12,
                    n20,n21,n22);

        z01 = i64dilatation9_red( n01,n02,n03,
                    n11,n12,n13,
                    n21,n22,n23);

        z02 = i64dilatation9_red( n02,n03,n04,
                    n12,n13,n14,
                    n22,n23,n24);

        max = i64erosion4(z00,z01,z02,factor);

        store2(Y,i,j,max);

//-------------------0-1-------------------

        z20 = i64dilatation9_red( n30,n31,n32,
                    n40,n41,n42,
                    n50,n51,n52);

        z21 = i64dilatation9_red( n31,n32,n33,
                    n41,n42,n43,
                    n51,n52,n53);

        z22 = i64dilatation9_red( n32,n33,n34,
                    n42,n43,n44,
                    n52,n53,n54);

        max = i64erosion4(z20,z21,z22,factor);

        store2(Y,i+1,j,max);

        break;

        case 2:
        n04 = load2(X,i-2,j+2);
        n14 = load2(X,i-1,j+2);
        n24 = load2(X,i  ,j+2);
        n34 = load2(X,i+1,j+2);
        n44 = load2(X,i+2,j+2);
        n54 = load2(X,i+3,j+2);

//-----------------calc factor-------------------
        z10 = i64dilatation9_red( n10,n11,n12,
                    n20,n21,n22,
                    n30,n31,n32);

        z20 = i64dilatation9_red( n20,n21,n22,
                    n30,n31,n32,
                    n40,n41,n42);

        z11 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z21 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);

        z12 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z22 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);

        factor = i64erosion6_red(z10,z11,z12,z20,z21,z22);

//-----------------0-0-------------------

        z00 = i64dilatation9_red( n00,n01,n02,
                    n10,n11,n12,
                    n20,n21,n22);

        z01 = i64dilatation9_red( n01,n02,n03,
                    n11,n12,n13,
                    n21,n22,n23);

        z02 = i64dilatation9_red( n02,n03,n04,
                    n12,n13,n14,
                    n22,n23,n24);

        max = i64erosion4(z00,z01,z02,factor);

        store2(Y,i,j,max);

//-------------------0-1-------------------

        z20 = i64dilatation9_red( n30,n31,n32,
                    n40,n41,n42,
                    n50,n51,n52);

        z21 = i64dilatation9_red( n31,n32,n33,
                    n41,n42,n43,
                    n51,n52,n53);

        z22 = i64dilatation9_red( n32,n33,n34,
                    n42,n43,n44,
                    n52,n53,n54);

        max = i64erosion4(z20,z21,z22,factor);

        store2(Y,i+1,j,max);

//-----------------------1-0--------------------
        n00 = load2(X,i-2,j+3);
        n10 = load2(X,i-1,j+3);
        n20 = load2(X,i  ,j+3);
        n30 = load2(X,i+1,j+3);
        n40 = load2(X,i+2,j+3);
        n50 = load2(X,i+3,j+3);
        
//-----------------calc factor-------------------
        z10 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z20 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);
                    
        z11 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z21 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);

        z12 = i64dilatation9_red( n13,n14,n10,
                    n23,n24,n20,
                    n33,n34,n30);

        z22 = i64dilatation9_red( n23,n24,n20,
                    n33,n34,n30,
                    n43,n44,n40);

        factor = i64erosion6_red(z10,z11,z12,z20,z21,z22);
//--------------------------------------

        z00 = i64dilatation9_red( n01,n02,n03,
                    n11,n12,n13,
                    n21,n22,n23);

        z01 = i64dilatation9_red( n02,n03,n04,
                    n12,n13,n14,
                    n22,n23,n24);

        z02 = i64dilatation9_red( n03,n04,n00,
                    n13,n14,n10,
                    n23,n24,n20);

        max = i64erosion4(z00,z01,z02,factor);

        store2(Y,i,j+1,max);

//-------------------1-1------------------

        z20 = i64dilatation9_red( n31,n32,n33,
                    n41,n42,n43,
                    n51,n52,n53);


        z21 = i64dilatation9_red( n32,n33,n34,
                    n42,n43,n44,
                    n52,n53,n54);


        z22 = i64dilatation9_red( n33,n34,n30,
                    n43,n44,n40,
                    n53,n54,n50);


        max = i64erosion4(z20,z21,z22,factor);


        store2(Y,i+1,j+1,max);

        break;
    
        case 3:
        n04 = load2(X,i-2,j+2);
        n14 = load2(X,i-1,j+2);
        n24 = load2(X,i  ,j+2);
        n34 = load2(X,i+1,j+2);
        n44 = load2(X,i+2,j+2);
        n54 = load2(X,i+3,j+2);

//-----------------calc factor-------------------
        z10 = i64dilatation9_red( n10,n11,n12,
                    n20,n21,n22,
                    n30,n31,n32);

        z20 = i64dilatation9_red( n20,n21,n22,
                    n30,n31,n32,
                    n40,n41,n42);

        z11 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z21 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);

        z12 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z22 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);

        factor = i64erosion6_red(z10,z11,z12,z20,z21,z22);

//-----------------0-0-------------------

        z00 = i64dilatation9_red( n00,n01,n02,
                    n10,n11,n12,
                    n20,n21,n22);

        z01 = i64dilatation9_red( n01,n02,n03,
                    n11,n12,n13,
                    n21,n22,n23);

        z02 = i64dilatation9_red( n02,n03,n04,
                    n12,n13,n14,
                    n22,n23,n24);

        max = i64erosion4(z00,z01,z02,factor);

        store2(Y,i,j,max);

//-------------------0-1-------------------

        z20 = i64dilatation9_red( n30,n31,n32,
                    n40,n41,n42,
                    n50,n51,n52);

        z21 = i64dilatation9_red( n31,n32,n33,
                    n41,n42,n43,
                    n51,n52,n53);

        z22 = i64dilatation9_red( n32,n33,n34,
                    n42,n43,n44,
                    n52,n53,n54);

        max = i64erosion4(z20,z21,z22,factor);

        store2(Y,i+1,j,max);

//-----------------------1-0--------------------
        n00 = load2(X,i-2,j+3);
        n10 = load2(X,i-1,j+3);
        n20 = load2(X,i  ,j+3);
        n30 = load2(X,i+1,j+3);
        n40 = load2(X,i+2,j+3);
        n50 = load2(X,i+3,j+3);
        
//-----------------calc factor-------------------
        z10 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z20 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);
                    
        z11 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z21 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);

        z12 = i64dilatation9_red( n13,n14,n10,
                    n23,n24,n20,
                    n33,n34,n30);

        z22 = i64dilatation9_red( n23,n24,n20,
                    n33,n34,n30,
                    n43,n44,n40);

        factor = i64erosion6_red(z10,z11,z12,z20,z21,z22);
//--------------------------------------

        z00 = i64dilatation9_red( n01,n02,n03,
                    n11,n12,n13,
                    n21,n22,n23);

        z01 = i64dilatation9_red( n02,n03,n04,
                    n12,n13,n14,
                    n22,n23,n24);

        z02 = i64dilatation9_red( n03,n04,n00,
                    n13,n14,n10,
                    n23,n24,n20);

        max = i64erosion4(z00,z01,z02,factor);

        store2(Y,i,j+1,max);

//-------------------1-1------------------

        z20 = i64dilatation9_red( n31,n32,n33,
                    n41,n42,n43,
                    n51,n52,n53);


        z21 = i64dilatation9_red( n32,n33,n34,
                    n42,n43,n44,
                    n52,n53,n54);


        z22 = i64dilatation9_red( n33,n34,n30,
                    n43,n44,n40,
                    n53,n54,n50);


        max = i64erosion4(z20,z21,z22,factor);


        store2(Y,i+1,j+1,max);


        n01 = load2(X,i-2,j+4);
        n11 = load2(X,i-1,j+4);
        n21 = load2(X,i  ,j+4);
        n31 = load2(X,i+1,j+4);
        n41 = load2(X,i+2,j+4);
        n51 = load2(X,i+3,j+4);
        
//-----------------calc factor-------------------
        z10 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z20 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);
        
        z11 = i64dilatation9_red( n13,n14,n10,
                    n23,n24,n20,
                    n33,n34,n30);

        z21 = i64dilatation9_red( n23,n24,n20,
                    n33,n34,n30,
                    n43,n44,n40);

        z12 = i64dilatation9_red( n14,n10,n11,
                    n24,n20,n21,
                    n34,n30,n31);

        z22 = i64dilatation9_red( n24,n20,n21,
                    n34,n30,n31,
                    n44,n40,n41);
        
        factor = i64erosion6_red(z10,z11,z12,z20,z21,z22);
//---------------2-0----------------------

        z00 = i64dilatation9_red( n02,n03,n04,
                    n12,n13,n14,
                    n22,n23,n24);

        z01 = i64dilatation9_red( n03,n04,n00,
                    n13,n14,n10,
                    n23,n24,n20);

        z02 = i64dilatation9_red( n04,n00,n01,
                    n14,n10,n11,
                    n24,n20,n21);

        max = i64erosion4(z00,z01,z02,factor);
        
        store2(Y,i,j+2,max);

//---------------2-1----------------------

        z20 = i64dilatation9_red( n32,n33,n34,
                    n42,n43,n44,
                    n52,n53,n54);

        z21 = i64dilatation9_red( n33,n34,n30,
                    n43,n44,n40,
                    n53,n54,n50);

        z22 = i64dilatation9_red( n34,n30,n31,
                    n44,n40,n41,
                    n54,n50,n51);
//--------------------------------------

        max = i64erosion4(z20,z21,z22,factor);
        
        store2(Y,i+1,j+2,max);

        break;

        case 4:
        n04 = load2(X,i-2,j+2);
        n14 = load2(X,i-1,j+2);
        n24 = load2(X,i  ,j+2);
        n34 = load2(X,i+1,j+2);
        n44 = load2(X,i+2,j+2);
        n54 = load2(X,i+3,j+2);

//-----------------calc factor-------------------
        z10 = i64dilatation9_red( n10,n11,n12,
                    n20,n21,n22,
                    n30,n31,n32);

        z20 = i64dilatation9_red( n20,n21,n22,
                    n30,n31,n32,
                    n40,n41,n42);

        z11 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z21 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);

        z12 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z22 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);

        factor = i64erosion6_red(z10,z11,z12,z20,z21,z22);

//-----------------0-0-------------------

        z00 = i64dilatation9_red( n00,n01,n02,
                    n10,n11,n12,
                    n20,n21,n22);

        z01 = i64dilatation9_red( n01,n02,n03,
                    n11,n12,n13,
                    n21,n22,n23);

        z02 = i64dilatation9_red( n02,n03,n04,
                    n12,n13,n14,
                    n22,n23,n24);

        max = i64erosion4(z00,z01,z02,factor);

        store2(Y,i,j,max);

//-------------------0-1-------------------

        z20 = i64dilatation9_red( n30,n31,n32,
                    n40,n41,n42,
                    n50,n51,n52);

        z21 = i64dilatation9_red( n31,n32,n33,
                    n41,n42,n43,
                    n51,n52,n53);

        z22 = i64dilatation9_red( n32,n33,n34,
                    n42,n43,n44,
                    n52,n53,n54);

        max = i64erosion4(z20,z21,z22,factor);

        store2(Y,i+1,j,max);

//-----------------------1-0--------------------
        n00 = load2(X,i-2,j+3);
        n10 = load2(X,i-1,j+3);
        n20 = load2(X,i  ,j+3);
        n30 = load2(X,i+1,j+3);
        n40 = load2(X,i+2,j+3);
        n50 = load2(X,i+3,j+3);
        
//-----------------calc factor-------------------
        z10 = i64dilatation9_red( n11,n12,n13,
                    n21,n22,n23,
                    n31,n32,n33);

        z20 = i64dilatation9_red( n21,n22,n23,
                    n31,n32,n33,
                    n41,n42,n43);
                    
        z11 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z21 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);

        z12 = i64dilatation9_red( n13,n14,n10,
                    n23,n24,n20,
                    n33,n34,n30);

        z22 = i64dilatation9_red( n23,n24,n20,
                    n33,n34,n30,
                    n43,n44,n40);

        factor = i64erosion6_red(z10,z11,z12,z20,z21,z22);
//--------------------------------------

        z00 = i64dilatation9_red( n01,n02,n03,
                    n11,n12,n13,
                    n21,n22,n23);

        z01 = i64dilatation9_red( n02,n03,n04,
                    n12,n13,n14,
                    n22,n23,n24);

        z02 = i64dilatation9_red( n03,n04,n00,
                    n13,n14,n10,
                    n23,n24,n20);

        max = i64erosion4(z00,z01,z02,factor);

        store2(Y,i,j+1,max);

//-------------------1-1------------------

        z20 = i64dilatation9_red( n31,n32,n33,
                    n41,n42,n43,
                    n51,n52,n53);


        z21 = i64dilatation9_red( n32,n33,n34,
                    n42,n43,n44,
                    n52,n53,n54);


        z22 = i64dilatation9_red( n33,n34,n30,
                    n43,n44,n40,
                    n53,n54,n50);


        max = i64erosion4(z20,z21,z22,factor);


        store2(Y,i+1,j+1,max);


        n01 = load2(X,i-2,j+4);
        n11 = load2(X,i-1,j+4);
        n21 = load2(X,i  ,j+4);
        n31 = load2(X,i+1,j+4);
        n41 = load2(X,i+2,j+4);
        n51 = load2(X,i+3,j+4);
        
//-----------------calc factor-------------------
        z10 = i64dilatation9_red( n12,n13,n14,
                    n22,n23,n24,
                    n32,n33,n34);

        z20 = i64dilatation9_red( n22,n23,n24,
                    n32,n33,n34,
                    n42,n43,n44);
        
        z11 = i64dilatation9_red( n13,n14,n10,
                    n23,n24,n20,
                    n33,n34,n30);

        z21 = i64dilatation9_red( n23,n24,n20,
                    n33,n34,n30,
                    n43,n44,n40);

        z12 = i64dilatation9_red( n14,n10,n11,
                    n24,n20,n21,
                    n34,n30,n31);

        z22 = i64dilatation9_red( n24,n20,n21,
                    n34,n30,n31,
                    n44,n40,n41);
        
        factor = i64erosion6_red(z10,z11,z12,z20,z21,z22);
//---------------2-0----------------------

        z00 = i64dilatation9_red( n02,n03,n04,
                    n12,n13,n14,
                    n22,n23,n24);

        z01 = i64dilatation9_red( n03,n04,n00,
                    n13,n14,n10,
                    n23,n24,n20);

        z02 = i64dilatation9_red( n04,n00,n01,
                    n14,n10,n11,
                    n24,n20,n21);

        max = i64erosion4(z00,z01,z02,factor);
        
        store2(Y,i,j+2,max);

//---------------2-1----------------------

        z20 = i64dilatation9_red( n32,n33,n34,
                    n42,n43,n44,
                    n52,n53,n54);

        z21 = i64dilatation9_red( n33,n34,n30,
                    n43,n44,n40,
                    n53,n54,n50);

        z22 = i64dilatation9_red( n34,n30,n31,
                    n44,n40,n41,
                    n54,n50,n51);
//--------------------------------------

        max = i64erosion4(z20,z21,z22,factor);
        
        store2(Y,i+1,j+2,max);

        n02 = load2(X,i-2,j+5);
        n12 = load2(X,i-1,j+5);
        n22 = load2(X,i  ,j+5);
        n32 = load2(X,i+1,j+5);
        n42 = load2(X,i+2,j+5);
        n52 = load2(X,i+3,j+5);

//-----------------calc factor-------------------
        z10 = i64dilatation9_red(  n13,n14,n10,
                                n23,n24,n20,
                                n33,n34,n30);

        z20 = i64dilatation9_red(  n23,n24,n20,
                                n33,n34,n30,
                                n43,n44,n40);

        z11 = i64dilatation9_red(  n14,n10,n11,
                                n24,n20,n21,
                                n34,n30,n31);

        z21 = i64dilatation9_red(  n24,n20,n21,
                                n34,n30,n31,
                                n44,n40,n41);

        z12 = i64dilatation9_red(  n10,n11,n12,
                                n20,n21,n22,
                                n30,n31,n32);

        z22 = i64dilatation9_red(  n20,n21,n22,
                                n30,n31,n32,
                                n40,n41,n42);
        factor = i64erosion6_red(z10,z11,z12,z20,z21,z22);
//------------------3-0--------------------

        z00 = i64dilatation9_red(  n03,n04,n00,
                                n13,n14,n10,
                                n23,n24,n20);

        z01 = i64dilatation9_red(  n04,n00,n01,
                                n14,n10,n11,
                                n24,n20,n21);

        z02 = i64dilatation9_red(  n00,n01,n02,
                                n10,n11,n12,
                                n20,n21,n22);

        max = i64erosion4(z00,z01,z02,factor);
       
        store2(Y,i,j+3,max);

//------------------3-1--------------------


        z20 = i64dilatation9_red(  n33,n34,n30,
                                n43,n44,n40,
                                n53,n54,n50);

        z21 = i64dilatation9_red(  n34,n30,n31,
                                n44,n40,n41,
                                n54,n50,n51);

        z22 = i64dilatation9_red(  n30,n31,n32,
                                n40,n41,n42,
                                n50,n51,n52);
//--------------------------------------

        max = i64erosion4(z20,z21,z22,factor);
       
        store2(Y,i+1,j+3,max);

        break;
    }
}
// -------------------------------------------------------------------------------------------------------
void fermeture3_ui64matrix_swp64_basic(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y, uint64 **Z)
// -------------------------------------------------------------------------------------------------------
{
    dilatation3_ui64matrix_swp64_basic(X, i0-1, i1+1, j0-1, j1+1, Y);
    erosion3_ui64matrix_swp64_basic(Y, i0,   i1,   j0,   j1,   Z);
}
// --------------------------------------------------------------------------------------------------
void fermeture3_ui64matrix_swp64_fusion_basic(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y)
// --------------------------------------------------------------------------------------------------
{
    for(int i = i0; i <= i1 ; i++ ){
        line_fermeture3_ui64matrix_swp64_fusion_basic(X, i, j0, j1, Y);
    }
}
// ------------------------------------------------------------------------------------------------
void fermeture3_ui64matrix_swp64_fusion_red(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y)
// ------------------------------------------------------------------------------------------------
{
    for(int i = i0; i <= i1 ; i++ ){
        line_fermeture3_ui64matrix_swp64_fusion_red(X, i, j0, j1, Y);
    }
}
// -----------------------------------------------------------------------------------------------------
void fermeture3_ui64matrix_swp64_fusion_ilu3_red(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y)
// -----------------------------------------------------------------------------------------------------
{
    for(int i = i0; i <= i1 ; i++ ){
        line_fermeture3_ui64matrix_swp64_fusion_ilu3_red(X, i, j0, j1, Y);
    }
}
// -----------------------------------------------------------------------------------------------------
void fermeture3_ui64matrix_swp64_fusion_elu2_red(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y)
// -----------------------------------------------------------------------------------------------------
{
    int i, r = (i1 - i0 + 1)%2;
    for(int i = i0; i <= i1-r ; i+=2 ){
        line_fermeture3_ui64matrix_swp64_fusion_elu2_red(X, i, j0, j1, Y);
    }
    if(r){
        line_fermeture3_ui64matrix_swp64_fusion_basic(X, i, j0, j1, Y);
    }
}
// ----------------------------------------------------------------------------------------------------------
void fermeture3_ui64matrix_swp64_fusion_ilu3_elu2_red(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y)
// ----------------------------------------------------------------------------------------------------------
{
    int i, r = (i1 - i0 + 1)%2;
    for(int i = i0; i <= i1-r ; i+=2 ){
        line_fermeture3_ui64matrix_swp64_fusion_ilu3_elu2_red(X, i, j0, j1, Y);
    }
    if(r){
        line_fermeture3_ui64matrix_swp64_fusion_basic(X, i, j0, j1, Y);
    }
}
// -----------------------------------------------------------------------------------------------------------------
void fermeture3_ui64matrix_swp64_fusion_ilu3_elu2_red_factor(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y)
// -----------------------------------------------------------------------------------------------------------------
{
    int i, r = (i1 - i0 + 1)%2;
    for(int i = i0; i <= i1-r ; i+=2 ){
        line_fermeture3_ui64matrix_swp64_fusion_ilu3_elu2_red_factor(X, i, j0, j1, Y);
    }
    if(r){
        line_fermeture3_ui64matrix_swp64_fusion_basic(X, i, j0, j1, Y);
    }
}
// ----------------------------------------------------------------------------------------------------------------
void fermeture3_ui64matrix_swp64_pipeline_basic(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y, uint64 **Z)
// ----------------------------------------------------------------------------------------------------------------
{
    //prologue pipeline
    line_dilatation3_ui64matrix_swp64_basic(X,i0-1,j0-1,j1+1,Y);
    line_dilatation3_ui64matrix_swp64_basic(X,i0  ,j0-1,j1+1,Y);

    //pipeline
    for(int i = i0; i <= i1-1; i++){
        line_dilatation3_ui64matrix_swp64_basic(X,i+1,j0-1,j1+1,Y);
        line_erosion3_ui64matrix_swp64_basic(Y,i  ,j0  ,j1  ,Z);
    }

    //epilogue pipeline
    line_erosion3_ui64matrix_swp64_basic(Y,i1  ,j0  ,j1  ,Z);
}
// --------------------------------------------------------------------------------------------------------------
void fermeture3_ui64matrix_swp64_pipeline_red(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y, uint64 **Z)
// --------------------------------------------------------------------------------------------------------------
{
    //prologue pipeline
    line_dilatation3_ui64matrix_swp64_red(X,i0-1,j0-1,j1+1,Y);
    line_dilatation3_ui64matrix_swp64_red(X,i0  ,j0-1,j1+1,Y);

    //pipeline
    for(int i = i0; i <= i1-1; i++){
        line_dilatation3_ui64matrix_swp64_red(X,i+1,j0-1,j1+1,Y);
        line_erosion3_ui64matrix_swp64_red(Y,i  ,j0  ,j1  ,Z);
    }

    //epilogue pipeline
    line_erosion3_ui64matrix_swp64_red(Y,i1  ,j0  ,j1  ,Z);
}
// -------------------------------------------------------------------------------------------------------------------
void fermeture3_ui64matrix_swp64_pipeline_ilu3_red(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y, uint64 **Z)
// -------------------------------------------------------------------------------------------------------------------
{
    //prologue pipeline
    line_dilatation3_ui64matrix_swp64_ilu3_red(X,i0-1,j0-1,j1+1,Y);
    line_dilatation3_ui64matrix_swp64_ilu3_red(X,i0  ,j0-1,j1+1,Y);

    //pipeline
    for(int i = i0; i <= i1-1; i++){
        line_dilatation3_ui64matrix_swp64_ilu3_red(X,i+1,j0-1,j1+1,Y);
        line_erosion3_ui64matrix_swp64_ilu3_red(Y,i  ,j0  ,j1  ,Z);
    }

    //epilogue pipeline
    line_erosion3_ui64matrix_swp64_ilu3_red(Y,i1  ,j0  ,j1  ,Z);
}
// -------------------------------------------------------------------------------------------------------------------
void fermeture3_ui64matrix_swp64_pipeline_elu2_red(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y, uint64 **Z)
// -------------------------------------------------------------------------------------------------------------------
{
    //prologue pipeline

    line_dilatation3_ui64matrix_swp64_elu2_red(X,i0-1,j0-1,j1+1,Y);

    //pipeline
    int r = (i1-i0+1)%2;
    for(int i = i0; i <= i1-r; i+=2){

        line_dilatation3_ui64matrix_swp64_elu2_red(X,i+1,j0-1,j1+1,Y);

        line_erosion3_ui64matrix_swp64_elu2_red(Y,i  ,j0  ,j1  ,Z);

    }

    //epilogue pipeline
        if(r){

            line_erosion3_ui64matrix_swp64_ilu3_red(Y,i1,j0  ,j1  ,Z);
        }
}
// --------------------------------------------------------------------------------------------------------------------------
void fermeture3_ui64matrix_swp64_pipeline_elu2_red_factor(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y, uint64 **Z)
// --------------------------------------------------------------------------------------------------------------------------
{
    //prologue pipeline
    line_dilatation3_ui64matrix_swp64_elu2_red_factor(X,i0-1,j0-1,j1+1,Y);

    //pipeline
    int r = (i1-i0+1)%2;
    for(int i = i0; i <= i1-r; i+=2){

        line_dilatation3_ui64matrix_swp64_elu2_red_factor(X,i+1,j0-1,j1+1,Y);
        line_erosion3_ui64matrix_swp64_elu2_red_factor(Y,i  ,j0  ,j1  ,Z);
    }

    //epilogue pipeline
        if(r){
            line_erosion3_ui64matrix_swp64_ilu3_red(Y,i1,j0  ,j1  ,Z);
        }
}
// ------------------------------------------------------------------------------------------------------------------------
void fermeture3_ui64matrix_swp64_pipeline_ilu3_elu2_red(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y, uint64 **Z)
// ------------------------------------------------------------------------------------------------------------------------
{
    //prologue pipeline
    line_dilatation3_ui64matrix_swp64_ilu3_elu2_red(X,i0-1,j0-1,j1+1,Y);

    //pipeline
    int r = (i1-i0+1)%2;
    for(int i = i0; i <= i1-r; i+=2){

        line_dilatation3_ui64matrix_swp64_ilu3_elu2_red(X,i+1,j0-1,j1+1,Y);
        line_erosion3_ui64matrix_swp64_ilu3_elu2_red(Y,i  ,j0  ,j1  ,Z);
    }

    //epilogue pipeline
        if(r){
            line_erosion3_ui64matrix_swp64_ilu3_red(Y,i1,j0  ,j1  ,Z);
        }
}
// -------------------------------------------------------------------------------------------------------------------------------
void fermeture3_ui64matrix_swp64_pipeline_ilu3_elu2_red_factor(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y, uint64 **Z)
// -------------------------------------------------------------------------------------------------------------------------------
{
    line_dilatation3_ui64matrix_swp64_ilu3_elu2_red_factor(X,i0-1,j0-1,j1+1,Y);

    //pipeline
    int r = (i1-i0+1)%2;
    for(int i = i0; i <= i1-r; i+=2){

        line_dilatation3_ui64matrix_swp64_ilu3_elu2_red_factor(X,i+1,j0-1,j1+1,Y);

        line_erosion3_ui64matrix_swp64_ilu3_elu2_red_factor(Y,i  ,j0  ,j1  ,Z);

    }

    //epilogue pipeline
        if(r){
            line_erosion3_ui64matrix_swp64_ilu3_red(Y,i1,j0  ,j1  ,Z);
        }
}
