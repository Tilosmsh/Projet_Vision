/* --------------------------------- */
/* --- morpho_dilatation_swp64.c --- */
/* --------------------------------- */

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
#include "morpho_dilatation_swp64.h"
/*
a1 a2 a3 a4  a5 a6 a7 a8  a9 a10 a11 a12
b1 b2 b3 b4  b5 b6 b7 b8  b9 b10 b11 b12
c1 c2 c3 c4  c5 c6 c7 c8  c9 c10 c11 c12
*/
// ----------------------------------------------------------------------------
void line_dilatation3_ui64matrix_swp64_basic(uint64 **X, int i, int j0, int j1, uint64 **Y)
// ----------------------------------------------------------------------------
{

    for(int j=j0; j<=j1; j++) {
        Y[i][j] = or9  (((X[i-1][j]<<1)|((X[i-1][j-1])>>63)), X[i-1][j], (X[i-1][j]>>1)|((X[i-1][j+1])<<63),
                        ((X[i  ][j]<<1)|((X[i  ][j-1])>>63)), X[i  ][j], (X[i  ][j]>>1)|((X[i  ][j+1])<<63),
                        ((X[i+1][j]<<1)|((X[i+1][j-1])>>63)), X[i+1][j], (X[i+1][j]>>1)|((X[i+1][j+1])<<63));
    }

}
// --------------------------------------------------------------------------
void line_dilatation3_ui64matrix_swp64_reg(uint64 **X, int i, int j0, int j1, uint64 **Y)
// --------------------------------------------------------------------------
{
    for(int j=j0; j<=j1; j++) {

        uint64   x0=load2(X,i-1,j-1), x1=load2(X,i-1,j), x2=load2(X,i-1,j+1),
                x3=load2(X,i  ,j-1), x4=load2(X,i  ,j), x5=load2(X,i  ,j+1),
                x6=load2(X,i+1,j-1), x7=load2(X,i+1,j), x8=load2(X,i+1,j+1);

        Y[i][j] = i64dilatation9(x0, x1, x2,
                                x3, x4, x5,
                                x6, x7, x8);
    }
}
// --------------------------------------------------------------------------
void line_dilatation3_ui64matrix_swp64_rot(uint64 **X, int i, int j0, int j1, uint64 **Y)
// --------------------------------------------------------------------------
{
    uint64   x0=load2(X,i-1,j0-1), x1=load2(X,i  ,j0-1), x2=load2(X,i+1,j0-1),
            res0 = i64dilatation3(x0, x1, x2),

            x3=load2(X,i-1,j0  ), x4=load2(X,i  ,j0  ), x5=load2(X,i+1,j0  ),
            res1 = i64dilatation3(x3, x4, x5),

            x6                  , x7                  , x8                  ,
            res2,

            res;

    for (int j=j0; j<=j1; j++) {
        x6=load2(X,i-1,j+1), x7=load2(X,i  ,j+1), x8=load2(X,i+1,j+1),
        res2 = i64dilatation3(x6, x7, x8);//load
        res=i64dilatation3(res0,res1,res2);//calc
        store2(Y, i, j, res);//store
        res0=res1;//rot res1->res0
        res1=res2;//rot res2->res1
    }
}
// --------------------------------------------------------------------------
void line_dilatation3_ui64matrix_swp64_red(uint64 **X, int i, int j0, int j1, uint64 **Y)
// --------------------------------------------------------------------------
{
    // uint64   x0=load2(X,i-1,j0-1), x1=load2(X,i  ,j0-1), x2=load2(X,i+1,j0-1),
    //         res0 = i64dilatation3(x0, x1, x2),

    //         x3=load2(X,i-1,j0  ), x4=load2(X,i  ,j0  ), x5=load2(X,i+1,j0  ),
    //         res1 = i64dilatation3(x3, x4, x5),

    //         x6                  , x7                  , x8                  ,
    //         res2,

    //         res;

    // for (int j=j0; j<=j1; j++) {
    //     x6=load2(X,i-1,j+1), x7=load2(X,i  ,j+1), x8=load2(X,i+1,j+1),
    //     res2 = i64dilatation3(x6, x7, x8);//load
    //     res=i64dilatation3(res0,res1,res2);//calc
    //     store2(Y, i, j, res);//store
    //     res0=res1;//rot res1->res0
    //     res1=res2;//rot res2->res1
    // }
    //v2
    // uint64  x00,x01,x02,
    //         x10,x11,x12,
    //         x20,x21,x22;
    // uint64  res0,res1,res2,res;

    // x00 = load2(X,i-1,j0-1);
    // x10 = load2(X,i  ,j0-1);
    // x20 = load2(X,i+1,j0-1);

    // x01 = load2(X,i-1,j0  );
    // x11 = load2(X,i  ,j0  );
    // x21 = load2(X,i+1,j0  );

    // res0 = or3(x00,x01,x02);
    // res1 = or3(x10,x11,x12);

    // for (int j=j0; j<=j1; j++) {
    //     x02=load2(X,i-1,j+1);
    //     x12=load2(X,i  ,j+1);
    //     x22=load2(X,i+1,j+1);
    //     res2 = or3(x20, x21, x22);//load
    //     res= i64dilatation3(res0,res1,res2);//calc
    //     store2(Y, i, j, res);//store
    //     res0=res1;//rot res1->res0
    //     res1=res2;//rot res2->res1
    // }

    //v3
    uint64  x00,x01,x02,
            x10,x11,x12,
            x20,x21,x22;

    uint64  res0,res1,res2,res;

    x00 = load2(X,i-1,j0-1); x01 = load2(X,i-1,j0  );
    x10 = load2(X,i  ,j0-1); x11 = load2(X,i  ,j0  );
    x20 = load2(X,i+1,j0-1); x21 = load2(X,i+1,j0  );
    
    res0 = or3(x00,x10,x20);
    res1 = or3(x01,x11,x21);

    for (int j=j0; j<=j1; j++) {

        x02=load2(X,i-1,j+1);
        x12=load2(X,i  ,j+1);
        x22=load2(X,i+1,j+1);

        res2 = or3(x02,x12,x22);
        res  = i64dilatation3(res0,res1,res2);//calc
        
        store2(Y, i, j, res);//store
        res0=res1;//rot res1->res0
        res1=res2;//rot res2->res1
    }
}
// ---------------------------------------------------------------------------
void line_dilatation3_ui64matrix_swp64_ilu3(uint64 **X, int i, int j0, int j1, uint64 **Y)
// ---------------------------------------------------------------------------
{
    uint64  a0,a1,a2,
            b0,b1,b2,
            c0,c1,c2,
            y;

    a0=load2(X, i-1, j0-1);b0=load2(X, i-1, j0);
    a1=load2(X, i  , j0-1);b1=load2(X, i  , j0);
    a2=load2(X, i+1, j0-1);b2=load2(X, i+1, j0);
    //9 regs pour la val cru, 3 regs pour le and colonne, 1 reg pour resultat

    //boucle
    int j,n=j1-j0+1,r=n%3;
    for(j=j0;j<j1;j=j+3){
        
        c0=load2(X,i-1,j+1);
        c1=load2(X,i  ,j+1);
        c2=load2(X,i+1,j+1);
        y=i64dilatation9_red(   a0,b0,c0,
                            a1,b1,c1,
                            a2,b2,c2);
        store2(Y,i,j  ,y);

        a0=load2(X,i-1,j+2);
        a1=load2(X,i  ,j+2);
        a2=load2(X,i+1,j+2);
        y=i64dilatation9_red(   b0,c0,a0,
                            b1,c1,a1,
                            b2,c2,a2);
        store2(Y,i,j+1,y);
        
        b0=load2(X,i-1,j+3);
        b1=load2(X,i  ,j+3);
        b2=load2(X,i+1,j+3);
        y=i64dilatation9_red(   c0,a0,b0,
                            c1,a1,b1,
                            c2,a2,b2);
        store2(Y,i,j+2,y);
    }

    //epilogue
    if(r==1){
        c0=load2(X,i-1,j+1);
        c1=load2(X,i  ,j+1);
        c2=load2(X,i+1,j+1);
        y=i64dilatation9_red(   a0,b0,c0,
                            a1,b1,c1,
                            a2,b2,c2);
        store2(Y,i,j  ,y);
    }
    else if(r==2){
        c0=load2(X,i-1,j+1);
        c1=load2(X,i  ,j+1);
        c2=load2(X,i+1,j+1);
        y=i64dilatation9_red(   a0,b0,c0,
                            a1,b1,c1,
                            a2,b2,c2);
        store2(Y,i,j  ,y);

        a0=load2(X,i-1,j+2);
        a1=load2(X,i  ,j+2);
        a2=load2(X,i+1,j+2);
        y=i64dilatation9_red(   b0,c0,a0,
                            b1,c1,a1,
                            b2,c2,a2);
        store2(Y,i,j+1,y);
    }
}
// -------------------------------------------------------------------------------
void line_dilatation3_ui64matrix_swp64_ilu3_red(uint64 **X, int i, int j0, int j1, uint64 **Y)
// -------------------------------------------------------------------------------
{
    uint64  a0,a1,a2,
            b0,b1,b2,
            c0,c1,c2,
            y;

    a0=load2(X, i-1, j0-1);b0=load2(X, i-1, j0);
    a1=load2(X, i  , j0-1);b1=load2(X, i  , j0);
    a2=load2(X, i+1, j0-1);b2=load2(X, i+1, j0);
    //9 regs pour la val cru, 3 regs pour le and colonne, 1 reg pour resultat

    //boucle
    int j,n=j1-j0+1,r=n%3;
    for(j=j0;j<j1;j=j+3){
        
        c0=load2(X,i-1,j+1);
        c1=load2(X,i  ,j+1);
        c2=load2(X,i+1,j+1);
        y=i64dilatation9_red(   a0,b0,c0,
                            a1,b1,c1,
                            a2,b2,c2);
        store2(Y,i,j  ,y);

        a0=load2(X,i-1,j+2);
        a1=load2(X,i  ,j+2);
        a2=load2(X,i+1,j+2);
        y=i64dilatation9_red(   b0,c0,a0,
                            b1,c1,a1,
                            b2,c2,a2);
        store2(Y,i,j+1,y);
        
        b0=load2(X,i-1,j+3);
        b1=load2(X,i  ,j+3);
        b2=load2(X,i+1,j+3);
        y=i64dilatation9_red(   c0,a0,b0,
                            c1,a1,b1,
                            c2,a2,b2);
        store2(Y,i,j+2,y);
    }

    //epilogue
    if(r==1){
        c0=load2(X,i-1,j+1);
        c1=load2(X,i  ,j+1);
        c2=load2(X,i+1,j+1);
        y=i64dilatation9_red(   a0,b0,c0,
                            a1,b1,c1,
                            a2,b2,c2);
        store2(Y,i,j  ,y);
    }
    else if(r==2){
        c0=load2(X,i-1,j+1);
        c1=load2(X,i  ,j+1);
        c2=load2(X,i+1,j+1);
        y=i64dilatation9_red(   a0,b0,c0,
                            a1,b1,c1,
                            a2,b2,c2);
        store2(Y,i,j  ,y);

        a0=load2(X,i-1,j+2);
        a1=load2(X,i  ,j+2);
        a2=load2(X,i+1,j+2);
        y=i64dilatation9_red(   b0,c0,a0,
                            b1,c1,a1,
                            b2,c2,a2);
        store2(Y,i,j+1,y);
    }
}
// -------------------------------------------------------------------------------
void line_dilatation3_ui64matrix_swp64_elu2_red(uint64 **X, int i, int j0, int j1, uint64 **Y)
// -------------------------------------------------------------------------------
{
    uint64 x0, x1 , x2,
          x3, x4 , x5,
          x6, x7 , x8,
          x9, x10, x11;
    uint64 res;

        x0=load2(X,i-1,j0-1), x1 =load2(X,i-1,j0),
        x3=load2(X,i  ,j0-1), x4 =load2(X,i  ,j0), 
        x6=load2(X,i+1,j0-1), x7 =load2(X,i+1,j0),
        x9=load2(X,i+2,j0-1), x10=load2(X,i+2,j0);

        
    for (int j=j0; j<=j1; j++) {

        x2=load2(X,i-1,j+1), x5=load2(X,i  ,j+1), x8=load2(X,i+1,j+1), x11=load2(X,i+2,j+1);
        
        res = i64dilatation9_red(   x0,x1,x2,
                                x3,x4,x5,
                                x6,x7,x8);
        store2(Y,i,j,res);

        res = i64dilatation9_red(  x3,x4 , x5,
                               x6,x7 , x8,
                               x9,x10,x11);
        store2(Y,i+1,j,res);

        x0=x1 ;x1=x2;
        x3=x4 ;x4=x5;
        x6=x7 ;x7=x8;
        x9=x10;x10=x11;
    }
}
// --------------------------------------------------------------------------------------
void line_dilatation3_ui64matrix_swp64_elu2_red_factor(uint64 **X, int i, int j0, int j1, uint64 **Y)
// --------------------------------------------------------------------------------------
{
    uint64 x0, x1 , x2,
          x3, x4 , x5,
          x6, x7 , x8,
          x9, x10, x11;
    uint64 res;

    uint64 cse;

        x0=load2(X,i-1,j0-1), x1 =load2(X,i-1,j0),
        x3=load2(X,i  ,j0-1), x4 =load2(X,i  ,j0), 
        x6=load2(X,i+1,j0-1), x7 =load2(X,i+1,j0),
        x9=load2(X,i+2,j0-1), x10=load2(X,i+2,j0);

        
    for (int j=j0; j<=j1; j++) {
        x2=load2(X,i-1,j+1), x5=load2(X,i  ,j+1), x8=load2(X,i+1,j+1), x11=load2(X,i+2,j+1);

        cse=or(i64dilatation3(x3,x4,x5),i64dilatation3(x6,x7,x8));

        res = or(i64dilatation3(x0,x1,x2),cse);
        store2(Y,i,j,res);

        res = or(i64dilatation3(x9,x10,x11),cse);
        store2(Y,i+1,j,res);

        x0=x1 ;x1=x2;
        x3=x4 ;x4=x5;
        x6=x7 ;x7=x8;
        x9=x10;x10=x11;
    }
}
// ------------------------------------------------------------------------------------
void line_dilatation3_ui64matrix_swp64_ilu3_elu2_red(uint64 **X, int i, int j0, int j1, uint64 **Y)
// ------------------------------------------------------------------------------------
{
    uint64  a0,b0,c0,
            a1,b1,c1,
            a2,b2,c2,
            a3,b3,c3,// 12 regs
            y0,y1;// 2 regs pour les resultats

    //init des regs pour le 1re tour de boucle
    a0=load2(X, i-1, j0-1);
    a1=load2(X, i  , j0-1);
    a2=load2(X, i+1, j0-1);
    a3=load2(X, i+2, j0-1);

    b0=load2(X, i-1, j0);
    b1=load2(X, i  , j0);
    b2=load2(X, i+1, j0);
    b3=load2(X, i+2, j0);


    //boucle
    int j,n=j1-j0+1,r=n%3;
    for(j=j0;j<j1-r;j=j+3){
        //load colonne c
        c0=load2(X,i-1,j+1);
        c1=load2(X,i  ,j+1);
        c2=load2(X,i+1,j+1);
        c3=load2(X,i+2,j+1);
        //calc & store(i,j)
        y0=i64dilatation9_red(  a0,b0,c0,
                                a1,b1,c1,
                                a2,b2,c2);
        store2(Y,i,j,y0);  
        //calc & store (i+1,j)
        y1=i64dilatation9_red(  a1,b1,c1,
                                a2,b2,c2,
                                a3,b3,c3);
        store2(Y,i+1,j,y1);

        //load colonne a
        a0=load2(X,i-1,j+2);
        a1=load2(X,i  ,j+2);
        a2=load2(X,i+1,j+2);
        a3=load2(X,i+2,j+2);
        //calc & store (i,j+1)
        y0=i64dilatation9_red(  b0,c0,a0,
                                b1,c1,a1,
                                b2,c2,a2);
        store2(Y,i,j+1,y0);
        //calc & store (i+1,j+1)
        y1=i64dilatation9_red(  b1,c1,a1,
                                b2,c2,a2,
                                b3,c3,a3);
        store2(Y,i+1,j+1,y1);
        
        //load colonne b
        b0=load2(X,i-1,j+3);
        b1=load2(X,i  ,j+3);
        b2=load2(X,i+1,j+3);
        b3=load2(X,i+2,j+3);
        //calc & store (i,j+2)
        y0=i64dilatation9_red(  c0,a0,b0,
                                c1,a1,b1,
                                c2,a2,b2);
        store2(Y,i,j+2,y0); 
        //calc & store (i+1,j+2)
        y1=i64dilatation9_red(  c1,a1,b1,
                                c2,a2,b2,
                                c3,a3,b3);
        store2(Y,i+1,j+2,y1);  
    }

    //epilogue
    if(r==1){
        //load colonne c
        c0=load2(X,i-1,j+1);
        c1=load2(X,i  ,j+1);
        c2=load2(X,i+1,j+1);
        c3=load2(X,i+2,j+1);
        //calc & store(i,j)
        y0=i64dilatation9_red(  a0,b0,c0,
                                a1,b1,c1,
                                a2,b2,c2);
        store2(Y,i,j,y0);  
        //calc & store (i+1,j)
        y1=i64dilatation9_red(  a1,b1,c1,
                                a2,b2,c2,
                                a3,b3,c3);
        store2(Y,i+1,j,y1);
    }
    else if(r==2){
        //load colonne c
        c0=load2(X,i-1,j+1);
        c1=load2(X,i  ,j+1);
        c2=load2(X,i+1,j+1);
        c3=load2(X,i+2,j+1);
        //calc & store(i,j)
        y0=i64dilatation9_red(  a0,b0,c0,
                                a1,b1,c1,
                                a2,b2,c2);
        store2(Y,i,j,y0);  
        //calc & store (i+1,j)
        y1=i64dilatation9_red(  a1,b1,c1,
                                a2,b2,c2,
                                a3,b3,c3);
        store2(Y,i+1,j,y1);

        //load colonne a
        a0=load2(X,i-1,j+2);
        a1=load2(X,i  ,j+2);
        a2=load2(X,i+1,j+2);
        a3=load2(X,i+2,j+2);
        //calc & store (i,j+1)
        y0=i64dilatation9_red(  b0,c0,a0,
                                b1,c1,a1,
                                b2,c2,a2);
        store2(Y,i,j+1,y0);
        //calc & store (i+1,j+1)
        y1=i64dilatation9_red(  b1,c1,a1,
                                b2,c2,a2,
                                b3,c3,a3);
        store2(Y,i+1,j+1,y1);
    }
}
// -------------------------------------------------------------------------------------------
void line_dilatation3_ui64matrix_swp64_ilu3_elu2_red_factor(uint64 **X, int i, int j0, int j1, uint64 **Y)
// -------------------------------------------------------------------------------------------
{
    uint64  a0,b0,c0,
            a1,b1,c1,
            a2,b2,c2,
            a3,b3,c3,// 12 regs
            factor,
            y0,y1;// 2 regs pour les resultats

    a0=load2(X, i-1, j0-1);
    a1=load2(X, i  , j0-1);
    a2=load2(X, i+1, j0-1);
    a3=load2(X, i+2, j0-1);

    b0=load2(X, i-1, j0);
    b1=load2(X, i  , j0);
    b2=load2(X, i+1, j0);
    b3=load2(X, i+2, j0);

    //boucle
    int j,n=j1-j0+1,r=n%3;
    for(j=j0;j<j1-r;j=j+3){
        
        //load colonne c
        c0=load2(X,i-1,j+1);
        c1=load2(X,i  ,j+1);
        c2=load2(X,i+1,j+1);
        c3=load2(X,i+2,j+1);
        //calc factor
        factor = i64dilatation3(or(a1,a2),or(b1,b2),or(c1,c2));
        //calc & store(i,j)
        y0=or(factor,i64dilatation3(a0,b0,c0));
        store2(Y,i,j,y0);  
        //calc & store (i+1,j)
        y1=or(factor,i64dilatation3(a3,b3,c3));
        store2(Y,i+1,j,y1); 

        //load colonne a
        a0=load2(X,i-1,j+2);
        a1=load2(X,i  ,j+2);
        a2=load2(X,i+1,j+2);
        a3=load2(X,i+2,j+2);
        //calc factor
        factor = i64dilatation3(or(b1,b2),or(c1,c2),or(a1,a2));
        //calc & store (i,j+1)
        y0=or(factor,i64dilatation3(b0,c0,a0));
        store2(Y,i,j+1,y0); 
        //calc & store (i+1,j+1)
        y1=or(factor,i64dilatation3(b3,c3,a3));
        store2(Y,i+1,j+1,y1);
        
        //load colonne b
        b0=load2(X,i-1,j+3);
        b1=load2(X,i  ,j+3);
        b2=load2(X,i+1,j+3);
        b3=load2(X,i+2,j+3);
        //calc factor
        factor = i64dilatation3(or(c1,c2),or(a1,a2),or(b1,b2));
        //calc & store (i,j+2)
        y0=or(factor,i64dilatation3(c0,a0,b0));
        store2(Y,i,j+2,y0); 
        //calc & store (i+1,j+2)
        y1=or(factor,i64dilatation3(c3,a3,b3));
        store2(Y,i+1,j+2,y1);
    }

    //epilogue
    if(r==1){
        //load colonne c
        c0=load2(X,i-1,j+1);
        c1=load2(X,i  ,j+1);
        c2=load2(X,i+1,j+1);
        c3=load2(X,i+2,j+1);
        //calc factor
        factor = i64dilatation3(or(a1,a2),or(b1,b2),or(c1,c2));
        //calc & store(i,j)
        y0=or(factor,i64dilatation3(a0,b0,c0));
        store2(Y,i,j,y0);  
        //calc & store (i+1,j)
        y1=or(factor,i64dilatation3(a3,b3,c3));
        store2(Y,i+1,j,y1); 
    }
    else if(r==2){
        //load colonne c
        c0=load2(X,i-1,j+1);
        c1=load2(X,i  ,j+1);
        c2=load2(X,i+1,j+1);
        c3=load2(X,i+2,j+1);
        //calc factor
        factor = i64dilatation3(or(a1,a2),or(b1,b2),or(c1,c2));
        //calc & store(i,j)
        y0=or(factor,i64dilatation3(a0,b0,c0));
        store2(Y,i,j,y0);  
        //calc & store (i+1,j)
        y1=or(factor,i64dilatation3(a3,b3,c3));
        store2(Y,i+1,j,y1); 

        //load colonne a
        a0=load2(X,i-1,j+2);
        a1=load2(X,i  ,j+2);
        a2=load2(X,i+1,j+2);
        a3=load2(X,i+2,j+2);
        //calc factor
        factor = i64dilatation3(or(b1,b2),or(c1,c2),or(a1,a2));
        //calc & store (i,j+1)
        y0=or(factor,i64dilatation3(b0,c0,a0));
        store2(Y,i,j+1,y0); 
        //calc & store (i+1,j+1)
        y1=or(factor,i64dilatation3(b3,c3,a3));
        store2(Y,i+1,j+1,y1);
    }
}
// --------------------------------------------------------------------------------
void dilatation3_ui64matrix_swp64_basic(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y)
// --------------------------------------------------------------------------------
{
    for(int i=i0;i<=i1;i++){
        line_dilatation3_ui64matrix_swp64_basic(X,i,j0,j1,Y);
    }
}
// ------------------------------------------------------------------------------
void dilatation3_ui64matrix_swp64_reg(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y)
// ------------------------------------------------------------------------------
{
    for(int i=i0;i<=i1;i++){
        line_dilatation3_ui64matrix_swp64_reg(X,i,j0,j1,Y);
    }
}
// ------------------------------------------------------------------------------
void dilatation3_ui64matrix_swp64_rot(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y)
// ------------------------------------------------------------------------------
{
    for(int i=i0;i<=i1;i++){
        line_dilatation3_ui64matrix_swp64_rot(X,i,j0,j1,Y);
    }
}
// ------------------------------------------------------------------------------
void dilatation3_ui64matrix_swp64_red(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y)
// ------------------------------------------------------------------------------
{
    for(int i=i0;i<=i1;i++){
        line_dilatation3_ui64matrix_swp64_red(X,i,j0,j1,Y);
    }
}
// -------------------------------------------------------------------------------
void dilatation3_ui64matrix_swp64_ilu3(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y)
// -------------------------------------------------------------------------------
{
    for(int i=i0;i<=i1;i++){
        line_dilatation3_ui64matrix_swp64_ilu3(X,i,j0,j1,Y);
    }
}
// -----------------------------------------------------------------------------------
void dilatation3_ui64matrix_swp64_ilu3_red(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y)
// -----------------------------------------------------------------------------------
{
    for(int i=i0;i<=i1;i++){
        line_dilatation3_ui64matrix_swp64_ilu3_red(X,i,j0,j1,Y);
    }
}
// -----------------------------------------------------------------------------------
void dilatation3_ui64matrix_swp64_elu2_red(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y)
// -----------------------------------------------------------------------------------
{
    int r = (i1-i0+1)%2,i;
    for(i=i0;i<i1;i+=2){
        line_dilatation3_ui64matrix_swp64_elu2_red(X,i,j0,j1,Y);
    }
    if(r==1){
        line_dilatation3_ui64matrix_swp64_basic(X,i,j0,j1,Y);
    }
}
// ------------------------------------------------------------------------------------------
void dilatation3_ui64matrix_swp64_elu2_red_factor(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y)
// ------------------------------------------------------------------------------------------
{
    int r = (i1-i0+1)%2,i;
    for(i=i0;i<i1;i+=2){
        line_dilatation3_ui64matrix_swp64_elu2_red_factor(X,i,j0,j1,Y);
    }
    if(r==1){
        line_dilatation3_ui64matrix_swp64_basic(X,i,j0,j1,Y);
    }
}
// ----------------------------------------------------------------------------------------
void dilatation3_ui64matrix_swp64_ilu3_elu2_red(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y)
// ----------------------------------------------------------------------------------------
{
    int r = (i1-i0+1)%2,i;
    for(i=i0;i<i1;i+=2){
        line_dilatation3_ui64matrix_swp64_ilu3_elu2_red(X,i  ,j0,j1,Y);
    }
    if(r==1){
        line_dilatation3_ui64matrix_swp64_basic(X,i,j0,j1,Y);
    }
}
// -----------------------------------------------------------------------------------------------
void dilatation3_ui64matrix_swp64_ilu3_elu2_red_factor(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y)
// -----------------------------------------------------------------------------------------------
{
    int r = (i1-i0+1)%2,i;
    for(i=i0;i<i1;i+=2){
        line_dilatation3_ui64matrix_swp64_ilu3_elu2_red_factor(X,i  ,j0,j1,Y);
    }
    if(r==1){
        line_dilatation3_ui64matrix_swp64_basic(X,i,j0,j1,Y);
    }
}