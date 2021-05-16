/* --------------------- */
/* --- macro.h --- */
/* --------------------- */

#ifndef __MACRO_H__
#define __MACRO_H__

#ifdef __cplusplus
extern "C" {
#endif

#define load2(T, i, j) T[i][j]

#define store(T,x) T=x
#define store2(T, i, j, x) T[i][j]=x

#define max(a, b) (a>b?a:b)
#define max3(a, b, c) max(max(a,b),c)
#define max9(a,b,c,d,e,f,g,h,i) max(max(max(max(max(max(max(max(a,b),c),d),e),f),g),h),i)
#define max_red(x, y)  (x^((x^y)&-(x<y)))
#define max3_red(a, b, c) max_red(max_red(a,b),c)
#define max9_red(a,b,c,d,e,f,g,h,i) max_red(max_red(max_red(max_red(max_red(max_red(max_red(max_red(a,b),c),d),e),f),g),h),i)
#define max9_red_mat33(x,i0,j0) max_red(max_red(max_red(max_red(max_red(max_red(max_red(max_red(x[i0][j0],x[i0][j0+1]),x[i0][j0+2]),x[i0+1][j0]),x[i0+1][j0+1]),x[i0+1][j0+2]),x[i0+2][j0]),x[i0+2][j0+1]),x[i0+2][j0+2])

#define min(a, b) (a<b?a:b)
#define min3(a, b, c) min(min(a,b),c)
#define min9(a,b,c,d,e,f,g,h,i) min(min(min(min(min(min(min(min(a,b),c),d),e),f),g),h),i)
#define min_red(x, y)  (y^((x^y)&-(x<y)))// reference
#define min3_red(a, b, c) min_red(min_red(a,b),c)
#define min9_red(a,b,c,d,e,f,g,h,i) min_red(min_red(min_red(min_red(min_red(min_red(min_red(min_red(a,b),c),d),e),f),g),h),i)
#define min9_red_mat33(x,i0,j0) min_red(min_red(min_red(min_red(min_red(min_red(min_red(min_red(x[i0][j0],x[i0][j0+1]),x[i0][j0+2]),x[i0+1][j0]),x[i0+1][j0+1]),x[i0+1][j0+2]),x[i0+2][j0]),x[i0+2][j0+1]),x[i0+2][j0+2])

#define add_red(x, y)  x+y
#define add3_red(a, b, c) a+b+c
#define add9_red(a,b,c,d,e,f,g,h,i) a+b+c+d+e+f+g+h+i
#define add9_mat33(x,i0,j0) x[i0][j0]+x[i0][j0+1]+x[i0][j0+2]+x[i0+1][j0]+x[i0+1][j0+1]+x[i0+1][j0+2]+x[i0+2][j0]+x[i0+2][j0+1]+x[i0+2][j0+2]

#define or(x0, x1) (x0)|(x1)
#define or3(x0, x1, x2) (x0)|(x1)|(x2)
#define or4(x0, x1, x2, x3) (x0)|(x1)|(x2)|(x3)
#define or5(x0, x1, x2, x3, x4) (x0)|(x1)|(x2)|(x3)|(x4)
#define or6(a,b,c,d,e,f) (a)|(b)|(c)|(d)|(e)|(f)
#define or9(a,b,c,d,e,f,g,h,i) (a)|(b)|(c)|(d)|(e)|(f)|(g)|(h)|(i)
#define or9_mat33(x,i0,j0) x[i0][j0]|x[i0][j0+1]|x[i0][j0+2]|x[i0+1][j0]|x[i0+1][j0+1]|x[i0+1][j0+2]|x[i0+2][j0]|x[i0+2][j0+1]|x[i0+2][j0+2]

#define and(x0, x1) (x0)&(x1)
#define and3(x0, x1, x2) (x0)&(x1)&(x2)
#define and5(x0, x1, x2, x3, x4) (x0)&(x1)&(x2)&(x3)&(x4)
#define and6(x0, x1, x2, x3, x4,x5) (x0)&(x1)&(x2)&(x3)&(x4)&(x5)
#define and9(a,b,c,d,e,f,g,h,i) (a)&(b)&(c)&(d)&(e)&(f)&(g)&(h)&(i)
#define and9_mat33(x,i0,j0) x[i0][j0]&x[i0][j0+1]&x[i0][j0+2]&x[i0+1][j0]&x[i0+1][j0+1]&x[i0+1][j0+2]&x[i0+2][j0]&x[i0+2][j0+1]&x[i0+2][j0+2]


#define i64right(a,b,n) (((a)>>n)|((b)<<(64-n)))
#define i64left(b,c,n) (((c)<<n)|((b)>>(64-n)))

#define i64erosion_mat33(X,i,j) and9(i64left(X[i-1][j-1],X[i-1][j],1),X[i-1][j],i64right(X[i-1][j],X[i-1][j+1],1),i64left(X[i  ][j-1],X[i  ][j],1),X[i  ][j],i64right(X[i  ][j],X[i  ][j+1],1),i64left(X[i+1][j-1],X[i+1][j],1),X[i+1][j],i64right(X[i+1][j],X[i+1][j+1],1))
#define i64dilatation_mat33(X,i,j) or9(i64left(X[i-1][j-1],X[i-1][j],1),X[i-1][j],i64right(X[i-1][j],X[i-1][j+1],1),i64left(X[i  ][j-1],X[i  ][j],1),X[i  ][j],i64right(X[i  ][j],X[i  ][j+1],1),i64left(X[i+1][j-1],X[i+1][j],1),X[i+1][j],i64right(X[i+1][j],X[i+1][j+1],1))


#define i64erosion3(x0,x1,x2) and3(i64left(x0,x1,1),x1,i64right(x1,x2,1))
#define i64dilatation3(x0,x1,x2) or3(i64left(x0,x1,1),x1,i64right(x1,x2,1))

#define i64erosion3c_left(x00,x01,x10,x11,x20,x21) and3(i64left(x00,x01,1),i64left(x10,x11,1),i64left(x20,x21,1))
#define i64dilatation3c_left(x00,x01,x10,x11,x20,x21) or3(i64left(x00,x01,1),i64left(x10,x11,1),i64left(x20,x21,1))

#define i64erosion3c_right(x00,x01,x10,x11,x20,x21) and3(i64right(x00,x01,1),i64right(x10,x11,1),i64right(x20,x21,1))
#define i64dilatation3c_right(x00,x01,x10,x11,x20,x21) or3(i64right(x00,x01,1),i64right(x10,x11,1),i64right(x20,x21,1))

#define i64erosion4(x0,x1,x2,factor) and(i64erosion3(x0,x1,x2),factor)
#define i64dilatation4(x0,x1,x2,factor) or(i64dilatation3(x0,x1,x2),factor)

#define i64erosion6(x0,x1,x2,x3,x4,x5) and(i64erosion3(x0,x1,x2),i64erosion3(x3,x4,x5))
#define i64dilatation6(x0,x1,x2,x3,x4,x5) or(i64dilatation3(x0,x1,x2),i64dilatation3(x3,x4,x5))

#define i64erosion6_red(x0,x1,x2,x3,x4,x5) i64erosion3(and(x0,x3),and(x1,x4),and(x2,x5))
#define i64dilatation6_red(x0,x1,x2,x3,x4,x5) i64dilatation3(or(x0,x3),or(x1,x4),or(x2,x5))

#define i64erosion9(x0,x1,x2,x3,x4,x5,x6,x7,x8) and9(i64left(x0,x1,1),x1,i64right(x1,x2,1),i64left(x3,x4,1),x4,i64right(x4,x5,1),i64left(x6,x7,1),x7,i64right(x7,x8,1))
#define i64dilatation9(x0,x1,x2,x3,x4,x5,x6,x7,x8) or9(i64left(x0,x1,1),x1,i64right(x1,x2,1),i64left(x3,x4,1),x4,i64right(x4,x5,1),i64left(x6,x7,1),x7,i64right(x7,x8,1))

#define i64erosion9_red(x0,x1,x2,x3,x4,x5,x6,x7,x8) i64erosion3(and3(x0,x3,x6),and3(x1,x4,x7),and3(x2,x5,x8))
#define i64dilatation9_red(x0,x1,x2,x3,x4,x5,x6,x7,x8) i64dilatation3(or3(x0,x3,x6),or3(x1,x4,x7),or3(x2,x5,x8))

#ifdef __cplusplus
}
#endif

#endif // __MACRO_H__
