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

#ifdef __cplusplus
}
#endif

#endif // __MACRO_H__
