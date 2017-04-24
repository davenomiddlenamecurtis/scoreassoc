#include "dvector.hpp"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

extern "C" {
extern void dsvdc( double **, int, int,
            double *, double *, double **, double **,
            double *, int, int *);

};

#ifndef MAX
#define MAX(a,b)  ( (a) > (b) ? (a) : (b) )
#define MIN(a,b)  ( (a) < (b) ? (a) : (b) )
#endif

int dsvd::dcmp(dv2d &A)
{
// C++ wrapper for dsvdc single value decomposition routine, which
// is modified from clinpack csvdc routine - Dave Curtis 1998-2002

    int info,i,j,high,wide;
    high=A.get_height();
    wide=A.get_width();
	double ctemp, **a, **x, **ut, **vt, *s, *e, *work;
	doublevec wwork(high), ss(MIN(wide,high+1)),ee(wide);
	dv2d aa(wide,MAX(high,wide)),tempUT(0,0);
//	aa is A in column-row format, zero-padded if necessary
    for (i=0;i<MAX(high,wide);++i)
      for (j=0;j<wide;++j)
        { aa[j][i]=i<high?A[i][j]:0; }
	a=aa.vec();
	work=wwork.vec(); s=ss.vec(); e=ee.vec();
	if (wide<=high)
	  ut=UT.vec();
	else
	  {
	  tempUT.resize(wide,wide);
	  ut=tempUT.vec();
	  // extra rows, have to copy into UT later
	  }
	vt=VT.vec();

dsvdc( a, high, wide, s, e, ut, vt, work, 21, &info );

if (info!=0) return 0; // failed
if (wide>high)
  for (i=0;i<high;++i)
    for (j=0;j<wide;++j)
        UT[j][i]=j<high?tempUT[j][i]:0;
/* set columns >= high to zero so as not to confuse testing routine
   they should not be used as undefined and all W[i][i], i>=high seem
   to be zero, though below in fact only does this for i>=high+1
*/
for (i=0;i<wide;++i)
    for (j=0;j<wide;++j)
        W[i][j]=(j==i && i<MIN(wide,high+1))?s[i]:0;
done=1;
return 1;
}

void dsvd::normal(int code)      //1 to subtract means, 2 to divide by sd
 {
 int i,j;
 for (i=0;i<U.wide;++i)
   mean[i]=sd[i]=0.0;
 for (i=0;i<U.high;++i)
    for (j=0;j<U.wide;++j)
      {
      double x;
      x=U[i][j];
      mean[j]+=x;
      sd[j]+=x*x;
      }
 for (i=0;i<U.wide;++i)
    {
    mean[i]/=U.high;
    sd[i]=sqrt(sd[i]/U.high-mean[i]*mean[i]);
    }
 if (code==1||code==2)
  for (i=0;i<U.high;++i)
    for (int j=0;j<U.wide;++j)
       {
       U[i][j]-=mean[j];
       if (code==2) U[i][j]/=sd[j];
       }
 }





