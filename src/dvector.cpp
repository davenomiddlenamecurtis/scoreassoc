#include "dvector.hpp"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

// double precision vector, matrix arithmetic
// Dave Curtis, 1998

int default_dvec_error(char *s)
 {
 fprintf(stderr,"%s\n",s);
 return 1;
 }

int  (*dvec_error)(char*)=default_dvec_error;

dv2d dv2d::operator+(dv2d const & mat2)
// Matrix addition: SUM=A+B
  {
  dv2d sum(high,wide);
  for (int i=0;i<high;++i)
    for (int j=0;j<wide;++j)
      sum.rv[i][j]=rv[i][j]+mat2.rv[i][j];
  return sum;
  }

dv2d dv2d::operator*(dv2d const & mat2)
// Matrix multiplication: SUM=A*B
  {
  dv2d sum(high,mat2.wide);
  for (int i=0;i<sum.high;++i)
    {
    for (int j=0;j<sum.wide;++j)
     {
     double temp=0;
       {
       for (int k=0;k<wide;++k)
         temp+=rv[i][k]*mat2.rv[k][j];
       }
     sum.rv[i][j]=temp;
     }
    }
  return sum;
  }

doublevec::doublevec(int s)
// Constructor, s is number of elements
 {
 if ((v=new double[s])==0)
  {
  dvec_error("Out of memory in doublevec(2)");
  sz=0;
  }
 else sz=s;
 }

dv2d& dv2d::operator=(dv2d const & othermat)
// Matrix assignment: A=B
 {
 if (wide!=othermat.wide || high!=othermat.high)
    resize(othermat.high,othermat.wide);
 if (size())
    for (int i=0;i<high;++i)
      memmove(rv[i],othermat.rv[i],wide*sizeof(double));
 return *this;
 }

doublevec::~doublevec()
// Destructor
 {
 if (v) delete v;
 }

dv2d::dv2d(int h,int w)
// Constructor, h rows and w columns
 {
 high=wide=0;
 rv=0;
 resize(h,w);
 }

dv2d::dv2d(dv2d const &othermat)
// Copy constructor
 {
 high=wide=0;
 rv=0;
 if (resize(othermat.high,othermat.wide))
   operator=(othermat);
 }

dv2d::~dv2d()
// Destructor
 {
 if (rv)
   {
   for (int i=0;i<high;++i)
    delete rv[i];
   delete rv;
   }
 }

int dv2d::resize(int h,int w)
// Carefully allocate and deallocate memory as required
// It might be more efficient not to resize if more
// than enough memory is already allocated, but that
// would involve additional house-keeping
 {
 int i;
 if (h==high && w==wide) return 1;
 if (rv)
   {
   for (i=0;i<high;++i)
     delete rv[i];
   delete rv;
   }
 max_h=h;
 high=h;
 wide=w;
 if (h==0 && w==0) return 1; // a mistake for only one to be 0
 rv=new DBLPOINTER[max_h];
 if (rv==NULL)
  {
  dvec_error("Out of memory in dv2d (1)");
  high=wide=0;
  return 0;
  }
 for (i=0;i<h;++i)
  if ((rv[i]=new double[w])==NULL)
  {
  dvec_error("Out of memory in dv2d (2)");
  for (;i;--i) delete rv[i-1];
  delete rv;
  high=wide=max_h=0;
  rv=0;
  return 0;
  }
 return 1;
 }

int dsvd::OK(void)
// Check that memory is allocated appropriately by constructor etc.
{ return UT.size()
         && W.size()
         && W.size()==VT.size(); }

dsvd::dsvd(int hi,int wi)
 : UT(wi,hi), W(wi,wi), VT(wi,wi), mean(wi), sd(wi)
// Constructor, only allocate memory for transposed matrices for now
 { 
 done=0; 
 gotU=gotV=gotWinv=0; 
 }

int dsvd::setWinv(double lim)
// Set up inverse matrix of weights, zeroing inverses of small weights
// (as in Numerical Recipes in C, Press et al)
// Helps deal with ill-conditioned matrices, it says here
{
int i;
if (!done) return 0;
if (gotWinv) return 1;
if (!Winv.resize(W.wide,W.high)) return 0;
Winv=W;
double highest=0.0;
for (i=0;i<Winv.wide;++i)
  if (Winv[i][i]>highest) highest=Winv[i][i];
lim *=highest;
for (i=0;i<Winv.wide;++i)
  if (Winv[i][i]<lim) Winv[i][i]=0;
  else Winv[i][i]=1/Winv[i][i];
gotWinv=1;
return 1;
}

int dsvd::setU()
// If U is needed, get from UT (transposed U)
{
if (!done) return 0;
if (gotU) return 1;
if (!U.resize(UT.wide,UT.high)) return 0;
for (int i=0;i<U.high;++i)
  for (int j=0;j<U.wide;++j)
    U[i][j]=UT[j][i];
return 1;
}

int dsvd::setV()
// If V is needed, get from VT (transposed V)
{
if (!done) return 0;
if (gotV) return 1;
if (!V.resize(VT.wide,VT.high)) return 0; // (actually VT is square)
for (int i=0;i<V.high;++i)
  for (int j=0;j<V.wide;++j)
    V[i][j]=VT[j][i];
return 1;
}

void dsvd::zero(double lim)
// Set to zero all weights less than lim*highest weight
 {
 double highest=0.0;
 int i;
 for (i=0;i<W.wide;++i)
  if (W[i][i]>highest) highest=W[i][i];
 lim *=highest;
 for (i=0;i<W.wide;++i)
  if (W[i][i]<lim) W[i][i]=0;
 }

void dv2d::normal(int code)
// Normalise columns
// Might want to do this before principal components analysis
// 1 to subtract means, 2 to divide by sd
 {
 int i,j;
 doublevec mean(wide),sd(wide);
 for (i=0;i<wide;++i)
   mean[i]=sd[i]=0.0;
 for (i=0;i<high;++i)
    for (j=0;j<wide;++j)
      {
      double x;
      x=rv[i][j];
      mean[j]+=x;
      sd[j]+=x*x;
      }
 for (i=0;i<wide;++i)
    {
    mean[i]/=high;
    sd[i]=sqrt(sd[i]/high-mean[i]*mean[i]);
    }
 if (code==1||code==2)
  for (i=0;i<high;++i)
    for (j=0;j<wide;++j)
       {
       rv[i][j]-=mean[j];
       if (code==2) rv[i][j]/=sd[j];
       }
 }

dv2d dsvd::solve(dv2d const & b)
// Solve for X in AX=b, where A is original matrix before decomposition
 {
 dv2d temp(W.wide,b.wide);
 setV();
 setWinv();
 temp=V*Winv*UT*b;
 return temp;
 }

#if 0
dv2d dsvd::svbksb(dv2d& b)
// Solve for X in AX=b where A is original matrix before decomposition
// Adapted from Numerical recipes in C, Press et al
   {
   dv2d x(W.wide,1);  // as high as W is wide
   int jj,j,i;
   double s;
   setU(); setV(); // could rewrite this later so won't need them
   int n=U.wide,m=U.high;
   doublevec tmp(n+1);
   for (j=0;j<n;j++) {
      s=0.0;
      if (W[j][j]) {
         for (i=0;i<m;i++) s += U[i][j]*b[i][0];
         s /= W[j][j];
      }
      tmp[j]=s;
   }
   for (j=0;j<n;j++) {
      s=0.0;
      for (jj=0;jj<n;jj++) s += V[j][jj]*tmp[jj];
      x[j][0]=s;
   }
return x;
}
#endif

dv2d dv2d::inv()
// Return inverse of matrix, SUM=A.inv()
{
dv2d sum(1,1);
if (high==1 && wide==1) { sum[0][0]=1/rv[0][0]; return sum; }
dsvd s(high,wide);
if (s.dcmp(*this)==0) dvec_error("Matrix decomposition failed");
s.setWinv();
s.setV();
sum=s.V*s.Winv*s.UT;
return sum;
}

doublevec::doublevec(doublevec const &othervec)
// Copy constructor
 {
 if ((v=new double[othervec.sz])==0)
  {
  dvec_error("Out of memory in doublevec(1)");
  sz=0;
  }
 else
  {
  sz=othervec.sz;
  for (int i=0;i<sz;++i) v[i]=othervec.v[i];
  }
 }

doublevec::doublevec(void)
// Constructor of zero-size vector, to be resized later
 {
 sz=0;
 v=NULL;
 }

doublevec doublevec::operator+(doublevec const & v2)
// Vector addition: SUM=A+B
 {
 doublevec sum(sz);
 for (int i=0;i<sz;++i)
    sum.v[i]=v[i]+v2.v[i];
 return sum;
 }

doublevec doublevec::operator-(doublevec const & v2)
// Vector subtraction: SUM=A-B
 {
 doublevec sum(sz);
 for (int i=0;i<sz;++i)
    sum.v[i]=v[i]-v2.v[i];
 return sum;
 }


dv2d dv2d::transpose()
// Matrix transposition SUM=A.transpose()
 {
 dv2d trans(wide,high);
 if (trans.size())
  for (int i=0;i<high;++i)
   for (int j=0;j<wide;++j)
     trans[j][i]=rv[i][j];
 return trans;
 }

dv2d& dv2d::operator+=(dv2d const & mat2)
// Matrix addition: A=A+B, written A+=B
  {
  for (int i=0;i<high;++i)
    for (int j=0;j<wide;++j)
        rv[i][j]+=mat2.rv[i][j];
  return *this;
  }

dv2d dv2d::operator-(dv2d const & mat2)
// Matrix subtraction: A=A-B, written A-=B
  {
  dv2d sum(high,wide);
  for (int i=0;i<high;++i)
    for (int j=0;j<wide;++j)
        sum[i][j]=rv[i][j]-mat2.rv[i][j];
  return sum;
  }

dv2d& dv2d::operator*=(dv2d const & mat2)
// Matrix multiplication: A=A*B, written A*=B
{
*this=(*this)*mat2;  // no shortcut to multiplying
return *this;
}

dv2d dv2d::operator*(double f)
// Scale whole matrix by a constant: SUM=A*f
  {
  dv2d sum(high,wide);
  for (int i=0;i<high;++i)
    for (int j=0;j<wide;++j)
       sum.rv[i][j]=rv[i][j]*f;
  return sum;
  }

doublevec& doublevec::operator=(doublevec const & othervec)
// Vector assignment: A=B
 {
 if (sz!=othervec.sz)
   {
   if (v) delete v;
   if ((v=new double[othervec.sz])==0)
     {
     dvec_error("Out of memory in doublevec(3)");
     sz=0;
     return *this;
     }
   else sz=othervec.sz;
   }
 for (int i=0;i<sz;++i)
    v[i]=othervec.v[i];
 return *this;
 }

dv2d& dv2d::operator*=(double f)
// Scale whole matrix by a constant: A=A*f, written A*=f
{
int i,j;
for (i=0;i<high;++i)
   for (j=0;j<wide;++i)
        rv[i][j]*=f;
return *this;
}

void doublevec::fprintf(FILE *fp,char *form)
// Print out whole vector
// Each element formatted according to form, which can be omitted
 {
 for (int i=0;i<sz;++i)
   ::fprintf(fp,form,v[i]);
 }

void dv2d::fprintf(FILE *fp,char *form)
// Print out whole matrix
// Each element formatted according to form, which can be omitted
 {
 for (int i=0;i<high;++i)
   {
   for (int k=0;k<wide;++k) ::fprintf(fp,form,rv[i][k]);
   ::fprintf(fp,"\n");
   }
 }

void dv2d::fscanf(FILE *fp)
// Read in matrix from file
 {
 for (int i=0;i<high;++i)
   for (int k=0;k<wide;++k) ::fscanf(fp,"%lf",&rv[i][k]);
 }

dv2d dsvd::operator*(dv2d const & mat)
// Multiply by reconstituted matrix, just for checking really
// SUM=A*B where A was original matrix
 {
 setU();
 dv2d temp(U.high,mat.wide);
 temp=U*W*VT*mat;
 return temp;
 }

dsvd& dsvd::operator=(dsvd const & other)
// Assignment of dsvd object, doubt it's necessary
 {
 done=other.done;
 gotU=other.gotU;
 gotV=other.gotV;
 gotWinv=other.gotWinv;
 UT=other.UT;
 VT=other.VT;
 Winv=other.Winv;
 W=other.W;
 U=other.U;
 V=other.V;
 mean=other.mean;
 sd=other.sd;
 return *this;
 }

dsvd::dsvd(dsvd const & other)
 : UT(other.UT), VT(other.VT), W(other.W), U(other.U), V(other.V),
   sd(other.sd), mean(other.mean)
// Copy constructor, likewise probably unnecessary
 { 
 done=other.done; 
 gotU=gotV=gotWinv=0; 
 }
