#include "dvector.hpp" 
#include <math.h>
#include <stdlib.h>

int error_flag;
#define error(x,y) (error_flag=x,fprintf(stderr,"Error %d: %s\n",x,y),exit(x))

#define LOWLIMIT -1.0e100

#define ITMAX 200
static float sqrarg;
#define SQR(a) (sqrarg=(a),sqrarg*sqrarg)

int linmin(double **pvar,doublevec &xi,int n,float *fret,float (*func)(void));
int powell(double **pvar,int n,float ftol,float *fret,float (*func)(void));
void mnbrak(float *ax,float *bx,float *cx,float *fa,float *fb,float *fc,float (*func)(float));
float brent(float ax,float bx,float cx,float (*f)(float),float tol,float *xmin);

int powell(double **pvar,int n,float ftol,float *fret,float (*func)(void))
// pvar points to variables
{
	int i,ibig,j,iter;
	float t,fptt,fp,del;

	dv2d xi(n,n);
	doublevec ptt(n),psave(n),xit(n);

	for (i=0;i<n;++i)
	{
		psave[i]=*pvar[i];
		for (j=0;j<n;++j)
			xi[i][j]=i==j; // unit matrix
	}

	*fret=func();
	for (iter=1;;iter++) {
		fp=(*fret);
		ibig=0;
		del=0.0;
		for (i=0;i<n;i++) {
			for (j=0;j<n;j++) xit[j]=xi[j][i];
			fptt=(*fret);
			if (!linmin(pvar,xit,n,fret,func)||error_flag) // error
				return 0;
			if (fabs(fptt-(*fret)) > del) {
				del=fabs(fptt-(*fret));
				ibig=i;
			}
		}
		if (2.0*fabs(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))) {
			return iter;
		}
		if (iter == ITMAX)
		{
			error(11,"routine POWELL");
			return 0;
		}
		for (j=0;j<n;j++) {
			ptt[j]=2.0* *pvar[j]-psave[j];
			xit[j]=*pvar[j]-psave[j];
			psave[j]=*pvar[j];
			*pvar[j]=ptt[j];
		}
		fptt=func();
		for (j=0;j<n;++j)
			*pvar[j]=psave[j];
		if (fptt < fp) {
			t=2.0*(fp-2.0*(*fret)+fptt)*SQR(fp-(*fret)-del)-del*SQR(fp-fptt);
			if (t < 0.0) {
				if (!linmin(pvar,xit,n,fret,func)||error_flag) //error
					return 0;
				for (j=0;j<n;j++) xi[j][ibig]=xit[j];
			}
		}
	}
}

#undef ITMAX
#undef SQR

#define TOL 2.0e-4

int ncom=0;	/* defining declarations */
double **pcom=0;
doublevec *xicom=0;
float (*nrfunc)();

int linmin(double **pvar,doublevec &xi,int n,float *fret,float (*func)(void))
{
	int j;
	float xx,xmin,fx,fb,fa,bx,ax;
	float f1dim(float);

	ncom=n;
	doublevec psave(n);
	xicom=new doublevec(n);
	pcom=new DBLPOINTER[n];
	if (!pcom || !xicom->size() || !psave.size())
	{
		error(9,"in linmin()");
		return 0;
	}
	nrfunc=func;
	for (j=0;j<n;j++) {
		psave[j]=*pvar[j];
		pcom[j]=pvar[j];
		(*xicom)[j]=xi[j];
	}
	ax=0.0;
	xx=1.0;
	bx=2.0;
	mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
	if (error_flag) return 0;
	*fret=brent(ax,xx,bx,f1dim,TOL,&xmin);
	if (error_flag) return 0;
	for (j=0;j<n;j++) {
		xi[j] *= xmin;
		*pvar[j]=psave[j] + xi[j];
	}
	delete xicom;
	delete pcom;
	return 1;
}

#undef TOL

#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

float brent(float ax,float bx,float cx,float (*f)(float),float tol,float *xmin)
{
	int iter;
	float a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	float e=0.0;

	a=((ax < cx) ? ax : cx);
	b=((ax > cx) ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=(*f)(x);
	for (iter=1;iter<=ITMAX;iter++) {
		if (error_flag) return 0.0;
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			*xmin=x;
			return fx;
		}
		if (fabs(e) > tol1) {
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			else {
				d=p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			}
		}
		else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu=(*f)(u);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			SHFT(v,w,x,u)
				SHFT(fv,fw,fx,fu)
		}
		else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			}
			else if (fu <= fv || v == x || v == w) {
				v=u;
				fv=fu;
			}
		}
	}
	error(11,"numerical routines. Too many iterations in BRENT");
	*xmin=x;
	return fx;
}

#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SIGN

extern int ncom;	/* defined in LINMIN */
extern double **pcom;
extern doublevec *xicom;

float f1dim(float x)
{
	int j;
	float f,*xt;

	doublevec pcomsave(ncom);
	for (j=0;j<ncom;j++)
	{
		pcomsave[j]=*pcom[j];
		*pcom[j]+=x*(*xicom)[j];
	}
	f=nrfunc();
	for (j=0;j<ncom;j++)
		*pcom[j]=pcomsave[j];
	return f;
}

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

void mnbrak(float *ax,float *bx,float *cx,float *fa,float *fb,float *fc,float (*func)(float))
{
	float ulim,u,r,q,fu,dum;

	*fa=(*func)(*ax);
	*fb=(*func)(*bx);
	if (*fb > *fa) {
		SHFT(dum,*ax,*bx,dum)
			SHFT(dum,*fb,*fa,dum)
	}
	*cx=(*bx)+GOLD*(*bx-*ax);
	*fc=(*func)(*cx);
	while (*fb > *fc) {
		if (error_flag) return;
		else if (*fb<LOWLIMIT ||
			*fa<LOWLIMIT ||
			*fc<LOWLIMIT)
		{
			error(11,"numerical routines. Underflow in MNBRAK");
			return;
		}

		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
			(2.0*SIGN(MAX(fabs(q-r),TINY),q-r));
		ulim=(*bx)+GLIMIT*(*cx-*bx);
		if ((*bx-u)*(u-*cx) > 0.0) {
			fu=(*func)(u);
			if (fu < *fc) {
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				return;
			}
			else if (fu > *fb) {
				*cx=u;
				*fc=fu;
				return;
			}
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
		}
		else if ((*cx-u)*(u-ulim) > 0.0) {
			fu=(*func)(u);
			if (fu < *fc) {
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
					SHFT(*fb,*fc,fu,(*func)(u))
			}
		}
		else if ((u-ulim)*(ulim-*cx) >= 0.0) {
			u=ulim;
			fu=(*func)(u);
		}
		else {
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
		}
		SHFT(*ax,*bx,*cx,u)
			SHFT(*fa,*fb,*fc,fu)
	}
}

#undef GOLD
#undef GLIMIT
#undef TINY
#undef MAX
#undef SIGN
#undef SHFT

