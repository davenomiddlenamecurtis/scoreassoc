#include <stdlib.h>
#include <math.h>
#include <dlib/optimization.h>
#include "dcerror.hpp"
#include "lrModel.hpp"

using namespace dlib;
typedef matrix<double,0,1> column_vector;
double derivative_eps = 1e-7; // used to get gradient of lnL by beta
double minimumP = 1e-8;
double betaLimit = 20,tLimit=20;

class lrModelMaximiser
{
public:
	double operator() (const column_vector& arg) const;
	double (*func)(void);
	lrModelMaximiser(double (*f)(void)) { func = f; }
};

// looks like I will have to use globals and assume that only one object is trying to use minimisation function

lrModel *modelToFit;
double *betasToFit;
int nBetasToFit;

double lrModelMaximiser::operator() (const column_vector& c) const
{
	int i;
	double f;
	for (i = 0; i<c.nr(); ++i)
		betasToFit[i] = c(i, 0);
	f = func();
	return f;
}


double getMinusModelLnL()
{
	int c,cc;
	double lnL;
	for (c=cc=0;c<modelToFit->nCol+1;++c)
		if (modelToFit->toFit[c])
		{
			modelToFit->beta[c]=betasToFit[cc];
			++cc;
		}
	lnL = modelToFit->getLnL();
	return -lnL;	
}

double lrModel::maximiseLnL()
{
	int c,cc;
	double d;
	float ftol; // need to find an appropriate value for this
	ftol = 0.01;
	if (betasToFit != 0)
		free(betasToFit);
	betasToFit = (double *)calloc(nCol + 1,sizeof(double));
	nBetasToFit = 0;
	for (c = 0; c < nCol + 1; ++c)
		if (toFit[c])
			++nBetasToFit;
	column_vector starting_point(nBetasToFit);
	for (c = cc = 0; c<nCol + 1; ++c)
		if (toFit[c])
		{
			starting_point(cc,0)=beta[c] ;
			++cc;
		}
	modelToFit = this;
	d = find_min_using_approximate_derivatives(cg_search_strategy(),
		objective_delta_stop_strategy(ftol),
		lrModelMaximiser(getMinusModelLnL), starting_point, -20, derivative_eps); // 0.001 is derivative_eps
	// if this fails an exception gets thrown
	for (c = cc = 0; c<nCol + 1; ++c)
		if (toFit[c])
		{
			beta[c]=starting_point(cc,0);
			++cc;
		}
	return -d;
}

double lrModel::getLnL()
{
	double lnL,eT,p;
	int r, c;
	lnL = 0;
#if 0
	for (c = 0; c < nCol; ++c)
		if (toUse[c])
		{
			if (beta[c] > betaLimit)
				beta[c] = betaLimit;
			if (beta[c] < -betaLimit)
				beta[c] = -betaLimit;
		}
#endif
	for (r = 0; r < nRow; ++r)
	{
		t[r] = 0;
		for (c = 0; c < nCol; ++c)
			if (toUse[c])
				t[r] += beta[c] * X[r][c];
		if (toUse[nCol])
			t[r] += beta[nCol];
		if (t[r] > tLimit)
			t[r] = tLimit;
		if (t[r] < -tLimit)
			t[r] = -tLimit;
		eT = exp(t[r]); // maximisation routine can set very large negative beta, so this can be 0
		p = eT / (eT + 1);
		if (p < minimumP)
			p = minimumP;
		if (p > 1 - minimumP)
			p = 1-minimumP;
		lnL += log(Y[r] == 1 ? p : 1 - p);
	}
	return lnL;
}

void lrModel::freeAll()
{
	int r;
	if (beta)
	{
		free(beta);
		beta = 0;
	}
	if (toFit)
	{
		free(toFit);
		toFit = 0;
	}
	if (toUse)
	{
		free(toUse);
		toFit = 0;
	}
	if (t)
	{
		free(t);
		t = 0;
	}
	if (sigmaT)
	{
		free(sigmaT);
		sigmaT = 0;
	}
	if (Y)
	{
		free(Y);
		Y = 0;
	}
	if (X)
	{
		for (r = 0; r < nRow; ++r)
			if (X[r])
				free(X[r]);
		free(X);
		X = 0;
	}
	nRow = nCol = 0;
}

int lrModel::init(int r, int c)
{
	int rr,rrr;
	freeAll();
	X = (float **)calloc(r, sizeof(float*));
	if (X == 0)
	{
		dcerror(1, "Memory allocation error in lrModel::init(), r=%d c=%d", r, c);
		return 0;
	}
	for (rr = 0; rr < r; ++rr)
	{
		X[rr] = (float *)calloc(c, sizeof(float));
		if (X[rr] == 0)
		{
			for (rrr = 0; rrr < rr; ++rrr)
				free(X[rrr]);
			free(X);
			X = 0;
			dcerror(1, "Memory allocation error in lrModel::init(), r=%d c=%d", r, c);
			return 0;
		}
	}
	nRow = r;
	nCol = c;
	beta = (double *)calloc(c + 1, sizeof(double));
	toFit = (int *)calloc(c + 1, sizeof(int));
	toUse = (int *)calloc(c + 1, sizeof(int));
	t = (double *)calloc(r, sizeof(double));
	sigmaT = (double *)calloc(r, sizeof(double));
	Y = (float *)calloc(r, sizeof(float));
	if (beta == 0 || toFit == 0 || toUse == 0 || t == 0 || sigmaT == 0 || Y == 0)
	{
		dcerror(1, "Memory allocation error in lrModel::init(), r=%d c=%d", r, c);
		return 0;
	}
	return 1;
}
