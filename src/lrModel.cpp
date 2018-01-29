#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <dlib/optimization.h>
#include "dcerror.hpp"
#include "lrModel.hpp"

using namespace dlib;
typedef matrix<double, 0, 1> column_vector;
double derivative_eps = 1e-7; // used to get gradient of lnL by beta
double second_derivative_eps = 1e-5; // I may be wrong but I think rounding errors otherwise
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
	ftol = 1E-7; // 0.001;
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
		lrModelMaximiser(getMinusModelLnL), starting_point, -20, derivative_eps);
	// if this fails an exception gets thrown
	for (c = cc = 0; c<nCol + 1; ++c)
		if (toFit[c])
		{
			beta[c]=starting_point(cc,0);
			++cc;
		}
	return -d;
}

void lrModel::getSEs()
{
	int i,ii,j,jj;
	double fittedLike, LPlus, LMinus, dPlus, dMinus,keptBetaI,keptBetaJ,d2;
	fittedLike = getLnL();
	// make hessian matrix then invert it and take square roots of diagonal elements
	matrix<double> minusHessian(nBetasToFit,nBetasToFit);
	for(i = 0,ii=0; i < nCol + 1; ++i)
	{
		if(toFit[i])
		{
			keptBetaI=beta[i];
			for(j=0,jj=0;j<nCol+1;++j)
			{
				if(toFit[j])
				{
					if(i==j)
					{
						beta[i] = keptBetaI - second_derivative_eps*2;
						LMinus = getLnL();
						dMinus=(fittedLike-LMinus)/ (second_derivative_eps*2);
						beta[i] = keptBetaI + second_derivative_eps*2;
						LPlus = getLnL();
						dPlus = (LPlus - fittedLike) / (second_derivative_eps*2);
						d2=(dMinus- dPlus)/ (second_derivative_eps*2);
						// SE[i] = sqrt(1 / ((dMinus- dPlus)/ (second_derivative_eps*2)));
						// SE[i] = (second_derivative_eps*2)/sqrt(2 * fittedLike - LMinus - LPlus);
						// https://math.stackexchange.com/questions/302160/correct-way-to-calculate-numeric-derivative-in-discrete-time

					}
					else
					{
						keptBetaJ=beta[j];
						beta[j] = keptBetaJ - second_derivative_eps;
						beta[i] = keptBetaI - second_derivative_eps;
						LMinus = getLnL();
						beta[i] = keptBetaI + second_derivative_eps;
						LPlus = getLnL();
						dMinus=(LPlus-LMinus)/ (second_derivative_eps*2);
						beta[j] = keptBetaJ + second_derivative_eps;
						beta[i] = keptBetaI - second_derivative_eps;
						LMinus = getLnL();
						beta[i] = keptBetaI + second_derivative_eps;
						LPlus = getLnL();
						dPlus=(LPlus-LMinus)/ (second_derivative_eps*2);
						d2=(dMinus- dPlus)/ (second_derivative_eps*2);
						beta[j]=keptBetaJ;
					}
					minusHessian(ii,jj)=minusHessian(jj,ii)=d2;
					++jj;
				}
			}
			beta[i]=keptBetaI;
			++ii;
		}
	}
	matrix<double> inverseMinusHessian=inv(minusHessian);
	for (i = 0,ii=0; i < nCol + 1; ++i)
	{
		if (toFit[i])
		{
			SE[i]=sqrt(inverseMinusHessian(ii,ii));
			++ii;
		}
		else
			SE[i] = 0;
	}
}

double lrModel::getLnL()
{
	double lnL,eT,p;
	int r, c;
	lnL = 0;
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
	if (name)
	{
		free(name);
		name = 0;
	}
	if (beta)
	{
		free(beta);
		beta = 0;
	}
	if (SE)
	{
		free(SE);
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
	name = (char **)calloc(c + 1, sizeof(char*));
	beta = (double *)calloc(c + 1, sizeof(double));
	SE = (double *)calloc(c + 1, sizeof(double));
	toFit = (int *)calloc(c + 1, sizeof(int));
	toUse = (int *)calloc(c + 1, sizeof(int));
	t = (double *)calloc(r, sizeof(double));
	sigmaT = (double *)calloc(r, sizeof(double));
	Y = (float *)calloc(r, sizeof(float));
	if (name==0 || beta == 0 || SE == 0 || toFit == 0 || toUse == 0 || t == 0 || sigmaT == 0 || Y == 0)
	{
		freeAll();
		dcerror(1, "Memory allocation error in lrModel::init(), r=%d c=%d", r, c);
		return 0;
	}
	return 1;
}
