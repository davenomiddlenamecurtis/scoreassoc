
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
double tLimit=8;
double stop_limit_increment = 1e-3; // 1e-7; 

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

void lrModel::getMeans()
{
	int c,r;
	if(gotMeans)
		return;
	for (c = 0; c < nCol; ++c)
		mean[c] = 0;
	for(c = 0; c < nCol; ++c)
		for(r = 0; r < nRow; ++r)
			mean[c] += X[r][c];
	for(c = 0; c < nCol; ++c)
		mean[c] /= nRow;
	gotMeans = 1;
}

#if 0

Minimisation can be difficult if some sets of values are much higher than others, so this code can normalise them.
Instead of fitting B with X, fit B' with Z.
Need to adjust B0 for intercept accordingly.
Note that the mean and SD always refer to the original values, not the normalised ones which replace them.

Xij=(Zij*Si+Mi)
Zij=(Xij-Mi)/Si
B' is for the normalised values
Ti=S(Bj*Xij) + B0
Ti=S(B'j*Zij)+B'0
=S(B'j*Xij/Sj)-S(B'j*Mj/Sj)+B'0
So Bj=B'j/Sj
And B0=B'0-S(B'j*Mj/Sj)

B'j=Bj*Sj
B'0=B0+S(Bj*Mj)

#endif

void lrModel::normalise()
{
	int c,r;
	double *X2,incBeta0,incSEBeta0;
	if(isNormalised)
		return;
	assert((X2=(double*)calloc(nCol,sizeof(double)))!=0);
	for (c = 0; c < nCol; ++c)
		mean[c] = 0;
	for(c = 0; c < nCol; ++c)
		for(r = 0; r < nRow; ++r)
		{
			mean[c] += X[r][c];
			X2[c]+=X[r][c]*X[r][c];
		}
	for(c = 0; c < nCol; ++c)
	{
		mean[c] /= nRow;
		SD[c]=sqrt(X2[c]/nRow-mean[c]*mean[c]);
	}
	incSEBeta0=incBeta0 = 0;
	for (c = 0; c < nCol; ++c)
	{
		for (r = 0; r < nRow; ++r)
		{
			X[r][c] -= mean[c];
			if (SD[c])
				X[r][c] /= SD[c];
		}
		if (SD[c])
		{
			incBeta0 += beta[c] * SD[c]; 
			incSEBeta0 += SE[c] * SD[c];
			beta[c] *= SD[c];
			SE[c] *= SD[c];
		}
	}
	beta[nCol] += incBeta0;
	SE[nCol]+=incSEBeta0;
	free(X2);
	isNormalised=gotMeans = 1;
}

void lrModel::deNormalise()
{

	int c,r;
	double incBeta0,incSEBeta0;
	if(!isNormalised)
		return;
	incSEBeta0=incBeta0 = 0;
	for(c = 0; c < nCol; ++c)
	{
		for(r = 0; r < nRow; ++r)
		{
			if(SD[c])
				X[r][c] *= SD[c];
			X[r][c] += mean[c];
		}
		if(SD[c])
		{
			incBeta0 -= beta[c] *mean[c]/ SD[c];
			incSEBeta0-=SE[c] *mean[c]/ SD[c];
			beta[c] /= SD[c];
			SE[c] /= SD[c];
		}
	}
	beta[nCol] += incBeta0;
	SE[nCol]+=incSEBeta0;
	isNormalised = 0;
}

double getRegularisedMinusModelLnL()
{
	// penalise betas, see here: http://openclassroom.stanford.edu/MainFolder/DocumentPage.php?course=MachineLearning&doc=exercises/ex5/ex5.html
	// main purpose is to stop beta hitting very large values and failing to fit properly, not to prevent over-fitting
	int c, cc;
	double lnL,halfLamda,penalty;
	halfLamda = modelToFit->nRow*0.0;
	if (!modelToFit->gotMeans)
		modelToFit->getMeans();
	for (c = cc = 0; c<modelToFit->nCol + 1; ++c)
		if (modelToFit->toFit[c])
		{
			modelToFit->beta[c] = betasToFit[cc];
			++cc;
		}
	lnL = modelToFit->getLnL();
	for (c = 0, penalty = 0; c<modelToFit->nCol; ++c) // do not include intercept
		if (modelToFit->toFit[c])
			penalty += halfLamda * modelToFit->beta[c] * modelToFit->beta[c] * (modelToFit->isNormalised?1: modelToFit->mean[c] * modelToFit->mean[c]);
	// parameters with large values (like weights) should not have large betas
	return -lnL+penalty;
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
	int c,cc,changed;
	double d,savedBeta,dd,bestD;
	double (*func)();
	func=getRegularisedMinusModelLnL;
	if (betasToFit != 0)
		free(betasToFit);
	betasToFit = (double *)calloc(nCol + 1,sizeof(double));
	nBetasToFit = 0;
	for (c = 0; c < nCol + 1; ++c)
		if (toFit[c])
			++nBetasToFit;
	column_vector starting_point(nBetasToFit);
	modelToFit = this;
	for (c = cc = 0; c < nCol + 1; ++c)
		if (toFit[c])
		{
			starting_point(cc, 0)= beta[c];
			++cc;
		}
	d = find_min_using_approximate_derivatives(cg_search_strategy(),
			objective_delta_stop_strategy(stop_limit_increment),
			lrModelMaximiser(func), starting_point, -20, derivative_eps);
		// if this fails an exception gets thrown
	for (c = cc = 0; c < nCol + 1; ++c)
			if (toFit[c])
			{
				beta[c] = starting_point(cc, 0);
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
	double lnL,eT,p,ll;
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
		{
			if (Y[r] == 1)
				ll = 0;
			else
				ll = -t[r];
		}
		else if (t[r] < -tLimit)
		{ 
			if (Y[r] == 0)
				ll = 0;
			else
				ll = t[r];
		}
		else
		{
			eT = exp(t[r]);
			p = eT / (eT + 1);
			ll = log(Y[r] == 1 ? p : 1 - p);
		}
		lnL += ll;
	}
	return lnL;
}

void lrModel::freeAll()
{
	int r;
	isNormalised=gotMeans = 0;
	if(mean)
	{
		free(mean);
		mean = 0;
	}
	if(SD)
	{
		free(SD);
		SD = 0;
	}
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
	mean = (double *)calloc(c + 1,sizeof(double));
	SD = (double *)calloc(c + 1,sizeof(double));
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