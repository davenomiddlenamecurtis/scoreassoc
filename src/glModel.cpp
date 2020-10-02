#if 0
Copyright 2018 David Curtis

This file is part of the scoreassoc package.

scoreassoc is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

scoreassoc is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with scoreassoc.If not, see <http://www.gnu.org/licenses/>.
#endif


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <dlib/optimization.h>
#include "dcerror.hpp"
#include "glModel.hpp"

using namespace dlib;
typedef matrix<double, 0, 1> column_vector;
double derivative_eps = 1e-7; // used to get gradient of lnL by beta
double second_derivative_eps = 1e-3; // I may be wrong but I think rounding errors otherwise
double tLimit=8;
double stop_limit_increment = 1e-7; // 1e-7; 

class glModelMaximiser
{
public:
	double operator() (const column_vector& arg) const;
	double (*func)(void);
	glModelMaximiser(double (*f)(void)) { func = f; }
};

// looks like I will have to use globals and assume that only one object is trying to use minimisation function

glModel *modelToFit;
double *betasToFit;
int nBetasToFit;

double glModelMaximiser::operator() (const column_vector& c) const
{
	int i;
	double f;
	for (i = 0; i<c.nr(); ++i)
		betasToFit[i] = c(i, 0);
	f = func();
	return f;
}

void glModel::getMeans()
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
Instead of fitting B with X, fit Bprime with Z.
Need to adjust B0 for intercept accordingly.
Note that the mean and SD always refer to the original values, not the normalised ones which replace them.

Xij=(Zij*Si+Mi)
Zij=(Xij-Mi)/Si
Bprime is for the normalised values
Ti=S(Bj*Xij) + B0
Ti=S(Bprimej*Zij)+Bprime0
=S(Bprimej*Xij/Sj)-S(Bprimej*Mj/Sj)+Bprime0
So Bj=Bprimej/Sj
And B0=Bprime0-S(Bprimej*Mj/Sj)

Bprimej=Bj*Sj
Bprime0=B0+S(Bj*Mj)

#endif

void glModel::normalise()
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
		if (SD[c] && toUse[c])
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

void glModel::deNormalise()
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
		if(SD[c] && toUse[c])
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

double glRidgePenaltyModel::penaltyFunction()
{
	double penalty;
	int c;
	penalty = 0;
	for (c = 0; c<nCol; ++c) // do not include intercept
		if (toFit[c])
			penalty += lamda * beta[c] * beta[c] * (isNormalised ? 1 : mean[c] * mean[c]);
	// parameters with large values (like weights) should not have large betas
	return penalty;
}

double getRegularisedMinusModelLnL()
{
	// penalise betas, see here: http://openclassroom.stanford.edu/MainFolder/DocumentPage.php?course=MachineLearning&doc=exercises/ex5/ex5.html
	// main purpose is to stop beta hitting very large values and failing to fit properly, not to prevent over-fitting
	int c, cc;
	double lnL,penalty;
	if (!modelToFit->isNormalised)
		modelToFit->normalise();
	for (c = cc = 0; c<modelToFit->nCol + 1; ++c)
		if (modelToFit->toFit[c])
		{
			modelToFit->beta[c] = betasToFit[cc];
			++cc;
		}
	lnL = modelToFit->getLnL();
	penalty = modelToFit->penaltyFunction();
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

double glModel::maximiseLnL()
{
	int c,cc;
	double d;
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
//	d = find_min_using_approximate_derivatives(cg_search_strategy(),
	d = find_min_using_approximate_derivatives(lbfgs_search_strategy(10),
			objective_delta_stop_strategy(stop_limit_increment),
			glModelMaximiser(func), starting_point, -20, derivative_eps);
		// if this fails an exception gets thrown
	for (c = cc = 0; c < nCol + 1; ++c)
			if (toFit[c])
			{
				beta[c] = starting_point(cc, 0);
				++cc;
			}
	return -d;
}

// incorporate penalty along with LL for overall function value
void glModel::getSEs()
{
	int i,ii,j,jj;
	double fittedLike, LPlus, LMinus, dPlus, dMinus,keptBetaI,keptBetaJ,d2,h;
	fittedLike = getLnL() - penaltyFunction();
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
						LMinus = getLnL()-penaltyFunction();
						dMinus=(fittedLike-LMinus)/ (second_derivative_eps*2);
						beta[i] = keptBetaI + second_derivative_eps*2;
						LPlus = getLnL() - penaltyFunction();
						dPlus = (LPlus - fittedLike) / (second_derivative_eps*2);
						d2=(dMinus- dPlus)/ (second_derivative_eps*2);
						// SE[i] = sqrt(1 / ((dMinus- dPlus)/ (second_derivative_eps*2)));
						// SE[i] = (second_derivative_eps*2)/sqrt(2 * fittedLike - LMinus - LPlus);
						// https://math.stackexchange.com/questions/302160/correct-way-to-calculate-numeric-derivative-in-discrete-time
						if (d2 < 0)
							fprintf(stderr, "Warning, d2 for beta[%d] is negative, local minimum not found\n", i);
					}
					else
					{
						keptBetaJ=beta[j];
						beta[j] = keptBetaJ - second_derivative_eps;
						beta[i] = keptBetaI - second_derivative_eps;
						LMinus = getLnL() - penaltyFunction();
						beta[i] = keptBetaI + second_derivative_eps;
						LPlus = getLnL() - penaltyFunction();
						dMinus=(LPlus-LMinus)/ (second_derivative_eps*2);
						beta[j] = keptBetaJ + second_derivative_eps;
						beta[i] = keptBetaI - second_derivative_eps;
						LMinus = getLnL() - penaltyFunction();
						beta[i] = keptBetaI + second_derivative_eps;
						LPlus = getLnL() - penaltyFunction();
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
	// matrix<double> inverseMinusHessian=inv(minusHessian);
	matrix<double> inverseMinusHessian = pinv(minusHessian);
	for (i = 0,ii=0; i < nCol + 1; ++i)
	{
		if (toFit[i])
		{
			h = inverseMinusHessian(ii, ii);
			if (h < 0)
				fprintf(stderr, "Warning, element %d of inverse negative Hessian diagonal is negative, SE will be nAn\n", ii);
			SE[i]=sqrt(h);
			++ii;
		}
		else
			SE[i] = 0;
	}
}

double glModel::getLnL()
{
	double lnL,eT,p,ll;
	int r, c;
	if (thing == LRLNLIKE)
	{
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
	}
	else
	{
		double sigmaX = 0, sigmaX2 = 0,diff,m,s2;
		for (r = 0; r < nRow; ++r)
		{
			t[r] = 0;
			for (c = 0; c < nCol; ++c)
				if (toUse[c])
					t[r] += beta[c] * X[r][c];
			if (toUse[nCol])
				t[r] += beta[nCol];
			diff = Y[r] - t[r];
			sigmaX += diff;
			sigmaX2 += diff * diff;
		}
		m = sigmaX / nRow;
		s2 = sigmaX2 / nRow - m * m; // do not use either of these
		// lnL = -nRow / 2 * (LN2PI + log(s2) + 1);
		lnL = -nRow / 2 * (LN2PI + log(sigmaX2/nRow) + 1);
		// we are using the MLE of the variance given mean is t[r], which is sigmaX2/nRow, hence the last term is 1
		// https://en.wikipedia.org/wiki/Maximum_likelihood_estimation
		// lnL = -n/2*log(2pi*s2)-1/(2*s2)*sigma((xi-m)^2)
	}
	return lnL;
}

void glModel::freeAll()
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

int glModel::init(int r, int c)
{
	int rr,rrr;
	freeAll();
	X = (double **)calloc(r, sizeof(double*));
	if (X == 0)
	{
		dcerror(1, "Memory allocation error in glModel::init(), r=%d c=%d", r, c);
		return 0;
	}
	for (rr = 0; rr < r; ++rr)
	{
		X[rr] = (double *)calloc(c, sizeof(double));
		if (X[rr] == 0)
		{
			for (rrr = 0; rrr < rr; ++rrr)
				free(X[rrr]);
			free(X);
			X = 0;
			dcerror(1, "Memory allocation error in glModel::init(), r=%d c=%d", r, c);
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
	Y = (double *)calloc(r, sizeof(double));
	if (name==0 || beta == 0 || SE == 0 || toFit == 0 || toUse == 0 || t == 0 || sigmaT == 0 || Y == 0)
	{
		freeAll();
		dcerror(1, "Memory allocation error in glModel::init(), r=%d c=%d", r, c);
		return 0;
	}
	return 1;
}
