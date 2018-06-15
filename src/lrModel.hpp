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

#ifndef LRMODELHPP
#define LRMODELHPP
class lrModel {
public:
	double **X, *Y;
	double *sigmaT,*t,*beta,*SE,*mean,*SD;
	int *toFit,*toUse;
	int nRow, nCol,gotMeans,isNormalised;
	char **name;
	void getMeans();
	void normalise();
	void deNormalise();
	double getLnL();
	void getSEs();
	lrModel() { toFit = 0; X = 0; beta = 0; mean = 0; nRow = nCol = 0; name = 0; gotMeans = isNormalised=0;  }
	~lrModel() { freeAll(); }
	int init(int r, int c);
	void freeAll();
	double maximiseLnL();
	virtual double penaltyFunction() { return 0;  }
};

class lrRidgePenaltyModel : public lrModel {
public:
	lrRidgePenaltyModel() { lamda = 0; }
	virtual double penaltyFunction();
	float lamda; // for ridge regression
};

#endif