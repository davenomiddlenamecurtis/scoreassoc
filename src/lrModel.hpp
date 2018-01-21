#ifndef LRMODELHPP
#define LRMODELHPP
class lrModel {
public:
	float **X, *Y;
	double *sigmaT,*t,*beta;
	int *toFit,*toUse;
	int nRow, nCol;
	double getLnL();
	lrModel() { toFit =0; X = 0; beta = 0; nRow = nCol = 0; }
	~lrModel() { freeAll(); }
	int init(int r, int c);
	void freeAll();
	double maximiseLnL();
};

#endif