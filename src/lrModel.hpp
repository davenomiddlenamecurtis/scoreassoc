#ifndef LRMODELHPP
#define LRMODELHPP
class lrModel {
public:
	float **X, *Y;
	double *sigmaT,*t,*beta,*SE;
	int *toFit,*toUse;
	int nRow, nCol;
	char **name;
	double getLnL();
	void getSEs();
	lrModel() { toFit = 0; X = 0; beta = 0; nRow = nCol = 0; name = 0; }
	~lrModel() { freeAll(); }
	int init(int r, int c);
	void freeAll();
	double maximiseLnL();
};

#endif