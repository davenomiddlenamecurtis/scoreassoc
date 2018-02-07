#ifndef LRMODELHPP
#define LRMODELHPP
class lrModel {
public:
	float **X, *Y;
	double *sigmaT,*t,*beta,*SE,*mean,*SD;
	int *toFit,*toUse;
	int nRow, nCol,gotMeans,isNormalised;
	char **name;
	void getMeans();
	void normalise();
	void deNormalise();
	double getLnL();
	void getSEs();
	lrModel() { toFit = 0; X = 0; beta = 0; mean = 0; nRow = nCol = 0; name = 0; gotMeans = isNormalised=0; }
	~lrModel() { freeAll(); }
	int init(int r, int c);
	void freeAll();
	double maximiseLnL();
};

#endif