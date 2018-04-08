#include <ctype.h>
#include <string.h>
#include <string>
#include <math.h>
#include <assert.h>
#include "scoreassoc.hpp"
#include "lrModel.hpp"
#include "dcerror.hpp"

char *intercept="intercept";

float runTestFile(FILE *fo, char *fn, lrModel *m, par_info *pi, sa_par_info *spi)
{
	FILE *ft;
	int b, scoreVar, toUse[MAXLRVARIABLES + 1], toFit0[MAXLRVARIABLES + 1], toFit1[MAXLRVARIABLES + 1],t0,t1;
	float startBetas[MAXLRVARIABLES + 1],chisq,df,p,sB,SLP,L0,L1;
	char varName[MAXLRVARIABLENAMELENGTH];
	ft = fopen(fn, "r");
	if (ft == 0)
	{
		dcerror(1, "Could not open test file %s\n", fn);
		return 0;
	}
	for (b = 0; b < m->nCol; ++b)
		toUse[b] = toFit0[b] = toFit1[b] = 0;
	toUse[b] = toFit0[b] = toFit1[b] = 1;
	startBetas[b] = 0;
	scoreVar = -1;
	df = 0;
	while (fgets(long_line, 200, ft) && sscanf(long_line, "%s %f %d %d", varName, &sB,&t0, &t1) == 4)
	{
		for (b = 0; b < spi->numVars; ++b)
			if (!strcmp(varName, allVars[b].name))
				break;
		if (b == spi->numVars)
		{
			dcerror(1, "Test file %s contains unknown variable called %s\n", fn, varName);
			return 0;
		}
		toUse[b] = 1;
		startBetas[b] = sB;
		toFit0[b] = t0;
		toFit1[b] = t1;
		if (!strcmp(varName, "score"))
			scoreVar = b;
		df += t1 - t0;
	}
	L0 = evaluateModel(fo, m, toUse, startBetas, toFit0, "L0");
	L1 = evaluateModel(fo, m, toUse, startBetas, toFit1, "L1");
	chisq = 2 * (L1 - L0);
	p = chistat(chisq, df);
	if (fo)
		fprintf(fo,
			"chisq = %.2f, %.0f df, p = %10.8f\n", chisq,df, p);
	if (scoreVar != -1 && df == 1 && toFit0[scoreVar] == 0 && toFit1[scoreVar] == 1)
	{
		SLP = -log10(p);
		if (m->beta[scoreVar] < 0)
			SLP *= -1;
		if (fo)
			fprintf(fo, "tSLP = %8.2f (signed log10(p), positive if cases have more variants than controls)\n", SLP);
	}
	else
	{
		fprintf(fo,"tMLP = %8.2f (minus log10(p))\n",-log10(p));
	}
	return p;


}

void fillModelWithVars(lrModel *m, subject **sub, int nsub, par_info *pi, sa_par_info *spi)
{
	int s, b;
	m->init(nsub, spi->numVars);
	m->lamda=spi->lamda;
	for (s = 0; s < nsub; ++s)
	{
		for (b = 0; b<spi->numVars; ++b)
			m->X[s][b] = allVars[b].val[s];
		m->Y[s] = sub[s]->cc;
	}
	for (b = 0; b<spi->numVars; ++b)
		m->name[b] = allVars[b].name;
	m->name[b]=intercept;
}

float evaluateModel(FILE *fo, lrModel *m, int *toUse, float *startBetas, int *toFit,char *name)
{
	float lnL;
	int b;
	for (b = 0; b < m->nCol+1; ++b)
	{
		m->toUse[b] = toUse[b];
		m->beta[b] = startBetas[b];
		m->toFit[b] = toFit[b];
	}
	m->normalise();
 	lnL = m->maximiseLnL();
	m->getSEs();
	m->deNormalise();
	if (fo)
		printModel(fo, name,lnL,m);
	return lnL;
}

void printModel(FILE *fo, char *LLstr,double LL,lrModel *m)
{
	// change this to be row-wise and use names
	int b,bb;
	fprintf(fo, "\n%s = %.2f\n", LLstr, LL);
	fprintf(fo, "%-" LOCUS_NAME_LENGTH_STR "s %-10s %-10s %-10s\n", "beta", "value", "SE","z");
	for (b = 0; b < m->nCol+1; ++b)
	{
		bb = (b + m->nCol) % (m->nCol + 1); // print last first
		if (m->toUse[bb])
			fprintf(fo, "%-" LOCUS_NAME_LENGTH_STR "s %10.5f %10.5f %10.5f\n", m->name[bb], m->beta[bb], m->SE[bb], m->toFit[bb]?m->beta[bb] /m->SE[bb]:0.0);
	}
}

float do_onetailed_LRT(FILE *fo, lrModel *m,par_info *pi, sa_par_info *spi)
{
	double L0, L1, p, SLP,chisq;
	int b,scoreVar,toUse[MAXLRVARIABLES+1],toFit[MAXLRVARIABLES+1];
	float startBetas[MAXLRVARIABLES+1];
	for(b=0;b<spi->numVars;++b)
		if(!strcmp(allVars[b].name,"score"))
			break;
	if(b==spi->numVars)
	{
		dcerror(1,"Variable named score is not available\n");
		return 0;
	}
	else
		scoreVar=b;
	for(b=0;b<spi->numVars;++b)
	{
		toUse[b]=1;
		toFit[b]=b!=scoreVar;
		startBetas[b]=0;
	}
	toUse[b]=toFit[b]=1;
	startBetas[b]=0;
	L0=evaluateModel(fo,m,toUse,startBetas,toFit,"L0");
	toFit[scoreVar]=1;
	L1=evaluateModel(fo,m,toUse,startBetas,toFit,"L1");
	chisq=2 * (L1 - L0);
	p = chistat(chisq,1.0); 
	SLP = -log10(p);
	if (m->beta[scoreVar] < 0)
		SLP *= -1;
	if (fo)
		fprintf(fo,
			"chisq = %.2f, 1 df, p = %10.8f\n"
			"lrSLP = %8.2f (signed log10(p), positive if cases have more variants than controls)\n",chisq,p,SLP);
	return SLP;
}

void fillModel(lrModel *m, float *score, subject **sub, int nsub, par_info *pi, sa_par_info *spi)
{
	int s;
	m->init(nsub, 1);
	for (s = 0; s < nsub; ++s)
	{
		m->X[s][0] = score[s];
		m->Y[s] = sub[s]->cc;
	}

}

