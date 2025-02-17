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

#include <ctype.h>
#include <string.h>
#include <string>
#include <math.h>
#include <assert.h>
#include "scoreassoc.hpp"
#include "glModel.hpp"
#include "dcerror.hpp"

char *intercept="intercept";

float runLinTestFile(FILE* fo, char* fn, glModel* m, lr_test_par_info* spi)
{
	float f;
	m->useLinearRegression(1);
	f=runTestFile(fo, fn, m, spi);
	m->useLinearRegression(0);
	return f;
}

float runTestFile(FILE* fo, char* fn, glModel* m, lr_test_par_info* spi)
{
	FILE *ft;
	int b, scoreVar, toUse[MAXLRVARIABLES + 1], toFit0[MAXLRVARIABLES + 1], toFit1[MAXLRVARIABLES + 1],t0,t1;
	float startBetas[MAXLRVARIABLES + 1],chisq,df,p,sB,SLP,L0,L1;
	char varName[MAXLRVARIABLENAMELENGTH];
	ft = fopen(fn, "r");
	hereOK();
	if (ft == 0)
	{
		dcerror(1, "Could not open test file %s\n", fn);
		return 0;
	}
	hereOK();
	for (b = 0; b < m->nCol; ++b)
		toUse[b] = toFit0[b] = toFit1[b] = 0;
	toUse[b] = toFit0[b] = toFit1[b] = 1;
	startBetas[b] = 0;
	scoreVar = -1;
	df = 0;
	while (fgets(long_line, 200, ft) && sscanf(long_line, "%s %f %d %d", varName, &sB,&t0, &t1) == 4)
	{
		if (!strcmp(varName, "intercept"))
			b = m->nCol;
		else
		{
			for (b = 0; b < spi->numVars; ++b)
				if (!strcmp(varName, allVars[b].name))
					break;
			if (b == spi->numVars)
			{
				dcerror(1, "Test file %s contains unknown variable called %s\n", fn, varName);
				return 0;
			}
		}
		toUse[b] = 1;
		startBetas[b] = sB;
		toFit0[b] = t0;
		toFit1[b] = t1;
		if (!strcmp(varName, "score"))
			scoreVar = b;
		df += t1 - t0;
	}
	fclose(ft);
	hereOK();
	L0 = evaluateModel(fo, m, toUse, startBetas, toFit0, "L0");
	hereOK();
	if (spi->start_from_fitted)
		for (b = 0; b < m->nCol + 1; ++b)
			startBetas[b] = m->beta[b];
	L1 = evaluateModel(fo, m, toUse, startBetas, toFit1, "L1");
	chisq = 2 * (L1 - L0);
	p = chistat(chisq, df);
	if (fo)
		fprintf(fo,
			"\nchisq = %.2f, %.0f df, p = %10.8f\n", chisq,df, p);
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
		fprintf(fo,"tMLP = %8.2f (minus log10(p))\n",SLP=-log10(p)); // this is a MLP not SLP but this is the what we want to return
	}
	return SLP;
}

int fillModelWithVars(glModel *m, int nsub, lr_test_par_info *spi,int which)
{
	int s, b;
	if (which == -1)
	{
		if (!m->init(nsub, spi->numVars))
			return 0;
		for (s = 0; s < nsub; ++s)
		{
			for (b = 0; b < spi->numVars; ++b)
				m->X[s][b] = allVars[b].val[s];
		}
		for (b = 0; b < spi->numVars; ++b)
			m->name[b] = allVars[b].name;
		m->name[b] = intercept;
	}
	else
		for (s = 0; s < nsub; ++s)
			m->X[s][which] = allVars[which].val[s];
	return 1;
}

float evaluateModel(FILE *fo, glModel *m, int *toUse, float *startBetas, int *toFit,char *name)
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
//	m->deNormalise(); // avoid repeated normalising and deNormalising
	if (fo)
		printModel(fo, name,lnL,m);
	return lnL;
}

void printModel(FILE *fo, char *LLstr,double LL,glModel *m)
{
	// change this to be row-wise and use names
	int b,bb,c;
	double *realBeta, *realSE;
	double incBeta0, incSEBeta0;
	realBeta = (double*)calloc(m->nCol + 1, sizeof(double));
	realSE = (double*)calloc(m->nCol + 1, sizeof(double));
	fprintf(fo, "\n%s = %.2f\n", LLstr, LL);
	fprintf(fo, "%-" LOCUS_NAME_LENGTH_STR "s %-10s %-10s %-10s\n", "beta", "value", "SE","z");
	if (m->isNormalised)
	{ 
		incSEBeta0 = incBeta0 = 0;
		for (c = 0; c < m->nCol; ++c)
		{
			if (m->SD[c] && m->toUse[c])
			{
				incBeta0 -= m->beta[c] * m->mean[c] / m->SD[c];
				incSEBeta0 -= m->SE[c] * m->mean[c] / m->SD[c];
				realBeta[c] = m->beta[c] / m->SD[c];
				realSE[c] = m->SE[c] / m->SD[c];
			}
		}
		realBeta[m->nCol] = m->beta[m->nCol]+incBeta0;
		realSE[m->nCol] = m->SE[m->nCol]+incSEBeta0;
	}
	else 
		for (b = 0; b < m->nCol + 1; ++b)
		{
			realBeta[b] = m->beta[b];
			realSE[b] = m->SE[b];
		}
	for (b = 0; b < m->nCol+1; ++b)
	{
		bb = (b + m->nCol) % (m->nCol + 1); // print last first
		if (m->toUse[bb])
			fprintf(fo, "%-" LOCUS_NAME_LENGTH_STR "s %10.5f %10.5f %10.5f\n", m->name[bb], realBeta[bb], realSE[bb], m->toFit[bb]?realBeta[bb] /realSE[bb]:0.0);
	}
	free(realBeta);
	free(realSE);
}

float do_onetailed_LRT(FILE *fo, glModel *m,lr_test_par_info *spi,int do_linear)
{
	double L0, L1, p, SLP,chisq;
	int b,scoreVar,toUse[MAXLRVARIABLES+1],toFit[MAXLRVARIABLES+1];
	float startBetas[MAXLRVARIABLES+1];
	m->useLinearRegression(do_linear);
	for(b=0;b<spi->numVars;++b)
		if(!strcmp(allVars[b].name,"score")|| !strcmp(allVars[b].name, "recscore"))
			break;
	if(b==spi->numVars)
	{
		dcerror(1,"Variable named score or recscore is not available for default regression test\n");
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
			"\nchisq = %.2f, 1 df, p = %10.8f\n"
			"%s = %8.2f %s",
			chisq,p,
			do_linear ? "linrSLP" : "lrSLP",
			SLP,
			do_linear ? "(signed log10(p), positive if weighted burden of variants correlates positively with phenotype)\n" 
			: "(signed log10(p), positive if cases have more variants than controls)\n");
	m->useLinearRegression(0);
	return SLP;
}

#define MISSING -999
int readVarFiles(std::map<std::string, int> subIDs, int nSub, lr_test_par_info  *spi)
{
	int i, s, idCol, c, colIndex[MAXLRVARIABLES], nCol;
	char colValue[MAXLRVARIABLENAMELENGTH], *ptr, *sptr;
	for (i = 0; i<spi->numVarFiles; ++i)
	{
		spi->varFiles[i].fp = fopen(spi->varFiles[i].fn, "r");
		if (spi->varFiles[i].fp == 0)
		{
			dcerror(1, "Could not open variable file %s\n", spi->varFiles[i].fn);
			return 0;
		}
		fgets(long_line, LONG_LINE_LENGTH, spi->varFiles[i].fp);
		ptr = long_line;
		idCol = -1;
		while (*ptr && isspace(*ptr))
			++ptr;
		for (c = 0; 1; ++c)
		{
			sptr = colValue;
			while (*ptr && !isspace(*ptr))
				*sptr++ = *ptr++;
			*sptr = '\0';
			while (*ptr && isspace(*ptr))
				++ptr;
			if (!strcmp(colValue, "IID"))
			{
				idCol = c;
				break;
			}
			if (*ptr == '\0')
				break;
		}
		if (idCol == -1)
		{
			dcerror(1, "Variable file %s did not contain a column headed IID\n", spi->varFiles[i].fn);
			return 0;
		}
		for (c = idCol + 1; 1; ++c)
		{
			sptr = colValue;
			while (*ptr && !isspace(*ptr))
				*sptr++ = *ptr++;
			*sptr = '\0';
			while (*ptr && isspace(*ptr))
				++ptr;
			if (colValue[0])
			{
				if (!strcmp(colValue, "score"))
				{
					dcerror(1, "Variable file %s contains a column headed score, which is not allowed because program will calculate scores\n", spi->varFiles[i].fn);
					return 0;
				}
				std::map<std::string, lrVariable *>::const_iterator varIter = varMap.find(colValue);
				if (varIter != varMap.end())
				{
					dcerror(1, "Variable file %s contains a column headed %s but that variable name is already in use\n", spi->varFiles[i].fn, colValue);
					return 0;
				}
				if (spi->numVars >= MAXLRVARIABLES)
				{
					dcerror(1, "Variable file %s contains a column headed %s but the maximum number of variables as specified by MAXLRVARIABLES is %d and has already been reached.\nPlease increase MAXLRVARIABLES then recompile.\n", 
						spi->varFiles[i].fn, colValue,MAXLRVARIABLES);
					return 0;
				}
				strcpy(allVars[spi->numVars].name, colValue);
				assert((allVars[spi->numVars].val = (double*)calloc(nSub, sizeof(double))) != 0);
				for (s = 0; s<nSub; ++s)
					allVars[spi->numVars].val[s] = MISSING;
				colIndex[c] = spi->numVars;
				varMap[colValue] = &allVars[spi->numVars];
				++spi->numVars;
			}
			if (*ptr == '\0')
				break;
		}
		nCol = c + 1;
		while (fgets(long_line, LONG_LINE_LENGTH, spi->varFiles[i].fp))
		{
			ptr = long_line;
			while (*ptr && isspace(*ptr))
				*ptr++;
			for (c = 0; c <= idCol; ++c)
			{
				sptr = colValue;
				while (*ptr && !isspace(*ptr))
					*sptr++ = *ptr++;
				*sptr = '\0';
				while (*ptr && isspace(*ptr))
					++ptr;
				if (*ptr == '\0')
					break;
			}
			
			if (colValue[0] == '\0')
			{
				dcerror(1, "IID value missing from variable file %s in this line:\n%s\n", spi->varFiles[i].fn, long_line);
				return 0;
			}
			std::map<std::string, int>::const_iterator idIter = subIDs.find(colValue);
			if (idIter == subIDs.end())
			{
				continue; // this is not an error - maybe values for subjects not included in the analysis
						  // dcerror(1,"Unknown IID value in variable file %s in this line:\n%s\n",spi->varFiles[i].fn,long_line);
						  // return 0;
			}
			else
				s = idIter->second;
			for (; c<nCol; ++c)
			{
				sptr = colValue;
				while (*ptr && !isspace(*ptr))
					*sptr++ = *ptr++;
				*sptr = '\0';
				while (*ptr && isspace(*ptr))
					++ptr;
				if (sscanf(colValue, "%lf", &allVars[colIndex[c]].val[s]) != 1)
				{
					dcerror(1, "Not enough values in variable file %s in this line:\n%s\n", spi->varFiles[i].fn, long_line);
					return 0;
				}
			}
		}
		int wentWrong=0;
		for (std::map<std::string, int>::iterator i1=subIDs.begin(); i1!=subIDs.end(); ++i1)
			if (allVars[colIndex[idCol + 1]].val[s=i1->second] == MISSING) // only have to check one to know if the subject line was found
			{
				
//				dcerror(1, "Missing values in variable file %s for subject %s (s=%d, idCol+1=%d, colIndex[idCol+1]=%d)\n", 
				printf("Missing values in variable file %s for subject %s (s=%d, idCol+1=%d, colIndex[idCol+1]=%d)\n", 
					spi->varFiles[i].fn, (char*)(i1->first.c_str()),
					s, idCol + 1, colIndex[idCol + 1]);
//				return 0;
				wentWrong=1;
			}
		if (wentWrong)
		{
			dcerror(1,"Aborting because of missing values in variable file %s\n",spi->varFiles[i].fn);
			return 0;
		}
		fclose(spi->varFiles[i].fp);
		spi->varFiles[i].fp = 0;
	}
	return 1;
}
