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

extern "C"
{
#include "cdflib.h" 
};
#include "dcerror.hpp"
#include "scoreassoc.hpp"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string>
#include <map>

#define PROGRAM "testGeneSets"
#define PAVERSION "1.0"

#define MAXTESTS 10
#define MAXGENES 50000
#define MAXPATHWAYGENES 12000

class geneResult {
public:
	char geneName[60];
	float result[MAXTESTS];
};

class tgsParams {
public:
	char resultsFileName[1000],outputFileName[1000],pathwayFileName[1000], outputFilePrefix[1000], outputFileSuffix[1000];
	float geneLevelOutputThreshold,correctionFactor;
	int numTests,numGenes,ignoreUninformative;
	geneResult* results;
	std::map<std::string,geneResult*> geneIndex;
	int readParms(int argc,char *argv[]);
	int getNextArg(char *nextArg,int argc,char *argv[],FILE **fpp,int *argNum);
	FILE *summaryOutputFile,*resultsFile;
};

#define LONGLINELENGTH 40000
char long_line[LONGLINELENGTH+1],rest[LONGLINELENGTH+1];

#define isArgType(a) (a[0]=='-' && a[1]=='-')
#define FILLARG(str) (strcmp(arg,str) ? 0 : ((getNextArg(arg, argc, argv, &fp, &argNum) && !isArgType(arg)) ? 1 : (dcerror(1,"No value provided for argument: %s\n",str), 0)))

double tstat(double t,double df)
{
	int which,status;
	double p,q,f,dfn,dfd,bound;
	f=t*t;
	which=1;
	dfn=1;
	dfd=df;
	cdff(&which,&p,&q,&f,&dfn,&dfd,&status,&bound);
	return q;
}

double chistat(double x, double df)
{
	double p, q, d1 = 1.0, bound;
	int status, which = 1;
	if (x == 0.0) return 1.0;
	// if (x<0.0 && x>-1.0) return 1.0; 
	if (x<0.0) return 1.0;
	/* do not worry about negative lrt values */
	cdfchi(&which, &p, &q, &x, &df, &status, &bound);
	if (status != 0)
		dcerror(1,"cdfchi() failed in chistat(%f,%f)", x,df);
	return q;
}

int tgsParams::readParms(int argc, char *argv[])
{
	char arg[2000],*ptr;
	int argNum;
	FILE *fp;
	fp=0;
	argNum=1;
	pathwayFileName[0] = resultsFileName[0] = outputFileName[0] = '\0';
	strcpy(outputFilePrefix, "tgs.");
	strcpy(outputFileSuffix, ".tgso");
	resultsFile=summaryOutputFile=0;
	geneLevelOutputThreshold = 1000;
	correctionFactor = 1.0;
	ignoreUninformative = 1; // default
	numTests = 0;
	while (getNextArg(arg, argc, argv, &fp, &argNum))
	{
		if (!isArgType(arg))
		{
			dcerror(1, "Expected argument beginning -- but got this: %s\n", arg);
			return 0;
		}
		else if (FILLARG("--arg-file"))
		{
			if (fp!=NULL)
				fclose(fp);
			fp=fopen(arg,"r");
			if (fp == NULL)
			{
				dcerror(1,"Could not open arg file: %s\n",arg);
				return 0;
			}
		}
		else if (FILLARG("--output-file-spec"))
		{
			if ((ptr = strstr(arg, "PATHWAY")) == 0)
			{
				dcerror(1, "--output-file-spec does not contain the string PATHWAY");
				return 0;
			}
			else
			{
				*ptr = '\0';
				strcpy(outputFilePrefix, arg);
				strcpy(outputFileSuffix, ptr + 7);
			}
		}
		else if (FILLARG("--ignore-uninformative"))
		{
			ignoreUninformative = atoi(arg);
		}
		else if (FILLARG("--gene-level-output-threshold"))
		{
			geneLevelOutputThreshold = atof(arg);
		}
		else if (FILLARG("--correction-factor"))
		{
			correctionFactor = atof(arg);
		}
		else if (FILLARG("--pathway-file"))
		{
			strcpy(pathwayFileName, arg);
		}
		else if (FILLARG("--results-file"))
		{
			strcpy(resultsFileName, arg);
		}
		else if (FILLARG("--summary-output-file"))
		{
			strcpy(outputFileName, arg);
		}
		else
			dcerror(1,"Did not recognise argument specifier %s\n",arg);
	}
	return 1;
}

int tgsParams::getNextArg(char *nextArg, int argc,char *argv[], FILE **fpp, int *argNum)
{
	*nextArg='\0';
	if (*fpp)
	{
		if (fscanf(*fpp,"%s ",nextArg)==1)
			return 1;
		else
		{
			fclose(*fpp);
			*fpp = NULL;
		}
	}
	if (*argNum < argc)
	{
		strcpy(nextArg,argv[*argNum]);
		++ *argNum;
		return 1;
	}
	else
		return 0;
}

float runOnePathway(char *line, tgsParams *pp, int writeFile)
{
	char pathwayName[1000], pathwayURL[1000], gene[MAXPATHWAYGENES][50], scoreFileName[1000], outputFileName[1000], thisGene[50];
	FILE *fs, *fo;
	int nGene, missing[MAXPATHWAYGENES], g,nValidGenes,t,over;
	float sigmaChisq[MAXTESTS],MLP;

	fo = 0;
	if (sscanf(line, "%s %s %[^\n]", pathwayName, pathwayURL, rest) != 3)
		return 0;
	for (nGene = 0; strcpy(line, rest), *rest = '\0', sscanf(line, "%s %[^\n]", gene[nGene], rest) >= 1; ++nGene)
		;
	nValidGenes = 0;
	for (t = 0; t < pp->numTests; ++t)
		sigmaChisq[t] = 0;
	if (writeFile)
	{
		if (pp->summaryOutputFile != 0)
			fprintf(pp->summaryOutputFile, "%s\t", pathwayName);
		sprintf(line, "%s%s%s", pp->outputFilePrefix, pathwayName, pp->outputFileSuffix);
		fo = fopen(line, "w");
		if (fo == 0)
		{
			dcerror(2, "Could not open output file %s\n", line);
		}
		else
			fprintf(fo, "%s\n%s\n\n", pathwayName, pathwayURL);
	}

	for (g = 0; g < nGene; ++g)
	{
		std::map<std::string, geneResult*>::const_iterator geneIter = pp->geneIndex.find(gene[g]);
		if (geneIter == pp->geneIndex.end())
			missing[g] = 1;
		else
		{
			++nValidGenes;
			for (t = 0; t < pp->numTests; ++t)
			{
				sigmaChisq[t] += 2*log(10)*fabs(geneIter->second->result[t])/pp->correctionFactor; // convert SLP to chisq
			}
		}
	}
	if (fo)
	{
		fprintf(fo, "%d genes, %d with results\n", nGene, nValidGenes);
		fprintf(fo, "Summed chisq:\t");
		for (t = 0; t < pp->numTests; ++t)
			fprintf(fo, "%.2f\t",sigmaChisq[t]);
		fprintf(fo, "\n");
		fprintf(fo, "Set MLP:\t");
		for (t = 0; t < pp->numTests; ++t)
		{
			MLP = -log10(chistat(sigmaChisq[t], 2 *nValidGenes));
			fprintf(fo, "%.2f ", MLP);
			if (pp->summaryOutputFile != 0)
				fprintf(pp->summaryOutputFile, "%.2f\t", MLP);
		}
		if (pp->summaryOutputFile != 0)
			fprintf(pp->summaryOutputFile, "\n");
		fprintf(fo, "\n\nGenes with result exceeding threshold of %.2f:\n", pp->geneLevelOutputThreshold);
		for (g = 0; g < nGene; ++g)
		{
			over = 0;
			std::map<std::string, geneResult*>::const_iterator geneIter = pp->geneIndex.find(gene[g]);
			if (geneIter != pp->geneIndex.end())
			{
				for (t = 0; t < pp->numTests; ++t)
				{
					if (fabs(geneIter->second->result[t]) > pp->geneLevelOutputThreshold)
						over = 1;
				}
				if (over == 1)
				{
					fprintf(fo, "%s\t", gene[g]);
					for (t = 0; t < pp->numTests; ++t)
						fprintf(fo, "%.2f\t", geneIter->second->result[t]);
					fprintf(fo, "\n");
				}
			}
		}
	}
	if (fo)
		fclose(fo);
	return MLP; // just return the value for the last test
}

int readResultsFile(tgsParams *pp)
{
	char gene[50],rest[2000],t,informative;
	pp->geneIndex.clear();
	pp->resultsFile = fopen(pp->resultsFileName, "r");
	if (pp->resultsFile == 0)
	{
		dcerror(1, "Could not open results file: %s\n", pp->resultsFileName);
		return 0;
	}
	fgets(long_line, LONGLINELENGTH, pp->resultsFile);
	sscanf(long_line, "%s %[^\n]", gene, rest);
	if (pp->summaryOutputFile)
		fprintf(pp->summaryOutputFile, "GeneSet\t");
	pp->numTests = 0;
	while (strcpy(long_line, rest), *rest = '\0', sscanf(long_line, "%s %[^\n]", gene, rest) > 0)
	{
		++pp->numTests;
		if (pp->summaryOutputFile)
			fprintf(pp->summaryOutputFile, "%s\t", gene); // just filling out the column headings
	}
	if (pp->summaryOutputFile)
		fprintf(pp->summaryOutputFile, "\n");
	pp->numGenes = 0;
	while (fgets(long_line, LONGLINELENGTH, pp->resultsFile) && sscanf(long_line, "%s %[^\n]", pp->results[pp->numGenes].geneName, rest)==2)
	{
		t = 0;
		informative = 0;
		while (strcpy(long_line, rest), *rest = '\0', sscanf(long_line, "%f %[^\n]", &pp->results[pp->numGenes].result[t], rest) > 0)
			if (pp->results[pp->numGenes].result[t++]!=0)
				informative=1;
		if (pp->ignoreUninformative == 0 || informative)
		{
			pp->geneIndex[pp->results[pp->numGenes].geneName] = &pp->results[pp->numGenes];
			++pp->numGenes;
		}
	}
	return 1;
}

int main(int argc, char *argv[])
{
	tgsParams pp;
	FILE *fp;
	int s;
	printf("%s v%s\n",PROGRAM,PAVERSION);
	pp.ignoreUninformative = 1;
	if (!pp.readParms(argc,argv))
		exit(1);
	assert((pp.results = (geneResult * )calloc(MAXGENES, sizeof(geneResult))) != 0);
	if ((fp = fopen(pp.pathwayFileName,"r")) == 0)
	{
		dcerror(2,"Could not open gene set file %s\n",pp.pathwayFileName);
		return 1;
	}
	if (pp.outputFileName)
	{
		pp.summaryOutputFile = fopen(pp.outputFileName, "w");
		if (pp.summaryOutputFile==0)
		{
			dcerror(2, "Could not open summary output file %s\n", pp.outputFileName);
			return 1;
		}
	}
	if (!readResultsFile(&pp))
		exit(1);
	while (fgets(long_line, LONGLINELENGTH, fp))
	{
		runOnePathway(long_line,&pp,1);
	}
	if (pp.summaryOutputFile!=0)
		fclose(pp.summaryOutputFile);
	free(pp.results);
	return 0;
}
