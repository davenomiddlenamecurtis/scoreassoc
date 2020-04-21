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
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#define PROGRAM "combineCohorts"
#define CCVERSION "1.0"

class ccParams {
public:
	char scoreFileSpec[1000],outputFileName[1000],geneListFileName[1000],cohortListFileName[1000];
	int readParms(int argc, char *argv[]);
	int getNextArg(char *nextArg,int argc,char *argv[],FILE **fpp,int *argNum);
};

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

int ccParams::readParms(int argc,char *argv[])
{
	char arg[2000],*ptr;
	int argNum;
	FILE *fp;
	fp=0;
	argNum=1;
	scoreFileSpec[0]=outputFileName[0]=geneListFileName[0]=cohortListFileName[0]='\0';
	while (getNextArg(arg,argc,argv,&fp,&argNum))
	{
		if (!isArgType(arg))
		{
			dcerror(1,"Expected argument beginning -- but got this: %s\n",arg);
			return 0;
		}
		else if (FILLARG("--arg-file"))
		{
			if (fp != NULL)
				fclose(fp);
			fp=fopen(arg,"r");
			if (fp == NULL)
			{
				dcerror(1,"Could not open arg file: %s\n",arg);
				return 0;
			}
		}
		else if (FILLARG("--score-file-spec"))
		{
			if ((ptr = strstr(arg,"GENE")) == 0)
			{
				dcerror(1,"--score-file-spec does not contain the string GENE\n");
				return 0;
			}
			else if ((ptr = strstr(arg,"COHORT")) == 0)
			{
				dcerror(1,"--score-file-spec does not contain the string GENE\n");
				return 0;
			}
			else
			{
				strcpy(scoreFileSpec,arg);
			}
		}
		else if (FILLARG("--output-file"))
		{
			strcpy(outputFileName,arg);
		}
		else if (FILLARG("--gene-list-file"))
		{
			strcpy(geneListFileName,arg);
		}
		else if (FILLARG("--cohort-list-file"))
		{
			strcpy(cohortListFileName,arg);
		}
		else
			dcerror(1,"Did not recognise argument specifier %s\n",arg);
	}
	// do checks
	return 1;
}

int ccParams::getNextArg(char *nextArg,int argc,char *argv[],FILE **fpp,int *argNum)
{
	*nextArg='\0';
	if (*fpp)
	{
		if (fscanf(*fpp,"%s ",nextArg) == 1)
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

#define MAX_COHORTS 100
#define MAX_COHORT_LENGTH 30

class SLPcalculator {
public:
	int nSub[2];
	float mean[2],var[2],sd[2],sigmaX[2],sigmaX2[2],SSQ[2],tval,SLP;
	double p;
	void addScore(int cc,float score);
	void clear();
	float getSLP(int useSSQs=0);
	void augment(SLPcalculator &source);
	SLPcalculator() { clear(); }
};

void SLPcalculator::addScore(int cc,float score)
{
	sigmaX[cc]+=score;
	sigmaX2[cc]+=score*score;
	++nSub[cc];
}

void SLPcalculator::clear()
{
	int cc;
	for (cc=0;cc<2;++cc)
		SSQ[cc]=sigmaX[cc]=sigmaX2[cc]=nSub[cc]=0;
}

void SLPcalculator::augment(SLPcalculator &source)
{
	int cc;
	float sourceMean;
	sourceMean=(source.sigmaX[0]+source.sigmaX[1])/(source.nSub[0]+source.nSub[1]);
	for (cc=0;cc<2;++cc)
	{
		sigmaX[cc]+=source.sigmaX[cc]-sourceMean*source.nSub[cc];
		SSQ[cc]+=source.SSQ[cc];
		nSub[cc]+=source.nSub[cc];
	}
}

float SLPcalculator::getSLP(int useSSQs)
{
	int cc;
	float SE,s2;
	for (cc=0;cc<2;++cc)
	{
		mean[cc]=sigmaX[cc]/nSub[cc];
		if (!useSSQs)
			SSQ[cc]=(sigmaX2[cc]-sigmaX[cc]*sigmaX[cc]/nSub[cc]);
		var[cc]=SSQ[cc]/(nSub[cc]-1);
		sd[cc]=sqrt(var[cc]);
	}
	s2=((nSub[0]-1)*var[0]+(nSub[1]-1)*var[1])/(nSub[0]+nSub[1]-2);
	if (s2==0)
	{
		SE=0;
		tval=0;
		p=0.5;
	}
	else
	{
		SE=sqrt(s2*(1/(float)nSub[0]+1/(float)nSub[1]));
		tval=(mean[1]-mean[0])/SE;
		p=tstat(tval,nSub[0]+nSub[1]-2.0)/2; // one-tailed
	}
	SLP=log10(2*p)*(mean[0]>=mean[1]?1:-1);
	return SLP;
}

int writeSLPs(char *scoreFileSpec,char cohorts[MAX_COHORTS][MAX_COHORT_LENGTH],int nCohorts,char *geneName,FILE *fo)
{
	SLPcalculator totalSC,cohortSC;
	int c,cc;
	float score;
	FILE *fs;
	char fn[1000],specBuff1[1000],specBuff2[1000],*ptr,*sptr,line[1001];
	strcpy(specBuff1,scoreFileSpec);
	for (sptr=specBuff1,specBuff2[0]='\0';ptr=strstr(sptr,"GENE");sptr=ptr+4)
	{
		*ptr='\0';
		strcat(specBuff2,sptr);
		strcat(specBuff2,geneName);
	}
	strcat(specBuff2,sptr);
	strcpy(specBuff1,specBuff2);
	fprintf(fo,"%s\t",geneName);
	for (c=0;c<nCohorts;++c)
	{
		strcpy(specBuff2,specBuff1);
		for (sptr=specBuff2,fn[0]='\0';ptr=strstr(sptr,"COHORT");sptr=ptr+6)
		{
			*ptr='\0';
			strcat(fn,sptr);
			strcat(fn,cohorts[c]);
		}
		strcat(fn,sptr);
		fs=fopen(fn,"r"); // not an error if this fails
		if (fs)
		{
			cohortSC.clear();
			while (fgets(line,1000,fs) && sscanf(line,"%*s %d %f",&cc,&score)==2)
			{
				cohortSC.addScore(cc,score);
			}
			fprintf(fo,"%.2f\t",cohortSC.getSLP());
			totalSC.augment(cohortSC); // must be after have called getSLP()
			fclose(fs);
		}
		else
			fprintf(fo,"0\t");
	}
	fprintf(fo,"%.2f\t",totalSC.getSLP(1));
	for (cc=0;cc<2;++cc)
		fprintf(fo,"%d\t",totalSC.nSub[cc]);
	for (cc=0;cc<2;++cc)
		fprintf(fo,"%.2f\t",totalSC.mean[cc]);
	for (cc=0;cc<2;++cc)
		fprintf(fo,"%.2f\t",totalSC.sd[cc]);
	fprintf(fo,"%.2f\n",totalSC.tval);
	return 1;
}

int main(int argc,char *argv[])
{
	ccParams cp;
	FILE *fo,*fc,*fg;
	char cohorts[MAX_COHORTS][MAX_COHORT_LENGTH],line[1001],geneName[1001];
	int nCohorts,c;
	printf("%s v%s\n",PROGRAM,CCVERSION);
	dcerror.kill();
	if (!cp.readParms(argc,argv))
		exit(1);
	if (cp.cohortListFileName[0]=='\0')
		dcerror(1,"Need to specify --cohort-list-file");
	if (cp.geneListFileName[0]=='\0')
		dcerror(1,"Need to specify --gene-list-file");
	if (cp.scoreFileSpec[0]=='\0')
		dcerror(1,"Need to specify --score-file-spec");
	if (cp.outputFileName[0]=='\0')
		dcerror(1,"Need to specify --output-file");
	fo=fopen(cp.outputFileName,"w");
	if (fo==0)
		dcerror(1,"Could not open output file %s",cp.outputFileName);
	fc=fopen(cp.cohortListFileName,"r");
	if (fc==0)
		dcerror(1,"Could not open cohort list file %s",cp.cohortListFileName);
	fg=fopen(cp.geneListFileName,"r");
	if (fg==0)
		dcerror(1,"Could not open gene list file %s",cp.geneListFileName);
	nCohorts=0;
	while (fgets(line,1000,fc) && sscanf(line,"%s",cohorts[nCohorts])==1)
		++nCohorts;
	fclose(fc);
	fprintf(fo, "Gene\t");
	for (c=0;c<nCohorts;++c)
		fprintf(fo,"%s\t",cohorts[c]);
	fprintf(fo,"SLP\tnCont\tnCase\tmCont\tmCase\tsdCont\tsdCase\tt\n");
	while (fgets(line,1000,fg) && sscanf(line,"%s",geneName)==1)
	{
		writeSLPs(cp.scoreFileSpec,cohorts,nCohorts,geneName,fo);
	}
	fclose(fo);
	return 0;
}
