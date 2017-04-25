#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "dcerror.hpp"

#define PROGRAM "fitScores"
#define FSVERSION "1.0"
#define MAXSUB 15000
#define MAXSCORESTOREAD 20
#define MAXPARAMS 100

struct param_t { float val; char name[200]; int toFit; };
typedef struct param_t param;


struct subject_t { char ID[50]; int cc; float **varScore; float totScore[MAXPARAMS*2+1]; };
typedef struct subject_t subject;

float getTStat();
extern int powell(double **pvar,int n,float ftol,float *fret,float (*func)(void));

// globals
subject *sub;
param par[MAXPARAMS];
int testTrain[2][MAXSUB],nVarType,nGeneSet,nSub,whichTotScore,paramToFit[MAXPARAMS],nParamToFit,tt;
double fittedPar[MAXPARAMS];
#define LONGLINELENGTH 20000
char line[LONGLINELENGTH+1],rest[LONGLINELENGTH+1];
FILE *resultsFile;

class fsParams { 
public:
	int fit,nReadScoresFiles;
	char readParamsFileName[200],writeParamsFileName[200],readScoresFileName[MAXSCORESTOREAD][200],writeScoresFileName[200],varScoresFileName[200],ttestFileName[200],testTrainFileName[200];
	int readArgs(int argc,char *argv[]);
	int getNextArg(char *nextArg,int argc,char *argv[],FILE **fpp,int *argNum);
	int check();
};

#define isArgType(a) (a[0]=='-' && a[1]=='-')
#define FILLARG(str) (strcmp(arg,str) ? 0 : ((getNextArg(arg, argc, argv, &fp, &argNum) && !isArgType(arg)) ? 1 : (dcerror(1,"No value provided for argument: %s\n",str), 0)))

int fsParams::readArgs(int argc,char *argv[])
{
	char arg[2000],*ptr;
	int argNum;
	FILE *fp;
	fp=0;
	argNum=1;
	readParamsFileName[0]=writeParamsFileName[0]=writeScoresFileName[0]=varScoresFileName[0]=ttestFileName[0]=testTrainFileName[0]='\0';
	fit=nReadScoresFiles=nVarType=nGeneSet=nSub=0;
	while (getNextArg(arg,argc,argv,&fp,&argNum))
	{
		if (!isArgType(arg))
		{
			dcerror(1,"Expected argument beginning -- but got this: %s\n",arg);
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
		else if (FILLARG("--var-scores-file"))
		{
			strcpy(varScoresFileName,arg);
		}
		else if (FILLARG("--read-params-file"))
		{
			strcpy(readParamsFileName,arg);
		}
		else if (FILLARG("--write-params-file"))
		{
			strcpy(writeParamsFileName,arg);
		}
		else if (FILLARG("--write-scores-file"))
		{
			strcpy(writeScoresFileName,arg);
		}
		else if (FILLARG("--read-scores-file"))
		{
			assert(nReadScoresFiles<MAXSCORESTOREAD);
			strcpy(readScoresFileName[nReadScoresFiles++],arg);
		}
		else if (FILLARG("--ttest-file"))
		{
			strcpy(ttestFileName,arg);
		}
		else if (FILLARG("--test-train-file"))
		{
			strcpy(testTrainFileName,arg);
		}
		else if (FILLARG("--fit"))
		{
			fit=atoi(arg);
		}
		else if (FILLARG("--n-var-type"))
		{
			nVarType=atoi(arg);
		}
		else if (FILLARG("--n-gene-set"))
		{
			nGeneSet=atoi(arg);
		}
		else
			dcerror(1,"Did not recognise argument specifier %s\n",arg);
	}
	// do checks
	return 1;
}

int fsParams::getNextArg(char *nextArg,int argc,char *argv[],FILE **fpp,int *argNum)
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

int fsParams::check()
{
	if (readParamsFileName[0]=='\0' || testTrainFileName[0]=='\0' || varScoresFileName[0]=='\0')
	{
		dcerror(1,"Must provide values for --read-params-file, --test-train-file and --var-scores-file"); return 0;
	}
	if (nVarType==0 || nGeneSet==0)
	{
		dcerror(1,"Must provide values for --n-gene-set and --n-var-type"); return 0;
	}
	return 1;
}

int readParams(fsParams *fs,FILE *readParamsFile,param *par,int nGeneSet,int nVarType)
{
	int p;
	for (p=0;p<nGeneSet;++p)
		if (!fgets(line,200,readParamsFile) || sscanf(line,"%f %s %d",&par[p].val,par[p].name,&par[p].toFit)!=3)
		{
			dcerror(1,"Could not read all geneList values\n"); return 0;
		}
	for (p=0;p<nVarType;++p)
		if (!fgets(line,200,readParamsFile) || sscanf(line,"%f %s %d",&par[nGeneSet+p].val,par[nGeneSet+p].name,&par[nGeneSet+p].toFit)!=3)
		{
			dcerror(1,"Could not read all varType values\n"); return 0;
		}
	return 1;
}

int writeParams(fsParams *fs,FILE *writeParamsFile,param *par,int nGeneSet,int nVarType)
{
	int p,savedNParamToFit;
	float tStat,savedVal;
	savedNParamToFit=nParamToFit;
	nParamToFit=0;
	for (p=0;p<nGeneSet+nVarType;++p)
	{
		fprintf(writeParamsFile,"%.2f\t%s\t%d\t",par[p].val,par[p].name,par[p].toFit);
		tStat=getTStat();
		fprintf(writeParamsFile,"%.4f\t",-tStat);
		savedVal=par[p].val;
		if (fabs(savedVal>1.0))
		{\
			par[p].val=savedVal*0.5;
			tStat=getTStat();
			fprintf(writeParamsFile,"%.4f\t",-tStat);
			par[p].val=savedVal*1.5;
			tStat=getTStat();
			fprintf(writeParamsFile,"%.4f\t",-tStat);
		}
		else
		{
			par[p].val=savedVal-0.5;
			tStat=getTStat();
			fprintf(writeParamsFile,"%.4f\t",-tStat);
			par[p].val=savedVal+0.5;
			tStat=getTStat();
			fprintf(writeParamsFile,"%.4f\t",-tStat);
		}
		par[p].val=savedVal;
		fprintf(writeParamsFile,"\n");
	}
	nParamToFit=savedNParamToFit;
	return 1;
}

int readTestTrain(fsParams *fs,FILE *testTrainFile,int testTrain[2][MAXSUB])
{
	int s;
	for (s=0;fgets(line,200,testTrainFile) && sscanf(line,"%d %d",&testTrain[0][s],&testTrain[1][s])==2;++s)
		;
	return s;
}

int readVarScores(fsParams *fs,FILE *varScoresFile,subject *sub,int nSub,int nGeneSet,int nVarType)
{
	int s,t,v;
	if (!fgets(line,LONGLINELENGTH,varScoresFile)) // dump first line
		{ dcerror(1,"Could not read enough lines from --var-scores-file\n"); return 1; }
#if 1
	for (s=0;s<nSub;++s)
	{
		if (fscanf(varScoresFile,"%s %d ",sub[s].ID,&sub[s].cc)!=2)
		{
			dcerror(1,"Incomplete data in --var-scores-file, got to s=%d t=%d v=%d sub[s].ID=%s\n",s,t,v,sub[s].ID); return 1;
		}
		for (t=0;t<nGeneSet;++t)
			for (v=0;v<nVarType;++v)
				if (fscanf(varScoresFile,"%f ",&sub[s].varScore[t][v])<1)
				{
					dcerror(1,"Incomplete data in --var-scores-file, got to s=%d t=%d v=%d sub[s].ID=%s\n",s,t,v,sub[s].ID); return 1;
				}
		
	}
#else
	for (s=0;s<nSub;++s)
	{
		if (!fgets(line,LONGLINELENGTH,varScoresFile))
			{ dcerror(1,"Could not read enough lines from --var-scores-file\n"); return 1; }
		if (sscanf(line,"%s %d %[^\n]",sub[s].ID,&sub[s].cc,rest)!=3)
			{ dcerror(1,"Incomplete line in --var-scores-file:\n%s\n",line); return 1; }
		strcpy(line,rest);
		for (t=0;t<nGeneSet;++t)
			for (v=0;v<nVarType;++v)
				if (sscanf(line,"%f %[^\n]",&sub[s].varScore[t][v],rest)<1)
					{ dcerror(1,"Incomplete line in --var-scores-file:\n%s\n",line); return 1; }
				else
					strcpy(line,rest);

	}
#endif
	return 1;
}

float getTStat()
{
	int p,s,t,v,i,cc,n[2];
	float score,sigma_x[2],sigma_x2[2],mean[2],var[2],SE,tval,s2;
	for (p=0;p<nParamToFit;++p)
		par[paramToFit[p]].val=fittedPar[p];
	for (i=0;i<2;++i)
		sigma_x[i]=sigma_x2[i]=n[i]=0;
	for (s=0;s<nSub;++s)
	{
		if (testTrain[tt][s]==0)
			continue;
		score=0;
		for (t=0;t<nGeneSet;++t)
			for (v=0;v<nVarType;++v)
				score+=sub[s].varScore[t][v]*par[t].val*par[nGeneSet+v].val;
		sub[s].totScore[whichTotScore]=score;
		cc=sub[s].cc;
		++n[cc];
		sigma_x[cc]+=score;
		sigma_x2[cc]+=score*score;
	}
	for (i = 0; i < 2; ++i)
	{
		var[i] = (sigma_x2[i] - sigma_x[i] * sigma_x[i] / n[i]) / (n[i] - 1);
		mean[i] = sigma_x[i] / n[i];
	}
	s2=((n[0]-1)*var[0]+(n[1]-1)*var[1])/(n[0]+n[1]-2);
	SE=sqrt(s2*(1/(float)n[0]+1/(float)n[1]));
	if (SE==0)
		tval=0;
	else
		tval=(mean[1]-mean[0])/SE;
	if (resultsFile)
	{
		fprintf(resultsFile,
				"             Controls  Cases     \n"
				"N            %9d %9d\n"
				"Mean score   %9.3f %9.3f\n"
				"t (%d df) = %6.3f\n",
				n[0],n[1],mean[0],mean[1],n[0] + n[1] - 2,tval);

	}
	return -tval;
	// plan to get normalised values of scores before writing them to a file so they can be combined later
}

int writeScores(fsParams *fs,FILE *writeScoresFile)
{
	int s,n;
	float sigma_x,sigma_x2,mean,var,score,sd;
	whichTotScore=0;
	tt=1;
	getTStat();
	sigma_x=sigma_x2=n=0;
	for (s=0;s<nSub;++s)
	{
		if (testTrain[tt][s]==0)
			continue;
		score=sub[s].totScore[0];
		sigma_x+=score;
		sigma_x2+=score*score;
		++n;
	}
	mean=sigma_x/n;
	var=(n*sigma_x2-sigma_x*sigma_x)/(n*(n-1));
	sd=sqrt(var);
	for (s=0;s<nSub;++s)
	{
		if (testTrain[tt][s]==0)
			continue;
		fprintf(writeScoresFile,"%s %d %.2f %f\n",sub[s].ID,sub[s].cc,sub[s].totScore[0],(sub[s].totScore[0]-mean)/sd);
	}
	return 1;
}

int processOption(fsParams &fs,char *option,char *value)
{
	FILE *readParamsFile,*testTrainFile,*varScoresFile,*writeParamsFile,*writeScoresFile;
	int s,t,p;
	float tval;
	double *toFitPtr[MAXPARAMS];
	if (strncmp(option,"--",2))
		{ dcerror(1,"Option should start with -- but had this:%s\n",option); return 0; }
	if (!strcmp(option,"-read-params-file"))
	{
		strcpy(fs.readParamsFileName,value);
		readParamsFile=fopen(fs.readParamsFileName,"r");
		if (!readParamsFile)
		{
			dcerror(1,"Could not open %s\n",fs.readParamsFileName); exit(1);
		}
		if (!readParams(&fs,readParamsFile,par,nGeneSet,nVarType))
		{
			dcerror(1,"Could not read parameter starting values from %s\n",fs.readParamsFileName); exit(1);
		}
		fclose(readParamsFile);
	}
	else if (!strcmp(option,"--test-train-file"))
	{
		strcpy(fs.testTrainFileName,value);
		testTrainFile=fopen(fs.testTrainFileName,"r");
		if (!testTrainFile)
		{
			dcerror(1,"Could not open %s\n",fs.testTrainFileName); exit(1);
		}
		if ((nSub=readTestTrain(&fs,testTrainFile,testTrain))==0)
		{
			dcerror(1,"Could not read parameter starting values from %s\n",fs.testTrainFileName); exit(1);
		}
		fclose(testTrainFile);
	}
	else if (!strcmp(option,"--var-scores-file"))
	{
		strcpy(fs.varScoresFileName,value);
		varScoresFile=fopen(fs.varScoresFileName,"r");
		if (!varScoresFile)
		{
			dcerror(1,"Could not open %s\n",fs.varScoresFileName); exit(1);
		}
		if (!readVarScores(&fs,varScoresFile,sub,nSub,nGeneSet,nVarType))
		{
			dcerror(1,"Problem reading %s\n",fs.varScoresFileName); exit(1);
		}
		fclose(varScoresFile);
	}
	else if (!strcmp(option,"--fit"))
	{ 
		fs.fit=atoi(value);
		if (fs.fit)
		{
			tt=0;
			for (p=0;p<nParamToFit;++p)
				toFitPtr[p]=&fittedPar[p];
			powell(toFitPtr,nParamToFit,0.00001,&tval,getTStat);
		}
	}
	else if (!strcmp(option,"--ttest-file"))
	{ 
		strcpy(fs.ttestFileName,value);
		if (fs.ttestFileName[0])
		{
			resultsFile=fopen(fs.ttestFileName,"w");
			if (!resultsFile)
			{
				dcerror(1,"Could not open %s\n",fs.ttestFileName); exit(1);
			}
			tt=1;
			tval=-getTStat();
			fclose(resultsFile);
			resultsFile=0;
		}
	}
	else if (!strcmp(option,"--write-params-file"))
	{
		strcpy(fs.writeParamsFileName,value);
		if (fs.writeParamsFileName[0])
		{
			tt=1;
			writeParamsFile=fopen(fs.writeParamsFileName,"w");
			if (!writeParamsFile)
			{
				dcerror(1,"Could not open %s\n",fs.writeParamsFileName); exit(1);
			}
			if (!writeParams(&fs,writeParamsFile,par,nGeneSet,nVarType))
			{
				dcerror(1,"Could not write parameter values to %s\n",fs.writeParamsFileName); exit(1);
			}
			fclose(writeParamsFile);
		}

	}
	else if (!strcmp(option,"--write-scores-file"))
	{
		strcpy(fs.writeScoresFileName,value);
		if (fs.writeScoresFileName[0])
		{
			writeScoresFile=fopen(fs.writeScoresFileName,"w");
			if (!writeScoresFile)
			{
				dcerror(1,"Could not open %s\n",fs.writeScoresFileName); exit(1);
			}
			if (!writeScores(&fs,writeScoresFile))
			{
				dcerror(1,"Could not write scores to %s\n",fs.writeScoresFileName); exit(1);
			}
			fclose(writeScoresFile);
		}
	}
	else
		{ dcerror(1,"Unrecognised option: %s\n",option); return 0; }
	return 1;
}

int main(int argc,char *argv[])
{
	fsParams fs;
	FILE *readParamsFile,*testTrainFile,*varScoresFile,*writeParamsFile,*writeScoresFile;
	int s,t,p;
	float tval;
	double *toFitPtr[MAXPARAMS];
	printf("%s v%s\n",PROGRAM,FSVERSION);
	if (!fs.readArgs(argc,argv))
		exit(1);
	if (!fs.check())
		exit(1);
	readParamsFile=fopen(fs.readParamsFileName,"r");
	if (!readParamsFile)
	{
		dcerror(1,"Could not open %s\n",fs.readParamsFileName); exit(1);
	}
	if (!readParams(&fs,readParamsFile,par,nGeneSet,nVarType))
	{
		dcerror(1,"Could not read parameter starting values from %s\n",fs.readParamsFileName); exit(1);
	}
	fclose(readParamsFile);
	testTrainFile=fopen(fs.testTrainFileName,"r");
	if (!testTrainFile)
	{
		dcerror(1,"Could not open %s\n",fs.testTrainFileName); exit(1);
	}
	if ((nSub=readTestTrain(&fs,testTrainFile,testTrain))==0)
	{
		dcerror(1,"Could not read parameter starting values from %s\n",fs.testTrainFileName); exit(1);
	}
	fclose(testTrainFile);
	assert((sub=(subject *)calloc(nSub,sizeof(subject)))!=0);
	for (s=0;s<nSub;++s)
	{
		assert((sub[s].varScore=(float**)calloc(nGeneSet,sizeof(float*)))!=0);
		for (t=0;t<nGeneSet;++t)
			assert((sub[s].varScore[t]=(float*)calloc(nVarType,sizeof(float)))!=0);
	}
	varScoresFile=fopen(fs.varScoresFileName,"r");
	if (!varScoresFile)
	{
		dcerror(1,"Could not open %s\n",fs.varScoresFileName); exit(1);
	}
	if (!readVarScores(&fs,varScoresFile,sub,nSub,nGeneSet,nVarType))
		{ dcerror(1,"Problem reading %s\n",fs.varScoresFileName); exit(1); }
	fclose(varScoresFile);
	for (p=0,nParamToFit=0;p<nGeneSet+nVarType;++p)
		if (par[p].toFit)
		{
			paramToFit[nParamToFit]=p;
			fittedPar[nParamToFit]=par[p].val;
			++nParamToFit;
		}
	if (fs.fit)
	{
		tt=0;
		for (p=0;p<nParamToFit;++p)
			toFitPtr[p]=&fittedPar[p];
		powell(toFitPtr,nParamToFit,0.00001,&tval,getTStat);
	}
	if (fs.ttestFileName[0])
	{
		resultsFile=fopen(fs.ttestFileName,"w");
		if (!resultsFile)
		{
			dcerror(1,"Could not open %s\n",fs.ttestFileName); exit(1);
		}
		tt=1;
		tval=-getTStat();
		fclose(resultsFile);
		resultsFile=0;
	}
	if (fs.writeParamsFileName[0])
	{
		tt=1;
		writeParamsFile=fopen(fs.writeParamsFileName,"w");
		if (!writeParamsFile)
		{
			dcerror(1,"Could not open %s\n",fs.writeParamsFileName); exit(1);
		}
		if (!writeParams(&fs,writeParamsFile,par,nGeneSet,nVarType))
		{
			dcerror(1,"Could not write parameter values to %s\n",fs.writeParamsFileName); exit(1);
		}
		fclose(writeParamsFile);
	}
	if (fs.writeScoresFileName[0])
	{
		writeScoresFile=fopen(fs.writeScoresFileName,"w");
		if (!writeScoresFile)
		{
			dcerror(1,"Could not open %s\n",fs.writeScoresFileName); exit(1);
		}
		if (!writeScores(&fs,writeScoresFile))
		{
			dcerror(1,"Could not write scores to %s\n",fs.writeScoresFileName); exit(1);
		}
		fclose(writeScoresFile);
	}
	return 0;
}