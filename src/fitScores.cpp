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

using namespace dlib;
typedef matrix<double,0,1> column_vector;

#include "dcerror.hpp"

#define PROGRAM "fitScores"
#define FSVERSION "1.0"
#define MAXSUB 15000
#define MAXSCORESTOREAD 20
#define MAXPARAMS 100

struct param_t { double val; char name[200]; int toFit; };
typedef struct param_t param;


struct subject_t { char ID[50]; int cc; float **varScore; float totScore[MAXPARAMS*2+1]; };
typedef struct subject_t subject;

double getTStat();

float float_getTStat() { return (float)getTStat(); }
#ifdef CANUSENRPOWELL
extern int NR_float_powell(double **pvar,int n,float ftol,float *fret,float (*func)(void));
#endif


int NR_powell(double **pvar,int n,float ftol,float *fret,double (*func)(void))
{
#ifdef CANUSENRPOWELL
	return NR_float_powell(pvar,n,ftol,fret,float_getTStat);
#else
	dcerror(1,"NR_powell() is not implemented - use only dlib_powell()");
	exit(1);
	return 0;
#endif
}

int NR_powell(double **pvar,int n,float ftol,float *fret,double (*func)(void));

int (*powell)(double **pvar,int n,float ftol,float *fret,double (*func)(void));

// globals
subject *sub;
param par[MAXPARAMS];
int testTrain[2][MAXSUB],nVarType,nGeneSet,nSub,whichTotScore,paramToFit[MAXPARAMS],nParamToFit,tt;
double fittedPar[MAXPARAMS];
#define LONGLINELENGTH 20000
char fsline[LONGLINELENGTH+1],rest[LONGLINELENGTH+1];
FILE *resultsFile,*logFile;

class powell_function
{
public:
	// I don't think I need to keep any state data
	double operator() (const column_vector& arg) const ;
	double (*func)(void);
	powell_function(double (*f)(void)) { func=f; }
};

// may have failed because function was returning float, not double
int dlib_powell(double **pvar,int n,float ftol,float *fret,double (*func)(void))
{
	int i;
	double d;
	column_vector starting_point(n);
	for (i=0;i<n;++i)
		starting_point(i,0)=*pvar[i];
	// I'm not sure if below will work
	d=find_min_using_approximate_derivatives(cg_search_strategy(),
		objective_delta_stop_strategy(ftol),
		powell_function(func),starting_point,-20);
	for (i=0;i<n;++i)
		*pvar[i]=starting_point(i,0);
	*fret=d;
	return 1;
}

double powell_function::operator() (const column_vector& c) const
{
	int i;
	double f;
	for (i=0;i<c.nr();++i)
		fittedPar[i]=c(i,0);
	f=func();
	return f;
}

class fsParams { 
public:
	int fit,nReadScoresFiles,doGrid,shouldSaveT;
	float stepwiseTol,powellTol,savedMaxT,scaleMean;
	char readParamsFileName[200],writeParamsFileName[200],readScoresFileName[MAXSCORESTOREAD][200],extraArgsFileName[200],
		writeScoresFileName[200],varScoresFileName[200],ttestFileName[200],testTrainFileName[200],logFileName[200];
	int readArgs(int argc,char *argv[]);
	int getNextArg(char *nextArg,int argc,char *argv[],FILE **fpp,int *argNum);
	int check();
};

int processOption(fsParams &fs,char *option,char *value);

#define isArgType(a) (a[0]=='-' && a[1]=='-')
#define FILLARG(str) (strcmp(arg,str) ? 0 : ((getNextArg(arg, argc, argv, &fp, &argNum) && !isArgType(arg)) ? 1 : (dcerror(1,"No value provided for argument: %s\n",str), 0)))

int fsParams::readArgs(int argc,char *argv[])
{
	char arg[2000],*ptr;
	int argNum;
	FILE *fp;
	fp=0;
	argNum=1;
	stepwiseTol=0;
	powellTol=0.00001; // 0.001 and 0.0001 did not fit well tval-=0.4
	readParamsFileName[0]=writeParamsFileName[0]=writeScoresFileName[0]=varScoresFileName[0]=ttestFileName[0]=testTrainFileName[0]=logFileName[0]=extraArgsFileName[0]='\0';
	fit=nReadScoresFiles=nVarType=nGeneSet=nSub=0;
	powell=dlib_powell;
	scaleMean=10;
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
			processOption(*this,"--var-scores-file",arg);
		else if (FILLARG("--test-train-file"))
			processOption(*this,"--test-train-file",arg);
		else if (FILLARG("--read-params-file"))
			processOption(*this,"--read-params-file",arg);
		else if (FILLARG("--write-params-file"))
			processOption(*this,"--write-params-file",arg);
		else if (FILLARG("--write-scores-file"))
			processOption(*this,"--write-scores-file",arg);
		else if (FILLARG("--read-scores-file"))
			processOption(*this,"--read-scores-file",arg);
		else if (FILLARG("--ttest-file"))
			processOption(*this,"--ttest-file",arg);
		else if (FILLARG("--extra-arg-file"))
			processOption(*this,"--extra-arg-file",arg);
		else if (FILLARG("--log-file"))
			processOption(*this,"--log-file",arg);
		else if (FILLARG("--fit"))
			processOption(*this,"--fit",arg);
		else if (FILLARG("--n-var-type"))
			processOption(*this,"--n-var-type",arg);
		else if (FILLARG("--n-gene-set"))
			processOption(*this,"--n-gene-set",arg);
		else if (FILLARG("--stepwise-tol"))
			processOption(*this,"--stepwise-tol",arg);
		else if (FILLARG("--do-grid"))
			processOption(*this,"--do-grid",arg);
		else if (FILLARG("--use-NR-powell"))
			processOption(*this,"--use-NR-powell",arg);
		else if (FILLARG("--test-train"))
			processOption(*this,"--save-max-t",arg);
		else if (FILLARG("--save-max-t"))
			processOption(*this,"--save-max-t",arg);
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
		if (!fgets(fsline,200,readParamsFile) || sscanf(fsline,"%lf %s %d",&par[p].val,par[p].name,&par[p].toFit)!=3)
		{
			dcerror(1,"Could not read all geneList values\n"); return 0;
		}
	for (p=0;p<nVarType;++p)
		if (!fgets(fsline,200,readParamsFile) || sscanf(fsline,"%lf %s %d",&par[nGeneSet+p].val,par[nGeneSet+p].name,&par[nGeneSet+p].toFit)!=3)
		{
			dcerror(1,"Could not read all varType values\n"); return 0;
		}
	return 1;
}

int writeParams(fsParams *fs,FILE *writeParamsFile,param *par,int nGeneSet,int nVarType)
{
	int p,savedNParamToFit,r,side,steps,nFitted,doNotScale;
	float tStat,savedVal,ratio,val,CL,savedTStat,step,scaleFactor,tot,diff;
	savedNParamToFit=nParamToFit;
	nParamToFit=0;
	if (fs->scaleMean!=0)
	{
		doNotScale=nFitted=tot=0;
		for (p=0;p<nGeneSet;++p)
			if (par[p].toFit)
			{
				tot+=fabs(par[p].val);
				++nFitted;
			}
			else
				if (par[p].val!=0)
					doNotScale=1;
		if (doNotScale==0)
		{
			scaleFactor=fs->scaleMean/(tot/nFitted);
			for (p=0;p<nGeneSet;++p)
				par[p].val*=scaleFactor;
		}
		doNotScale=nFitted=tot=0;
		for (p=nGeneSet;p<nGeneSet+nVarType;++p)
			if (par[p].toFit)
			{
				tot+=fabs(par[p].val);
				++nFitted;
			}
			else
				if (par[p].val!=0)
					doNotScale=1;
		if (doNotScale==0)
		{
			scaleFactor=fs->scaleMean/(tot/nFitted);
			for (p=nGeneSet;p<nGeneSet+nVarType;++p)
				par[p].val*=scaleFactor;
		}
	}
	fprintf(writeParamsFile,"value\tname\ttoFit\t");
	if (fs->doGrid)
	{
		fprintf(writeParamsFile,"tval\tlolim\thilim\t");
		for (val=-100;val<-0.05;val/=10)
			for (ratio=1.0;ratio>=0.05;ratio-=0.1)
				fprintf(writeParamsFile,"%.2f\t",val*ratio);
		fprintf(writeParamsFile,"0.0\t");
		for (val=0.1;val<500;val*=10)
			for (ratio=0.1;ratio<1.05;ratio+=0.1)
				fprintf(writeParamsFile,"%.1f\t",val*ratio);
	}
	fprintf(writeParamsFile,"\n");
	for (p=0;p<nGeneSet+nVarType;++p)
	{
		fprintf(writeParamsFile,"%.2f\t%s\t%d\t",par[p].val,par[p].name,par[p].toFit);
		if (fs->doGrid)
		{
			tStat=-getTStat();
			fprintf(writeParamsFile,"%.4f\t",tStat);
			savedVal=par[p].val;
			// savedTStat=tStat;
			savedTStat=fs->savedMaxT;
			diff=1;
			for (side=-1;side<=1;side+=2)
			{
				CL=0;
				for (step=100;step>=0.05;step/=10.0)
				{
					for (steps=0;steps<=10;++steps)
					{
						par[p].val=CL*side+savedVal;
						tStat=-getTStat();
						if (tStat<savedTStat-diff)
						{
							CL-=step;
						}
						else
							CL+=step;
					}
					if (CL==1000)
						break;
				}
				fprintf(writeParamsFile,"%.2f\t",CL*side+savedVal);
			}

			for (val=-100;val<-0.05;val/=10)
				for (ratio=1.0;ratio>=0.05;ratio-=0.1)
				{
					par[p].val=val*ratio;
					tStat=getTStat();
					fprintf(writeParamsFile,"%.4f\t",-tStat);
				}
			par[p].val=0.0;
			tStat=getTStat();
			fprintf(writeParamsFile,"%.4f\t",-tStat);
			for (val=0.1;val<500;val*=10)
				for (ratio=0.1;ratio<1.05;ratio+=0.1)
				{
					par[p].val=val*ratio;
					tStat=getTStat();
					fprintf(writeParamsFile,"%.4f\t",-tStat);
				}
			par[p].val=savedVal;
		}
		fprintf(writeParamsFile,"\n");
	}
	nParamToFit=savedNParamToFit;
	return 1;
}

int readTestTrain(fsParams *fs,FILE *testTrainFile,int testTrain[2][MAXSUB])
{
	int s;
	for (s=0;fgets(fsline,200,testTrainFile) && sscanf(fsline,"%d %d",&testTrain[0][s],&testTrain[1][s])==2;++s)
		;
	return s;
}

int readVarScores(fsParams *fs,FILE *varScoresFile,subject *sub,int nSub,int nGeneSet,int nVarType)
{
	int s,t,v;
	if (!fgets(fsline,LONGLINELENGTH,varScoresFile)) // dump first line
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
		if (!fgets(fsline,LONGLINELENGTH,varScoresFile))
			{ dcerror(1,"Could not read enough lines from --var-scores-file\n"); return 1; }
		if (sscanf(fsline,"%s %d %[^\n]",sub[s].ID,&sub[s].cc,rest)!=3)
			{ dcerror(1,"Incomplete line in --var-scores-file:\n%s\n",fsline); return 1; }
		strcpy(fsline,rest);
		for (t=0;t<nGeneSet;++t)
			for (v=0;v<nVarType;++v)
				if (sscanf(fsline,"%f %[^\n]",&sub[s].varScore[t][v],rest)<1)
					{ dcerror(1,"Incomplete line in --var-scores-file:\n%s\n",fsline); return 1; }
				else
					strcpy(fsline,rest);

	}
#endif
	return 1;
}

double getTStat()
{
	int p,s,t,v,i,cc,n[2];
	double score,sigma_x[2],sigma_x2[2],mean[2],var[2],SE,tval,s2;
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
	if (logFile)
	{
		fprintf(logFile,"%.3f\t",tval);
		for (p=0;p<nParamToFit;++p)
			fprintf(logFile,"%.2f\t",fittedPar[p]);
		fprintf(logFile,"\n");
	}
	return -tval;
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

int stepwise(fsParams *fs)
{
	float tval,highTval,savedTval,savedVal;
	double *toFitPtr[MAXPARAMS];
	int p,pp,savedToFit,temp,highestP,localToFit[MAXPARAMS],localNParamToFit;
	FILE *myLog;
	myLog=logFile;
	logFile=0;
	tt=0;
	for (p=0,nParamToFit=0;p<nGeneSet+nVarType;++p)
		if (par[p].toFit)
		{
			paramToFit[nParamToFit]=p;
			fittedPar[nParamToFit]=par[p].val;
			++nParamToFit;
		}
	for (p=0;p<nParamToFit;++p)
		toFitPtr[p]=&fittedPar[p];
	powell(toFitPtr,nParamToFit,fs->powellTol,&savedTval,getTStat);
	//for (p=0;p<nGeneSet+nVarType;++p)
		//savedToFit[p]=par[p].toFit;
	savedTval=savedTval*-1;
	if (fs->shouldSaveT)
		fs->savedMaxT=savedTval;
	if (myLog)
		fprintf(myLog,"Original t = %.4f\n",savedTval);
	do
	{
		highTval=-1000;
		highestP=-1;
		if (myLog)
		{
			for (p=0;p<nParamToFit;++p)
				fprintf(myLog,"%s\t",par[paramToFit[p]].name);
			fprintf(myLog,"\n");
			for (p=0;p<nParamToFit;++p)
				fprintf(myLog,"%.2f\t",par[paramToFit[p]].val);
			fprintf(myLog,"\n");
		}
		localNParamToFit=nParamToFit;
		for (pp=0;pp<localNParamToFit;++pp)
			localToFit[pp]=paramToFit[pp];
		for (pp=0;pp<localNParamToFit;++pp) 
		{
			if (pp!=0)
			{
				par[savedToFit].val=savedVal;
				par[savedToFit].toFit=1;
			}
			savedToFit=localToFit[pp];
			savedVal=par[savedToFit].val;
			par[savedToFit].val=0;
			par[savedToFit].toFit=0;
			for (p=0,nParamToFit=0;p<nGeneSet+nVarType;++p)
				if (par[p].toFit)
				{
					paramToFit[nParamToFit]=p;
					fittedPar[nParamToFit]=par[p].val;
					++nParamToFit;
				}
			for (p=0;p<nParamToFit;++p)
				toFitPtr[p]=&fittedPar[p];
			tval=getTStat();
			tval=tval*-1;
			if (tval>highTval)
			{
				highTval=tval;
				highestP=localToFit[pp];
			}
			if (myLog)
				fprintf(myLog,"%.2f\t",tval);
		}
		par[savedToFit].val=savedVal;
		par[savedToFit].toFit=1;
		if (savedTval-highTval>fs->stepwiseTol)
			break;
		else
		{
			par[highestP].val=0;
			par[highestP].toFit=0;
			for (p=0,nParamToFit=0;p<nGeneSet+nVarType;++p)
				if (par[p].toFit)
				{
					paramToFit[nParamToFit]=p;
					fittedPar[nParamToFit]=par[p].val;
					++nParamToFit;
				}
			for (p=0;p<nParamToFit;++p)
				toFitPtr[p]=&fittedPar[p];
			powell(toFitPtr,nParamToFit,fs->powellTol,&tval,getTStat);
		}

	} while (nParamToFit>0);
	for (p=0,nParamToFit=0;p<nGeneSet+nVarType;++p)
		if (par[p].toFit)
		{
			paramToFit[nParamToFit]=p;
			fittedPar[nParamToFit]=par[p].val;
			++nParamToFit;
		}
	for (p=0;p<nParamToFit;++p)
		toFitPtr[p]=&fittedPar[p];
	for (p=0;p<nParamToFit;++p)
		par[paramToFit[p]].val=1;
	powell(toFitPtr,nParamToFit,fs->powellTol,&tval,getTStat);
	if (myLog)
	{
		fprintf(myLog,"\n\nFinished stepwise exclusion, t = %.2f\n",-tval);
		for (p=0;p<nParamToFit;++p)
			fprintf(myLog,"%s\t",par[paramToFit[p]].name);
		fprintf(myLog,"\n");
		for (p=0;p<nParamToFit;++p)
			fprintf(myLog,"%.2f\t",par[paramToFit[p]].val);
		fprintf(myLog,"\n");
	}
	logFile=myLog;
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
	if (!strcmp(option,"--var-scores-file"))
	{
		strcpy(fs.varScoresFileName,value);
		varScoresFile=fopen(fs.varScoresFileName,"r");
		if (!varScoresFile)
		{
			dcerror(1,"Could not open %s\n",fs.varScoresFileName); exit(1);
		}
		assert((sub=(subject *)calloc(nSub,sizeof(subject)))!=0); // assume this will only be called once
		for (s=0;s<nSub;++s)
		{
			assert((sub[s].varScore=(float**)calloc(nGeneSet,sizeof(float*)))!=0);
			for (t=0;t<nGeneSet;++t)
				assert((sub[s].varScore[t]=(float*)calloc(nVarType,sizeof(float)))!=0);
		}
		if (!readVarScores(&fs,varScoresFile,sub,nSub,nGeneSet,nVarType))
		{
			dcerror(1,"Problem reading %s\n",fs.varScoresFileName); exit(1);
		}
		fclose(varScoresFile);
	}
	else if (!strcmp(option,"--read-params-file"))
	{
		strcpy(fs.readParamsFileName,value);
		readParamsFile=fopen(fs.readParamsFileName,"r");
		if (!readParamsFile)
		{
			dcerror(1,"Could not open params file %s\n",fs.readParamsFileName); exit(1);
		}
		if (!readParams(&fs,readParamsFile,par,nGeneSet,nVarType))
		{
			dcerror(1,"Could not read parameter starting values from %s\n",fs.readParamsFileName); exit(1);
		}
		fclose(readParamsFile);
		for (p=0,nParamToFit=0;p<nGeneSet+nVarType;++p)
			if (par[p].toFit)
			{
				paramToFit[nParamToFit]=p;
				fittedPar[nParamToFit]=par[p].val;
				++nParamToFit;
			}
	}
	else if (!strcmp(option,"--test-train-file"))
	{
		strcpy(fs.testTrainFileName,value);
		testTrainFile=fopen(fs.testTrainFileName,"r");
		if (!testTrainFile)
		{
			dcerror(1,"Could not open test/train file %s\n",fs.testTrainFileName); exit(1);
		}
		if ((nSub=readTestTrain(&fs,testTrainFile,testTrain))==0)
		{
			dcerror(1,"Could not read parameter starting values from %s\n",fs.testTrainFileName); exit(1);
		}
		fclose(testTrainFile);
	}
	else if (!strcmp(option,"--log-file"))
	{
		strcpy(fs.logFileName,value);
		if (!strcmp(fs.logFileName,"close"))
		{
			if (logFile)
			{
				fclose(logFile);
				logFile=0;
			}
		}
		else
		{
			logFile=fopen(fs.logFileName,"w");
			if (!logFile)
			{
				dcerror(1,"Could not open %s\n",fs.logFileName); exit(1);
			}
		}

	}
	else if (!strcmp(option,"--fit"))
	{
		fs.fit=atoi(value);
		if (fs.fit)
		{
			if (logFile)
			{
				fprintf(logFile,"tVal\t");
				for (p=0;p<nParamToFit;++p)
					fprintf(logFile,"%s\t",par[paramToFit[p]].name);
				fprintf(logFile,"\n");
			}
			tt=0;
			for (p=0;p<nParamToFit;++p)
				toFitPtr[p]=&fittedPar[p];
			powell(toFitPtr,nParamToFit,fs.powellTol,&tval,getTStat);
			if (fs.shouldSaveT)
				fs.savedMaxT=-tval;
		}
	}
	else if (!strcmp(option,"--n-var-type"))
	{
		nVarType=atoi(value);
	}
	else if (!strcmp(option,"--n-gene-set"))
	{
		nGeneSet=atoi(value);
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
			if (fs.shouldSaveT)
				fs.savedMaxT=tval;
			fclose(resultsFile);
			resultsFile=0;
		}
	}
	else if (!strcmp(option,"--write-params-file"))
	{
		strcpy(fs.writeParamsFileName,value);
		if (fs.writeParamsFileName[0])
		{
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
	else if (!strcmp(option,"--stepwise-tol"))
	{
		fs.stepwiseTol=atof(value);
		if (fs.stepwiseTol)
			stepwise(&fs);
	}
	else if (!strcmp(option,"--use-NR-powell"))
	{
		if (atoi(value))
			powell=NR_powell;
		else
			powell=dlib_powell;
	}
	else if (!strcmp(option,"--save-max-t"))
	{
		if (atoi(value))
			fs.shouldSaveT=1;
		else
			fs.shouldSaveT=0;
	}
	else if (!strcmp(option,"--test-train"))
	{
		if (atoi(value))
			tt=1;
		else
			tt=0;
	}
	else if (!strcmp(option,"--do-grid"))
		fs.doGrid=atoi(value);
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
	char arg[100];
	printf("%s v%s\n",PROGRAM,FSVERSION);
	if (!fs.readArgs(argc,argv))
		exit(1);
#if 0
	if (!fs.check())
		exit(1);
	if (!fs.readParamsFileName[0])
	{
		dcerror(1,"Must set value for --read-params-file"); exit(1);
	}
	processOption(fs,"--read-params-file",fs.readParamsFileName);
	if (!fs.testTrainFileName[0])
	{
		dcerror(1,"Must set value for --test-train-file"); exit(1);
	}
	processOption(fs,"--test-train-file",fs.testTrainFileName);
	if (!fs.varScoresFileName[0])
	{
		dcerror(1,"Must set value for --var-scores-file"); exit(1);
	}
	processOption(fs,"--var-scores-file",fs.varScoresFileName);
	if (fs.logFileName[0])
		processOption(fs,"--log-file",fs.logFileName);
	if (fs.fit)
	{
		sprintf(arg,"%d",fs.fit);
		processOption(fs,"--fit",arg);
	}
	if (fs.stepwiseTol)
	{
		sprintf(arg,"%f",fs.stepwiseTol);
		processOption(fs,"--stepwise-tol",arg);
	}
	if (fs.ttestFileName[0])
		processOption(fs,"--ttest-file",fs.ttestFileName);
	if (fs.writeParamsFileName[0])
		processOption(fs,"--write-params-file",fs.writeParamsFileName);
	if (fs.writeScoresFileName[0])
		processOption(fs,"--write-scores-file",fs.writeScoresFileName);
#endif
	return 0;
}
