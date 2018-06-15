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

#include <string>
#include <map>

typedef std::pair<std::string,int> TStrIntPair;
typedef std::map<std::string,int> TStrIntMap;

#define PROGRAM "permPathwayAssoc"
#define PPAVERSION "1.3"

#ifndef MAX_LOCI
#define MAX_LOCI 12000
#endif
#ifndef MAX_SUB
#define MAX_SUB 15000
#endif
#ifndef MAX_ID_LENGTH 
#define MAX_ID_LENGTH 35
#endif
#define LONGLINELENGTH 20000
#define MAX_PATHWAY_LENGTH 500
#define MAX_PATHWAYS 2000
#define MAX_GENES_IN_PATHWAY 2000
#define MAX_GENES 25000

#if 0
MSVC implements long with 32 bits but unix is 64 bits so use __int64
#endif
#include <inttypes.h>
// http://stackoverflow.com/questions/9225567/how-to-print-a-int64-t-type-in-c
#ifdef MSDOS
// i.e. MSVC
#define FILEOFFSET __int64
#define fseek _fseeki64
#define ftell _ftelli64
#else
#define FILEOFFSET long
#endif

//#define USEDCASSERT 1
#ifdef USEDCASSERT
// catch assertion failures so Windows does not do something clever with them
#ifdef  __cplusplus
extern "C" {
#endif

void __cdecl _assert(void *str, void *fn, unsigned line)
{
fprintf(stderr,"Assertion failed: %s, file %s, line %d\n",(char *)str,(char *)fn,line);
exit(1);
}

#ifdef  __cplusplus
}
#endif

#endif

class ppaParams {
public:
	char scoreFilePrefix[1000],scoreFileSuffix[1000],outputFilePrefix[1000],outputFileSuffix[1000],pathwayFileName[1000],scoreTableFileName[1000];
	std::map<std::string,FILEOFFSET> geneIndex;
	int nSub,nPathways,nGenes,nTopPathways,topPathways[MAX_PATHWAYS],nPerms,nOver[MAX_PATHWAYS];
	float geneLevelOutputThreshold,pathwayThreshold;
	int readParms(int argc,char *argv[]);
	int getNextArg(char *nextArg,int argc,char *argv[],FILE **fpp,int *argNum);
	FILE *summaryOutputFile,*SLPFile,*fullOutputFile,*scoreTableFile;
};

char line[LONGLINELENGTH+1],rest[LONGLINELENGTH+1];

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

int ppaParams::readParms(int argc, char *argv[])
{
	char arg[2000],*ptr;
	int argNum;
	FILE *fp;
	fp=0;
	argNum=1;
	scoreFilePrefix[0]=scoreFileSuffix[0]=pathwayFileName[0]='\0';
	summaryOutputFile=0;
	pathwayThreshold=geneLevelOutputThreshold=1000;
	nTopPathways=0;
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
		else if (FILLARG("--score-file-spec"))
		{
			if ((ptr = strstr(arg, "GENE")) == 0)
			{
				dcerror(1,"--score-file-spec does not contain the string GENE");
				return 0;
			}
			else
			{
				*ptr='\0';
				strcpy(scoreFilePrefix,arg);
				strcpy(scoreFileSuffix,ptr+4);
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
		else if (FILLARG("--summary-file"))
		{
			summaryOutputFile = fopen(arg, "w");
			if (summaryOutputFile == 0)
			{
				dcerror(1, "Could not open summary file %s\n", arg);
				return 0;
			}
		}
		else if (FILLARG("--output-file"))
		{
			fullOutputFile = fopen(arg, "w");
			if (fullOutputFile == 0)
			{
				dcerror(1, "Could not open output file %s\n", arg);
				return 0;
			}
		}
		else if (FILLARG("--SLP-file"))
		{
			SLPFile=fopen(arg,"w");
			if (SLPFile == 0)
			{
				dcerror(1,"Could not open SLP file %s\n",arg);
				return 0;
			}
		}
		else if (FILLARG("--pathway-file"))
		{
			strcpy(pathwayFileName,arg);
		}
		else if (FILLARG("--nPerms"))
		{
			nPerms=atoi(arg);
		}
		else if (FILLARG("--gene-level-output-threshold"))
		{
			geneLevelOutputThreshold=atof(arg);
		}
		else if (FILLARG("--pathway-threshold"))
		{
			pathwayThreshold=atof(arg);
		}
		else if (FILLARG("--score-table-file"))
		{
			strcpy(scoreTableFileName,arg);
		}
		else
			dcerror(1,"Did not recognise argument specifier %s\n",arg);
	}
	// need to do checks here at some point
	return 1;
}

int ppaParams::getNextArg(char *nextArg, int argc,char *argv[], FILE **fpp, int *argNum)
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

class pathway {
public:
	char name[MAX_PATHWAY_LENGTH],scoreFileName[1000];
	int geneList[MAX_GENES_IN_PATHWAY];
	int nGenes;
	int nOver;
	float SLP;
};

#define NTHRESHOLD 5
float threshold[NTHRESHOLD] = { 1,2,3,4,5 };

int getSLPs(pathway **allPathways, float **geneScores, float **pathwayScores, int *cc, ppaParams *pp, int realData, float analysisWiseResults[NTHRESHOLD + 1],int nOverAnalysisWise[NTHRESHOLD + 1])
{
	int nGene,missing[MAX_LOCI],s,nSub,g,gg,n[2],i,p,t,nOverThreshold[NTHRESHOLD];
	float sigma_x[2],sigma_x2[2],mean[2],var[2],SLP,SE,tval,s2,score[MAX_SUB],maxSLP;
	pathway *path;
	double pVal;
	char shortName[11];
	shortName[10]='\0';
	maxSLP=0;
	if (realData)
	{
		for (p = 0; p < pp->nPathways; ++p)
			fprintf(pp->SLPFile,"%s\t",allPathways[p]->name);
		fprintf(pp->SLPFile,"maxSLP\t");
		for (t = 0; t < NTHRESHOLD; ++t)
			fprintf(pp->SLPFile, "nSLP>%.1f\t", threshold[t]);
		fprintf(pp->SLPFile, "\n");
		for (p = 0; p < pp->nPathways; ++p)
			fprintf(pp->fullOutputFile, "%s\t", allPathways[p]->name);
		fprintf(pp->fullOutputFile,"maxSLP\t");
		for (t = 0; t < NTHRESHOLD; ++t)
			fprintf(pp->fullOutputFile, "nSLP>%.1f\t", threshold[t]);
		fprintf(pp->fullOutputFile, "\n");
	}
	for (t = 0; t < NTHRESHOLD; ++t)
		nOverThreshold[t] = 0;
	for (p = 0; p < pp->nPathways; ++p)
	{
		path = allPathways[p];
		if (realData)
		{
			assert((pathwayScores[p] = (float *)calloc(pp->nSub, sizeof(float))) != 0);
			for (s = 0; s<pp->nSub; ++s)
				score[s] = 0;
			for (g = 0; g < path->nGenes; ++g)
			{
				gg = path->geneList[g];
				for (s = 0; s < pp->nSub; ++s)
				{
					score[s] += geneScores[gg][s];
					pathwayScores[p][s] += geneScores[gg][s];
				}
			}
		}
		for (i=0;i<2;++i)
			sigma_x[i]=sigma_x2[i]=n[i]=0;
	for (s=0;s<pp->nSub;++s)
    {
		 ++n[cc[s]];
		 sigma_x[cc[s]]+=pathwayScores[p][s];
		 sigma_x2[cc[s]]+=pathwayScores[p][s]*pathwayScores[p][s];
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
    pVal=tstat(tval,n[0]+n[1]-2.0)/2; // one-tailed
	SLP=log10(2*pVal)*(mean[0]>=mean[1]?1:-1);
	fprintf(pp->SLPFile,"%10.2f\t",SLP);
	if (SLP>maxSLP)
		maxSLP=SLP;
	if (realData)
	{
		path->SLP=SLP;
		path->nOver = 0;
		if (SLP>pp->pathwayThreshold)
			pp->topPathways[pp->nTopPathways++]=p;
	}
	else
	{
		if (SLP >= path->SLP)
			++path->nOver;
	}
	for (t = 0; t < NTHRESHOLD; ++t)
		if (SLP>=threshold[t])
			++nOverThreshold[t];
	}
	if (realData)
	{
		analysisWiseResults[0] = maxSLP;
		for (t = 0; t < NTHRESHOLD; ++t)
			analysisWiseResults[t + 1] = nOverThreshold[t];
	}
	else
	{
		if (maxSLP >= analysisWiseResults[0])
			++nOverAnalysisWise[0];
		for (t = 0; t < NTHRESHOLD; ++t)
			if (nOverThreshold[t] >= analysisWiseResults[t + 1])
				++nOverAnalysisWise[t + 1];
	}
	fprintf(pp->SLPFile, "%.2f\t", maxSLP);
	for (t = 0; t < NTHRESHOLD; ++t)
		fprintf(pp->SLPFile, "%d\t", nOverThreshold[t]);
	fprintf(pp->SLPFile, "\n");
	return 1;
}


int fillGeneScores(char *line, pathway **allPathways, float **geneScores, int *cc, TStrIntMap &geneIndex, ppaParams *pp)
{
	pathway *p;
	char geneName[200],scoreFileName[1000],localLine[200],thisGene[50];
	TStrIntMap::iterator it;
	FILE *fs;
	int g,s,c;
	p=new pathway;
	if (sscanf(line,"%s %*s %[^\n]",p->name,rest)!=2)
		return 0;
	for (p->nGenes = 0; strcpy(line, rest), *rest = '\0', sscanf(line, "%s %[^\n]", geneName, rest) >= 1;)
	{
		it=geneIndex.find(geneName);
		if (it == geneIndex.end())
		{
			if (pp->scoreTableFile)
			{
				std::map<std::string,FILEOFFSET>::const_iterator geneIter=pp->geneIndex.find(geneName);
				if (geneIter==pp->geneIndex.end())
					geneIndex.insert(TStrIntPair(geneName,-1));
				else
				{
					if (pp->nSub==-1)
					{
						s=0;
						fseek(pp->scoreTableFile,0L,SEEK_SET);
						fscanf(pp->scoreTableFile,"%*s"); // ID
						while (1)
						{
							do {
								c=fgetc(pp->scoreTableFile);
							} while (c!=EOF && c!='\n' && isspace(c));
							if (c=='\n')
								break;
							ungetc(c,pp->scoreTableFile);
							fscanf(pp->scoreTableFile,"%*s");
						}
						fscanf(pp->scoreTableFile,"%*s"); // CC
						for (s=0;s<pp->nSub;++s)
							assert(fscanf(pp->scoreTableFile,"%d",cc[s])==1);
					}

					fseek(pp->scoreTableFile,geneIter->second,SEEK_SET);
					assert(fscanf(pp->scoreTableFile,"%s",thisGene)==1);
					if (strcmp(thisGene,geneName))
					{
#if 0
						dcerror(1,"Problem finding %s in %s. Tried to seek to position %" PRId64 " but got this: %s\n",
							pp->scoreTableFileName,geneName,geneIter->second,thisGene);
#else						
						dcerror(1,"Problem finding %s in %s. Tried to seek to position %uld but got this: %s\n",
							pp->scoreTableFileName,geneName,(unsigned long)geneIter->second,thisGene);
#endif
						exit(1);
					}
					for (s=0;s<pp->nSub;++s)
						assert(fscanf(pp->scoreTableFile,"%f",&geneScores[g][s])==1);
					geneIndex.insert(TStrIntPair(geneName,g));
					++pp->nGenes;
					p->geneList[p->nGenes++]=g;
				}
			}
			else
			{
				sprintf(scoreFileName,"%s%s%s",pp->scoreFilePrefix,geneName,pp->scoreFileSuffix);
				fs = fopen(scoreFileName,"r");
				if (fs == 0)
				{
					geneIndex.insert(TStrIntPair(geneName,-1));
				}
				else
				{
					if (pp->nSub == -1)
					{
						for (s=0;fgets(localLine,199,fs) && sscanf(localLine,"%*s %d",&cc[s])==1;++s)
							;
						pp->nSub=s;
						fseek(fs,0L,SEEK_SET);
					}
					g=pp->nGenes;
					assert((geneScores[g]=(float *)calloc(pp->nSub,sizeof(float)))!=0);
					for (s=0;s<pp->nSub;++s)
						assert(fgets(localLine,199,fs) && sscanf(localLine,"%*s %*d %f",&geneScores[g][s])==1);
					fclose(fs);
					geneIndex.insert(TStrIntPair(geneName,g));
					++pp->nGenes;
					p->geneList[p->nGenes++]=g;
				}
			}
		}
		else
		{
			g=it->second;
			if (g!=-1)
				p->geneList[p->nGenes++]=g;
		}
	}
	allPathways[pp->nPathways++]=p;
	return 1;
}

int ISTEMP;
#define INT_SWAP(x,y) (ISTEMP=x,x=y,y=ISTEMP)

int perm(int *cc, int nSub)
{
int i,k;
for (i=nSub-1;i>=1;--i)
  {
  k=rand()%(i+1);
  INT_SWAP(cc[i],cc[k]);
  }
return 1;
}


int indexGenes(ppaParams *pp)
{
	int ch,i;
	FILEOFFSET pos;
	char gene[50];
	pp->geneIndex.clear();
	fseek(pp->scoreTableFile,0L,SEEK_SET);
	for (i=0;i<2;++i)
		while ((ch=fgetc(pp->scoreTableFile))!='\n' && ch!=EOF)
			;
	while (1)
	{
		pos=ftell(pp->scoreTableFile);
		if (fscanf(pp->scoreTableFile,"%s",gene)!=1)
			break;
		pp->geneIndex[gene]=pos;
		while ((ch=fgetc(pp->scoreTableFile))!='\n' && ch!=EOF)
			;
		if (ch==EOF)
			break;
	}
	return 1;
}

int cc[MAX_SUB];

int main(int argc, char *argv[])
{
	ppaParams pp;
	FILE *fp;
	int s,i,nOverAnalysisWise[NTHRESHOLD + 1];
	float **geneScores,**pathwayScores,analysisWiseResults[NTHRESHOLD + 1];
	pathway **allPathways;
	TStrIntMap geneIndex;
	printf("%s v%s\n",PROGRAM,PPAVERSION);
	printf("MAX_LOCI=%d\nMAX_SUB=%d\n",MAX_LOCI,MAX_SUB);
	assert((geneScores=(float**)calloc(MAX_GENES,sizeof(float *)))!=0);
	assert((allPathways=(pathway **)calloc(MAX_PATHWAYS,sizeof(pathway *)))!=0);
	if (!pp.readParms(argc,argv))
		exit(1);
	if ((fp = fopen(pp.pathwayFileName, "r")) == 0)
	{
		dcerror(2,"Could not open pathway file %s\n",pp.pathwayFileName);
		return 1;
	}
	if (pp.scoreTableFileName[0])
	{
		pp.scoreTableFile=fopen(pp.scoreTableFileName,"rb"); //  see if fseek() works OK with binary mode
		if (pp.scoreTableFile==0)
		{
			dcerror(1,"Could not open score table file %s\n",pp.scoreTableFileName);
			return 1;
		}
		indexGenes(&pp);
	}
	pp.nSub=-1; // first time
	pp.nPathways=pp.nGenes=0;
	while (fgets(line, LONGLINELENGTH, fp))
	{
		fillGeneScores(line,allPathways,geneScores,cc,geneIndex,&pp);
	}
	assert((pathwayScores=(float**)calloc(pp.nPathways,sizeof(float *)))!=0);
	getSLPs(allPathways,geneScores,pathwayScores,cc,&pp,1,analysisWiseResults,nOverAnalysisWise);
	for (i = 0; i<NTHRESHOLD + 1; ++i)
		nOverAnalysisWise[i]=0;
	for (i = 0; i < pp.nPerms; ++i)
	{
		perm(cc,pp.nSub);
		getSLPs(allPathways, geneScores, pathwayScores, cc, &pp, 0,analysisWiseResults,nOverAnalysisWise);
	}
	for(i=0;i<pp.nPathways;++i)
		fprintf(pp.fullOutputFile,"%d\t",allPathways[i]->nOver);
	for(i=0;i<NTHRESHOLD+1;++i)
		fprintf(pp.fullOutputFile,"%d\t",nOverAnalysisWise[i]);
	fprintf(pp.fullOutputFile,"\n");
	for(i=0;i<pp.nPathways;++i)
		fprintf(pp.fullOutputFile,"%.6f\t",(allPathways[i]->nOver+1.0)/(pp.nPerms+1.0));
	for(i=0;i<NTHRESHOLD+1;++i)
		fprintf(pp.fullOutputFile,"%.6f\t",(nOverAnalysisWise[i]+1.0)/(pp.nPerms+1.0));
	fprintf(pp.fullOutputFile,"\n");
	fprintf(pp.summaryOutputFile,"Pathway\tSLP\tt_test_p\tempirical_p\n");
	for (i = 0; i < pp.nTopPathways; ++i)
		fprintf(pp.summaryOutputFile, "%s\t%.2f\t%.6f\t%.6f\n",
			allPathways[pp.topPathways[i]]->name,
			allPathways[pp.topPathways[i]]->SLP,
			(allPathways[pp.topPathways[i]]->SLP>0 ? pow(10.0,-allPathways[pp.topPathways[i]]->SLP)/2 : 1-pow(10.0,allPathways[pp.topPathways[i]]->SLP)/2),
			(allPathways[pp.topPathways[i]]->nOver+1.0)/(pp.nPerms+1.0));
	fclose(pp.summaryOutputFile);
	fclose(pp.SLPFile);
	return 0;
}