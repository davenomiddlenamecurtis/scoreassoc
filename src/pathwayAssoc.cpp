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

#define PROGRAM "pathwayAssoc"
#define PAVERSION "2.0"

#ifndef MAX_GENES
#define MAX_GENES 12000
#endif
#ifndef MAX_SUB
#define MAX_SUB 15000
#endif
#ifndef MAX_ID_LENGTH 
#define MAX_ID_LENGTH 100
#endif

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

class paParams : public lr_test_par_info {
public:
	char scoreFilePrefix[1000],scoreFileSuffix[1000],outputFilePrefix[1000],outputFileSuffix[1000],pathwayFileName[1000],scoreTableFileName[1000];
	std::map<std::string,FILEOFFSET> geneIndex;
	int nSub;
	float geneLevelOutputThreshold;
	int readParms(int argc,char *argv[]);
	int getNextArg(char *nextArg,int argc,char *argv[],FILE **fpp,int *argNum);
	FILE *summaryOutputFile,*scoreTableFile;
	int do_ttest, do_lrtest,do_linrtest,is_quantitative;
};

class pathwaySubject {
public:
	char id[MAX_ID_LENGTH+1]; 
	int cc;
	float pheno;
	float score[MAX_GENES];
};

float *totScore;

#define LONGLINELENGTH 40000
char long_line[LONGLINELENGTH+1],rest[LONGLINELENGTH+1];

lrVariable allVars[MAXLRVARIABLES];
std::map<std::string, lrVariable *> varMap;

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

int paParams::readParms(int argc, char *argv[])
{
	char arg[2000],*ptr;
	int argNum;
	FILE *fp;
	fp=0;
	argNum=1;
	scoreFilePrefix[0]=scoreFileSuffix[0]=pathwayFileName[0]=scoreTableFileName[0]='\0';
	summaryOutputFile=0;
	geneLevelOutputThreshold=1000;
	scoreTableFile=0;
	do_ttest = 1;
	do_lrtest = do_linrtest=is_quantitative=0;
	start_from_fitted = 1;
	numVars = numVarFiles = numTestFiles = numLinTestFiles = 0;
	lamda = DEFAULT_LAMDA;
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
			summaryOutputFile=fopen(arg,"w");
			if (summaryOutputFile == 0)
			{
				dcerror(1,"Could not open summary file %s\n",arg);
				return 0;
			}
		}
		else if (FILLARG("--pathway-file"))
		{
			strcpy(pathwayFileName,arg);
		}
		else if (FILLARG("--score-table-file"))
		{
			strcpy(scoreTableFileName,arg);
		}
		else if (FILLARG("--gene-level-output-threshold"))
		{
			geneLevelOutputThreshold = atof(arg);
		}
		else if (FILLARG("--lamda"))
		{
			lamda = atof(arg);
		}
		else if (FILLARG("--dottest"))
		{
			do_ttest = atoi(arg);
		}
		else if (FILLARG("--dolrtest"))
		{
			do_lrtest = atoi(arg);
		}
		else if (FILLARG("--isquantitative"))
		{
			is_quantitative = atoi(arg);
		}
		else if (FILLARG("--dolinrtest"))
		{
			do_linrtest = atoi(arg);
		}
		else if (FILLARG("--start-from-fitted"))
		{
			start_from_fitted = atoi(arg);
		}
		else if (FILLARG("--varfile"))
		{
			strcpy(varFiles[numVarFiles++].fn, arg);
		}
		else if (FILLARG("--testfile"))
		{
			strcpy(testFiles[numTestFiles++].fn, arg);
		}
		else if (FILLARG("--lintestfile"))
		{
			strcpy(linTestFiles[numLinTestFiles++].fn, arg);
		}
		else
			dcerror(1,"Did not recognise argument specifier %s\n",arg);
	}
	// do checks
	return 1;
}

int paParams::getNextArg(char *nextArg, int argc,char *argv[], FILE **fpp, int *argNum)
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

int readSubIds(pathwaySubject **sub,paParams *pp)
{
	FILE *fs,*fp;
	char gene[50],scoreFileName[200],pathwayName[1000],pathwayURL[1000];
	int s,c;
	if(pp->scoreTableFile)
	{
		s=0;
		fseek(pp->scoreTableFile,0L,SEEK_SET);
		fscanf(pp->scoreTableFile,"%*s"); // ID
		while(1)
		{
			do {
				c=fgetc(pp->scoreTableFile);
			} while(c!=EOF && c!='\n' && isspace(c));
			if(c=='\n')
				break;
			ungetc(c,pp->scoreTableFile);
			assert(fscanf(pp->scoreTableFile,"%s",sub[s++]->id)==1);
		}
		pp->nSub=s;
		fscanf(pp->scoreTableFile,"%*s"); // CC
		if (pp->is_quantitative)
			for (s = 0; s < pp->nSub; ++s)
				assert(fscanf(pp->scoreTableFile, "%f", &sub[s]->pheno) == 1);
		else 
			for(s=0;s<pp->nSub;++s)
				assert(fscanf(pp->scoreTableFile,"%d",&sub[s]->cc)==1);
	}
	else
	{
		if((fp = fopen(pp->pathwayFileName,"r")) == 0)
		{
			dcerror(2,"Could not open gene set file %s\n",pp->pathwayFileName);
			return 0;
		}
		fgets(long_line,LONGLINELENGTH,fp);
		fclose(fp);
		if(sscanf(long_line,"%s %s %[^\n]",pathwayName,pathwayURL,rest)!=3)
			{ dcerror(1,"Could not read any genes from first line in gene set file:\n%s\n",long_line); return 0; }
		do {
			strcpy(long_line,rest);
			*rest='\0';
			if (sscanf(long_line,"%s %[^\n]",gene,rest)<1)
			{
				dcerror(1,"Could not read any valid genes first pathway\n",long_line); return 0;
			}
			sprintf(scoreFileName,"%s%s%s",pp->scoreFilePrefix,gene,pp->scoreFileSuffix);
			fs=fopen(scoreFileName,"rb"); // because that's what I do for scoreTableFile
		} while(fs==0);
		if (pp->is_quantitative)
			for (s = 0; fgets(long_line, 1000, fs) && sscanf(long_line, "%s %f", sub[s]->id, &sub[s]->pheno) == 2; ++s)
			;
		else
			for (s = 0; fgets(long_line, 1000, fs) && sscanf(long_line, "%s %d", sub[s]->id, &sub[s]->cc) == 2; ++s)
			;
		pp->nSub=s;
		fclose(fs);
	}
	return pp->nSub;
}

float runOnePathway(char *line, pathwaySubject **sub, glModel *model,paParams *pp, int writeFile)
{
	char pathwayName[1000], pathwayURL[1000], gene[MAX_LOCI][50], scoreFileName[1000], outputFileName[1000], thisGene[50];
	FILE *fs, *fo;
	int nGene, missing[MAX_LOCI], s, g, n[2], i, cc, c,t,nValidGenes;
	float sigma_x[2], sigma_x2[2], mean[2], var[2], SLP, SE, tval, s2, score;
	double p, pathway_p;
	fo = 0;
	if (sscanf(line, "%s %s %[^\n]", pathwayName, pathwayURL, rest) != 3)
		return 0;
	for (nGene = 0; strcpy(line, rest), *rest = '\0', sscanf(line, "%s %[^\n]", gene[nGene], rest) >= 1; ++nGene)
		;
	for (s = 0; s < pp->nSub; ++s)
		totScore[s] = 0;
	nValidGenes = 0;
	for (g = 0; g < nGene; ++g)
	{
		if (pp->scoreTableFile)
		{
			std::map<std::string, FILEOFFSET>::const_iterator geneIter = pp->geneIndex.find(gene[g]);
			if (geneIter == pp->geneIndex.end())
				missing[g] = 1;
			else
			{
				fseek(pp->scoreTableFile, geneIter->second, SEEK_SET);
				assert(fscanf(pp->scoreTableFile, "%s", thisGene) == 1);
				if (strcmp(thisGene, gene[g]))
				{
#if 0
					dcerror(1, "Problem finding %s in %s. Tried to seek to position %" PRId64 " but got this: %s\n",
						pp->scoreTableFileName, gene[g], geneIter->second, thisGene);
#else
					dcerror(1, "Problem finding %s in %s. Tried to seek to position %uld but got this: %s\n",
						pp->scoreTableFileName, gene[g], (unsigned long)geneIter->second, thisGene);
#endif
					exit(1);
				}
				missing[g] = 0;
				++nValidGenes;
				for (s = 0; s < pp->nSub; ++s)
				{
					if (fscanf(pp->scoreTableFile, "%f", &sub[s]->score[g]) != 1)
					{
						dcerror(1, "Not enough entries in %s for gene %s\n", pp->scoreTableFileName, gene[g]); exit(1);
					}
					totScore[s] += sub[s]->score[g];
				}
			}
		}
		else
		{
			sprintf(scoreFileName, "%s%s%s", pp->scoreFilePrefix, gene[g], pp->scoreFileSuffix);
			fs = fopen(scoreFileName, "r");
			if (fs == 0)
				missing[g] = 1;
			else
			{
				missing[g] = 0;
				++nValidGenes;
				for (s = 0; s < pp->nSub; ++s)
				{
					if (!fgets(line, 1000, fs))
					{
						dcerror(1, "Not enough lines in scores file %s\n", scoreFileName); exit(1);
					}
					if (pp->is_quantitative) {
						if (sscanf(line, "%s %f %f", sub[s]->id, &sub[s]->pheno, &sub[s]->score[g]) != 3)
						{
							dcerror(1, "Not enough entries on this line in scores file %s:\n%s\n", scoreFileName, line); exit(1);
						}
					}
					else {
						if (sscanf(line, "%s %d %f", sub[s]->id, &sub[s]->cc, &sub[s]->score[g]) != 3)
						{
							dcerror(1, "Not enough entries on this line in scores file %s:\n%s\n", scoreFileName, line); exit(1);
						}
					}
					totScore[s] += sub[s]->score[g];
				}
				fclose(fs);
			}

		}
	}
	for (s = 0; s < pp->nSub; ++s)
		totScore[s] /= nValidGenes; // use average score not total score
	if (writeFile)
	{
		if (pp->summaryOutputFile != 0)
			fprintf(pp->summaryOutputFile, "%s\t%s\t", pathwayName, pathwayURL);
		sprintf(line, "%s%s%s", pp->outputFilePrefix, pathwayName, pp->outputFileSuffix);
		fo = fopen(line, "w");
		if (fo == 0)
		{
			dcerror(2, "Could not open output file %s\n", line);
		}
		else
			fprintf(fo, "%s\n%s\n\n", pathwayName, pathwayURL);
	}

	if (pp->do_ttest)
	{
	for (i = 0; i < 2; ++i)
		sigma_x[i] = sigma_x2[i] = n[i] = 0;
	for (s = 0; s < pp->nSub; ++s)
	{
		cc = sub[s]->cc;
		++n[cc];
		score = totScore[s];
		sigma_x[cc] += score;
		sigma_x2[cc] += score * score;
	}
	for (i = 0; i < 2; ++i)
	{
		var[i] = (sigma_x2[i] - sigma_x[i] * sigma_x[i] / n[i]) / (n[i] - 1);
		mean[i] = sigma_x[i] / n[i];
	}

	s2 = ((n[0] - 1)*var[0] + (n[1] - 1)*var[1]) / (n[0] + n[1] - 2);
	SE = sqrt(s2*(1 / (float)n[0] + 1 / (float)n[1]));
	if (SE == 0)
		tval = 0;
	else
		tval = (mean[1] - mean[0]) / SE;
	pathway_p = tstat(tval, n[0] + n[1] - 2.0) / 2; // one-tailed
	SLP = log10(2 * pathway_p)*(mean[0] >= mean[1] ? 1 : -1);
	if (fo != NULL)
	{
		fprintf(fo, "             Controls  Cases     \n"
			"N            %9d %9d\n"
			"Mean score   %9.3f %9.3f\n"
			"t (%d df) = %6.3f\n"
			"p = %10.8f\n"
			"SLP = %8.2f (signed log10(p), positive if cases score higher than controls)\n",
			n[0], n[1], mean[0], mean[1], n[0] + n[1] - 2, tval, 2 * pathway_p, SLP);
		if (pp->summaryOutputFile != 0)
			fprintf(pp->summaryOutputFile, "%f\t", SLP);
		if (SLP > pp->geneLevelOutputThreshold)
		{
			fprintf(fo, "\n\nSLPs for individual genes above threshold:\n");
			for (g = 0; g < nGene; ++g)
			{
				if (missing[g] == 1)
					continue;
				for (i = 0; i < 2; ++i)
					sigma_x[i] = sigma_x2[i] = n[i] = 0;
				for (s = 0; s < pp->nSub; ++s)
				{
					cc = sub[s]->cc;
					++n[cc];
					score = sub[s]->score[g];
					sigma_x[cc] += score;
					sigma_x2[cc] += score * score;
				}
				for (i = 0; i < 2; ++i)
				{
					var[i] = (sigma_x2[i] - sigma_x[i] * sigma_x[i] / n[i]) / (n[i] - 1);
					mean[i] = sigma_x[i] / n[i];
				}

				s2 = ((n[0] - 1)*var[0] + (n[1] - 1)*var[1]) / (n[0] + n[1] - 2);
				SE = sqrt(s2*(1 / (float)n[0] + 1 / (float)n[1]));
				if (SE == 0)
					tval = 0;
				else
					tval = (mean[1] - mean[0]) / SE;
				p = tstat(tval, n[0] + n[1] - 2.0) / 2; // one-tailed
				SLP = log10(2 * p)*(mean[0] >= mean[1] ? 1 : -1);

				if (SLP > pp->geneLevelOutputThreshold)
					fprintf(fo, "%s %8.2f\n", gene[g], SLP);
			}
			fprintf(fo, "\n\nList of genes for which no score file was found:\n");
			for (g = 0; g < nGene; ++g)
				if (missing[g] == 1)
					fprintf(fo, "%s\n", gene[g]);
		}

	}
	if (mean[0] > mean[1])
		pathway_p = 1.0 - pathway_p;
	}
	if(pp->do_lrtest || pp->do_linrtest || pp->numTestFiles > 0 || pp->numLinTestFiles > 0)
		fillModelWithVars(model,pp->nSub,pp,pp->scoreCol);
	if (pp->do_lrtest)
	{
		SLP = do_onetailed_LRT(fo, model, pp,0);
		if (writeFile && pp->summaryOutputFile != 0)
			fprintf(pp->summaryOutputFile, "%f\t", SLP);
	}
	if (pp->do_linrtest)
	{
		SLP = do_onetailed_LRT(fo, model, pp,1);
		if (writeFile && pp->summaryOutputFile != 0)
			fprintf(pp->summaryOutputFile, "%f\t", SLP);
	}
	if (pp->numTestFiles > 0)
	{
		for (t = 0; t < pp->numTestFiles; ++t)
		{
			SLP = runTestFile(fo, pp->testFiles[t].fn, model, pp);
			if (writeFile && pp->summaryOutputFile != 0)
				fprintf(pp->summaryOutputFile, "%f\t", SLP);
		}
	}
	if (pp->numLinTestFiles > 0)
	{
		for (t = 0; t < pp->numLinTestFiles; ++t)
		{
			SLP = runLinTestFile(fo, pp->linTestFiles[t].fn, model, pp);
			if (writeFile && pp->summaryOutputFile != 0)
				fprintf(pp->summaryOutputFile, "%f\t", SLP);
		}
	}
	if (fo)
		fclose(fo);
	if (writeFile && pp->summaryOutputFile != 0)
		fprintf(pp->summaryOutputFile, "\n");
	return SLP;
}

int indexGenes(paParams *pp)
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

int main(int argc, char *argv[])
{
	paParams pp;
	FILE *fp;
	pathwaySubject **sub;
	int s,filledModel;
	glRidgePenaltyModel model;
	printf("%s v%s\n",PROGRAM,PAVERSION);
	printf("MAX_LOCI=%d\nMAX_SUB=%d\n",MAX_LOCI,MAX_SUB);
	if (!pp.readParms(argc,argv))
		exit(1);
	assert(sub=(pathwaySubject **)calloc(MAX_SUB,sizeof(pathwaySubject*)));
	for (s=0;s<MAX_SUB;++s)
		assert(sub[s]=(pathwaySubject *)calloc(1,sizeof(pathwaySubject)));
	if ((fp = fopen(pp.pathwayFileName,"r")) == 0)
	{
		dcerror(2,"Could not open gene set file %s\n",pp.pathwayFileName);
		return 1;
	}
	if(pp.scoreTableFileName[0])
	{
		pp.scoreTableFile=fopen(pp.scoreTableFileName,"rb"); //  see if fseek() works OK with binary mode
		if(pp.scoreTableFile==0)
		{
			dcerror(1,"Could not open score table file %s\n",pp.scoreTableFileName);
			return 1;
		}
		indexGenes(&pp);
	}
	readSubIds(sub,&pp);
	if(pp.numVarFiles)
	{
		std::map<std::string,int> subIDs;
		for(s = 0; s<pp.nSub; ++s)
			subIDs[sub[s]->id] = s;
		if(!readVarFiles(subIDs,pp.nSub,&pp))
			exit(1);
	}
	totScore = (float*)calloc(pp.nSub,sizeof(float));
	assert(totScore != 0);
	pp.scoreCol = pp.numVars;
	strcpy(allVars[pp.scoreCol].name,"score");
	allVars[pp.scoreCol].val = totScore;
	varMap["score"] = &allVars[pp.scoreCol];
	++pp.numVars;
	if (pp.do_lrtest || pp.do_linrtest || pp.numTestFiles > 0 || pp.numLinTestFiles > 0)
	{
		model.lamda = pp.lamda;
		fillModelWithVars(&model, pp.nSub, &pp);
		for (s = 0; s < pp.nSub; ++s)
			model.Y[s] = pp.is_quantitative?sub[s]->pheno:sub[s]->cc;
		filledModel = 1;
	}

	while (fgets(long_line, LONGLINELENGTH, fp))
	{
		runOnePathway(long_line,sub,&model,&pp,1);
	}
	if (pp.summaryOutputFile!=0)
		fclose(pp.summaryOutputFile);
//	free(totScore); will be done by the~lrVariable()
	return 0;
}
