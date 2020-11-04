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
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <string>
#include <map>
#include <assert.h>
#include "scoreassoc.hpp"

#include "dcexpr.hpp"
#include "safilterfuncs.hpp"

#define PROGRAM "getVarScores"
#define GVSVERSION "1.0"

#define MAXVARTYPES 100
#define MAXVARNAMELENGTH 100
struct varType_t { char name[MAXVARNAMELENGTH]; int flag; };
typedef struct varType_t varType;
varType varFlagTable[MAXVARTYPES];
float **varScore;

option opt[]=
{
	{ "flagfile", FLAGFILE },
	{ "psdatafile", PSDATAFILE },
	{ "gcdatafile", GCDATAFILE },
	{ "gendatafile", GENDATAFILE },
	{ "weightfile", WEIGHTFILE },
	{ "filterfile", FILTERFILE },
	{ "annotfile", ANNOTFILE },
	{ "locusfilterfile", LOCUSFILTERFILE },
	{ "locusnamefile", LOCUSNAMEFILE },
	{ "locusweightfile", LOCUSWEIGHTFILE },
	{ "triofile", TRIOFILE },
	{ "outfile", OUTFILE },
	{ "scorefile", SCOREFILE },
	{ "nostringtomatchthis", NUMDATAFILETYPES },
	{ "numloci", NUMLOCI },
// need this if using gc files, unless locusfilterfile provided
	{ "usehaps", USEHAPS },
	{ "weightfactor", WEIGHTFACTOR },
	{ "argfile", ARGFILE },
	{"", NUMOPTS}
};
// readable files must be listed before writable files

void usage()
{
	printf("getVarScores --psdatafile file || --gendatafile file || --gcdatafile file     [options]\n\nOptions:\n"
		"--flagfile (binary flags for variant types)"
"--flagfile file (attribute flags with binary codes)\n"
"--weightfactor x (weight for very rare variants, default 10)\n"
"--outfile file\n"
"--scorefile file "
"--weightfile file (specify functional weights)\n"
"--annotfile file (annotations from plink/seq)\n"
"--filterfile file\n"
"--locusweightfile file (specify functional weight for each locus)\n"
"--locusnamefile file (provide name for each locus)\n"
"--locusfilterfile file (exclude specific loci)\n"
"--triofile file\n"
"--usehaps\n"
"--numloci x (needed with --gcdatafile or --gendatafile)\n"
"--argfile file (additional arguments)\n"
);
}

#define MAXDEPTH 5

void getVarScores(varType *varFlagTable,int nVarTypes,float **varScore,float *weight,float *func_weight,float *missing,int *rarer,subject **sub,int nsub,par_info *pi,sa_par_info *spi)
{
	int l,ll,s,a,p,v,flags;
	float increment;
	for (s=0;s<nsub;++s)
	{
		for (v=0;v<nVarTypes;++v)
			varScore[s][v]=0;
		for (l=0;l<pi->n_loci_to_use;++l)
		{
			ll=pi->loci_to_use[l];
			if (spi->use_probs)
			{
				if (sub[s]->prob[ll][0]+sub[s]->prob[ll][1]+sub[s]->prob[ll][2]==0)
					increment=missing[l]*2;
				else if (rarer[l]==2)
					increment=(sub[s]->prob[ll][1]+sub[s]->prob[ll][2]*2)*weight[l];
				else
					increment=(sub[s]->prob[ll][1]+sub[s]->prob[ll][0]*2)*weight[l];
			}
			else
			{
				increment=0;
				for (a = 0; a < 2; ++a)
					if (sub[s]->all[ll][a] == rarer[l])
						increment+=weight[l];
					else if (sub[s]->all[ll][a] == 0)
						increment+=missing[l];
			}
			flags=(int)func_weight[ll];
			for (v=0;v<nVarTypes;++v)
			{
				if (varFlagTable[v].flag&flags)
					varScore[s][v]+=increment;
			}
		}
	}
}

void writeVarScores(FILE *fs,subject **sub,int nsub,int nVarTypes,float **varScore)
{
	int s,v;
	for (s=0;s<nsub;++s) 
	{
		fprintf(fs,"%20s %d ",sub[s]->id,sub[s]->cc);
		for (v=0;v<nVarTypes;++v)
			fprintf(fs,"%8.4f ",varScore[s][v]);
		fprintf(fs,"\n");
	}
}

int readFlagTable(varType *vt,sa_par_info *spi)
{
	int nVarTypes;
	char line[1000];
	for (nVarTypes=0;fgets(line,999,spi->df[FLAGFILE].fp)&&sscanf(line,"%s %d",vt[nVarTypes].name,&vt[nVarTypes].flag)==2;++nVarTypes)
		;
	return nVarTypes;
}

int getNextArg(char *nextArg, int argc,char *argv[], FILE *fp[MAXDEPTH],int *depth, int *argNum)
{
	*nextArg='\0';
	while (*depth>-1)
	{
		if (fscanf(fp[*depth],"%s ",nextArg)==1)
			return 1;
		else
		{
			fclose(fp[*depth]);
			--*depth;
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

int read_all_args(char *argv[],int argc, par_info *pi, sa_par_info *spi)
{
	char arg[2000];
	int a,arg_depth,arg_num;
	FILE *fp[MAXDEPTH];
	arg_depth=-1;
	arg_num=1;
	pi->nloci=0;
	spi->use_locus_names=spi->use_comments=1;
	spi->wfactor=10;
	spi->do_recessive_test=0;
	spi->LD_threshold=0.9;
	spi->use_haplotypes=spi->use_trios=0;
	spi->use_cc_freqs[0]=spi->use_cc_freqs[1]=0;
	spi->use_probs=0;
	while (getNextArg(arg, argc, argv, fp,&arg_depth, &arg_num))
	{
		if (strncmp(arg, "--", 2))
		{
			printf("Need to specify options with --\n%s\n\n", arg);
			usage();
			exit(1);
		}
		for (a=0;a<NUMOPTS-1;++a)
			if (!strncmp(arg+2,opt[a].str,strlen(opt[a].str)))
				break;
		if (a == NUMOPTS - 1)
			{
			printf("Unrecognised option: \n%s\n\n", arg);
			usage();
			exit(1);
			}
		int error=0;
		switch (opt[a].o)
		{
		case ARGFILE:
			if (++arg_depth >= MAXDEPTH)
			{
				dcerror(1, "Attempting to recurse too deeply into arg-files with this one: %s\n", arg);
				return 0;
			}
			else if (getNextArg(arg, argc, argv, fp,&arg_depth, &arg_num) == 0 || arg[0]=='-')
				error=1;
			else
			{
				fp[arg_depth] = fopen(arg, "r");
				if (fp[arg_depth] == NULL)
				{
					dcerror(1, "Could not open arg file: %s\n", arg);
					return 0;
				}
			}
			break;
		case NUMLOCI:
			if (getNextArg(arg, argc, argv, fp,&arg_depth, &arg_num) == 0 || sscanf(arg, "%d", &pi->nloci) != 1)
				error=1;
			break;
		case WEIGHTFACTOR:
			if (getNextArg(arg, argc, argv, fp,&arg_depth, &arg_num) == 0 || sscanf(arg, "%f", &spi->wfactor) != 1)
				error=1;
			break;
		case USEHAPS:
			spi->use_haplotypes=1;
			break;
		case FLAGFILE:
		case PSDATAFILE:
		case GCDATAFILE:
		case GENDATAFILE:
		case WEIGHTFILE:
		case ANNOTFILE:
		case FILTERFILE:
		case LOCUSFILTERFILE:
		case LOCUSWEIGHTFILE:
		case LOCUSNAMEFILE:
		case OUTFILE:
		case SCOREFILE:
		case TRIOFILE:
			if (getNextArg(arg, argc, argv, fp,&arg_depth, &arg_num) == 0 || arg[0]=='-' || sscanf(arg, "%s",spi->df[opt[a].o].fn) != 1)
				error=1;
			break;
		}
		if (error)
			{
				dcerror(1, "Error reading %s\n %s", opt[a].str,arg); exit(1);
			}
	}
	return 1;
}

int process_options(par_info *pi, sa_par_info *spi)
{
	int a,l;
	if (!spi->df[SCOREFILE].fn[0])
	{
		dcerror(1,"Must specify --scorefile"); exit(1);
	}
	if (!spi->df[FLAGFILE].fn[0])
	{
		dcerror(1,"Must specify --flagfile"); exit(1);
	}
	if ((spi->df[PSDATAFILE].fn[0]==0)+(spi->df[GCDATAFILE].fn[0]==0)+(spi->df[GENDATAFILE].fn[0]==0) < 2)
	{
		dcerror(1,"Must specify only one of --psdatafile, --gcdatafile or --gendatafile"); exit(1);
	}
	// put below into scoreassoc
	if ((spi->df[GCDATAFILE].fn[0]!='\0' || spi->df[GENDATAFILE].fn[0]!='\0' ) && spi->df[LOCUSFILTERFILE].fn[0]=='\0'  && pi->nloci==0)
	{
		dcerror(1,"Need to specify --numloci or --locusfilterfile when using --gcdatafile or --gendatafile"); exit(1);
	}
	if (spi->df[GENDATAFILE].fn[0] != '\0')
	{
		if (spi->use_trios)
		{
			dcerror(1, "Cannot use trios with --gendatafile"); exit(1);
		}
		spi->use_probs = 1;
	}

	for (a=0;a<OUTFILE;++a)
		if (spi->df[a].fn[0] && (spi->df[a].fp = fopen(spi->df[a].fn, "r")) == 0)
		{
			dcerror(1, "Could not open %s %s\n", opt[a].str,spi->df[a].fn); exit(1);
		}
	for (a=OUTFILE;a<=SCOREFILE;++a)
		if (spi->df[a].fn[0] && (spi->df[a].fp = fopen(spi->df[a].fn, "w")) == 0)
		{
			dcerror(1, "Could not open %s %s\n", opt[a].str,spi->df[a].fn); exit(1);
		}
	if (spi->df[WEIGHTFILE].fp && spi->df[ANNOTFILE].fp)
		spi->use_func_weights=1;
	else if (spi->df[WEIGHTFILE].fp && !spi->df[ANNOTFILE].fp)
		printf("Warning: weightfile specified but not annotfile so will not assign weights according to function\n");
	else if (!spi->df[WEIGHTFILE].fp && spi->df[ANNOTFILE].fp)
		printf("Warning: annotfile specified but not weightfile so will not assign weights according to function\n");
	if (spi->df[LOCUSWEIGHTFILE].fp)
	{
		spi->use_func_weights=1;
		if (spi->df[WEIGHTFILE].fp)
			printf("Warning: weightfile specified but will read weights from locusweightfile instead\n");
	}
	if (spi->df[TRIOFILE].fp)
		spi->use_trios=1;
}

char *skip_word(char *ptr)
{
	while (*ptr && !isspace(*ptr))
		++ptr;
	while (*ptr && isspace(*ptr))
		++ptr;
	return ptr;
}

int read_all_gen_subjects(FILE *fi,subject **sub,int *nsub,par_info *pi)
{
char id[MAX_ID_LENGTH+1],*ptr;
int n_to_skip,s,i,l;
s=0;
while (fgets(long_line,LONG_LINE_LENGTH,fi) && sscanf(long_line,"%s",id)==1)
  {
  if (s>=MAX_SUB)
    { error("Number of subjects exceeds MAX_SUB",""); return 0; }
  if (sscanf(long_line,"%s %d",sub[s]->id,&sub[s]->cc)!=2)
    { error("Syntax error in subject line: ",long_line); return 0; }
  n_to_skip=2;
  ptr=long_line;
  for (i=0;i<n_to_skip;++i)
  {
	  while (*ptr && isspace(*ptr++)) ;
	  while (*ptr &&!isspace(*ptr++)) ;
  }
  for (l = 0; l < pi->nloci; ++l)
  {
	  if (sscanf(ptr,"%f %f %f",&sub[s]->prob[l][0],&sub[s]->prob[l][1],&sub[s]->prob[l][2])!=3)
	  { error("Not enough genotype probabilities in subject line: ",long_line); return 0; }
	  n_to_skip=3;
	  for (i=0;i<n_to_skip;++i)
	  {
		while (*ptr && isspace(*ptr++)) ;
		while (*ptr &&!isspace(*ptr++)) ;
	  }
  }
  ++s;
  }
*nsub=s;
return 1;
}

// need to consider that not all subjects may be typed for all loci
// code relies on fact that sub is allocated with calloc, which will have set all alleles to 0
// NB the annot file MUST have two sets of counts - for cases and controls

int read_ps_datafile(par_info *pi, sa_par_info *spi, subject **sub, int *nsubptr, char names[MAX_LOCI][20], char comments[MAX_LOCI][MAX_COMMENT_LENGTH], float func_weight[MAX_LOCI],	
	std::map<std::string,float> weightMap,	std::map<std::string,std::string> effectMap)
{
	int nsub,first,s,a,l,f;
	char dline[1000],aline[1000],pos[100],rsname[100],ref[100],alls[100],all[2][100],sname[100],ccstr[100],*ptr;
	float weight;
	std::map<std::string,int> subIDs;
	first=1;
	nsub=0;
	pi->nloci=0;
	while (fgets(dline, 999, spi->df[PSDATAFILE].fp))
	{
		if (sscanf(dline,"%s",pos)!=1)
			continue;
		if (!strncmp("chr", pos, 3) && strchr(pos, ':')) //looks like new locus, not subject ID, will break if subject does have ID like this
		{
			++pi->nloci;
			l=pi->nloci-1;
			pi->n_alls[l]=2;
			if (sscanf(dline, "%s %s %s", pos, rsname, alls) != 3)
			{
				dcerror(1, "Could not read this locus line:\n%s\n", dline);
				exit(1);
			}
			if (sscanf(alls,"%[^/]",ref)!=1)
			{
				dcerror(1, "Could not read reference allele in this locus line:\n%s\n", dline);
				exit(1);
			}
			if (spi->df[WEIGHTFILE].fp && spi->df[ANNOTFILE].fp)
			{
				std::map<std::string,std::string>::const_iterator effIter=effectMap.find(pos);
				if (effIter == effectMap.end())
				{
					dcerror(1,"%s not found in annotation file %s\n",pos,spi->df[ANNOTFILE].fn);
					exit(1);
				}
				sprintf(comments[l],"%s_%s_%s", pos, effIter->second.c_str(),alls);
				std::map<std::string,float>::const_iterator weightIter =weightMap.find(effIter->second);
				if (weightIter==weightMap.end())
					{
						dcerror(1,"weight for effect %s not found in weight file %s\n",effIter->second.c_str(),spi->df[WEIGHTFILE].fn);
						exit(1);
					}
				else 
					func_weight[l]=weightIter->second;
			}
			else
			{
				sprintf(comments[l],"%s_%s", pos, alls);
				func_weight[l]=1;
			}
		}
		else
		{
			if (sscanf(dline,"%s %*d %s %[^/]/%s", sname, ccstr, all[0], all[1]) != 4)
			{
				dcerror(1, "Could not read this genotype line:\n%s\n", dline); exit(1);
			}
			std::map<std::string,int>::const_iterator iter =subIDs.find(sname);
			if (iter == subIDs.end()) // found a subject we haven't seen before
			{
				subIDs[sname]=nsub;
				s=nsub;
				strcpy(sub[s]->id,sname);
				if (!strcmp(ccstr,"CONTROL"))
					sub[s]->cc=0;
				else if (!strcmp(ccstr,"CASE"))
					sub[s]->cc=1;
				else
				{
					dcerror(1, "Unrecognised phenotype in this line:\n%s\n", dline); exit(1);
				}
				++nsub;
			}
			else
				s=iter->second;
			if (all[0][0]=='.')
				continue;
			else for (a=0;a<2;++a)
				if (!strcmp(all[a],ref))
					sub[s]->all[l][a]=1;
				else
					sub[s]->all[l][a]=2;
		}
	}
	*nsubptr=nsub;
	return 1;
}

int main(int argc, char *argv[])
{
	char arg_string[2000];
	int nsub,n_new_sub,real_nsub,nVarTypes;
	float **varScore,p,*score;
	int s,n_non_mendelian;
	non_mendelian *non_mendelians;
	char *non_mendelian_report;
	par_info pi;
	sa_par_info spi;
	subject **sub,**new_sub,**real_sub;
	pi.use_cc=1;
	printf("%s v%s\n",PROGRAM,GVSVERSION);
	printf("MAX_LOCI=%d\nMAX_SUB=%d\n",MAX_LOCI,MAX_SUB);

	assert(sub=(subject **)calloc(MAX_SUB,sizeof(subject*)));
	for (s=0;s<MAX_SUB;++s)
		assert(sub[s]=(subject *)calloc(1,sizeof(subject)));
	max_cc[0]=max_cc[1]=0;
	read_all_args(argv,argc, &pi, &spi);
	// make_arg_string(arg_string,argc,argv);
	// parse_arg_string(arg_string,&pi,&spi,&pspi);
	process_options(&pi,&spi);
	if (spi.df[FILTERFILE].fp)
		initExclusions(spi.df[FILTERFILE].fp);
	assert(score = (float*)calloc(MAX_SUB, sizeof(float))); // this is only here because read_all_data() may use it
	read_all_data(&pi,&spi,sub,&nsub,names,comments,func_weight,score);
if (spi.use_trios)
{
	if (atoi(comments[0])>22 || toupper(comments[0][0]) == 'X' || toupper(comments[0][0]) == 'Y' ||
		toupper(comments[0][0]) == 'C'&&toupper(comments[0][1]) == 'H'&&toupper(comments[0][2]) == 'R' && (atoi(comments[0] + 3) > 22 || toupper(comments[0][3]) == 'X' || toupper(comments[0][3]) == 'Y'))
	{
		error("Cannot at present use trios for genes on X or Y chromosome", "");
		return 1;
	}
}

if (spi.use_trios)
	{
		assert(new_sub=(subject **)calloc(nsub,sizeof(subject*)));
		for (s=0;s<nsub;++s)
			assert(new_sub[s]=(subject *)calloc(1,sizeof(subject)));
		assert(non_mendelians=(non_mendelian *)calloc(MAX_SUB,sizeof(non_mendelian)));
		if ((n_new_sub=sort_trios(sub,nsub,&pi,&spi,new_sub,non_mendelians,&n_non_mendelian,long_line))==0)
			exit(1);
		real_sub=sub;
		sub=new_sub;
		real_nsub=nsub;
		nsub=n_new_sub;
		assert(non_mendelian_report=(char*)malloc(strlen(long_line)+1));
		strcpy(non_mendelian_report,long_line);
	}
else
	non_mendelian_report=0;
	// not used, compiler error otherwise
get_freqs(sub,nsub,&pi,&spi,cc_freq,cc_count,cc_genocount);
applyExclusions(sub, nsub ,&pi);
spi.use_func_weights=0; // this is a trick to prevent the rarity weight being multiplied by the functional weight
set_weights(0,weight,missing_score,rarer,sub,nsub,&pi,&spi,func_weight,cc_freq,cc_count,max_cc,names,comments);
nVarTypes=readFlagTable(varFlagTable,&spi);
assert(varScore=(float **)malloc(nsub*sizeof(float*)));
for (s=0;s<nsub;++s)
	assert(varScore[s]=(float*)malloc(nVarTypes*sizeof(float)));

// allocate a table to hold the subject scores for each variant type, then output fill it and it
getVarScores(varFlagTable,nVarTypes,varScore,weight,func_weight,missing_score,rarer,sub,nsub,&pi,&spi);
// beware, func_weight is indexed differently
// weight[l]*=func_weight[pi->loci_to_use[l]];
// weight and missing score are only calcuated for loci to be used
writeVarScores(spi.df[SCOREFILE].fp,sub,nsub,nVarTypes,varScore);
fclose(spi.df[SCOREFILE].fp);
spi.df[SCOREFILE].fp=0;

fprintf(spi.df[OUTFILE].fp,"geneVarAssoc output\n");
fprintf(spi.df[OUTFILE].fp,"Used %d valid variants\n",pi.n_loci_to_use);
stateExclusions(spi.df[OUTFILE].fp);
fclose(spi.df[OUTFILE].fp);
spi.df[OUTFILE].fp=0; //  because otherwise the destructor will try to fclose it
printf("\nProgram run completed OK\n");
return 0;

}

int readVarFiles(std::map<std::string, int> subIDs, int nSub, lr_test_par_info* spi)
{
	dcerror(5, "readVarFiles() has been called in getVarScores but this function is not supported.");
	return 0;
}