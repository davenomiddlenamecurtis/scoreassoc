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
#include <assert.h>
#include <ctype.h>
#include <string>
#include <map>
#include "scoreassoc.hpp"
#include "dcerror.hpp"

option opt[]=
{
	{ "psdatafile", PSDATAFILE },
	{ "gcdatafile", GCDATAFILE },
	{ "gendatafile", GENDATAFILE },
	{ "weightfile", WEIGHTFILE },
	{ "annotfile", ANNOTFILE },
	{"filterfile",FILTERFILE},
	{ "locusfilterfile", LOCUSFILTERFILE },
	{ "locusnamefile", LOCUSNAMEFILE },
	{ "locusweightfile", LOCUSWEIGHTFILE },
	{"triofile",TRIOFILE},
	{"samplefile",SAMPLEFILE},
	{ "casefreqfile", CASEFREQFILE },
	{ "contfreqfile", CONTFREQFILE },
	{ "outfile", OUTFILE },
	{ "scorefile", SCOREFILE },
	{ "nostringtomatchthis", NUMDATAFILETYPES }, // data files above must be in same order as enum in header
	{ "numloci", NUMLOCI },
// need this if using gc files
	{ "ldthreshold", LDTHRESHOLD },
	{ "minweight", WEIGHTTHRESHOLD },
	{ "dorecessive", DORECESSIVE },
	{ "dottest",DOTTEST },
	{ "dolrtest",DOLRTEST },
	{ "start-from-fitted",STARTFROMFITTED },
	{ "usehaps", USEHAPS },
	{ "showhaplocusnames", SHOWHAPLOCUSNAMES },
	{ "weightfactor", WEIGHTFACTOR },
	{"lamda",LAMDA},
	{"varfile",VARFILE},
	{"testfile",TESTFILE},
	{"argfile",ARGFILE},
{"", NUMOPTS}
};
// readable files must be listed before writable files

void usage()
{
	printf("scoreassoc --psdatafile file || --gendatafile file || --gcdatafile file     [options]\n\nOptions:\n"
"--samplefile file (sample file to match IMPUTE2 datafile)\n"
"--weightfactor x (weight for very rare variants, default 10)\n"
"--outfile file\n"
"--scorefile file (optionallly output scores for each subject)"
"--weightfile file (specify functional weights)\n"
"--annotfile file (annotations from plink/seq)\n"
"--filterfile file\n"
"--locusweightfile file (specify functional weight for each locus)\n"
"--locusnamefile file (provide name for each locus)\n"
"--locusfilterfile file (exclude specific loci)\n"
"--casefreqfile file (provide allele frequency for each locus in cases)\n"
"--contfreqfile file (provide allele frequency for each locus in controls)\n"
"--triofile file\n"
"--ldthreshold x (to discard variants in LD for recessive analysis, default 0.9)\n"
"--minweight x (to include in recessive analysis, default 0, i.e. all variants)\n"
"--lamda x\n"
"--dorecessive x\n"
"--dottest x\n"
"--dolrtest x\n"
"--start-from-fitted x\n"
"--varfile file\n"
"--testfile file\n"
"--usehaps x\n"
"--showhaplocinames x\n"
"--numloci x (needed with --gcdatafile or --gendatafile)\n"
"--argfile file (additional arguments)\n"
);
}

#define MAXDEPTH 5


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
	spi->show_hap_locus_names=spi->use_haplotypes=spi->use_trios=0;
	spi->use_cc_freqs[0]=spi->use_cc_freqs[1]=0;
	spi->use_probs=0;
	spi->do_ttest = 1;
	spi->do_lrtest = 0;
	spi->start_from_fitted = 1;
	spi->numVars=spi->numVarFiles=spi->numTestFiles=0;
	spi->lamda=DEFAULT_LAMDA;
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
			if(getNextArg(arg,argc,argv,fp,&arg_depth,&arg_num) == 0 || sscanf(arg,"%f",&spi->wfactor) != 1)
				error=1;
			break;
		case LAMDA:
			if(getNextArg(arg,argc,argv,fp,&arg_depth,&arg_num) == 0 || sscanf(arg,"%f",&spi->lamda) != 1)
				error=1;
			break;
		case LDTHRESHOLD:
			if (getNextArg(arg, argc, argv, fp,&arg_depth, &arg_num) == 0 || sscanf(arg,"%f",&spi->LD_threshold)!=1)
				error=1;
			break;
		case WEIGHTTHRESHOLD:
			if (getNextArg(arg, argc, argv, fp,&arg_depth, &arg_num) == 0 || sscanf(arg,"%f",&spi->weight_threshold)!=1)
				error=1;
			break;
		case DORECESSIVE:
			if(getNextArg(arg,argc,argv,fp,&arg_depth,&arg_num) == 0 || sscanf(arg,"%d",&spi->do_recessive_test)!=1)
				error=1;
			break;
		case DOTTEST:
			if(getNextArg(arg,argc,argv,fp,&arg_depth,&arg_num) == 0 || sscanf(arg,"%d",&spi->do_ttest)!=1)
				error=1;
			break;
		case DOLRTEST:
			if (getNextArg(arg, argc, argv, fp, &arg_depth, &arg_num) == 0 || sscanf(arg, "%d", &spi->do_lrtest) != 1)
				error = 1;
			break;
		case STARTFROMFITTED:
			if (getNextArg(arg, argc, argv, fp, &arg_depth, &arg_num) == 0 || sscanf(arg, "%d", &spi->start_from_fitted) != 1)
				error = 1;
			break; 
		case USEHAPS:
			if (getNextArg(arg, argc, argv, fp, &arg_depth, &arg_num) == 0 || sscanf(arg, "%d", &spi->use_haplotypes) != 1)
				error = 1;
			break;
		case SHOWHAPLOCUSNAMES:
			if (getNextArg(arg, argc, argv, fp, &arg_depth, &arg_num) == 0 || sscanf(arg, "%d", &spi->show_hap_locus_names) != 1)
				error = 1;
			break;
		case VARFILE:
			if(getNextArg(arg,argc,argv,fp,&arg_depth,&arg_num) == 0 || sscanf(arg,"%s",spi->varFiles[spi->numVarFiles++].fn)!=1)
				error=1;
			break;
		case TESTFILE:
			if(getNextArg(arg,argc,argv,fp,&arg_depth,&arg_num) == 0 || sscanf(arg,"%s",spi->testFiles[spi->numTestFiles++].fn)!=1)
				error=1;
			break;
		case PSDATAFILE:
		case GCDATAFILE:
		case GENDATAFILE:
		case WEIGHTFILE:
		case ANNOTFILE:
		case FILTERFILE:
		case LOCUSFILTERFILE:
		case LOCUSWEIGHTFILE:
		case LOCUSNAMEFILE:
		case SAMPLEFILE:
		case CASEFREQFILE:
		case CONTFREQFILE:
		case OUTFILE:
		case SCOREFILE:
		case TRIOFILE:
			if (getNextArg(arg, argc, argv, fp,&arg_depth, &arg_num) == 0 || arg[0]=='-' || sscanf(arg, "%s",spi->df[a].fn) != 1)
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
	if ((spi->df[PSDATAFILE].fn[0]==0)+(spi->df[GCDATAFILE].fn[0]==0)+(spi->df[GENDATAFILE].fn[0]==0) < 2)
	{
		dcerror(1,"Must specify only one of --psdatafile, --gcdatafile or --gendatafile"); exit(1);
	}
	if ((spi->df[GCDATAFILE].fn[0]!='\0' || spi->df[GENDATAFILE].fn[0]!='\0' ) && pi->nloci==0)
	{
		dcerror(1,"Need to specify --nloci when using --gcdatafile or --gendatafile"); exit(1);
	}
	if (spi->df[GENDATAFILE].fn[0] != '\0')
	{
		if (spi->do_recessive_test)
		{
			dcerror(1, "Cannot do recessive test with --gendatafile"); exit(1);
		}
		if (spi->use_trios)
		{
			dcerror(1, "Cannot use trios with --gendatafile"); exit(1);
		}
		spi->use_probs = 1;
	}

	for (l = 0; l < pi->nloci; ++l)
		pi->n_alls[l]=2;
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

