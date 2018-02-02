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
	{ "filterfile", FILTERFILE },
	{ "annotfile", ANNOTFILE },
	{ "locusfilterfile", LOCUSFILTERFILE },
	{ "locusnamefile", LOCUSNAMEFILE },
	{ "locusweightfile", LOCUSWEIGHTFILE },
	{ "casefreqfile", CASEFREQFILE },
	{ "contfreqfile", CONTFREQFILE },
	{ "samplefile", SAMPLEFILE },
	{ "triofile", TRIOFILE },
	{ "outfile", OUTFILE },
	{ "scorefile", SCOREFILE },
	{ "nostringtomatchthis", NUMDATAFILETYPES },
	{ "numloci", NUMLOCI },
// need this if using gc files
	{ "ldthreshold", LDTHRESHOLD },
	{ "minweight", WEIGHTTHRESHOLD },
	{ "dorecessive", DORECESSIVE },
	{ "dottest",DOTTEST },
	{ "dolrtest",DOLRTEST },
	{ "usehaps", USEHAPS },
	{ "weightfactor", WEIGHTFACTOR },
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
"--dorecessive x\n"
"--dottest x\n"
"--dolrtest x\n"
"--varfile file\n"
"--testfile file\n"
"--usehaps x\n"
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

#define MISSING -999
int readVarFiles(subject **sub,int nSub,sa_par_info *spi)
{
	int i,s,idCol,c,colIndex[MAXLRVARIABLES],nCol;
	char colValue[MAXLRVARIABLENAMELENGTH],*ptr,*sptr;
	std::map<std::string,int> subIDs;
	if(spi->numVarFiles>0)
	{
		for(s=0;s<nSub;++s)
			subIDs[sub[s]->id]=s;
	}
	for(i=0;i<spi->numVarFiles;++i)
	{
		spi->varFiles[i].fp=fopen(spi->varFiles[i].fn,"r");
		if(spi->varFiles[i].fn==0)
		{
			dcerror(1,"Could not open variable file %s\n",spi->varFiles[i].fn);
			return 0;
		}
		fgets(long_line,LONG_LINE_LENGTH,spi->varFiles[i].fp);
		ptr=long_line;
		idCol=-1;
		while(*ptr && isspace(*ptr))
			++ptr;
		for (c=0;1;++c)
		{
			sptr=colValue;
			while(*ptr && isspace(*ptr))
				++ptr;
			*sptr='\0';
			while(*ptr && !isspace(*ptr))
				*sptr++=*ptr++;
			if(!strcmp(colValue,"IID"))
			{
				idCol=c;
				break;
			}
			if(*ptr=='\0')
				break;
		} 
		if(idCol==-1)
		{
			dcerror(1,"Variable file %s did not contain a column headed IID\n",spi->varFiles[i].fn);
			return 0;
		}
		for(c=idCol;1;++c)
		{
			sptr=colValue;
			while(*ptr && !isspace(*ptr))
				*sptr++=*ptr++;
			*sptr='\0';
			while(*ptr && isspace(*ptr))
				++ptr;
			if(colValue[0])
			{
				if (!strcmp(colValue,"score"))
				{
					dcerror(1, "Variable file %s contains a column headed score, which is not allowed because program will calculate scores\n", spi->varFiles[i].fn);
					return 0;
				}
				std::map<std::string,lrVariable *>::const_iterator varIter=varMap.find(colValue);
				if(varIter!=varMap.end())
				{
					dcerror(1,"Variable file %s contains a column headed %s but that variable name is already in use\n",spi->varFiles[i].fn,colValue);
					return 0;
				}
				strcpy(allVars[spi->numVars].name,colValue);
				assert((allVars[spi->numVars].val=(float*)calloc(nSub,sizeof(float)))!=0);
				for(s=0;s<nSub;++s)
					allVars[spi->numVars].val[s]=MISSING;
				colIndex[c]=spi->numVars;
				varMap[colValue]=allVars+spi->numVars;
				++spi->numVars;
			}
			if(*ptr=='\0')
				break;
		}
		nCol=c;
		while(fgets(long_line,LONG_LINE_LENGTH,spi->varFiles[i].fp))
		{
			while(*ptr && !isspace(*ptr))
				*sptr++=*ptr++;
			for(c=0;c<idCol;++c)
			{
				sptr=colValue;
				while(*ptr && !isspace(*ptr))
					*sptr++=*ptr++;
				*sptr='\0';
				while(*ptr && isspace(*ptr))
					++ptr;
				if(*ptr=='\0')
					break;
			}
			if(colValue[0]=='\0')
			{
				dcerror(1,"IID value missing from variable file %s in this line:\n%s\n",spi->varFiles[i].fn,long_line);
				return 0;
			}
			std::map<std::string,int>::const_iterator idIter=subIDs.find(colValue);
			if(idIter==subIDs.end())
			{
				dcerror(1,"Unknown IID value in variable file %s in this line:\n%s\n",spi->varFiles[i].fn,long_line);
				return 0;
			}
			else
				s=idIter->second;
			for(;c<nCol;++c)
			{
				sptr=colValue;
				while(*ptr && !isspace(*ptr))
					*sptr++=*ptr++;
				*sptr='\0';
				while(*ptr && isspace(*ptr))
					++ptr;
				if (sscanf(colValue,"%f",&allVars[colIndex[c]].val[s])!=1)
				{
					dcerror(1,"Not enough values in variable file %s in this line:\n%s\n",spi->varFiles[i].fn,long_line);
					return 0;
				}
			}
		}
		for(s=0;s<nSub;++s)
			if(allVars[idCol+1].val[s]==MISSING)
			{
				dcerror(1,"Missing values in variable file %s for subject %s\n",spi->varFiles[i].fn,sub[s]->id);
				return 0;
			}
		fclose(spi->varFiles[i].fp);
		spi->varFiles[i].fp=0;
	}
	return 1;
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
	spi->do_ttest = 1;
	spi->do_lrtest = 0;
	spi->numVars=spi->numVarFiles=spi->numTestFiles=0;
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
			if(getNextArg(arg,argc,argv,fp,&arg_depth,&arg_num) == 0 || sscanf(arg,"%d",&spi->do_lrtest)!=1)
				error=1;
			break;
		case USEHAPS:
			if(getNextArg(arg,argc,argv,fp,&arg_depth,&arg_num) == 0 || sscanf(arg,"%d",&spi->use_haplotypes)!=1)
				error=1;
			break;
		case VARFILE:
			if(getNextArg(arg,argc,argv,fp,&arg_depth,&arg_num) == 0 || sscanf(arg,"%s",spi->varFiles[spi->numVarFiles++])!=1)
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

int read_ps_datafile(par_info *pi, sa_par_info *spi, subject **sub, int *nsubptr, char names[MAX_LOCI][LOCUS_NAME_LENGTH], char comments[MAX_LOCI][MAX_COMMENT_LENGTH], float func_weight[MAX_LOCI],
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

int read_freqs_datafile(par_info *pi,sa_par_info *spi,int cc,float cc_freq[2][MAX_LOCI],float cc_count[2][MAX_LOCI],int max_cc[2])
{
	int f,l;
	FILE *fp;
	char *fn,*ptr;
	f=cc?CASEFREQFILE:CONTFREQFILE;
	fn=spi->df[f].fn;
	fp=spi->df[f].fp;
	spi->use_cc_freqs[cc]=1;
	if (!fgets(long_line,LONG_LINE_LENGTH,fp))
	{
		error("Could not read supplied frequencies in ",fn);
		return 0;
	}
	for (l=0, ptr=long_line;l<pi->nloci;++l)
	{
		if (sscanf(ptr,"%f",&cc_freq[cc][l])<1)
		{
			error("Could not read all supplied frequencies:\n",long_line);
			return 0;
		}
		while(isspace(*ptr++)) ;
		while(!isspace(*ptr++)) ;
	}
	if(!fgets(long_line,LONG_LINE_LENGTH,fp))
	{
		error("Could not read subject counts in ",fn);
		return 0;
	}
	max_cc[cc]=0;
	for(l=0,ptr=long_line;l<pi->nloci;++l)
	{
		if(sscanf(ptr,"%f",&cc_count[cc][l])<1)
		{
			error("Could not read all subject counts:\n",long_line);
			return 0;
		}
		if(max_cc[cc]<cc_count[cc][l])
			max_cc[cc]=cc_count[cc][l];
		while(isspace(*ptr++)) ;
		while(!isspace(*ptr++)) ;
	}
	return 1;
}

int read_all_data(par_info *pi,sa_par_info *spi,subject **sub,int *nsubptr,char names[MAX_LOCI][LOCUS_NAME_LENGTH],char comments[MAX_LOCI][MAX_COMMENT_LENGTH],float func_weight[MAX_LOCI])
{
	char aline[1000],pos[100],effect[100],*ptr;
	int func_pos,l,f,use;
	float wt;
	std::map<std::string,float> weightMap;
	std::map<std::string,std::string> effectMap;
	if (spi->df[WEIGHTFILE].fp)
	{
		while (fgets(aline, 999, spi->df[WEIGHTFILE].fp))
		{
			sscanf(aline,"%s %f",effect,&wt);
			weightMap[effect]=wt;
		}
	}
	if (spi->df[ANNOTFILE].fp)
	{
		fgets(aline, 999, spi->df[ANNOTFILE].fp); // ignore first line (though could use it to see how many cohorts there are)
		for (func_pos=0,ptr=aline;strncmp(ptr,"FUNC",4);ptr=skip_word(ptr))
			if (!*ptr)
			{
				dcerror(1, "Could not find FUNC in annotation file %s:\n%s\n", spi->df[ANNOTFILE].fn, aline); exit(1);
			}
			else
				++func_pos;
		while (fgets(aline, 999, spi->df[ANNOTFILE].fp))
		{
			sscanf(aline,"%s",pos);
			ptr=aline;
			for (f=0;f<func_pos;++f)
				ptr=skip_word(ptr);
			sscanf(ptr,"%s",effect);
			effectMap[pos]=effect;
		}
	}
	if (spi->df[PSDATAFILE].fp)
		read_ps_datafile(pi,spi,sub,nsubptr,names,comments, func_weight,weightMap,effectMap);
	else if (spi->df[GCDATAFILE].fp)
		read_all_subjects(spi->df[GCDATAFILE].fp,sub,nsubptr,pi);
	else if (spi->df[GENDATAFILE].fp)
		read_all_gen_subjects(spi->df[GENDATAFILE].fp,sub,nsubptr,pi);
	if (spi->df[LOCUSFILTERFILE].fp)
	{
		pi->n_loci_to_use=0;
		for (l = 0; l < pi->nloci; ++l)
		{
			if (fscanf(spi->df[LOCUSFILTERFILE].fp, " %d", &use) != 1)
			{
				dcerror(1, "Not enough values in locusfilterfile %s\n", spi->df[LOCUSFILTERFILE].fn); exit(1);
			}
			else if (use==1)
				pi->loci_to_use[pi->n_loci_to_use++]=l;
		}
	}
	else
	{
		pi->n_loci_to_use = pi->nloci;
		for (l = 0; l < pi->nloci; ++l)
			pi->loci_to_use[l] = l;
	}
	if (spi->df[LOCUSWEIGHTFILE].fp)
	{
		for (l = 0; l < pi->nloci; ++l)
			if (fscanf(spi->df[LOCUSWEIGHTFILE].fp,"%f ",&func_weight[l])!=1)
			{
				dcerror(1, "Not enough values in locusweightfile %s\n", spi->df[LOCUSWEIGHTFILE].fn); exit(1);
			}
	}
	else if (spi->df[WEIGHTFILE].fp ==0)
		for (l = 0; l < pi->nloci; ++l)
			func_weight[l]=1;
	if (spi->df[LOCUSNAMEFILE].fp)
	{
		for (l = 0; l < pi->nloci; ++l)
			if (fscanf(spi->df[LOCUSNAMEFILE].fp,"%s ",&comments[l])!=1)
			{
				dcerror(1, "Not enough values in locusnamefile %s\n", spi->df[LOCUSNAMEFILE].fn); exit(1);
			}
	}
	if(spi->df[CONTFREQFILE].fp)
	{
		read_freqs_datafile(pi,spi,0,cc_freq,cc_count,max_cc);
	}
	if(spi->df[CASEFREQFILE].fp)
	{
		read_freqs_datafile(pi,spi,1,cc_freq,cc_count,max_cc);
	}
	for (l = 0; l < pi->nloci; ++l)
	{
		strncpy(names[l],comments[l],NAME_LENGTH-1);
		names[l][NAME_LENGTH-1]='\0';
	}
	for(f=0;f<spi->numVarFiles;++f)
		if(!readVarFiles(sub,*nsubptr,spi))
			exit(1);
	return 1;
}

