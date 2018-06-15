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
#include "scoreassoc.hpp"
#include "dcexpr.hpp"
#include "safilterfuncs.hpp"

#define PROGRAM "scoreassoc"
#define SAVERSION "5.0"

int main(int argc, char *argv[])
{
	char arg_string[2000];
	int nsub,n_new_sub,real_nsub,filledModel;
	float *score,SLP,p;
	int s,n_non_mendelian,t;
	non_mendelian *non_mendelians;
	char *non_mendelian_report;
	par_info pi;
	sa_par_info spi;
	subject **sub,**new_sub,**real_sub;
	lrRidgePenaltyModel model;
	pi.use_cc=1;
	printf("%s v%s\n",PROGRAM,SAVERSION);
	printf("MAX_LOCI=%d\nMAX_SUB=%d\n",MAX_LOCI,MAX_SUB);

	assert(sub=(subject **)calloc(MAX_SUB,sizeof(subject*)));
	for (s=0;s<MAX_SUB;++s)
		assert(sub[s]=(subject *)calloc(1,sizeof(subject)));
	assert(score=(float *)calloc(MAX_SUB,sizeof(float)));
	max_cc[0]=max_cc[1]=0;
	read_all_args(argv,argc, &pi, &spi);
	// make_arg_string(arg_string,argc,argv);
	// parse_arg_string(arg_string,&pi,&spi,&pspi);
	process_options(&pi,&spi);
	if (spi.df[FILTERFILE].fp)
		initExclusions(spi.df[FILTERFILE].fp);
	read_all_data(&pi,&spi,sub,&nsub,names,comments,func_weight);
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
fprintf(spi.df[OUTFILE].fp,"scoreassoc output\n"
"Locus                                             contAA  contAB  contBB  contFreq  caseAA  caseAB  caseBB  caseFreq  MAF       rarer  weight   %s\n",
	spi.use_comments ? "comment" : "");
get_freqs(sub,nsub,&pi,&spi,cc_freq,cc_count,cc_genocount);
applyExclusions(&pi);
set_weights(spi.df[OUTFILE].fp,weight,missing_score,rarer,sub,nsub,&pi,&spi,func_weight,cc_freq,cc_count,max_cc,names,comments);
get_scores(score,weight,missing_score,rarer,sub,nsub,&pi,&spi);
filledModel=0;
strcpy(allVars[spi.numVars].name, "score");
spi.scoreCol=spi.numVars;
allVars[spi.scoreCol].val = score;
varMap["score"] = &allVars[spi.numVars];
++spi.numVars;

if (spi.do_ttest)
	SLP=do_score_onetailed_ttest(spi.df[OUTFILE].fp,score,sub,nsub,&pi,&spi,cc_freq,cc_count,max_cc,weight,missing_score,rarer);

model.lamda = spi.lamda;
if(spi.do_lrtest)
{
	if (!filledModel)
	{
		fillModelWithVars(&model, nsub, &spi);
		for (s = 0; s < nsub; ++s)
			model.Y[s] = sub[s]->cc;
		filledModel = 1;
	}
	SLP = do_onetailed_LRT(spi.df[OUTFILE].fp,&model,&spi);
}
if (spi.numTestFiles>0)
{
	if (!filledModel)
	{
		fillModelWithVars(&model, nsub, &spi);
		for (s = 0; s < nsub; ++s)
			model.Y[s] = sub[s]->cc;
		filledModel = 1;
	}
	for (t = 0; t < spi.numTestFiles;++t)
		p= runTestFile(spi.df[OUTFILE].fp, spi.testFiles[t].fn,&model, &spi);

}

if (spi.use_trios)
	output_nm_report(spi.df[OUTFILE].fp,&pi,non_mendelian_report);

if (spi.df[SCOREFILE].fp)
{
	write_scores(spi.df[SCOREFILE].fp,sub,nsub,score);
	fclose(spi.df[SCOREFILE].fp);
	spi.df[SCOREFILE].fp=0;
}
if (spi.do_recessive_test)
{
if (atoi(comments[0])>22 || toupper(comments[0][0]) == 'X' || toupper(comments[0][0]) == 'Y' ||
		toupper(comments[0][0]) == 'C'&&toupper(comments[0][1]) == 'H'&&toupper(comments[0][2]) == 'R' && (atoi(comments[0] + 3) > 22 || toupper(comments[0][3]) == 'X' || toupper(comments[0][3]) == 'Y'))
		// simple trick for now to avoid X and Y genes
	fprintf(spi.df[OUTFILE].fp,"\nCannot do recessive test for genes on X or Y chromosome.\n");
else if (spi.use_haplotypes)
	do_recessive_HWE_test_with_haplotypes(spi.df[OUTFILE].fp,score,sub,nsub,&pi,&spi,cc_freq,cc_count,max_cc,weight,missing_score,rarer,names);
else 
	do_recessive_HWE_test(spi.df[OUTFILE].fp,score,sub,nsub,&pi,&spi,cc_freq,cc_count,max_cc,weight,missing_score,rarer,names);
}

stateExclusions(spi.df[OUTFILE].fp);
fclose(spi.df[OUTFILE].fp);
spi.df[OUTFILE].fp=0; //  because otherwise the destructor will try to fclose it
printf("\nProgram run completed OK\n");
return 0;

}


