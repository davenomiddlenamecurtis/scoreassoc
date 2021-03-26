#if 0
Copyright 2020 David Curtis

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
#define SAVERSION "7.0"

int main(int argc, char *argv[])
{
	char arg_string[2000];
	int nsub,n_new_sub,real_nsub,filledModel;
	double** score;
	float SLP, p;
	int s,n_non_mendelian,t,sc;
	non_mendelian *non_mendelians;
	par_info pi;
	sa_par_info spi;
	subject **sub,**new_sub,**real_sub;
	glRidgePenaltyModel model;
	printf("%s v%s\n",PROGRAM,SAVERSION);
	printf("MAX_LOCI=%d\nMAX_SUB=%d\n",MAX_LOCI,MAX_SUB);

	assert(sub=(subject **)calloc(MAX_SUB,sizeof(subject*)));
	max_cc[0]=max_cc[1]=0;
	read_all_args(argv,argc, &pi, &spi);
	// make_arg_string(arg_string,argc,argv);
	// parse_arg_string(arg_string,&pi,&spi,&pspi);
	process_options(&pi,&spi); // use this to get number of weight files and set numScore here
	assert(score = (double**)calloc(numScores, sizeof(double*)));
	assert(weight = (double**)calloc(numScores, sizeof(double*)));
	assert(func_weight = (double**)calloc(numScores, sizeof(double*)));
	assert(missing_score = (float**)calloc(numScores, sizeof(double*)));
	for (sc = 0; sc < numScores; ++sc)
	{
		assert(score[sc] = (double*)calloc(MAX_SUB, sizeof(double)));
		assert(weight[sc] = (double*)calloc(pi.nloci? pi.nloci:MAX_LOCI, sizeof(double)));
		assert(func_weight[sc] = (double*)calloc(pi.nloci ? pi.nloci : MAX_LOCI, sizeof(double)));
		assert(missing_score[sc] = (float*)calloc(pi.nloci ? pi.nloci : MAX_LOCI, sizeof(double)));
	}
	// at this point, have allocated func_weight, weight and score - need to free them later
	try
	{
		for (s = 0; s < MAX_SUB; ++s)
		{
			if (spi.df[PSDATAFILE].fn[0]) {
				sub[s] = new subject(); // use default MAX_LOCI
			}
			else if (spi.use_probs) {
				sub[s] = new subject(pi.nloci, 1);
			}
			else if (spi.df[GENDATAFILE].fn[0]) {
				sub[s] = new subject(pi.nloci, 0);
			}
			else if (spi.useFlatFile) { // this has not been implemented
				sub[s] = new subject(0, 0);
			}
			else
				sub[s] = new subject(pi.nloci, 0); // could be new (std::nothrow) subject(pi.nloci, 0), then check for NULL pointer
		}
	}
	catch (std::bad_alloc& ba) // even though sizeof(subject) is small, call to new can sometimes throw exception
	{
		error("There is not enough memory for MAXSUB subjects and the requested number of loci", "");
		return 1;
	}

	if (spi.df[FILTERFILE].fp)
		initExclusions(spi.df[FILTERFILE].fp);
	read_all_data(&pi,&spi,sub,&nsub,names,comments,func_weight,score,numScores);
	if (nsub == 0)
	{
		error("There were zero subjects to input","");
		return 1;
	}
if (spi.use_trios)
{
	if (atoi(comments[0])>22 || toupper(comments[0][0]) == 'X' || toupper(comments[0][0]) == 'Y' ||
		toupper(comments[0][0]) == 'C'&&toupper(comments[0][1]) == 'H'&&toupper(comments[0][2]) == 'R' && (atoi(comments[0] + 3) > 22 || toupper(comments[0][3]) == 'X' || toupper(comments[0][3]) == 'Y'))
	{
		fprintf(spi.df[OUTFILE].fp, "Cannot at present use trios for genes on X or Y chromosome, but have this comment: %s", comments[0]);
		fclose(spi.df[OUTFILE].fp);
		error("Cannot at present use trios for genes on X or Y chromosome, but have this comment: ", comments[0]);
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
	}
else
	non_mendelians=0;
	// not used, compiler error otherwise
fprintf(spi.df[OUTFILE].fp, "scoreassoc output\n");
if (spi.df[INPUTSCOREFILE].fp==0)
{
	if (pi.is_quantitative)
		fprintf(spi.df[OUTFILE].fp, 
			"Locus                                         contAA : contAB : contBB  contFreq  meanAA   : meanAB   : meanBB    rarer  ");
	else
		fprintf(spi.df[OUTFILE].fp, 
			"Locus                                         contAA : contAB : contBB  contFreq  caseAA : caseAB : caseBB  caseFreq  MAF       rarer  ");
	if (spi.numScores == 1)
		fprintf(spi.df[OUTFILE].fp, "weight   ");
	else
		for (sc = 0; sc < spi.numScores; ++sc)
			fprintf(spi.df[OUTFILE].fp, "weight%-2d ", sc);
	fprintf(spi.df[OUTFILE].fp,"%s\n",spi.use_comments ? "comment" : "");
	get_freqs(sub, nsub, &pi, &spi, cc_freq, cc_count, cc_genocount);
	applyExclusions(sub, nsub, &pi);
	set_weights(spi.df[OUTFILE].fp, weight, missing_score, rarer, sub, nsub, &pi, &spi, func_weight, cc_freq, cc_count, max_cc, names, comments);
	get_scores(score, weight, missing_score, rarer, sub, nsub, &pi, &spi);
}
filledModel=0;
spi.scoreCol = spi.numVars; // first of possibly more than one score
for (sc = 0; sc < spi.numScores; ++sc)
{
	strcpy(allVars[spi.numVars].name, weightNames[sc]);
	allVars[spi.numVars].val = score[sc];
	varMap[allVars[spi.numVars].name] = &allVars[spi.numVars];
	++spi.numVars;
}
if (spi.do_lrtest || spi.do_linrtest || spi.numTestFiles > 0 || spi.numLinTestFiles > 0)
{
	if (!filledModel)
	{
		if (!fillModelWithVars(&model, nsub, &spi))
		{
			error("Failed to perform fillModelWithVars(&model, nsub, &spi), probably not enough memory","");
			return 1;
		}
	fillModelWithVars(&model, nsub, &spi);
	if (pi.is_quantitative)
		for (s = 0; s < nsub; ++s)
			model.Y[s] = sub[s]->pheno;
	else
		for (s = 0; s < nsub; ++s)
			model.Y[s] = sub[s]->cc;
	filledModel = 1;
	}
}
// above is here because can fail to allocate memory and I want to exit before producing t test output

if (spi.do_ttest)
	SLP=do_score_onetailed_ttest(spi.df[OUTFILE].fp,score,sub,nsub,&pi,&spi,cc_freq,cc_count,max_cc,weight,missing_score,rarer);

model.lamda = spi.lamda;
if (spi.do_lrtest)
{
	SLP = do_onetailed_LRT(spi.df[OUTFILE].fp, &model, &spi, 0);
}
if (spi.do_linrtest)
{
	SLP = do_onetailed_LRT(spi.df[OUTFILE].fp, &model, &spi, 1);
}

if (spi.numTestFiles > 0)
{
	for (t = 0; t < spi.numTestFiles; ++t)
		p = runTestFile(spi.df[OUTFILE].fp, spi.testFiles[t].fn, &model, &spi);
}

if (spi.numLinTestFiles > 0)
{
	for (t = 0; t < spi.numLinTestFiles; ++t)
		p = runLinTestFile(spi.df[OUTFILE].fp, spi.linTestFiles[t].fn, &model, &spi);
}

if (spi.use_trios)
	output_nm_report(spi.df[OUTFILE].fp, &pi, n_non_mendelian, non_mendelians);

if (spi.df[SCOREFILE].fp)
{
	write_scores(spi.df[SCOREFILE].fp,sub,nsub,score,spi.numScores,&pi);
	fclose(spi.df[SCOREFILE].fp);
	spi.df[SCOREFILE].fp=0;
}
#if 0
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
#endif

stateExclusions(spi.df[OUTFILE].fp);
fclose(spi.df[OUTFILE].fp);
spi.df[OUTFILE].fp=0; //  because otherwise the destructor will try to fclose it
for (s = 0; s < MAX_SUB; ++s)
	delete *(sub+s);
free(sub);
printf("\nProgram run completed OK\n");
return 0;

}


