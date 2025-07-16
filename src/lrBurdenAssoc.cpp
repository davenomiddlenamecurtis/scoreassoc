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
#include <ctype.h>
#include "lrBurdenAssoc.hpp"
#include "lrBAFilterFuncs.hpp"
#include "lrModel.hpp"
#include "dcerror.hpp"

#define error(s,t) error_func(__LINE__,__FILE__,s,t)
#define PROGRAM "lrBurdenAssoc"
#define LRBAVERSION "1.0"

int main(int argc, char *argv[])
{
	char arg_string[2000];
	int nsub, n_new_sub, real_nsub;
	float *score, SLP;
	int s, n_non_mendelian;
	non_mendelian *non_mendelians;
	char *non_mendelian_report;
	par_info pi;
	lrba_par_info spi;
	lrModel model;
	subject **sub, **new_sub, **real_sub;
	pi.use_cc = 1;
	printf("%s v%s\n", PROGRAM, LRBAVERSION);
	printf("MAX_LOCI=%d\nMAX_SUB=%d\n", MAX_LOCI, MAX_SUB);

	assert(sub = (subject **)calloc(MAX_SUB, sizeof(subject*)));
	for (s = 0; s<MAX_SUB; ++s)
		assert(sub[s] = (subject *)calloc(1, sizeof(subject)));
	assert(score = (float *)calloc(MAX_SUB, sizeof(float)));
	max_cc[0] = max_cc[1] = 0;
	read_all_args(argv, argc, &pi, &spi);
	// make_arg_string(arg_string,argc,argv);
	// parse_arg_string(arg_string,&pi,&spi,&pspi);
	process_options(&pi, &spi);
	if (spi.df[FILTERFILE].fp)
		initExclusions(spi.df[FILTERFILE].fp);
	read_all_data(&pi, &spi, sub, &nsub, names, comments, func_weight);
	fprintf(spi.df[OUTFILE].fp,
	"Locus                                   contAA  contAB  contBB  contFreq  caseAA  caseAB  caseBB  caseFreq  MAF       rarer  weight   %s\n",
		spi.use_comments ? "comment" : "");
	get_freqs(sub, nsub, &pi, &spi, cc_freq, cc_count, cc_genocount);
	applyExclusions(&pi);
	set_weights(spi.df[OUTFILE].fp, weight, missing_score, rarer, sub, nsub, &pi, &spi, func_weight, cc_freq, cc_count, max_cc, names, comments);
	set_weights(spi.df[OUTFILE].fp, weight, missing_score, rarer, sub, nsub, &pi, &spi, func_weight, cc_freq, cc_count, max_cc, names, comments);
	get_scores(score, weight, missing_score, rarer, sub, nsub, &pi, &spi);
	fillModel(&model,score,sub,nsub, &pi, &spi);
	SLP = do_onetailed_LRT(spi.df[OUTFILE].fp, &model);
	stateExclusions(spi.df[OUTFILE].fp);
	fclose(spi.df[OUTFILE].fp);
	spi.df[OUTFILE].fp = 0; //  because otherwise the destructor will try to fclose it
	printf("\nProgram run completed OK\n");
	return 0;
}