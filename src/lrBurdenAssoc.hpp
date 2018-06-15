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

#ifndef LRBURDENASSOCHPP
#define LRBURDENASSOCHPP 1
#include <string>
#include <map>

extern "C" {
#include "sagcutils.h"
};
#include "lrModel.hpp"

#define MAX_COMMENT_LENGTH 100
#define NAME_LENGTH 20

enum OPT {
	FLAGFILE = 0, PSDATAFILE, GCDATAFILE, GENDATAFILE, WEIGHTFILE, ANNOTFILE, FILTERFILE, LOCUSFILTERFILE, LOCUSNAMEFILE, LOCUSWEIGHTFILE, LOCUSGROUPFILE, SAMPLEFILE, TRIOFILE, OUTFILE, SCOREFILE, NUMDATAFILETYPES, NUMLOCI, WEIGHTFACTOR, ARGFILE, NUMOPTS
};

// FLAGFILE is only useed by getVarScores but handy to have it here

struct option_t {
	char *str;
	OPT o;
};

typedef struct option_t option;
extern option opt[];

class lrba_data_file_type {
public:
	FILE *fp;
	char fn[200];
	lrba_data_file_type() { fp = 0; fn[0] = '\0'; }
	~lrba_data_file_type() { if (fp != 0 && fp != stdout) fclose(fp); }
};

#define MAXLRVARIABLENAMELENGTH 100
struct lrVariable_t { char name[MAXLRVARIABLENAMELENGTH]; float *val; };
typedef struct lrVariable_t lrVariable;
#define MAXLRVARIABLES 50
extern lrVariable allVars[MAXLRVARIABLES];
extern std::map<std::string, lrVariable *> varMap;

struct lrba_par_info_t {
	lrba_data_file_type df[NUMDATAFILETYPES];
	float wfactor;
	int use_func_weights, use_locus_names, use_comments, use_trios, use_probs;
	int do_ttest, do_lrtest, num_vars,num_varFiles,num_testFiles;
	lrba_data_file_type varFiles[MAXLRVARIABLES], testFiles[MAXLRVARIABLES];
};

typedef struct lrba_par_info_t lrba_par_info;

enum { DE_NOVO = 0, NON_MENDELIAN };
struct non_mendelian_t {
	int loc, sub, nd;
};
typedef struct non_mendelian_t non_mendelian;

extern float weight[MAX_LOCI], missing_score[MAX_LOCI], func_weight[MAX_LOCI], cc_freq[2][MAX_LOCI], cc_count[2][MAX_LOCI], cc_genocount[2][3][MAX_LOCI];
extern int rarer[MAX_LOCI], group[MAX_LOCI], max_cc[2]; // may be able to group by e.g variant type, separae beta for each
extern char names[MAX_LOCI][NAME_LENGTH], comments[MAX_LOCI][MAX_COMMENT_LENGTH], trios_fn[500];
extern char long_line[LONG_LINE_LENGTH + 1];

void usage();
int read_all_args(char *argv[], int argc, par_info *pi, lrba_par_info *spi);
int process_options(par_info *pi, lrba_par_info *spi);
int read_all_data(par_info *pi, lrba_par_info *spi, subject **sub, int *nsubptr, char names[MAX_LOCI][20], char comments[MAX_LOCI][MAX_COMMENT_LENGTH], float func_weight[MAX_LOCI]);
float get_quadratic_weight(float freq, float wfactor);
void set_weights(FILE *f, float *weight, float *missing_score, int *rarer, subject **sub, int nsub, par_info *pi, lrba_par_info *spi, float *func_weight, float cc_freq[2][MAX_LOCI], float cc_count[2][MAX_LOCI], int max_cc[2], char names[MAX_LOCI][20], char comments[MAX_LOCI][MAX_COMMENT_LENGTH]);
void get_freqs(subject **sub, int nsub, par_info *pi, lrba_par_info *spi, float cc_freq[2][MAX_LOCI], float cc_count[2][MAX_LOCI], float cc_genocount[2][3][MAX_LOCI]);
void get_scores(float *score, float *weight, float *missing, int *rarer, subject **sub, int nsub, par_info *pi, lrba_par_info *spi);
void fillModel(lrModel *m, float *score, subject **sub, int nsub, par_info *pi, lrba_par_info *spi);
double do_onetailed_LRT(FILE *fo, lrModel *m);
double getlnLikeForH(int h, lrModel *m);

#endif