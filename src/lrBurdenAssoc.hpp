#ifndef LRBURDENASSOCHPP
#define LRBURDENASSOCHPP 1

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

struct lrba_par_info_t {
	lrba_data_file_type df[NUMDATAFILETYPES];
	float wfactor;
	int use_func_weights, use_locus_names, use_comments, use_trios, use_probs;
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