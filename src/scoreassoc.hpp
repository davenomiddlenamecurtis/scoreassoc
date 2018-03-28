#ifndef SCOREASSOCHPP
#define SCOREASSOCHPP

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <string>
#include <map>
#include "lrModel.hpp"

extern "C" {
#include "sagcutils.h"
};

#define MAX_COMMENT_LENGTH 1000
#define NAME_LENGTH 20
#define MAXLRVARIABLES 50
#define MAXLRVARIABLENAMELENGTH 100
#define LOCUS_NAME_LENGTH 40
#define LOCUS_NAME_LENGTH_STR "40"

enum OPT {
	FLAGFILE=0,PSDATAFILE,GCDATAFILE,GENDATAFILE,WEIGHTFILE,ANNOTFILE,FILTERFILE,LOCUSFILTERFILE,LOCUSNAMEFILE,LOCUSWEIGHTFILE,TRIOFILE,SAMPLEFILE,CASEFREQFILE,CONTFREQFILE,OUTFILE,SCOREFILE,NUMDATAFILETYPES,NUMLOCI,LDTHRESHOLD,WEIGHTTHRESHOLD,DORECESSIVE,DOTTEST,DOLRTEST,USEHAPS,WEIGHTFACTOR,VARFILE,TESTFILE,ARGFILE,NUMOPTS
};

// FLAGFILE is only used by getVarScores but handy to have it here

struct option_t {
	char *str;
	OPT o;
};

typedef struct option_t option;
extern option opt[];

class sa_data_file_type {
public:
	FILE *fp;
	char fn[200];
	sa_data_file_type() { fp = 0; fn[0] = '\0'; }
	~sa_data_file_type() { if (fp != 0 && fp != stdout) fclose(fp); }
};

struct sa_par_info_t {
sa_data_file_type df[NUMDATAFILETYPES];
float wfactor,LD_threshold,weight_threshold;
int use_func_weights,use_cc_freqs[2],use_locus_names,use_comments,do_recessive_test,use_haplotypes,use_trios,use_probs;
int do_ttest, do_lrtest,numVars,numVarFiles,numTestFiles;
sa_data_file_type varFiles[MAXLRVARIABLES],testFiles[MAXLRVARIABLES];
};

typedef struct sa_par_info_t sa_par_info;

enum { DE_NOVO=0, NON_MENDELIAN };
struct non_mendelian_t {
	int loc,sub,nd;
};
typedef struct non_mendelian_t non_mendelian;

class lrVariable {
public:
	char name[MAXLRVARIABLENAMELENGTH];
	float *val;
	lrVariable() { name[0]='\0'; val=0; }
	void clear() { name[0]='\0'; if(val) free(val); val=0; }
	~lrVariable() { if(val) free(val); }
};
extern lrVariable allVars[MAXLRVARIABLES];
extern std::map<std::string,lrVariable *> varMap;


void usage();
int read_all_args(char *argv[], int argc, par_info *pi, sa_par_info *spi);
int process_options(par_info *pi, sa_par_info *spi);
int read_all_data(par_info *pi, sa_par_info *spi, subject **sub, int *nsubptr, char names[MAX_LOCI][LOCUS_NAME_LENGTH], char comments[MAX_LOCI][MAX_COMMENT_LENGTH], float func_weight[MAX_LOCI]);

double do_score_onetailed_ttest(FILE *fo,float *score,subject **sub,int nsub,par_info *pi,sa_par_info *spi,float cc_freq[2][MAX_LOCI],float cc_count[2][MAX_LOCI],int max_cc[2],float *weight,float *missing,int *rarer);
void get_scores(float *score,float *weight,float *missing_score,int *rarer,subject **sub,int nsub,par_info *pi,sa_par_info *spi);
void set_weights(FILE *fo,float *weight,float *missing_score,int *rarer,subject **sub,int nsub,par_info *pi,sa_par_info *spi,float *func_weight,float cc_freq[2][MAX_LOCI],float cc_count[2][MAX_LOCI],int max_cc[2],char names[MAX_LOCI][LOCUS_NAME_LENGTH],char comments[MAX_LOCI][MAX_COMMENT_LENGTH]);
int read_score_assoc_par(FILE *fp,par_info *pi,float *wfactor,int *use_func_weights,int use_cc_freqs[2],float *func_weight,float cc_freq[2][MAX_LOCI],float cc_count[2][MAX_LOCI],int max_cc[2],int *use_locus_names,char names[MAX_LOCI][LOCUS_NAME_LENGTH],int *use_comments,char comments[MAX_LOCI][MAX_COMMENT_LENGTH], int *do_recessive_test,float *weight_threshold,float *LD_threshold,int *use_haplotypes);
int read_sa_par(FILE *fp,par_info *pi,sa_par_info *spi,float *func_weight,float cc_freq[2][MAX_LOCI],float cc_count[2][MAX_LOCI],int max_cc[2],char names[MAX_LOCI][LOCUS_NAME_LENGTH],char comments[MAX_LOCI][MAX_COMMENT_LENGTH]);
void get_freqs(subject **sub,int nsub,par_info *pi,sa_par_info *spi,float cc_freq[2][MAX_LOCI],float cc_count[2][MAX_LOCI],float cc_gencount[2][3][MAX_LOCI]);
void write_scores(FILE *fs,subject **sub,int nsub,float *score);
void do_recessive_HWE_test(FILE *fo,float *score,subject **sub,int nsub,par_info *pi,sa_par_info *spi,float cc_freq[2][MAX_LOCI],float cc_count[2][MAX_LOCI],int max_cc[2],float *weight,float *missing,int *old_rarer,char names[MAX_LOCI][LOCUS_NAME_LENGTH]);
void do_recessive_HWE_test_with_haplotypes(FILE *fo, float *score, subject **sub, int nsub, par_info *pi,sa_par_info *spi, float cc_freq[2][MAX_LOCI], float cc_count[2][MAX_LOCI], int max_cc[2], float *weight, float *missing, int *old_rarer, char names[MAX_LOCI][LOCUS_NAME_LENGTH]);
extern int sort_trios(subject **sub,int nsub,par_info *pi, sa_par_info *spi, subject **new_sub,non_mendelian *nm,int *n_non_mendelian,char *non_mendelian_report);
extern int output_nm_report(FILE *fp, par_info *pi, char *non_mendelian_report);
int read_all_gen_subjects(FILE *fi,subject **s,int *nsub,par_info *pi);
int read_freqs_datafile(par_info *pi,sa_par_info *spi,int cc,float cc_freq[2][MAX_LOCI],float cc_count[2][MAX_LOCI],int max_cc[2]);
float do_onetailed_LRT(FILE *fo,lrModel *m,par_info *pi,sa_par_info *spi);
void fillModelWithVars(lrModel *m,subject **sub,int nsub,par_info *pi,sa_par_info *spi);
void printModel(FILE *fo, char *LLstr, double LL, lrModel *m);
extern double cumulBinom(int N,int k,double p);
float evaluateModel(FILE *fo, lrModel *m, int *toUse, float *startBetas, int *toFit, char *name);
float runTestFile(FILE *fo, char *fn, lrModel *m, par_info *pi, sa_par_info *spi);
char *skip_word(char *ptr);
int read_ps_datafile(par_info *pi, sa_par_info *spi, subject **sub, int *nsubptr, char names[MAX_LOCI][LOCUS_NAME_LENGTH], char comments[MAX_LOCI][MAX_COMMENT_LENGTH], float func_weight[MAX_LOCI],
	std::map<std::string, float> weightMap, std::map<std::string, std::string> effectMap);
int readVarFiles(subject **sub, int nSub, sa_par_info *spi);

extern float weight[MAX_LOCI],missing_score[MAX_LOCI],func_weight[MAX_LOCI],cc_freq[2][MAX_LOCI],cc_count[2][MAX_LOCI],cc_genocount[2][3][MAX_LOCI];
extern int rarer[MAX_LOCI],max_cc[2];
extern char names[MAX_LOCI][LOCUS_NAME_LENGTH],comments[MAX_LOCI][MAX_COMMENT_LENGTH],trios_fn[500];

#define USEFILTERS 1
#endif