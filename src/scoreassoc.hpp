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

#ifndef SCOREASSOCHPP
#define SCOREASSOCHPP

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <string>
#include <map>
#include "glModel.hpp"

#include "sagcutils.hpp"

#define MAX_COMMENT_LENGTH 15000 // this was set to 1500 but this would miss off some of the VEP annotation
#define NAME_LENGTH 20
#ifndef MAXLRVARIABLES
#define MAXLRVARIABLES 50
#endif
#define MAXLRVARIABLENAMELENGTH 100
#define LOCUS_NAME_LENGTH 45
#define LOCUS_NAME_LENGTH_STR "45"
#define DEFAULT_LAMDA 1

// NB these must match the declared list of options or NUMOPTS will be incorrect
enum OPT {
	PSDATAFILE=0,
	GCDATAFILE,
	GENDATAFILE,
	WEIGHTFILE,
	ANNOTFILE,
	FILTERFILE,
	LOCUSFILTERFILE,
	LOCUSNAMEFILE,
	LOCUSWEIGHTNAMEFILE,
	LOCUSWEIGHTFILE,
	LOCUSRECWEIGHTNAMEFILE,
	LOCUSRECWEIGHTFILE,
	TRIOFILE,
	CASEFREQFILE,
	CONTFREQFILE,
	INPUTSCOREFILE,
	IDPHENOTYPEFILE,
	OUTFILE,
	SCOREFILE,
	RECSCOREFILE,
	NUMDATAFILETYPES,
	NUMLOCI,
	// LDTHRESHOLD,
	// WEIGHTTHRESHOLD,
	ISQUANTITATIVE,
	// DORECESSIVE,
	MISSINGZERO,
	DOTTEST,
	DOLRTEST,
	DOLINRTEST,
	STARTFROMFITTED,
	// USEHAPS,
	// SHOWHAPLOCUSNAMES,
	WEIGHTFACTOR,
	MAXMAF,
	LAMDA,
	VARFILE, 
	TESTFILE, 
	LINTESTFILE, 
	TRANSPOSEDATA,
	NUMSCORES,
	ARGFILE,
	NUMOPTS, 
	FLAGFILE
};
// must begin with datafiles followed by output files

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

class lr_test_par_info {
public:
	float lamda;
	int numVars, numVarFiles, numTestFiles, numLinTestFiles, start_from_fitted,scoreCol;
	sa_data_file_type varFiles[MAXLRVARIABLES], testFiles[MAXLRVARIABLES], linTestFiles[MAXLRVARIABLES];
};

#define MAXWEIGHTFILES 100
class sa_par_info : public lr_test_par_info {
public:
sa_data_file_type df[NUMDATAFILETYPES];
sa_data_file_type locusWeightFile[MAXWEIGHTFILES], locusRecWeightFile[MAXWEIGHTFILES];
int numLocusWeightFiles=0,numLocusRecWeightFiles = 0;
float wfactor,LD_threshold,weight_threshold,max_MAF;
int use_func_weights,use_cc_freqs[2],use_locus_names,use_comments,do_recessive_test,use_haplotypes,show_hap_locus_names,use_trios,use_probs;
int do_ttest, do_lrtest,do_linrtest;
int missingZero = 0;
int useFlatFile = 0;
int useTransposedFile = 0;
float LDThreshold2022 = 1;
int numScores = -1, numAddScores = -1, numRecScores = -1;
};

enum { DE_NOVO=0, NON_MENDELIAN };
struct non_mendelian_t {
	int loc,sub,nd;
	char *report;
};
typedef struct non_mendelian_t non_mendelian;

class lrVariable {
public:
	char name[MAXLRVARIABLENAMELENGTH];
	double *val;
	lrVariable() { name[0]='\0'; val=0; }
	void clear() { name[0]='\0'; if(val) free(val); val=0; }
	~lrVariable() { if(val) free(val); }
};
extern lrVariable allVars[MAXLRVARIABLES];
extern std::map<std::string,lrVariable *> varMap;

typedef struct recPair_t { int l[2]; } recPair; //  loci contributing to hightest compound heterozygote score


void usage();
int read_all_args(char *argv[], int argc, par_info *pi, sa_par_info *spi);
int process_options(par_info *pi, sa_par_info *spi);
int read_all_data(par_info *pi, sa_par_info *spi, subject **sub, int *nsubptr, char names[MAX_LOCI][LOCUS_NAME_LENGTH], char comments[MAX_LOCI][MAX_COMMENT_LENGTH], double **func_weight,double **score,int numScores);

double do_score_onetailed_ttest(FILE *fo,double **score,subject **sub,int nsub,par_info *pi,sa_par_info *spi,float cc_freq[2][MAX_LOCI],float cc_count[2][MAX_LOCI],int max_cc[2],double **weight,float **missing,int *rarer);
void get_scores(double** score, double** weight, float** missing_score, int* rarer, subject** sub, int nsub, par_info* pi, sa_par_info* spi);
void get_rec_scores(double** score, recPair ** rP,double** weight, float** missing_score, int* rarer, subject** sub, int nsub, par_info* pi, sa_par_info* spi);
void set_weights(FILE *fo,double **weight,float **missing_score,int *rarer,subject **sub,int nsub,par_info *pi,sa_par_info *spi,double **func_weight,float cc_freq[2][MAX_LOCI],float cc_count[2][MAX_LOCI],int max_cc[2],char names[MAX_LOCI][LOCUS_NAME_LENGTH],char comments[MAX_LOCI][MAX_COMMENT_LENGTH]);
void get_freqs(subject **sub,int nsub,par_info *pi,sa_par_info *spi,float cc_freq[2][MAX_LOCI],float cc_count[2][MAX_LOCI],float cc_gencount[2][3][MAX_LOCI]);
void write_scores(FILE* fs, subject** sub, int nsub, double** score, int numScores, par_info* pi);
void write_rec_scores(FILE* fs, subject** sub, int nsub, double** score, recPair** rP, int numRecScores, double** weight, char names[MAX_LOCI][LOCUS_NAME_LENGTH], par_info* pi, sa_par_info* spi);
extern int sort_trios(subject **sub,int nsub,par_info *pi, sa_par_info *spi, subject **new_sub,non_mendelian *nm,int *n_non_mendelian,char *non_mendelian_report);
extern int output_nm_report(FILE *fp, par_info *pi, int n_non_mendelians, non_mendelian *non_mendelians);
int read_all_gen_subjects(FILE *fi,subject **s,int *nsub,par_info *pi);
int read_freqs_datafile(par_info *pi,sa_par_info *spi,int cc,float cc_freq[2][MAX_LOCI],float cc_count[2][MAX_LOCI],int max_cc[2]);
float do_onetailed_LRT(FILE *fo,glModel *m,lr_test_par_info *spi,int do_linear);
int fillModelWithVars(glModel *m,int nsub, lr_test_par_info *spi,int which=-1);
void printModel(FILE *fo, char *LLstr, double LL, glModel *m);
extern double cumulBinom(int N,int k,double p);
float evaluateModel(FILE *fo, glModel *m, int *toUse, float *startBetas, int *toFit, char *name);
float runTestFile(FILE* fo, char* fn, glModel* m, lr_test_par_info* spi);
float runLinTestFile(FILE* fo, char* fn, glModel* m, lr_test_par_info* spi);
char *skip_word(char *ptr);
int read_ps_datafile(par_info *pi, sa_par_info *spi, subject **sub, int *nsubptr, char names[MAX_LOCI][LOCUS_NAME_LENGTH], char comments[MAX_LOCI][MAX_COMMENT_LENGTH], double *func_weight,
	std::map<std::string, double> weightMap, std::map<std::string, std::string> effectMap);
int readVarFiles(std::map<std::string, int> subIDs, int nSub, lr_test_par_info *spi);
int read_all_subject_scores(FILE* fi, subject** s, int* nsub, double **score,int numScores);
int read_phenotypes(FILE* fi, subject** s, int *nsub, double **scores,int numScores,int isquantitative);

extern float **missing_score,cc_freq[2][MAX_LOCI],cc_count[2][MAX_LOCI],cc_genocount[2][3][MAX_LOCI];
extern double **weight, **func_weight;
extern int rarer[MAX_LOCI],max_cc[2];
extern char names[MAX_LOCI][LOCUS_NAME_LENGTH],weightNames[MAX_LOCI][LOCUS_NAME_LENGTH], comments[MAX_LOCI][MAX_COMMENT_LENGTH],trios_fn[500];

#define USEFILTERS 1
#endif