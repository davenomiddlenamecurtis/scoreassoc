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

#include <ctype.h>
#include "dcerror.hpp"
#include "scoreassoc.hpp"

void write_scores(FILE *fs,subject **sub,int nsub,double **score,int numScores,par_info *pi)
{
	int s,sc;
	for (s = 0; s < nsub; ++s) {
		fprintf(fs, "%-20s ", sub[s]->id);
		if (pi->is_quantitative)
			fprintf(fs, "%f ", sub[s]->pheno);
		else
			fprintf(fs, "%d ", sub[s]->cc);
		for (sc = 0; sc < numScores; ++sc)
			fprintf(fs, "%8.4f ", score[sc][s]);
		fprintf(fs, "\n");
	}
}

double do_score_onetailed_ttest(FILE *fo, double **score, subject **sub, int nsub, par_info *pi, sa_par_info *spi, float cc_freq[2][MAX_LOCI], float cc_count[2][MAX_LOCI], int max_cc[2], double **weight, float **missing, int *rarer)
{
	int s, i, n[2], cc, pl, l, j,sc;
	float sigma_x[2], sigma_x2[2], mean[2], var[2], tval, SE, s2, rfreq, fscore, z, total_score;
	double p, pz;
	char SLPstr[10];
#ifndef USEMLP
	float SLP;
#endif
	for (sc = 0; sc < spi->numScores; ++sc)
	{
		for (i = 0; i < 2; ++i)
			sigma_x[i] = sigma_x2[i] = n[i] = 0;
		for (s = 0; s < nsub; ++s)
		{
			cc = sub[s]->cc;
			++n[cc];
			sigma_x[cc] += score[sc][s];
			sigma_x2[cc] += score[sc][s] * score[sc][s];
		}
		for (i = 0; i < 2; ++i)
		{
			if (spi->use_cc_freqs[i])
			{
				total_score = 0;
				n[i] = max_cc[i];
				var[i] = 0;
				for (pl = 0; pl < pi->n_loci_to_use; ++pl)
				{
					float tempvar = 0;
					float ngen[3]; // number of typed subjects who are AA,AB,BB
					l = pi->loci_to_use[pl];
					sigma_x2[i] = 0;
					sigma_x[i] = 0;
					if (rarer[pl] == 2)
						rfreq = cc_freq[i][l];
					// assume supplied frequency is of alt allele, i.e. allele 2
					// but score will be added to by rarer allele, even if this is allele 1
					else
						rfreq = 1 - cc_freq[i][l];
					// as actual counts missing, assume HWE
					ngen[0] = (1 - rfreq) * (1 - rfreq) * cc_count[i][l];
					ngen[1] = 2 * rfreq * (1 - rfreq) * cc_count[i][l];
					ngen[2] = rfreq * rfreq * cc_count[i][l];
					fscore = weight[sc][pl]; // average score for each subject for this locus
					sigma_x[i] += ngen[1] * fscore; // AB
					sigma_x2[i] += ngen[1] * fscore * fscore;
					sigma_x[i] += ngen[2] * 2 * fscore; // BB
					sigma_x2[i] += ngen[2] * 4 * fscore * fscore;
					fscore = missing[sc][l]; // as if missing scores were also independent
					sigma_x[i] += (max_cc[i] - cc_count[i][l]) * fscore;
					sigma_x2[i] += (max_cc[i] - cc_count[i][l]) * fscore * fscore;
					tempvar = (sigma_x2[i] - sigma_x[i] * sigma_x[i] / n[i]) / (n[i] - 1);
					var[i] += tempvar; // add variances due to each marker
					total_score += sigma_x[i];
				}
				mean[i] = total_score / n[i];
			}
			else
			{
				if (n[i] < 2)
				{
					mean[i] = var[i] = 0;
				}
				else
				{
					var[i] = (sigma_x2[i] - sigma_x[i] * sigma_x[i] / n[i]) / (n[i] - 1);
					mean[i] = sigma_x[i] / n[i];
				}
			}
		}
		if (n[0] + n[1] < 2)
			s2 = SE = 0;
		else
		{
			s2 = ((n[0] - 1) * var[0] + (n[1] - 1) * var[1]) / (n[0] + n[1] - 2);
			SE = sqrt(s2 * (1 / (float)n[0] + 1 / (float)n[1]));
		}
		if (SE == 0)
			tval = 0;
		else
			tval = (mean[1] - mean[0]) / SE;
		p = tstat(tval, n[0] + n[1] - 2.0) / 2; // one-tailed
		SLP = log10(2 * p) * (mean[0] >= mean[1] ? 1 : -1);
		if (spi->numScores == 1)
			strcpy(SLPstr, "SLP");
		else
			sprintf(SLPstr, "SLP%d", sc);
		if (fo != NULL)
		{
			if (spi->numScores > 1)
				fprintf(fo, "\n%s", weightNames[sc]);
			fprintf(fo, "\n             Controls  Cases     \n"
				"N            %9d %9d\n"
				"Mean score   %9.3f %9.3f\n"
				"SD           %9.3f %9.3f\n"
				"t (%d df) = %6.3f\n"
				"p = %10.8f\n"
				"%s = %8.2f (signed log10(p), positive if cases score higher than controls)\n",
				n[0], n[1], mean[0], mean[1], sqrt(var[0]), sqrt(var[1]), n[0] + n[1] - 2, tval, 2 * p, SLPstr,SLP);
		// I am writing SD because it will allow me to combine statistics later
		}
	}
	return SLP;
}


/* treats male subjects as homozygous females for X loci*/
void get_scores(double **score,double **weight,float **missing,int *rarer,subject **sub,int nsub,par_info *pi,sa_par_info *spi)
{
	int l,ll,s,a,p,sc;
	for (s=0;s<nsub;++s)
	{
		for (sc = 0; sc < spi->numScores; ++sc)
		{
			score[sc][s] = 0;
			for (l = 0; l < pi->n_loci_to_use; ++l)
			{
				ll = pi->loci_to_use[l];
				if (spi->use_probs)
				{
					if (sub[s]->prob[ll][0] + sub[s]->prob[ll][1] + sub[s]->prob[ll][2] == 0)
						score[sc][s] += missing[sc][l] * 2;
					else if (rarer[l] == 2)
						score[sc][s] += (sub[s]->prob[ll][1] + sub[s]->prob[ll][2] * 2) * weight[sc][l];
					else
						score[sc][s] += (sub[s]->prob[ll][1] + sub[s]->prob[ll][0] * 2) * weight[sc][l];
				}
				else
				{
					for (a = 0; a < 2; ++a)
						if (sub[s]->all[ll][a] == rarer[l])
							score[sc][s] += weight[sc][l];
						else if (sub[s]->all[ll][a] == 0)
							score[sc][s] += missing[sc][l];
				}
			}
		}
	}
}

// treats male subjects as homozygous females for X loci
// this function is here to allow us to exclude loci with wildly different allele frequencies, etc.
void get_freqs(subject **sub,int nsub,par_info *pi,sa_par_info *spi,float cc_freq[2][MAX_LOCI],float cc_count[2][MAX_LOCI],float cc_genocount[2][3][MAX_LOCI])
{
	int l,s,nh[2],cc,i,g;
	float ccfreq[2],gencount[2][3],vcount[2];
	for (l=0;l<pi->n_loci_to_use;++l)
	{
		nh[0]=nh[1]=vcount[0]=vcount[1]=0;
		for (i=0;i<3;++i)
			for (cc=0;cc<2;++cc)
				gencount[cc][i]=0;
		if (spi->use_probs)
			for (s = 0; s < nsub; ++s)
			{
			for (g=0;g<3;++g)
				gencount[sub[s]->cc][g]+=sub[s]->prob[pi->loci_to_use[l]][g];
			vcount[sub[s]->cc]+=sub[s]->prob[pi->loci_to_use[l]][1]+sub[s]->prob[pi->loci_to_use[l]][2]*2;
			// allele count
			nh[sub[s]->cc]+=(sub[s]->prob[pi->loci_to_use[l]][0]+sub[s]->prob[pi->loci_to_use[l]][1]+sub[s]->prob[pi->loci_to_use[l]][2])*2;
			// allow for possibility that unknowns could be coded as 0 0 0
			}
		else
		for (s=0;s<nsub;++s)
			if (sub[s]->all[pi->loci_to_use[l]][0]!=0)
			{
			i=(sub[s]->all[pi->loci_to_use[l]][0]==2)+(sub[s]->all[pi->loci_to_use[l]][1]==2);
			++gencount[sub[s]->cc][i];
			vcount[sub[s]->cc]+=(sub[s]->all[pi->loci_to_use[l]][0]==2)+(sub[s]->all[pi->loci_to_use[l]][1]==2);
			nh[sub[s]->cc]+=2;
			}
		for (cc=0;cc<2;++cc)
			if (!spi->use_cc_freqs[cc])
			{
				if (nh[cc]==0)
					ccfreq[cc]=0;
				else
					ccfreq[cc]=vcount[cc]/nh[cc];
				cc_freq[cc][pi->loci_to_use[l]]=ccfreq[cc];
				// this line is here because I will use this for filtering bad loci
				cc_count[cc][pi->loci_to_use[l]]=nh[cc]/2;
				for (g=0;g<3;++g)
					cc_genocount[cc][g][pi->loci_to_use[l]]=gencount[cc][g];
			}
	}
}

float get_quadratic_weight(float freq,float wfactor)
{
	float wt;
	wt=(4*wfactor-4)*freq*freq-(4*wfactor-4)*freq+wfactor;
	return wt;
}

float get_quartic_weight(float freq,float wfactor)
{
	static float a,b,kept_wfactor=-1;
	float wt;
	if (wfactor!=kept_wfactor)
	{
		a=(1/pow(0.4,2)-wfactor/pow(0.5,2))/(pow(0.4,2)-pow(0.5,2));
		b=(1/pow(0.4,4)-wfactor/pow(0.5,4))/(1/pow(0.4,2)-1/pow(0.5,2));
		kept_wfactor=wfactor;
	}
	wt=a*pow(freq-0.5,4)+b*pow(freq-0.5,2);
	return wt;
}

float get_zero_based_quadratic_weight(float freq,float wfactor)
/* wfactor is gradient when freq==0 or freq==1 - a=wfactor */
{
	float wt;
	wt=wfactor*(freq-0.5)*(freq-0.5);
	return wt;
}

void set_weights(FILE *f,double **weight,float **missing_score,int *rarer,subject **sub,int nsub,par_info *pi,sa_par_info *spi,double **func_weight,float cc_freq[2][MAX_LOCI],float cc_count[2][MAX_LOCI],int max_cc[2],char names[MAX_LOCI][LOCUS_NAME_LENGTH],char comments[MAX_LOCI][MAX_COMMENT_LENGTH])
{
	int l,ll,s,nh[2],cc,i,g,sc;
	float freq,ccfreq[2],vcount[2],gencount[2][3],qtot[3],fixedfreq;

	for (l=0;l<pi->n_loci_to_use;++l)
	{		
		nh[0]=nh[1]=vcount[0]=vcount[1]=0;
		for (i = 0; i < 3; ++i)
		{
			qtot[i] = 0;
			for (cc = 0; cc < 2; ++cc)
				gencount[cc][i] = 0;
		}
		ll=pi->loci_to_use[l];
		for (s = 0; s < nsub; ++s)
		{
			if (spi->use_probs)
			{
				for (g = 0; g < 3; ++g)
				{
					gencount[sub[s]->cc][g] += sub[s]->prob[ll][g];
					if (pi->is_quantitative)
						qtot[g] += sub[s]->prob[ll][g] * sub[s]->pheno;
				}
				vcount[sub[s]->cc] += sub[s]->prob[ll][1]+sub[s]->prob[ll][2]*2;
				nh[sub[s]->cc]+=(sub[s]->prob[ll][0]+sub[s]->prob[ll][1]+sub[s]->prob[ll][2])*2;
			}
			else
			{
				if (sub[s]->all[ll][0] != 0)
				{
					i = (sub[s]->all[ll][0] == 2) + (sub[s]->all[ll][1] == 2);
					++gencount[sub[s]->cc][i];
					if (pi->is_quantitative)
						qtot[i] += sub[s]->pheno;
					vcount[sub[s]->cc] += (sub[s]->all[ll][0] == 2) + (sub[s]->all[ll][1] == 2);
					nh[sub[s]->cc] += 2;
				}
			}
		}
		for (cc=0;cc<2;++cc)
			if (spi->use_cc_freqs[cc])
			{
				ccfreq[cc]=cc_freq[cc][ll];
				nh[cc]=2*cc_count[cc][ll];
			}
			else
			{
				if (nh[cc]==0)
					ccfreq[cc]=0;
				else
					ccfreq[cc]=vcount[cc]/nh[cc];
				cc_freq[cc][ll]=ccfreq[cc];
				// this line is here because I will use this for filtering bad loci
			}
		if (nh[0]+nh[1]==0)
			freq=0.0;
		else
			freq=(ccfreq[0]*nh[0]+ccfreq[1]*nh[1])/(nh[0]+nh[1]);
		if (freq>0.5)
		{
			freq=1.0-freq;
			rarer[l]=1;
		}
		else
			rarer[l]=2;
		if (freq==0.0) /* monomorphic */
			for (sc = 0; sc < spi->numScores; ++sc)
				weight[sc][l]=0;
		else
		{
			if (spi->max_MAF >= 0.5)
				fixedfreq = freq;
			else if (freq > spi->max_MAF)
				fixedfreq = 0.5;
			else
				fixedfreq = freq * 0.5/spi->max_MAF;
			for (sc = 0; sc < spi->numScores; ++sc)
				weight[sc][l] = get_quadratic_weight(fixedfreq, spi->wfactor);
		}
		if (spi->use_func_weights)
			for (sc = 0; sc < spi->numScores; ++sc) {
				weight[sc][l] *= func_weight[sc][pi->loci_to_use[l]];
				if (spi->missingZero)
					missing_score[sc][l] = 0;
				else
					missing_score[sc][l] = weight[sc][l] * freq; // I think
			}
		if (f!=0)
	{
		if (!spi->use_locus_names)
			sprintf(names[ll],"LOC%05d",ll+1);
		fprintf(f,"%-" LOCUS_NAME_LENGTH_STR "s ",names[ll]);
		if (pi->is_quantitative)
			fprintf(f,
				spi->use_probs ?
				"%6.2f : %6.2f : %6.2f  %8.6f  %8.3f : %8.3f : %8.3f  %d     " :
				"%6.0f : %6.0f : %6.0f  %8.6f  %8.3f : %8.3f : %8.3f  %d     ",
				gencount[0][0], gencount[0][1], gencount[0][2],
				ccfreq[0],
				gencount[0][0] ? qtot[0] / gencount[0][0]:0, 
				gencount[0][1] ? qtot[1] / gencount[0][1]:0,
				gencount[0][2] ? qtot[2] / gencount[0][2]:0,
				rarer[l]);
		else
			fprintf(f,
			spi->use_probs ?
			"%6.2f : %6.2f : %6.2f  %8.6f  %6.2f : %6.2f : %6.2f  %8.6f  %8.6f  %d     " :
			"%6.0f : %6.0f : %6.0f  %8.6f  %6.0f : %6.0f : %6.0f  %8.6f  %8.6f  %d     ",
			gencount[0][0],gencount[0][1],gencount[0][2],
			ccfreq[0],
			gencount[1][0],gencount[1][1],gencount[1][2],
			ccfreq[1],
			freq,rarer[l]);
		for (sc = 0; sc < spi->numScores; ++sc)
			fprintf(f, "%7.2f  ", weight[sc][l]);
		fprintf(f,"%s\n", spi->use_comments ? comments[ll] : "");

	}
	}
}

int read_sa_par(FILE *fp, par_info *pi, sa_par_info *spi, float *func_weight, float cc_freq[2][MAX_LOCI], float cc_count[2][MAX_LOCI], int max_cc[2], char names[MAX_LOCI][20], char comments[MAX_LOCI][MAX_COMMENT_LENGTH])
{
char *ptr;
int l,cc,c;
if (!fgets(long_line,LONG_LINE_LENGTH,fp) || sscanf(long_line,"%f",&spi->wfactor)<1)
  {  
  error("Could not read weight factor\n",long_line);
  return 0;
  }
spi->use_func_weights=spi->use_cc_freqs[0]=spi->use_cc_freqs[1]=spi->use_locus_names=spi->use_comments=spi->do_recessive_test=0;
spi->weight_threshold=0;
spi->LD_threshold=1;
spi->use_haplotypes=0;
spi->show_hap_locus_names = 0;
if (fgets(long_line,LONG_LINE_LENGTH,fp)) // extra info is optional
{
	sscanf(long_line,"%d %d %d %d %d %d %f %f %d",
		&spi->use_func_weights,&spi->use_cc_freqs[0],&spi->use_cc_freqs[1],&spi->use_locus_names,&spi->use_comments,&spi->do_recessive_test,&spi->weight_threshold,&spi->LD_threshold,&spi->use_haplotypes);
}
if (spi->use_func_weights)
{
	if (!fgets(long_line,LONG_LINE_LENGTH,fp))
	{  
		error("Could not read functional weights\n","");
		return 0;
	}
	for (l=0,ptr=long_line;l<pi->nloci;++l)
	{
		if (sscanf(ptr,"%f",&func_weight[l])<1)
		{
			error("Could not read all functional weights\n",long_line);
			return 0;
		}
		while(isspace(*ptr++)) ;
		while(!isspace(*ptr++)) ;

	}
}
else
	for (l=0;l<pi->nloci;++l)
		func_weight[l]=1.0;
for (cc=0;cc<2;++cc) if (spi->use_cc_freqs[cc])
{
	if (!fgets(long_line,LONG_LINE_LENGTH,fp))
	{  
		error("Could not read supplied frequencies\n","");
		return 0;
	}
	for (l=0,ptr=long_line;l<pi->nloci;++l)
	{
		if (sscanf(ptr,"%f",&cc_freq[cc][l])<1)
		{
			error("Could not read all supplied frequencies\n",long_line);
			return 0;
		}
		while(isspace(*ptr++)) ;
		while(!isspace(*ptr++)) ;
	}
	if (!fgets(long_line,LONG_LINE_LENGTH,fp))
	{  
		error("Could not read subject counts\n","");
		return 0;
	}
	max_cc[cc]=0;
	for (l=0,ptr=long_line;l<pi->nloci;++l)
	{
		if (sscanf(ptr,"%f",&cc_count[cc][l])<1)
		{
			error("Could not read all subject counts\n",long_line);
			return 0;
		}
		if (max_cc[cc]<cc_count[cc][l])
			max_cc[cc]=cc_count[cc][l];
		while(isspace(*ptr++)) ;
		while(!isspace(*ptr++)) ;
	}
}
if (spi->use_locus_names)
{
	if (!fgets(long_line,LONG_LINE_LENGTH,fp))
	{  
		error("Could not read locus names\n","");
		return 0;
	}
	for (l=0,ptr=long_line;l<pi->nloci;++l)
	{
		if (sscanf(ptr,"%s",names[l])<1)
		{
			error("Could not read all locus names\n",long_line);
			return 0;
		}
		while(isspace(*ptr++)) ;
		while(!isspace(*ptr++)) ;

	}
}
if (spi->use_comments)
{
	if (!fgets(long_line,LONG_LINE_LENGTH,fp))
	{  
		error("Could not read locus comments\n","");
		return 0;
	}
	for (l=0,ptr=long_line;l<pi->nloci;++l)
	{
		while(isspace(*ptr)) 
			++ptr;
		if (*ptr=='\0')
		{
			error("Could not read all locus comments\n",long_line);
			return 0;
		}
		c=0;
		while(!isspace(*ptr))
		{
			if (c == MAX_COMMENT_LENGTH - 1)
			{
				comments[l][c] = '\0';
				while (*ptr && !isspace(*ptr))
					++ptr;
				break;
			}
			else
			{
				comments[l][c++] = *ptr;
				++ptr;
			}
		}
		comments[l][c]='\0';
	}
}
return 1;
}

int read_score_assoc_par(FILE *fp,par_info *pi,float *wfactor,int *use_func_weights,int use_cc_freqs[2],float *func_weight,float cc_freq[2][MAX_LOCI],float cc_count[2][MAX_LOCI],int max_cc[2],int *use_locus_names,char names[MAX_LOCI][20],int *use_comments,char comments[MAX_LOCI][MAX_COMMENT_LENGTH],int *do_recessive_test,float *weight_threshold,float *LD_threshold,int *use_haplotypes)
{
char *ptr;
int l,cc,c;
if (!fgets(long_line,LONG_LINE_LENGTH,fp) || sscanf(long_line,"%f",wfactor)<1)
  {  
  error("Could not read weight factor\n",long_line);
  return 0;
  }
*use_func_weights=use_cc_freqs[0]=use_cc_freqs[1]=*use_locus_names=*use_comments=*do_recessive_test=0;
*weight_threshold=0;
*LD_threshold=1;
*use_haplotypes=0;
if (fgets(long_line,LONG_LINE_LENGTH,fp)) // extra info is optional
{
	sscanf(long_line,"%d %d %d %d %d %d %f %f %d",use_func_weights,&use_cc_freqs[0],&use_cc_freqs[1],use_locus_names,use_comments,do_recessive_test,weight_threshold,LD_threshold,use_haplotypes);
}
if (*use_func_weights)
{
	if (!fgets(long_line,LONG_LINE_LENGTH,fp))
	{  
		error("Could not read starting weights\n","");
		return 0;
	}
	for (l=0,ptr=long_line;l<pi->nloci;++l)
	{
		if (sscanf(ptr,"%f",&func_weight[l])<1)
		{
			error("Could not read all starting weights\n",long_line);
			return 0;
		}
		while(isspace(*ptr++)) ;
		while(!isspace(*ptr++)) ;

	}
}
else
	for (l=0;l<pi->nloci;++l)
		func_weight[l]=1.0;
for (cc=0;cc<2;++cc) if (use_cc_freqs[cc])
{
	if (!fgets(long_line,LONG_LINE_LENGTH,fp))
	{  
		error("Could not read supplied frequencies\n","");
		return 0;
	}
	for (l=0,ptr=long_line;l<pi->nloci;++l)
	{
		if (sscanf(ptr,"%f",&cc_freq[cc][l])<1)
		{
			error("Could not read all supplied frequencies\n",long_line);
			return 0;
		}
		while(isspace(*ptr++)) ;
		while(!isspace(*ptr++)) ;
	}
	if (!fgets(long_line,LONG_LINE_LENGTH,fp))
	{  
		error("Could not read subject counts\n","");
		return 0;
	}
	max_cc[cc]=0;
	for (l=0,ptr=long_line;l<pi->nloci;++l)
	{
		if (sscanf(ptr,"%f",&cc_count[cc][l])<1)
		{
			error("Could not read all subject counts\n",long_line);
			return 0;
		}
		if (max_cc[cc]<cc_count[cc][l])
			max_cc[cc]=cc_count[cc][l];
		while(isspace(*ptr++)) ;
		while(!isspace(*ptr++)) ;
	}
}
if (*use_locus_names)
{
	if (!fgets(long_line,LONG_LINE_LENGTH,fp))
	{  
		error("Could not read locus names\n","");
		return 0;
	}
	for (l=0,ptr=long_line;l<pi->nloci;++l)
	{
		if (sscanf(ptr,"%s",names[l])<1)
		{
			error("Could not read all locus names\n",long_line);
			return 0;
		}
		while(isspace(*ptr++)) ;
		while(!isspace(*ptr++)) ;

	}
}
if (*use_comments)
{
	if (!fgets(long_line,LONG_LINE_LENGTH,fp))
	{  
		error("Could not read locus comments\n","");
		return 0;
	}
	for (l=0,ptr=long_line;l<pi->nloci;++l)
	{
		while(isspace(*ptr)) 
			++ptr;
		if (*ptr=='\0')
		{
			error("Could not read all locus comments\n",long_line);
			return 0;
		}
		c=0;
		while(!isspace(*ptr))
		{
			if (c==MAX_COMMENT_LENGTH-1)
				comments[l][c]='\0';
			else
				comments[l][c++]=*ptr;
			++ptr;
		}
		comments[l][c]='\0';
	}
}
return 1;
}

int read_all_data(par_info *pi,sa_par_info *spi,subject **sub,int *nsubptr,char names[MAX_LOCI][LOCUS_NAME_LENGTH],char comments[MAX_LOCI][MAX_COMMENT_LENGTH],double **func_weight,double **score,int numScores)
{
	char aline[1000],pos[100],effect[100],*ptr;
	int func_pos,l,f,use,s;
	double wt;
	std::map<std::string,double> weightMap;
	std::map<std::string,std::string> effectMap;
	if (spi->df[WEIGHTFILE].fp)
	{
		while (fgets(aline, 999, spi->df[WEIGHTFILE].fp))
		{
			sscanf(aline,"%s %lf",effect,&wt);
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
		read_ps_datafile(pi,spi,sub,nsubptr,names,comments, func_weight[0],weightMap,effectMap);
	else if (spi->df[GCDATAFILE].fp && !spi->useTransposedFile)
		read_all_subjects(spi->df[GCDATAFILE].fp, sub, nsubptr, pi);
	else if (spi->df[GCDATAFILE].fp && spi->useTransposedFile)
		read_all_subjects_transposed(spi->df[GCDATAFILE].fp, sub, nsubptr, pi);
	else if (spi->df[INPUTSCOREFILE].fp)
		read_all_subject_scores(spi->df[INPUTSCOREFILE].fp, sub, nsubptr, score,spi->numScores);
	if (spi->df[IDPHENOTYPEFILE].fp)
		read_phenotypes(spi->df[IDPHENOTYPEFILE].fp, sub, nsubptr,score, spi->numScores, pi->is_quantitative);
	// overwrite previously defined phenotypes
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
	if (spi->numLocusWeightFiles == 0 && spi->numLocusRecWeightFiles == 0 && spi->df[WEIGHTFILE].fp == 0)
		for (l = 0; l < pi->nloci; ++l)
			func_weight[0][l] = 1;
	if (spi->numLocusWeightFiles>0)
	{
		for (s = 0; s < spi->numLocusWeightFiles; ++s)
		{
			for (l = 0; l < pi->nloci; ++l)
				if (fscanf(spi->locusWeightFile[s].fp, "%lf ", &func_weight[s][l]) != 1)
				{
					dcerror(1, "Not enough values in locusweightfile %s\n", spi->locusWeightFile[s].fn); exit(1);
				}
		}
	}
	if (spi->numLocusRecWeightFiles > 0)
	{
		for (s = 0; s < spi->numLocusRecWeightFiles; ++s)
		{
			for (l = 0; l < pi->nloci; ++l)
				if (fscanf(spi->locusRecWeightFile[s].fp, "%lf ", &func_weight[spi->numAddScores + s][l]) != 1)
				{
					dcerror(1, "Not enough values in locusrecweightfile %s\n", spi->locusWeightFile[s].fn); exit(1);
				}
		}
	}
	if (spi->df[LOCUSWEIGHTNAMEFILE].fp)
	{
		for (s = 0; s < spi->numAddScores; ++s)
			if (fscanf(spi->df[LOCUSWEIGHTNAMEFILE].fp, "%s", weightNames[s]) != 1)
			{
				dcerror(1, "Not enough names for weights in %s\n", spi->df[LOCUSWEIGHTNAMEFILE].fn); exit(1);
			}
	}
	else if (spi->numAddScores == 1)
		strcpy(weightNames[0], "score");
	else
		for (s = 0; s < spi->numAddScores; ++s)
			sprintf(weightNames[s], "score%d", s);
	if (spi->df[LOCUSRECWEIGHTNAMEFILE].fp)
	{
		for (s = 0; s < spi->numRecScores; ++s)
			if (fscanf(spi->df[LOCUSRECWEIGHTNAMEFILE].fp, "%s", weightNames[spi->numAddScores+s]) != 1)
			{
				dcerror(1, "Not enough names for recessive weights in %s\n", spi->df[LOCUSRECWEIGHTNAMEFILE].fn); exit(1);
			}
	}
	else if (spi->numRecScores == 1)
		strcpy(weightNames[spi->numAddScores], "recscore");
	else
		for (s = 0; s < spi->numAddScores; ++s)
			sprintf(weightNames[spi->numAddScores+s], "recscore%d", s);
	if (spi->df[LOCUSNAMEFILE].fp)
		{
		for (l = 0; l < pi->nloci; ++l)
#if 1
		{
			int ch,c;
			while ((ch = fgetc(spi->df[LOCUSNAMEFILE].fp)) != EOF && isspace(ch))
				;
			if (ch == EOF)
			{
				error("Could not read all locus names/comments from ", spi->df[LOCUSNAMEFILE].fn);
				return 0;
			}
			c = 0;
			while (ch!=EOF && !isspace(ch))
			{
				if (c == MAX_COMMENT_LENGTH - 1)
				{
					comments[l][c] = '\0';
					while ((ch = fgetc(spi->df[LOCUSNAMEFILE].fp)) != EOF && !isspace(ch))
						;
					break;
				}
				else
				{
					comments[l][c++] = ch;
					ch = fgetc(spi->df[LOCUSNAMEFILE].fp);
				}
			}
			comments[l][c] = '\0';
		}
#else
			if (fscanf(spi->df[LOCUSNAMEFILE].fp,"%s ",&comments[l])!=1)
			{
				dcerror(1, "Not enough values in locusnamefile %s\n", spi->df[LOCUSNAMEFILE].fn); exit(1);
			}
#endif
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
		strncpy(names[l],comments[l],LOCUS_NAME_LENGTH-1);
		names[l][LOCUS_NAME_LENGTH-1]='\0';
	}
	if (spi->numVarFiles)
	{
		std::map<std::string, int> subIDs;
		if (spi->numVarFiles>0)
		{
			for (s = 0; s<*nsubptr; ++s)
				subIDs[sub[s]->id] = s;
		}
		if (!readVarFiles(subIDs, *nsubptr, spi))
		exit(1);
	}
	return 1;
}

int read_freqs_datafile(par_info *pi, sa_par_info *spi, int cc, float cc_freq[2][MAX_LOCI], float cc_count[2][MAX_LOCI], int max_cc[2])
{
	int f, l;
	FILE *fp;
	char *fn, *ptr;
	f = cc ? CASEFREQFILE : CONTFREQFILE;
	fn = spi->df[f].fn;
	fp = spi->df[f].fp;
	spi->use_cc_freqs[cc] = 1;
	if (!fgets(long_line, LONG_LINE_LENGTH, fp))
	{
		error("Could not read supplied frequencies in ", fn);
		return 0;
	}
	for (l = 0, ptr = long_line; l<pi->nloci; ++l)
	{
		if (sscanf(ptr, "%f", &cc_freq[cc][l])<1)
		{
			error("Could not read all supplied frequencies:\n", long_line);
			return 0;
		}
		while (isspace(*ptr++));
		while (!isspace(*ptr++));
	}
	if (!fgets(long_line, LONG_LINE_LENGTH, fp))
	{
		error("Could not read subject counts in ", fn);
		return 0;
	}
	max_cc[cc] = 0;
	for (l = 0, ptr = long_line; l<pi->nloci; ++l)
	{
		if (sscanf(ptr, "%f", &cc_count[cc][l])<1)
		{
			error("Could not read all subject counts:\n", long_line);
			return 0;
		}
		if (max_cc[cc]<cc_count[cc][l])
			max_cc[cc] = cc_count[cc][l];
		while (isspace(*ptr++));
		while (!isspace(*ptr++));
	}
	return 1;
}

int read_ps_datafile(par_info *pi, sa_par_info *spi, subject **sub, int *nsubptr, char names[MAX_LOCI][LOCUS_NAME_LENGTH], char comments[MAX_LOCI][MAX_COMMENT_LENGTH], double *func_weight,
	std::map<std::string, double> weightMap, std::map<std::string, std::string> effectMap)
{
	int nsub, first, s, a, l, f;
	char dline[1000], aline[1000], pos[100], rsname[100], ref[100], alls[100], all[2][100], sname[100], ccstr[100], *ptr;
	float weight;
	std::map<std::string, int> subIDs;
	first = 1;
	nsub = 0;
	pi->nloci = 0;
	while (fgets(dline, 999, spi->df[PSDATAFILE].fp))
	{
		if (sscanf(dline, "%s", pos) != 1)
			continue;
		if (!strncmp("chr", pos, 3) && strchr(pos, ':')) //looks like new locus, not subject ID, will break if subject does have ID like this
		{
			++pi->nloci;
			l = pi->nloci - 1;
			pi->n_alls[l] = 2;
			if (sscanf(dline, "%s %s %s", pos, rsname, alls) != 3)
			{
				dcerror(1, "Could not read this locus line:\n%s\n", dline);
				exit(1);
			}
			if (sscanf(alls, "%[^/]", ref) != 1)
			{
				dcerror(1, "Could not read reference allele in this locus line:\n%s\n", dline);
				exit(1);
			}
			if (spi->df[WEIGHTFILE].fp && spi->df[ANNOTFILE].fp)
			{
				std::map<std::string, std::string>::const_iterator effIter = effectMap.find(pos);
				if (effIter == effectMap.end())
				{
					dcerror(1, "%s not found in annotation file %s\n", pos, spi->df[ANNOTFILE].fn);
					exit(1);
				}
				sprintf(comments[l], "%s_%s_%s", pos, effIter->second.c_str(), alls);
				std::map<std::string, double>::const_iterator weightIter = weightMap.find(effIter->second);
				if (weightIter == weightMap.end())
				{
					dcerror(1, "weight for effect %s not found in weight file %s\n", effIter->second.c_str(), spi->df[WEIGHTFILE].fn);
					exit(1);
				}
				else
					func_weight[l] = weightIter->second;
			}
			else
			{
				sprintf(comments[l], "%s_%s", pos, alls);
				func_weight[l] = 1;
			}
		}
		else
		{
			if (sscanf(dline, "%s %*d %s %[^/]/%s", sname, ccstr, all[0], all[1]) != 4)
			{
				dcerror(1, "Could not read this genotype line:\n%s\n", dline); exit(1);
			}
			std::map<std::string, int>::const_iterator iter = subIDs.find(sname);
			if (iter == subIDs.end()) // found a subject we haven't seen before
			{
				subIDs[sname] = nsub;
				s = nsub;
				strcpy(sub[s]->id, sname);
				if (!strcmp(ccstr, "CONTROL"))
					sub[s]->cc = 0;
				else if (!strcmp(ccstr, "CASE"))
					sub[s]->cc = 1;
				else
				{
					dcerror(1, "Unrecognised phenotype in this line:\n%s\n", dline); exit(1);
				}
				++nsub;
			}
			else
				s = iter->second;
			if (all[0][0] == '.')
				continue;
			else for (a = 0; a<2; ++a)
				if (!strcmp(all[a], ref))
					sub[s]->all[l][a] = 1;
				else
					sub[s]->all[l][a] = 2;
		}
	}
	*nsubptr = nsub;
	return 1;
}

int read_all_subject_scores(FILE* fi, subject** s, int* nsub, double **score,int numScores)
{
	char id[MAX_ID_LENGTH + 1],rest[2000];
	float pheno;
	int sc;
	long here;
	*nsub = 0;

	while (fgets(long_line, LONG_LINE_LENGTH, fi) && sscanf(long_line, "%s %f %[^\n]", id, &pheno, rest) == 3)
		if (*nsub == MAX_SUB)
		{
			error("Number of subjects exceeds MAX_SUB", ""); return 0;
		}
		else
		{
			here = ftell(fi);
			strcpy(s[*nsub]->id, id);
			s[*nsub]->pheno = pheno;
			for (sc=0;sc<numScores;++sc)
			{
				strcpy(long_line,rest);
				if (sscanf(long_line, "%f %[^\n]", &score[sc][*nsub], rest) < 1)
				{
					fseek(fi, here, SEEK_SET);
					fgets(long_line, LONG_LINE_LENGTH, fi);
					error("Not enough scores in this line:", long_line);
					return 0;
				}
			}
			++(*nsub);
		}
	return 1;
}

#define MISSINGPHENOTYPECODE -999 // this should never be used
int read_phenotypes(FILE* fi, subject** s, int *nsub, double **score, int numScores, int isquantitative)
{
	std::map<std::string, float> subPhenos;
	char id[MAX_ID_LENGTH + 1];
	float pheno;
	int ss,sss,sc;
	subject* tempSub; // now subject is a class can't have two pointers pointing to one subject
	while (pheno=MISSINGPHENOTYPECODE, fgets(long_line, LONG_LINE_LENGTH, fi) && sscanf(long_line, "%s %f", id, &pheno)>=1)
		subPhenos.insert(std::pair<std::string, float>(id, pheno));
	for (ss = 0; ss < *nsub; ++ss)
	{
		std::map<std::string, float>::const_iterator it = subPhenos.find(s[ss]->id);
		if (it == subPhenos.end())
		{
			// dcerror(2, "No phenotype found for subject %s\n", s[ss]->id);
			// return 0;
			s[ss]->pheno = MISSINGPHENOTYPECODE; // only analyse subjects in IDandPhenotype file, ignore others
		}
		else if (isquantitative)
			s[ss]->pheno= it->second;
		else
			s[ss]->cc = it->second;
	}
	for (ss = *nsub - 1; ss >= 0; --ss)
	{
		if (s[ss]->pheno == MISSINGPHENOTYPECODE || (isquantitative == 0 && s[ss]->cc != 0 && s[ss]->cc != 1))
		{
			tempSub = s[ss];
			for (sss = ss; sss < *nsub - 1; ++sss)
			{
				s[sss] = s[sss + 1];
				for (sc=0;sc<numScores;++sc)
				score[sc][sss] = score[sc][sss + 1];
			}
			s[sss] = tempSub; // pop it back at the end
			--*nsub;
		}
	}
	return 1;
}

int read_all_subjects_transposed(FILE* fi, subject** sub, int* nsub, par_info* pi)
{
	char id[MAX_ID_LENGTH + 1];
	int found_error, ch, finished, s, i;
	float pheno;
	found_error = 0;
	s = 0;
	finished = 0;
	while (fscanf(fi, "%s", id) == 1)
	{
		if (s >= MAX_SUB)
		{
			error("Number of subjects exceeds MAX_SUB", ""); return 0;
		}
		else
			strcpy(sub[s++]->id, id);
		do {
			ch = fgetc(fi);
			if (ch == '\n' || ch == EOF)
			{
				finished = 1;
				break;
			}
		} while (isspace(ch));
		if (finished)
			break;
		ungetc(ch, fi);
	}
	*nsub = s;
	finished = 0;
	s = 0;
	while (fscanf(fi, "%f", &pheno) == 1)
	{
		if (pi->is_quantitative)
		{
			sub[s]->pheno = pheno;
			sub[s]->cc = 0;
		}
		else
		{
			if (pheno != 0 && pheno != 1)
			{
				dcerror(1, "Bad case control status %f for subject %s\n", pheno, sub[s]->id);
				return 0;
			}
			else
				sub[s]->cc = pheno;
		}
		++s;
		do {
			ch = fgetc(fi);
			if (ch == '\n' || ch == EOF)
			{
				finished = 1;
				break;
			}
		} while (isspace(ch));
		if (finished)
			break;
		ungetc(ch, fi);
	}
	if (s != *nsub)
	{
		dcerror(1, "There are %d phenotype codes in second line of data file but %d subjects\n", s, *nsub);
		return 0;
	}
	for (i = 0; i < pi->nloci; ++i)
	{
		finished = 0;
		s = 0;
		while (fscanf(fi, "%d %d", &sub[s]->all[i][0], &sub[s]->all[i][1]) == 2)
		{
			if ((sub[s]->all[i][0] == 0) != (sub[s]->all[i][1] == 0))
			{
				printf("In subject %s there is only one zero allele for locus %d\n",
					sub[s]->id, i + 1);
				found_error = 1;
			}
			if (sub[s]->all[i][0]<0 || sub[s]->all[i][1]<0 || sub[s]->all[i][0]>pi->n_alls[i] || sub[s]->all[i][1]>pi->n_alls[i])
			{
				printf("In subject %s bad allele number for locus %d\n",
					sub[s]->id, i + 1);
				found_error = 1;
			}
			++s;
			do {
				ch = fgetc(fi);
				if (ch == '\n' || ch == EOF)
				{
					finished = 1;
					break;
				}
			} while (isspace(ch));
			if (finished)
				break;
			ungetc(ch, fi);
		}
		if (s != *nsub)
		{
			dcerror(1, "There are %d allele pairs but %d subjects for locus %d\n", s, *nsub, i);
			return 0;
		}
	}
	if (i != pi->nloci)
	{
		dcerror(1, "Could only read in data for %d out of %d loci\n", i, pi->nloci);
		return 0;
	}
	if (found_error)
		dcerror(1,"Error[s] found in data file\n");
	return !found_error;
}

