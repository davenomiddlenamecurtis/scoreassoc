#include <ctype.h>
#include <string.h>
#include <string>
#include <math.h>
#include <assert.h>
#include "lrBurdenAssoc.hpp"
#include "lrModel.hpp"
#include "dcerror.hpp"

double do_onetailed_LRT(FILE *fo, lrModel *m)
{
	double L0, L1, p, SLP;
	int b;
	L0 = getlnLikeForH(0, m);
	if (fo)
	{
		fprintf(fo,"\nL0 = %.2f, beta = %.5f",L0,m->beta[m->nCol]);
		for(b = 0; b < m->nCol;++b)
			fprintf(fo," %.5f",m->beta[b]);
		fprintf(fo,"\n");
	}
	L1 = getlnLikeForH(1, m);
	if (fo)
	{
		fprintf(fo, "\nL1 = %.2f, beta = %.5f", L1, m->beta[m->nCol]);
		for(b = 0; b < m->nCol;++b)
			fprintf(fo, " %.5f", m->beta[b]);
		fprintf(fo, "\n");
	}
	p = chistat(2 * (L1 - L0),1.0); // change for more complex models
	SLP = -log10(p);
	if (m->beta[0] < 0)
		SLP *= -1;
	if (fo)
		fprintf(fo,
			"p = %10.8f\n"
			"SLP = %8.2f (signed log10(p), positive if cases have more variants than controls)\n",p,SLP);
	return SLP;
}

double getlnLikeForH(int h,lrModel *m)
{
	int b;
	double lnL;
	assert(h == 0 || h == 1);
	for (b = 0; b < m->nCol + 1; ++b)
		m->toUse[b] = m->toFit[b] = 0;
	switch (h)
	{
	case 0:
		b= m->nCol;
		m->toUse[b] = m->toFit[b] = 1;
		break;
	case 1:
		for (b = 0; b < m->nCol + 1; ++b)
			m->toUse[b] = m->toFit[b] = 1;
		break;
	}
	lnL = m->maximiseLnL();
	return lnL;
}

void fillModel(lrModel *m, float *score, subject **sub, int nsub, par_info *pi, lrba_par_info *spi)
{
	int s;
	m->init(nsub, 1);
	for (s = 0; s < nsub; ++s)
	{
		m->X[s][0] = score[s];
		m->Y[s] = sub[s]->cc;
	}

}

void get_scores(float *score, float *weight, float *missing, int *rarer, subject **sub, int nsub, par_info *pi, lrba_par_info *spi)
{
	int l, ll, s, a, p;
	for (s = 0; s<nsub; ++s)
	{
		score[s] = 0;
		for (l = 0; l<pi->n_loci_to_use; ++l)
		{
			ll = pi->loci_to_use[l];
			if (spi->use_probs)
			{
				if (sub[s]->prob[ll][0] + sub[s]->prob[ll][1] + sub[s]->prob[ll][2] == 0)
					score[s] += missing[l] * 2;
				else if (rarer[l] == 2)
					score[s] += (sub[s]->prob[ll][1] + sub[s]->prob[ll][2] * 2)*weight[l];
				else
					score[s] += (sub[s]->prob[ll][1] + sub[s]->prob[ll][0] * 2)*weight[l];
			}
			else
			{
				for (a = 0; a < 2; ++a)
					if (sub[s]->all[ll][a] == rarer[l])
						score[s] += weight[l];
					else if (sub[s]->all[ll][a] == 0)
						score[s] += missing[l];
			}
		}
	}
}

// treats male subjects as homozygous females for X loci
// this function is here to allow us to exclude loci with wildly different allele frequencies, etc.
void get_freqs(subject **sub,int nsub,par_info *pi,lrba_par_info *spi, float cc_freq[2][MAX_LOCI], float cc_count[2][MAX_LOCI],float cc_genocount[2][3][MAX_LOCI])
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
			{
				if (nh[cc]==0)
					ccfreq[cc]=0;
				else
					ccfreq[cc]=vcount[cc]/nh[cc];
				cc_freq[cc][pi->loci_to_use[l]]=ccfreq[cc];
				cc_count[cc][pi->loci_to_use[l]]=nh[cc]/2;
				for (g=0;g<3;++g)
					cc_genocount[cc][g][pi->loci_to_use[l]]=gencount[cc][g];
			}
	}
}

void set_weights(FILE *f,float *weight,float *missing_score,int *rarer,subject **sub,int nsub,par_info *pi,lrba_par_info *spi,float *func_weight,float cc_freq[2][MAX_LOCI],float cc_count[2][MAX_LOCI],int max_cc[2],char names[MAX_LOCI][20],char comments[MAX_LOCI][MAX_COMMENT_LENGTH])
{
	int l,ll,s,nh[2],cc,i,g;
	float freq,ccfreq[2],vcount[2],gencount[2][3];

	for (l=0;l<pi->n_loci_to_use;++l)
	{		
		nh[0]=nh[1]=vcount[0]=vcount[1]=0;
		for (i=0;i<3;++i)
			for (cc=0;cc<2;++cc)
				gencount[cc][i]=0;
		ll=pi->loci_to_use[l];
		for (s = 0; s < nsub; ++s)
		{
			if (spi->use_probs)
			{
				for (g=0;g<3;++g)
					gencount[sub[s]->cc][g]+=sub[s]->prob[ll][g];
				vcount[sub[s]->cc] += sub[s]->prob[ll][1]+sub[s]->prob[ll][2]*2;
				nh[sub[s]->cc]+=(sub[s]->prob[ll][0]+sub[s]->prob[ll][1]+sub[s]->prob[ll][2])*2;
			}
			else
			{
				if (sub[s]->all[ll][0] != 0)
				{
					i = (sub[s]->all[ll][0] == 2) + (sub[s]->all[ll][1] == 2);
					++gencount[sub[s]->cc][i];
					vcount[sub[s]->cc] += (sub[s]->all[ll][0] == 2) + (sub[s]->all[ll][1] == 2);
					nh[sub[s]->cc] += 2;
				}
			}
		}
		for (cc=0;cc<2;++cc)
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
			weight[l]=0;
		else
			weight[l]=get_quadratic_weight(freq,spi->wfactor);
		if (spi->use_func_weights)
			weight[l]*=func_weight[pi->loci_to_use[l]];
		missing_score[l]=weight[l]*freq; // I think
	if (f!=0)
	{
		if (!spi->use_locus_names)
			sprintf(names[ll],"LOC%05d",ll+1);
		fprintf(f,"%-40s",names[ll]);
		fprintf(f,
			spi->use_probs ?
			"%6.2f  %6.2f  %6.2f  %8.6f  %6.2f  %6.2f  %6.2f  %8.6f  %8.6f  %d      %6.2f   %s\n" :
			"%6.0f  %6.0f  %6.0f  %8.6f  %6.0f  %6.0f  %6.0f  %8.6f  %8.6f  %d      %6.2f   %s\n",
			gencount[0][0],gencount[0][1],gencount[0][2],
			ccfreq[0],
			gencount[1][0],gencount[1][1],gencount[1][2],
			ccfreq[1],
			freq,rarer[l],weight[l],
			spi->use_comments?comments[ll]:"");

	}
	}
}

float get_quadratic_weight(float freq, float wfactor)
{
	float wt;
	wt = (4 * wfactor - 4)*freq*freq - (4 * wfactor - 4)*freq + wfactor;
	return wt;
}

