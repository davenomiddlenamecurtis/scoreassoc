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

int getValidPairs(int** usePair, int* rarer, subject * *sub, int nsub, par_info * pi, sa_par_info * spi)
{
	int l,ll1,l2,ll2,a,aa,n1,n2,i,j,s,t;
	float* called, * count, **pairCount,**pairCalled,expected;
	assert(count = (float*)calloc(pi->n_loci_to_use, sizeof(float)));
	assert(called = (float*)calloc(pi->n_loci_to_use, sizeof(float)));
	assert(pairCount = (float**)calloc(pi->n_loci_to_use, sizeof(float*)));
	assert(pairCalled = (float**)calloc(pi->n_loci_to_use, sizeof(float*)));
	for (i = 0; i < pi->n_loci_to_use; ++i)
	{
		assert(pairCount[i] = (float*)calloc(pi->n_loci_to_use, sizeof(float)));
		assert(pairCalled[i] = (float*)calloc(pi->n_loci_to_use, sizeof(float)));
	}
	for (s = 0; s < nsub; ++s)
	{
		for (l = 0; l < pi->n_loci_to_use; ++l)
		{
			ll1 = pi->loci_to_use[l];
			n1 = 0;
			for (a = 0; a < 2; ++a)
			{
				aa = sub[s]->all[ll1][a];
				if (aa != 0)
					++called[l];
				if (aa == rarer[l])
					++n1;
			}
			count[l] += n1;
			for (l2 = l + 1; l2 < pi->n_loci_to_use; ++l2)
			{
				ll2 = pi->loci_to_use[l2];
				n2 = 0;
				for (a = 0; a < 2; ++a)
				{
					aa = sub[s]->all[ll2][a];
					if (aa != 0)
						++pairCalled[l][l2];
					if (aa == rarer[l2])
						++n2;
				}
				t = n1 + n2;
				if (t > 1)
				{
					if (t == 4) // 2 and 2
						pairCount[l][l2] += 2;
					else if (t == 3) // 1 and 2
						pairCount[l][l2] += 1;
					else if (n1 == 1) // both 1
						pairCount[l][l2] += 0.5; // if 2 and 0 then add nothing
				}
			}
		}
	}
	for (l = 0; l < pi->n_loci_to_use; ++l)
	{
		usePair[l][l] = 1; // might not do this if there seems to be an excess of homozygotes
		for (l2 = 0; l2 < pi->n_loci_to_use; ++l2)
		{
			expected = count[l] / called[l] * count[l2] / called[l2] * pairCalled[l][l2];
			usePair[l][l2] = pairCount[l][l2] < expected * spi->LDThreshold2022;
		}
	}
	
	for (i = 0; i < pi->n_loci_to_use; ++i)
	{
		free(pairCount[i]);
		free(pairCalled[i]);
	}
	free(pairCount);
	free(pairCalled);
	free(count);
	free(called);
	return 1;
}

/* treats male subjects as homozygous females for X loci */
void get_rec_scores(double **score,recPair **rP,double **weight,float **missing,int *rarer,subject **sub,int nsub,par_info *pi,sa_par_info *spi)
{
	int l,ll1,l2,ll2,s,a,p,sc,s1,s2,n1,n2,t;
	float recAlls;
	int **usePair;
	float pairScore;
	if (spi->use_probs)
	{
		dcerror(1,"Recessive analyses with probabilistic genotypes are not supported yet");
		exit(1);
	}
	assert(usePair = (int **)calloc(pi->n_loci_to_use, sizeof(int*)));
	for (l = 0; l < pi->n_loci_to_use; ++l)
		assert(usePair[l] = (int*)calloc(pi->n_loci_to_use, sizeof(int)));
	getValidPairs(usePair,rarer, sub, nsub, pi, spi);
	for (s = 0; s < nsub; ++s)
	{
		for (sc = 0; sc < spi->numRecScores; ++sc)
		{
			score[spi->numAddScores + sc][s]=0;
			rP[sc][s].l[0] = rP[sc][s].l[1] = -1;
		}
		for (l = 0; l < pi->n_loci_to_use; ++l)
		{
			ll1 = pi->loci_to_use[l];
			if (sub[s]->all[ll1][0] == 0) // unknown
				n1 = -1;
			else
			{
				n1 = 0;
				for (a = 0; a < 2; ++a)
				{
					n1 += sub[s]->all[ll1][a] == rarer[l];
				}
			}
			if (n1 == 2)
			{
				for (sc = 0; sc < spi->numRecScores; ++sc)
				{
					pairScore = weight[sc][l] * 2;
					if (pairScore > score[spi->numAddScores + sc][s])
					{
						score[spi->numAddScores + sc][s] = pairScore;
						rP[sc][s].l[0] = rP[sc][s].l[1] = l;
					}
				}
			}
			for (l2 = l+1; l2 < pi->n_loci_to_use; ++l2)
			{
				if (usePair[l][l2] == 0) 
					continue;
				ll2 = pi->loci_to_use[l2];
				if (sub[s]->all[ll2][0] == 0)
					n2 = -1;
				else
				{
					n2 = 0;
					for (a = 0; a < 2; ++a)
					{
						n2 += sub[s]->all[ll2][a] == rarer[l2];
					}
				}
				if (n1 == -1 || n2 == -1)
				{
					; // handle missing
				}
				else 
				{
					t = n1 + n2;
					recAlls = 0;
					if (t > 1)
					{
						if (t == 4) // 2 and 2
							recAlls = 2; // as things stand this will not be used, because recessive score will be higher
						else if (t == 3) // 1 and 2
							recAlls = 1;
						else if (n1 == 1) // both 1
							recAlls = 0.5; // if 2 and 0 then add nothing
					}
					if (recAlls != 0)
					{
						for (sc = 0; sc < spi->numRecScores; ++sc)
						{
							pairScore = (weight[sc][l] + weight[sc][l2]) * recAlls;
							if (pairScore > score[spi->numAddScores + sc][s])
							{
								score[spi->numAddScores + sc][s] = pairScore;
								rP[sc][s].l[0] = l;
								rP[sc][s].l[1] = l2;
							}
						}
					}
					
				}
			}
		}
	}
	for (l = 0; l < pi->n_loci_to_use; ++l)
		free(usePair[l]);
	free(usePair);
}

void write_rec_scores(FILE* fs, subject** sub, int nsub, double** score, recPair** rP, int numRecScores, double** weight, char names[MAX_LOCI][LOCUS_NAME_LENGTH], par_info* pi, sa_par_info* spi)
{
	int s, sc,i,l,ll;
	for (s = 0; s < nsub; ++s) {
		fprintf(fs, "%-20s ", sub[s]->id);
		if (pi->is_quantitative)
			fprintf(fs, "%f ", sub[s]->pheno);
		else
			fprintf(fs, "%d ", sub[s]->cc);
		for (sc = 0; sc < numRecScores; ++sc)
		{
			fprintf(fs, "%8.4f ", score[spi->numAddScores + sc][s]); // add locus names here
			for (i = 0; i < 2; ++i)
			{
				l = rP[sc][s].l[i];
				ll = pi->loci_to_use[l];
				fprintf(fs, "%-" LOCUS_NAME_LENGTH_STR "s %8.4f ", names[ll], weight[spi->numAddScores + sc][l]);
			}
		}
		fprintf(fs, "\n");
	}
}
