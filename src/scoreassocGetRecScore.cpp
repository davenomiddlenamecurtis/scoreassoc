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

class locusIndex {
public:
	int nLoci,*loci;
	locusIndex() { nLoci = 0; loci = 0; }
	locusIndex(int l);
	~locusIndex() { if (loci) free(loci); }
	void fill(int* l, int n);
	void unpack(int* l) { for (int ll = 0; ll < nLoci; ++ll) l[loci[ll]] = 1; }
	void remove(int i);
};

locusIndex::locusIndex(int l) 
{
	nLoci = 0;
	if (l == 0)
		loci = 0;
	else
		assert((loci = (int*)calloc(l, sizeof(int))) != 0);
}

void locusIndex::remove(int i)
{
	--nLoci;
	for (; i < nLoci; ++i)
		loci[i] = loci[i + 1];
}

void locusIndex::fill(int* l, int n)
{
	if (loci) 
		free(loci); 
	nLoci = 0; 
	loci = (int*)calloc(n, sizeof(int));
	for (int ll = 0; ll < n; ++ll)
		if (l[ll])
			loci[nLoci++] = ll;
}

class variantContribution {
public:
	int index;
	float contribution;
};

int compContribution(const void* s1, const void* s2) 
{
	variantContribution* v1, * v2;
	v1 = (variantContribution*)s1;
	v2 = (variantContribution*)s2;
	return v2->contribution - v1->contribution;
}

int getPairsToUseIndex(locusIndex * usablePairs, double** weight, int* rarer, subject** sub, int nsub, par_info* pi, sa_par_info* spi)
{
	// assume usablePairs already zeroed
	// first get all genotype counts
	// then get most important loci and discard the rest
	// then search each one for usable pairs with it
	// discard any with no usable pairs
	int(* genoCounts)[3],s,i,c,a,l,l2,ll1,ll2,lll,n1,n2,*useThisRecLocus, * useSecondRecLocus,
		anyValidPair,thisPairValid,coOccur,d;
	variantContribution* contrib;
	assert((genoCounts = (int (*)[3])calloc(pi->n_loci_to_use, sizeof(int[3]))) != 0);
	assert((contrib = (variantContribution*)calloc(pi->nloci, sizeof(variantContribution))) != 0);
	assert((useThisRecLocus = (int*)calloc(pi->nloci, sizeof(int))) != 0);
	assert((useSecondRecLocus = (int*)calloc(pi->nloci, sizeof(int))) != 0);
	for (l = 0; l < pi->n_loci_to_use; ++l)
	{
		ll1 = pi->loci_to_use[l];
		for (s = 0; s < nsub; ++s)
		{
			for (c=0,a = 0; a < 2; ++a)
				c += sub[s]->all[ll1][a] == rarer[l];
			++genoCounts[l][c];
		}
		if (genoCounts[l][1] > 1)
		{
			contrib[ll1].index = ll1;
			contrib[ll1].contribution = weight[spi->numAddScores][l] * genoCounts[l][1];
			// homozygotes will be ignored
			// not interested in singleton heterozygotes
			// only the first set of recessive weights will be used to select most informative loci
			// if this is a problem, will either need to include all loci or run repeated analyses
		}
	}
	qsort(contrib, pi->nloci, sizeof(variantContribution), compContribution);
	for (ll1 = 0; ll1 < spi->maxRecLociToUse; ++ll1) // keep the best ones
	{
		if (contrib[ll1].contribution == 0 || ll1>=pi->nloci)
			break;
		else
			useThisRecLocus[contrib[ll1].index] = 1;
	}
	if (spi->df[DEBUGFILE].fp)
	{
		fprintf(spi->df[DEBUGFILE].fp, "Initial set of potentially informative loci for recessive analysis:\n");
		for (l = 0; l < pi->n_loci_to_use; ++l)
		{
			ll1 = pi->loci_to_use[l];
			if (useThisRecLocus[ll1] == 0)
				continue;
			fprintf(spi->df[DEBUGFILE].fp, "%-" LOCUS_NAME_LENGTH_STR "s\n", names[ll1]);
		}
	}
	for (l = 0; l < pi->n_loci_to_use; ++l)
	{
		ll1 = pi->loci_to_use[l];
		if (useThisRecLocus[ll1] == 0)
			continue;
		anyValidPair = 0;
		for (lll = 0; lll <= ll1; ++lll)
			useSecondRecLocus[lll] = 0;
		for (l2 = l + 1; l2 < pi->n_loci_to_use; ++l2)
		{
			ll2 = pi->loci_to_use[l2];
			if (useThisRecLocus[ll2] == 0)
				continue;
			thisPairValid = coOccur = 0;
			for (s = 0; s < nsub; ++s)
			{
					n1 = 0;
					n2 = 0;
					for (a = 0; a < 2; ++a)
					{
						n1 += sub[s]->all[ll1][a] == rarer[l];
						n2 += sub[s]->all[ll2][a] == rarer[l2];
					}
					if ((n1 == 2 && n2 > 0) || (n1 > 0 && n2 == 2))
					{
						thisPairValid = 0;
						if (spi->df[DEBUGFILE].fp)
						{
							fprintf(spi->df[DEBUGFILE].fp,
								"Excluding %-" LOCUS_NAME_LENGTH_STR "s and %-" LOCUS_NAME_LENGTH_STR "s because they have allele counts %d and %d in %s\n", 
								names[ll1], names[ll2], n1, n2, sub[s]->id);
						}
						break; // do not allow any where must be in same haplotype
					}
					else if (n1 == 1 && n2 == 1)
					{
						++coOccur;
						thisPairValid = 1;
					}
			}
			if (coOccur > spi->LDThreshold2022 * genoCounts[l][1] * genoCounts[l2][1] / nsub)
			{
				thisPairValid = 0;
				if (spi->df[DEBUGFILE].fp)
				{
					fprintf(spi->df[DEBUGFILE].fp, 
						"Excluding %-" LOCUS_NAME_LENGTH_STR "s and %-" LOCUS_NAME_LENGTH_STR "s because they co-occur too frequently, in %d subjects",
						names[ll1], names[ll2], coOccur);
				}
			}
			if (thisPairValid)
			{
				useSecondRecLocus[ll2] = 1;
				anyValidPair = 1;
			}
		}
		if (anyValidPair)
			usablePairs[ll1].fill(useSecondRecLocus, pi->nloci);
	}


	free(genoCounts);
	free(contrib);
	free(useThisRecLocus);
	free(useSecondRecLocus);
	return 1;
}

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
	float pairScore,**expandedWeights;
	locusIndex* usablePairs,usableLoci(pi->nloci);
	if (spi->use_probs)
	{
		dcerror(1,"Recessive analyses with probabilistic genotypes are not supported yet");
		exit(1);
	}
	usablePairs = new locusIndex[pi->nloci];
	getPairsToUseIndex(usablePairs,weight, rarer, sub, nsub, pi, spi);
	for (ll1 = 0; ll1 < pi->nloci; ++ll1)
		if (usablePairs[ll1].nLoci > 0)
			usableLoci.loci[usableLoci.nLoci++] = ll1;
	assert(expandedWeights = (float**)calloc(pi->nloci, sizeof(float*)));
	// note that expandedWeights is partial transpose of weight
	for (l = 0; l < pi->n_loci_to_use; ++l)
	{
		assert((expandedWeights[pi->loci_to_use[l]] = (float*)calloc(spi->numRecScores, sizeof(float))) != 0);
		for (sc = 0; sc < spi->numRecScores; ++sc)
			expandedWeights[pi->loci_to_use[l]][sc] = weight[spi->numAddScores + sc][l];
	}

	for (s = 0; s < nsub; ++s)
	{
		for (sc = 0; sc < spi->numRecScores; ++sc)
		{
			score[spi->numAddScores + sc][s]=0;
			rP[sc][s].l[0] = rP[sc][s].l[1] = -1;
		}
		for (l = 0; l < usableLoci.nLoci; ++l)
		{
			ll1 = usableLoci.loci[l];
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
			if (n1 !=1)
				continue; // ignoring homozygotes in all situations, worry about unknowns later
			for (l2 = 0; l2 < usablePairs[ll1].nLoci; ++l2)
			{
				ll2 = usablePairs[ll1].loci[l2];
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
				if (n2 != 1)
					continue;
				for (sc = 0; sc < spi->numRecScores; ++sc)
				{
					pairScore = (expandedWeights[ll1][sc] + expandedWeights[ll2][sc]);
					if (pairScore > score[spi->numAddScores + sc][s])
					{
						score[spi->numAddScores + sc][s] = pairScore;
						rP[sc][s].l[0] = ll1;
						rP[sc][s].l[1] = ll2;
					}
				}
			}
		}
	}
	delete[] usablePairs;
	for (l = 0; l < pi->n_loci_to_use; ++l)
		free(expandedWeights[pi->loci_to_use[l]]);
	free(expandedWeights);
}

void write_rec_scores(FILE* fs, subject** sub, int nsub, double** score, recPair** rP, int numRecScores, double** weight, char names[MAX_LOCI][LOCUS_NAME_LENGTH], par_info* pi, sa_par_info* spi)
{
	int s, sc,i,l,ll;
	float** expandedWeights;
	assert(expandedWeights = (float**)calloc(pi->nloci, sizeof(float*)));
	// note that expandedWeights is partial transpose of weight
	for (l = 0; l < pi->n_loci_to_use; ++l)
	{
		assert((expandedWeights[pi->loci_to_use[l]] = (float*)calloc(spi->numRecScores, sizeof(float))) != 0);
		for (sc = 0; sc < spi->numRecScores; ++sc)
			expandedWeights[pi->loci_to_use[l]][sc] = weight[spi->numAddScores + sc][l];
	}
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
				ll = rP[sc][s].l[i];
				if (ll != -1)
					fprintf(fs, "%-" LOCUS_NAME_LENGTH_STR "s %8.4f ", names[ll], expandedWeights[ll][spi->numAddScores + sc]);
				else
					fprintf(fs, "%-" LOCUS_NAME_LENGTH_STR "s %8.4f ", "NOLOCUS",0);
			}
		}
		fprintf(fs, "\n");
	}
	for (l = 0; l < pi->n_loci_to_use; ++l)
		free(expandedWeights[pi->loci_to_use[l]]);
	free(expandedWeights);
}
