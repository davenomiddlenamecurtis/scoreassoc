#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "dcerror.hpp"
#include <string>
#include <unordered_set>

#define PROGRAM "combineGeneScores"
#define MSTVERSION "1.1"

void usage()
{
	printf("Usage:\ncombineGeneScores geneSetTableName.txt geneSetList.txt nVarTypes *.sco\n");
}

struct geneSet_t { char name[100]; std::unordered_set<std::string> genes; };
#define MAXSETS 50
struct geneSet_t set[MAXSETS];
struct subject_t { char ID[50]; int cc; float **score; };


// for MSDOS need to link setargv.obj to expand wildcards
int main(int argc,char *argv[])
{
	FILE *fi,*fo,*fl,*fs;
	char geneName[100],*ptr,fn[400],ID[100],line[400],rest[400];
	int a,first,cc,nSet,nVarType,s,t,v,nSub,nMember,member[MAXSETS],m;
	float score,**scoreTable;
	struct subject_t *sub;
	printf("%s v%s\n",PROGRAM,MSTVERSION);
	if (argc<5)
		usage();
	assert((fl=fopen(argv[2],"r"))!=0);
	for (nSet=0;fgets(line,100,fl)&&sscanf(line,"%s",set[nSet].name)==1;++nSet)
	{
		assert((fs=fopen(set[nSet].name,"r"))!=0);
		while (fgets(line,100,fs)&&sscanf(line,"%s",geneName))
			set[nSet].genes.insert(geneName);
		fclose(fs);
	}
	fclose(fl);
	nVarType=atoi(argv[3]);
	assert((fi=fopen(argv[4],"r"))!=0);
	for (nSub=0;fgets(line,99,fi)&&sscanf(line,"%s",ID)==1;++nSub)
		;
	assert((sub=(struct subject_t *)calloc(nSub,sizeof(struct subject_t)))!=0);
	for (s=0;s<nSub;++s)
	{
		assert((sub[s].score=(float**)calloc(nSet,sizeof(float*)))!=0);
		for (t=0;t<nSet;++t)
			assert((sub[s].score[t]=(float*)calloc(nVarType,sizeof(float)))!=0);
	}
	fseek(fi,SEEK_SET,0L);
	for (s=0;s<nSub;++s)
		assert(fgets(line,99,fi)&&sscanf(line,"%s %d",sub[s].ID,&sub[s].cc)==2);
	fclose(fi);
	for (a=4;a<argc;++a)
	{
		strcpy(fn,argv[a]);
		fi=fopen(fn,"r");
		if (!fi)
		{
			dcerror(1,"Could not open file %s\n",fn);
			exit(1);
		}
		ptr=strstr(fn,".vsco");
		if (ptr==0)
		{
			dcerror(1,"Could not find .vsco in filename: %s\n",fn);
			exit(1);
		}
		*ptr='\0';
		while (*(--ptr)!='.')
			if (ptr<fn)
			{
				dcerror(1,"Could not find second . in filename: %s\n",argv[a]);
				exit(1);
			}
		strcpy(geneName,ptr+1); // relies on no gene containng .
		for (nMember=0,t=0;t<nSet;++t)
		{
			auto it=set[t].genes.find(geneName);
			if (it!=set[t].genes.end())
				member[nMember++]=t;
		}
		for (s=0;s<nSub;++s)
		{
			assert(fgets(line,99,fi)&&sscanf(line,"%*s %*d %[^\n]",rest)==1);
			for (v=0;v<nVarType;++v)
			{
				strcpy(line,rest);
				*rest='\0';
				assert(sscanf(line,"%f %[^\n]",&score,rest)>=1);
				for (m=0;m<nMember;++m)
					sub[s].score[member[m]][v]+=score;
			}
		}
		fclose(fi);
	}
	fo=fopen(argv[1],"w");
	if (!fo)
	{
		dcerror(1,"Could not open file %s\n",argv[1]);
		exit(1);
	}
	fprintf(fo,"ID\tCC\t");
	for (t=0;t<nSet;++t)
		for (v=0;v<nVarType;++v)
			fprintf(fo,"%s:%d\t",set[t].name,v);
	fprintf(fo,"\n");
	for (s=0;s<nSub;++s)
	{
		fprintf(fo,"%s\t%d\t",sub[s].ID,sub[s].cc);
		for (t=0;t<nSet;++t)
			for (v=0;v<nVarType;++v)
				fprintf(fo,"%8.4f\t",sub[s].score[t][v]);
		fprintf(fo,"\n");
	}
	fclose(fo);
	return 0;
}