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
struct subject_t { char ID[50]; float **score; };


// for MSDOS need to link setargv.obj to expand wildcards
int main(int argc,char *argv[])
{
	FILE *fi,*fo,*fl,*fs;
	char geneName[100],*ptr,fn[400],ID[100],line[400];
	int a,first,cc,nSet,nVarType,s,t,v,nSub;
	float score,**scoreTable;
	struct subject_t *sub;
	printf("%s v%s\n",PROGRAM,MSTVERSION);
	if (argc<6)
		usage();
	assert((fl=fopen(argv[3],"r"))!=0);
	for (nSet=0;fgets(line,100,fl)&&sscanf(line,"%s",set[nSet].name)==1;++nSet)
	{
		assert((fs=fopen(set[nSet].name,"r"))!=0);
		while (fgets(line,100,fs)&&sscanf(line,"%s",geneName))
			set[nSet].genes.insert(geneName);
		fclose(fs);
	}
	fclose(fl);
	nVarType=atoi(argv[4]);
	fo=fopen(argv[1],"w");
	if (!fo)
	{
		dcerror(1,"Could not open file %s\n",argv[1]);
		exit(1);
	}
	fprintf(fo,"ID\tCC\t");
	for (t=0;t<nSet;++t)
		for (v=0;v<nVarType;++v)
			fprintf(fo,"%s:%d",set[t].name,v);
	fprintf(fo,"\n");
	assert((fi=fopen(argv[5],"r"))!=0);
	
	for (a=5;a<argc;++a)
	{
		strcpy(fn,argv[a]);
		fi=fopen(fn,"r");
		if (!fo)
		{
			dcerror(1,"Could not open file %s\n",fn);
			exit(1);
		}
		ptr=strstr(fn,".sco");
		if (ptr==0)
		{
			dcerror(1,"Could not find .sco in filename: %s\n",fn);
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
		if (first)
		{
			first=0;
			fprintf(fo,"ID\t");
			while (fgets(line,399,fi) && sscanf(line,"%s",ID)==1)
				fprintf(fo,"%s\t",ID);
			fprintf(fo,"\n");
			fseek(fi,SEEK_SET,0L);
			fprintf(fo,"CC\t");
			while (fgets(line,399,fi) && sscanf(line,"%*s %d",&cc)==1)
				fprintf(fo,"%d\t",cc);
			fprintf(fo,"\n");
			fseek(fi,SEEK_SET,0L);
		}

		fprintf(fo,"%s\t",geneName);
		while (fgets(line,399,fi) && sscanf(line,"%*s %*d %f",&score)==1)
			fprintf(fo,"%f\t",score);
		fprintf(fo,"\n");
		fclose(fi);
	}
	fclose(fo);
	return 0;
}