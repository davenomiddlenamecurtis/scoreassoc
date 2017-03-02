#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "dcerror.hpp"

int main(int argc,char *argv[])
{
	FILE *fi,*fo;
	char geneName[100],*ptr,fn[400],ID[100],line[400];
	int a,first,cc;
	float score;
	fo=fopen(argv[1],"w");
	if (!fo)
	{
		dcerror(1,"Could not open file %s\n",argv[1]);
		exit(1);
	}
	first=1;
	for (a=2;a<argc;++a)
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