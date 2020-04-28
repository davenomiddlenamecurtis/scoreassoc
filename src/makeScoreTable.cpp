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

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "dcerror.hpp"

#define PROGRAM "makeScoreTable"
#define MSTVERSION "1.2"

void usage()
{
	printf("Usage:\nmakeScoreTable tableName.txt *.sco\n");
}

// for MSDOS need to link setargv.obj to expand wildcards
int main(int argc,char *argv[])
{
	FILE *fi,*fo;
	char geneName[100],*ptr,fn[400],ID[100],line[400],pheno[100];
	int a,first;
	float score;
	printf("%s v%s\n",PROGRAM,MSTVERSION);
	if (argc<3)
		usage();
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
			while (fgets(line,399,fi) && sscanf(line,"%*s %s",pheno)==1)
				fprintf(fo,"%s\t",pheno);
			fprintf(fo,"\n");
			fseek(fi,SEEK_SET,0L);
		}

		fprintf(fo,"%s\t",geneName);
		while (fgets(line,399,fi) && sscanf(line,"%*s %*s %f",&score)==1)
			fprintf(fo,"%f\t",score);
		fprintf(fo,"\n");
		fclose(fi);
	}
	fclose(fo);
	return 0;
}