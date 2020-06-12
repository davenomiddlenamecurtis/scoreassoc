#!/bin/bash
if [ .$2 == . ]
then
  echo Usage:  $0 folderContainingSAOFiles outputSummaryFileName
  exit
fi

workFolder=$1
summFile=$2

# The following line needs to be edited so that the column headings will match the tests contained in the SAO files
echo Gene$'\t'SLP$'\t'tSLP$'\t'tSLPPC> $summFile
 

getSLPs='
BEGIN { ORS=""; nSLP=0; } 
{
if ($1 == "SLP" || $1 == "tSLP" || $1 == "tMLP" || $1 == "linrSLP" || $1 == "lrSLP") 
	{
	nSLP=nSLP+1;
	SLPs[nSLP]=$3;
	}
}
END {
    if (nSLP > 0) {
	print gene "\t";
	for (i=1; i<=nSLP; ++i) {
	print SLPs[i] "\t";
	}
	print "\n";
	}
}
'


find  $resultsFolder -name '*.sao' | while read resultsFile
	do
	gene=${resultsFile%.sao}
	gene=${gene#*$model.}
	awk -v gene=$gene "$getSLPs" $resultsFile >> $summFile
	done
