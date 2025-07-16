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

#ifndef SAFILTERFUNCSHPP
#define SAFILTERFUNCSHPP

#include "scoreassoc.hpp"

#ifndef USEFILTERS
#define USEFILTERS
#endif

int initExclusions(FILE *fp,char *extras[]=0);
int applyExclusions(subject * *sub, int nsub, par_info *pi);
int stateExclusions(FILE *fp);

#endif