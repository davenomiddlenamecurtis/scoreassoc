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

#include "lrBurdenAssoc.hpp"

// lazily make these global to save allocating them and to avoid stack overflow
float weight[MAX_LOCI],missing_score[MAX_LOCI],func_weight[MAX_LOCI],cc_freq[2][MAX_LOCI],cc_count[2][MAX_LOCI],cc_genocount[2][3][MAX_LOCI];
int rarer[MAX_LOCI],max_cc[2];
char names[MAX_LOCI][20],comments[MAX_LOCI][MAX_COMMENT_LENGTH],trios_fn[500];
subject **global_sub;
lrVariable allVars[MAXLRVARIABLES];
std::map<std::string, lrVariable *> varMap;
