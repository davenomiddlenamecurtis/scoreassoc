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

										/* DCERROR.HPP */
/* Copyright Dave Curtis 1994 */ 
/* dcurtis@hgmp.mrc.ac.uk */ 
/* no warranty or liability of any kind is accepted, expressed or implied */ 
 
/* DCERROR.HPP */ 
/* Copyright Dave Curtis 1991 */ 
#ifndef DCERRORHPP
#define DCERRORHPP 1

#if defined wx_msw || defined wx_xview || defined wx_motif
#include "allwx.hpp"
#endif

#include <stdarg.h>
#include <errno.h>

extern void fg_error_display(char *s,va_list arg_ptr); // lives in bitmap.cpp

class error_object
{
friend void fg_error_display(char *s,va_list arg_ptr);
int fatal,showing,status;//,fg_inited;
//va_list arg_ptr;
void (*display_func)(char *,va_list);
//void show_message(char *format);
public:
int operator()(int e,char *format,...);
int operator()(int e) { return operator()(e,0); }
error_object();
~error_object();
void set_display(void (*df)(char *,va_list));
void reset_display() { display_func=0; }
int stat();
void clear();
void hide();
void show();
//void set_fg(int i);
void warn();
void kill();
};

extern error_object dcerror;

#ifndef hereOK
#define hereOK() fprintf(stderr,"Got to line %d in %s OK\n",__LINE__,__FILE__)
#endif

#endif

