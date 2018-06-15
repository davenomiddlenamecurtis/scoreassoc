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

#ifndef dvector_HPP
#define dvector_HPP 1

// header file for double precision dvector and matrix classes
// doublevec - double dvector indexed from 0
// dvector - double dvector indexed from offset, can't remember what it does
// dv2d - double matrix indexed from 0, supports matrix arithmetic
// dsvd - set of three matrices forming singular value decomposition of source matrix
// rank_table - not distributing this code yet


#include <stdio.h>

class doublevec{
protected:
double *v;
int sz;
public:
doublevec(doublevec const &);
doublevec(int);
doublevec(void);
doublevec& operator=(doublevec const &);
doublevec operator+(doublevec const &);
doublevec operator-(doublevec const &);
int size(void) {return sz;}
double* vec(void) {return v;}
~doublevec();
void fprintf(FILE *fp,char *form);
void print(char *form= "%3.2g  ") { fprintf(stdout,form); }
double& operator[](int i){return v[i];}
};

class dvector:public doublevec{
int offset; // if starts at non-zero index
public:
dvector(int l,int h):doublevec(h-l+1) {offset=l;}
~dvector() {;}
double& operator[](int i) {return v[i-offset];}
double* vec(void){return v;}
};

typedef doublevec *DBLVPOINTER;
typedef double *DBLPOINTER;

class dv2d{
double **rv;
int wide,high,max_h;
public:
dv2d(int h=0,int w=0);
dv2d(dv2d const &);
~dv2d();
double **vec() { return rv; }
int get_width(void) const {return wide;}
int get_height(void) const {return high;}
long size(void){return long(wide)*high;}
dv2d& operator=(dv2d const &);
dv2d operator+(dv2d const &);
dv2d operator*(dv2d const &);
dv2d operator*(double);
dv2d& operator*=(dv2d const &);
dv2d& operator+=(dv2d const &);
dv2d& operator*=(double);
dv2d operator-(dv2d const &);
DBLPOINTER operator[](int i) {return rv[i];}
// double determ();
dv2d inv();
// int dsvdcmp(dv2d&,dv2d&);
dv2d transpose();
int add_row(int rownum=-1);
int del_row(int rownum=-1);
int add_col(int colnum=-1);
int del_col(int colnum=-1);
void sort_rows(int,int);
int resize(int h,int w);
void fprintf(FILE *fp=stdout,char *form= "%3.2g  ");
void fscanf(FILE *fp=stdin);
void print(char *form= "%3.2g  ") { fprintf(stdout,form); }
friend class dsvd;
void normal(int);
};

class dsvd {
public:
dv2d U,W,V,UT,VT,Winv;
doublevec mean,sd;
// keep decomposition in UT,W,VT
// it may often be more efficient, if slightly less convenient
// to use UT and VT directly rather than U and V
int done,gotU,gotV,gotWinv;
int OK(void);
dsvd(int,int);
dsvd(dsvd const &);
//dsvd(dv2d&);
~dsvd(){};
dv2d solve(dv2d const &);
dv2d svbksb(dv2d const &);
int setU();
int setV();
int setWinv(double lim=1e-6);
void zero(double lim=1e-6);
dv2d operator*(dv2d const &);
int dcmp(dv2d &);
dsvd& operator=(dsvd const &);
void normal(int code);
};

class rank_table {
public:
double val_max_ks;
rank_table(int hi,int ng=2);
~rank_table() {;}
dv2d table;
doublevec index;
int num, ndiff, ngroups;
double median,mode,last_rank,quartile1,quartile2;
int add(double,int,int);
int do_ranks();
double get_ks();
void rank_sum(double *,int *,int *);
double get_rank(int);
double get_rank(double);
};

extern int default_dvec_error(char *s);
extern int (*dvec_error)(char*);
extern "C" {
extern int comp_row(void const *,void const *);
}

#endif

