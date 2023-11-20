#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* defines */

/* c functions */
//void assign_fnvars_fortran(char *pre, char *infile, char *outfile);

/* fortran functions */
/*
extern void  hgs_get_nnodes_(char *string, int *len_string, int *nn);
extern void  hgs_get_head_(char *string, int *len_string, double *hh, int *nn, char *mess_string);
extern void  hgs_get_coordinates_(char *string, int *len_string, double *x, double *y, double *z, int *nx, int *ny, int *nz, int *nn);
*/
extern void  read_coordinates_();
extern void  get_pm_heads_();
extern void  test_interface_();
extern void  set_pm_heads_();
extern void  assign_fnvars_(char *p, int *plen, char *in, char *out);


/* fortran variables */
extern double *__hgsdat_mp_head_;
extern int    __hgsdat_mp_nn_;
extern char   __hgsdat_mp_message_[80];
extern char   __hgsdat_mp_prefix_[100];
extern char   __hgsdat_mp_insuffix_[3];
extern char   __hgsdat_mp_outsuffix_[3];

/* global variables */
/*
char    prefix[]="test";
char    idin[3]="066";
char    idout[3]="xxx";
*/

