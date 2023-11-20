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
//extern void  read_coordinates();
extern void  lhgs_read_coordinates();
//extern void  get_pm_heads();
extern void  lhgs_get_heads();
//extern void  test_interface();
extern void  lhgs_test_interface();
//extern void  set_pm_heads();
extern void  lhgs_set_heads();
//extern void  assign_fnvars(char *p, int *plen, char *in, char *out);
extern void  lhgs_assign_hgsdat(char *p, int *plen, char *in, char *out);


/* fortran variables */
//extern double *__hgsdat_MOD_head;
extern double *lhgsvar_head;
extern double *lhgsvar_x;
extern double *lhgsvar_y;
extern double *lhgsvar_z;

extern int    hgsdat_MOD_nn;
extern char   hgsdat_MOD_message[80];
extern char   hgsdat_MOD_prefix[100];
extern char   hgsdat_MOD_insuffix[3];
extern char   hgsdat_MOD_outsuffix[3];

/* global variables */
/*
char    prefix[]="test";
char    idin[3]="066";
char    idout[3]="xxx";
*/

