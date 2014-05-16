#include <stdio.h>
#include "idl_export.h"

/*
DLM for provisional means algorithm, called from the object
class CPM__DEFINE in order to avoid IDL for-loop.
Compile using IDL procedure MAKE_DLL.
Under MS Windows OS, use MAKE_PROVMEANS.PRO,
assuming VC++ 6.0 or later compiler is installed.
For other OS, see IDL documentation for MAKE_DLL.
---------------------------------------------------------*/

IDL_VPTR IDL_CDECL provmeans(int argc, IDL_VPTR argv[]) 
{
   float  *Xs, *Ws, w; 
   double *mn, *cov, *sw, r, d[400];
   IDL_LONG i,j,k,N,n;
   IDL_VPTR ivReturn = IDL_Gettmp();
   ivReturn->type = IDL_TYP_INT;
   ivReturn->value.i = 0;
/* get the data fields */   
   Xs = (float *)  argv[0]->value.arr->data;
   Ws = (float *)  argv[1]->value.arr->data;
   N  = (IDL_LONG) argv[2]->value.l;
   n  = (IDL_LONG) argv[3]->value.l;
   sw = (double *) &(argv[4]->value.d);
   mn = (double *) argv[5]->value.arr->data;
   cov= (double *) argv[6]->value.arr->data;    
/* loop over observation vectors */   
   for (i=0; i<n; i++)
     {
      w = *Ws;
      *sw += *Ws++;
      r = w/(*sw); 
/* mean */
      for (j=0; j<N; j++)
         {
           d[j] = Xs[i*N+j]-mn[j];
           mn[j] += d[j]*r;
         }
/* covariance */
      for (j=0; j<N; j++)
        for (k=j; k<N; k++) cov[j*N+k] += d[j]*d[k]*(1-r)*w;	  
     };  
/* have to return an IDL_VPTR */                                     
   return ivReturn; 
}

/* the entry point, which loads the routine into IDL */
int IDL_Load(void) 
{ 
  static IDL_SYSFUN_DEF2 function_addr[] = { 
    { provmeans, "PROVMEANS", 0, 7, 0, 0 } 
  }; 
  return IDL_SysRtnAdd(function_addr, IDL_TRUE, 1); 
}  

