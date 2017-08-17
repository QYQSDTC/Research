//  Copyright (C) 2006,2007,2008,2009,2010,2011 George Hobbs, Russell Edwards

/*
 *    This file is part of TEMPO2. 
 * 
 *    TEMPO2 is free software: you can redistribute it and/or modify 
 *    it under the terms of the GNU General Public License as published by 
 *    the Free Software Foundation, either version 3 of the License, or 
 *    (at your option) any later version. 
 *    TEMPO2 is distributed in the hope that it will be useful, 
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of 
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 *    GNU General Public License for more details. 
 *    You should have received a copy of the GNU General Public License 
 *    along with TEMPO2.  If not, see <http://www.gnu.org/licenses/>. 
 */

/*
 *    If you use TEMPO2 then please acknowledge it by citing 
 *    Hobbs, Edwards & Manchester (2006) MNRAS, Vol 369, Issue 2, 
 *    pp. 655-672 (bibtex: 2006MNRAS.369..655H)
 *    or Edwards, Hobbs & Manchester (2006) MNRAS, VOl 372, Issue 4,
 *    pp. 1549-1574 (bibtex: 2006MNRAS.372.1549E) when discussing the
 *    timing model.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "tempo2.h"

void help() /* Display help */
{
}


extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
   int i,p;
   char parFile[MAX_PSR][MAX_FILELEN];
   char timFile[MAX_PSR][MAX_FILELEN];
   char fname[MAX_FILELEN];
   *npsr=0;
   printf("Graphical Interface: flattenErrors\n");
   printf("Author:              M. Keith\n");
   printf("CVS Version:         $Revision: 1.14 $\n");
   printf(" --- don't bother typing 'h' for help information\n");


   for (i=2;i<argc;i++)
   {
	  if (strcmp(argv[i],"-f")==0) {
		 strcpy(parFile[*npsr],argv[++i]); 
		 strcpy(timFile[*npsr],argv[++i]);
		 (*npsr)++;
	  }
   }


   readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
   readTimfile(psr,timFile,*npsr); /* Load the arrival times    */
   preProcess(psr,*npsr,argc,argv);

   for (p=0; p < *npsr; p++){
	  psr[p].nT2equad=0;
	  psr[p].nT2efac=0;
	  for (i=0; i < psr[p].nobs; i++){
		 psr[p].obsn[i].origErr=psr[p].obsn[i].toaErr;
		 psr[p].obsn[i].efac=1;
		 psr[p].obsn[i].equad=0;
	  }
	  sprintf(fname,"%s.flat.tim",psr[p].name);
	  writeTim(fname,psr+p,"tempo2");
   }
   return 0;
}

