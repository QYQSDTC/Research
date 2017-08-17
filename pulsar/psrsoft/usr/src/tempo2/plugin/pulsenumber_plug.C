//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russel Edwards

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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tempo2.h"


extern "C" int tempoOutput(int argc,char *argv[],pulsar *psr,int npsr) 
{  
   int i;
   for (i=0;i<psr[0].nobs;i++)
   {
	  if ((double)psr[0].obsn[i].sat > psr[0].param[param_start].val[0] && (double)psr[0].obsn[i].sat < psr[0].param[param_finish].val[0]){
		 printf("%g %g\n",(double)psr[0].obsn[i].sat,(double)(psr[0].obsn[i].pulseN - psr[0].obsn[0].pulseN));
		 strcpy(psr[0].obsn[i].flagID[psr[0].obsn[i].nFlags],"-pn");
		 sprintf(psr[0].obsn[i].flagVal[psr[0].obsn[i].nFlags],"%lld",psr[0].obsn[i].pulseN-psr[0].obsn[0].pulseN);
		 psr[0].obsn[i].nFlags++;
	  }
   }
   writeTim("pn.tim",psr,"tempo2");
}

