/* Copyright 2008 BORDEAUX I UNIVERSITY & INRIA 
**
** This file is part of the PaStiX parallel sparse matrix package.
**
** This software is governed by the CeCILL-C license under French law
** and abiding by the rules of distribution of free software. You can
** use, modify and/or redistribute the software under the terms of the
** CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
** URL: "http://www.cecill.info".
** 
** As a counterpart to the access to the source code and rights to copy,
** modify and redistribute granted by the license, users are provided
** only with a limited warranty and the software's author, the holder of
** the economic rights, and the successive licensors have only limited
** liability.
** 
** In this respect, the user's attention is drawn to the risks associated
** with loading, using, modifying and/or developing or reproducing the
** software by the user in light of its specific status of free software,
** that may mean that it is complicated to manipulate, and that also
** therefore means that it is reserved for developers and experienced
** professionals having in-depth computer knowledge. Users are therefore
** encouraged to load and test the software's suitability as regards
** their requirements in conditions enabling the security of their
** systems and/or data to be ensured and, more generally, to use and
** operate it in the same conditions as regards security.
** 
** The fact that you are presently reading this means that you have had
** knowledge of the CeCILL-C license and that you accept its terms.
*/
/*== Creation/Destruction de CSC ==*/

/******************************************************************************/
/* void CscOrder(CscMatrix *thecsc, char *Type, char *RhsType, FLOAT **rhs,  */
/*                FLOAT **rhs2, const Order *ord, ...)                        */
/*                                                                            */
/* Construction de la csc a partir de col/row/val et permutation a            */
/* partir du vecteur permutation fournit par Scotch                           */
/*                                                                            */
/* thecsc : La csc                                                            */
/* Type : type du HB (RUA, RSA ....)                                          */
/* RhsType : type pour les seconds membres                                    */
/* rhs : Vecteur second membre                                                */
/* rhs2 : Vecteur solution                                                    */
/* ord : la permutation                                                       */
/*                                                                            */
/* Type doit etre alloue avant l'appel, char Type[4]                          */
/* RhsType doit etre alloue avant l'appel, char RhsType[4]                    */
/* rhs et rhs2 sont alloue si besoin est.                                     */
/* RhsType[0] == '\0' si pas de second membre dans le fichier                 */
/******************************************************************************/

/* !!! FONCTION INUTILISEE !!! */
void CscOrder(CscMatrix *thecsc,
	      char *Type, char *RhsType, FLOAT **transcsc,
	      FLOAT **rhs, FLOAT **rhs2,
	      const Order *ord, 
	      INT Nrow, INT Ncol, INT Nnzero, 
	      INT *colptr, INT *rowind, FLOAT *val, INT forcetrans)
{
  char *crhs=NULL; /* 2nd member vector */
  char *crhs2=NULL; /* solution vector */ 

  INT index,itercol,iter,newcol,colidx, rowp1; 
  INT valnbr=0;
  INT *trscltb=NULL; /* coltab for transpose csc */
  INT *trowtab=NULL; /* rowtab for transpose csc */

#ifdef CSC_LOG
  fprintf(stdout, "-> CscOrder \n");
#endif

  /* Csc alloc */
  CSC_FNBR(thecsc) = 1;
  CSC_FTAB(thecsc) = (CscFormat *) memAlloc(1*sizeof(CscFormat));
  if (CSC_FTAB(thecsc) == NULL)
    errorPrint( "CscHbRead : Not enough memory for CSC_FTAB\n");
  CSC_COLNBR(thecsc,0) = Ncol;
  CSC_COLTAB(thecsc,0) = (INT *) memAlloc((Ncol+1)*sizeof(INT));
  if (CSC_COLTAB(thecsc,0) == NULL)
    errorPrint( "CscHbRead : Not enough memory for CSC_COLTAB\n");

  for (index=0; index < (Ncol+1); index++)
    CSC_COL(thecsc,0,index) = 0;

  if (Type[1] == 'S')
    {
      /* Symetric */
      CSC_ROWTAB(thecsc) = (INT*) memAlloc(2*Nnzero*sizeof(INT));
      if (CSC_ROWTAB(thecsc) == NULL)
	errorPrint( "CscHbRead : Not enough memory for CSC_ROWTAB\n");
      CSC_VALTAB(thecsc) = (FLOAT*) memAlloc(2*Nnzero*sizeof(FLOAT));
      if (CSC_VALTAB(thecsc) == NULL)
	errorPrint( "CscHbRead : Not enough memory for CSC_VALTAB\n");
      if (forcetrans)
	{
	  printf("Force transpose on symetric...\n");

	  (*transcsc) = (FLOAT *) memAlloc(2*Nnzero*sizeof(FLOAT));
	  if ((*transcsc) == NULL)
	    errorPrint( "CscHbRead : Not enough memory for (*transcsc)\n");
	  trowtab = (INT *) memAlloc(2*Nnzero*sizeof(INT));
	  if (trowtab == NULL)
	    errorPrint( "CscHbRead : Not enough memory for trowtab\n");
	}
    }
  else
    {
      /* Unsymmetric */
      CSC_ROWTAB(thecsc) = (INT*) memAlloc(Nnzero*sizeof(INT));
      if (CSC_ROWTAB(thecsc) == NULL)
	errorPrint( "CscHbRead : Not enough memory for CSC_ROWTAB\n");
      CSC_VALTAB(thecsc) = (FLOAT*) memAlloc(Nnzero*sizeof(FLOAT));
      if (CSC_VALTAB(thecsc) == NULL)
	errorPrint( "CscHbRead : Not enough memory for CSC_VALTAB\n");
      
      (*transcsc) = (FLOAT *) memAlloc(Nnzero*sizeof(FLOAT));
      if ((*transcsc) == NULL)
	errorPrint( "CscHbRead : Not enough memory for (*transcsc)\n");
      trowtab = (INT *) memAlloc(Nnzero*sizeof(INT));
      if (trowtab == NULL)
	errorPrint( "CscHbRead : Not enough memory for trowtab\n");
    }

  /* Computing good coltabs */
  for (itercol = 0; itercol < Ncol; itercol++)
    {
      newcol = ord->permtab[itercol];

      CSC_COL(thecsc,0,newcol) += colptr[itercol+1] - colptr[itercol];

      if (Type[1] == 'S')
	{
	  /* Symmetric */
	  for (iter=colptr[itercol]; iter<colptr[itercol+1]; iter++)
	    {
	      if ((rowind[iter-1]-1) != itercol)
		{
		  newcol = ord->permtab[rowind[iter-1]-1];
		  (CSC_COL(thecsc,0,newcol))++;
		}
	    }
	}
    }

  newcol = 0;
  for (index=0; index<(Ncol+1); index++)
    {
      colidx = CSC_COL(thecsc,0,index);
      CSC_COL(thecsc,0,index) = newcol;
      newcol += colidx;
    }

  if ((*transcsc) != NULL)
    {
      trscltb = (INT *) memAlloc((Ncol+1)*sizeof(INT));
      if (trscltb == NULL)
	errorPrint( "CscHbRead : not enough memory for trscltb\n");
      for (index=0; index<(Ncol+1); index++)
	{
	  trscltb[index] = CSC_COL(thecsc,0,index);
	}
    }

  /* Put the element */

  printf("Debut construction CSC\n");
  for (itercol = 0; itercol <Ncol; itercol++)
    {
      for (iter=colptr[itercol]; iter<colptr[itercol+1]; iter++)
	{
	  INT therow;

	  rowp1 = rowind[iter-1] -1;
	  newcol = ord->permtab[itercol];
	  colidx = CSC_COL(thecsc,0,newcol);
	  therow = ord->permtab[rowp1];
	  CSC_ROW(thecsc,colidx) = therow;
	  CSC_VAL(thecsc,colidx) = val[iter-1];	  
	  valnbr++;

	  CSC_COL(thecsc,0,newcol)++;

	  if ((*transcsc) != NULL)
	    {
	      (*transcsc)[trscltb[therow]] = val[iter-1];	      
	      trowtab[trscltb[therow]] = newcol;
	      trscltb[therow]++;
	    }


	  /* Symmetric */
	  if (Type[1] == 'S')
	    {
	      if (rowp1 != itercol)
		{
		  newcol = ord->permtab[rowp1];
		  colidx = CSC_COL(thecsc,0,newcol);
		  therow = ord->permtab[itercol];
		  CSC_ROW(thecsc,colidx) = therow;
		  CSC_VAL(thecsc,colidx) = val[iter-1];
		  valnbr++;

		  (CSC_COL(thecsc,0,newcol))++;
		  
		  if ((*transcsc) != NULL)
		    {
		      (*transcsc)[trscltb[therow]] = val[iter-1];
		      trowtab[trscltb[therow]] = newcol;
		      trscltb[therow]++;
		    }
		}
	    }
	  
	}
    }
  memFree_null(colptr);
  memFree_null(rowind);
  memFree_null(val);
  if (trscltb != NULL)
    memFree_null(trscltb);
  printf("Fin construction CSC\n");

  /* 2nd member */
  if (RhsType[0] != '\0')
    {
      printf("Not yet implemented\n");
      EXIT(MOD_SI,NOTIMPLEMENTED_ERR);

      (*rhs) = (FLOAT*) memAlloc(Ncol*sizeof(FLOAT));
      if ((*rhs) == NULL)
	errorPrint( "CscHbRead : Not enough memory for rhs\n");
      
      for (index=0; index<Ncol; index++)
	{
	  (*rhs)[ord->permtab[index]] = crhs[index];
	}

      if (RhsType[2] == 'X')
	{
	  /* Vector Solution */
	  (*rhs2) = (FLOAT*) memAlloc(Ncol*sizeof(FLOAT));
	  if ((*rhs2) == NULL)
	    errorPrint( "CscHbRead : not enough memory for rhs2\n");

	  for (index=0; index<Ncol; index++)
	    {
	      (*rhs2)[ord->permtab[index]] = crhs2[index];
	    }
	}
      
      memFree_null(crhs);
      memFree_null(crhs2);
    }
  else
    {
      (*rhs) = NULL;
      (*rhs2) = NULL;
    }
  printf("valnbr = %ld\n", (long)valnbr);

  /* good coltab */
  colidx = 0;
  for (index=0; index<Ncol; index++)
    {
      newcol = CSC_COL(thecsc,0,index);
      CSC_COL(thecsc,0,index) = colidx;
      colidx = newcol;
    }

  /* Sort on the row */
  for (index = 0; index < Ncol; index++)
    {
      /* bubble sort between coltab[index2] and coltab[index2+1] */
      INT *t = &(CSC_FROW(thecsc,0,index));
      FLOAT *v = &(CSC_FVAL(thecsc,0,index));
      INT n = CSC_COL(thecsc,0,index+1) - CSC_COL(thecsc,0,index);
	
      INT i,j;
      for (i=0; i<n; i++)
	for (j=0; j<n-i-1; j++)
	  if (t[j] > t[j+1])
	    {
	      INT tempt = t[j+1];
	      FLOAT tempv = v[j+1];
	      
	      t[j+1] = t[j];
	      v[j+1] = v[j];
	      t[j] = tempt;
	      v[j] = tempv;
	    }
    }

  if ((*transcsc) != NULL)
    {
      for (index=0; index<Ncol; index++)
	{
	  INT *t = &(trowtab[CSC_COL(thecsc,0,index)]);
	  FLOAT *v = &((*transcsc)[CSC_COL(thecsc,0,index)]);
	  
	  INT n = CSC_COL(thecsc,0,index+1) - CSC_COL(thecsc,0,index);
	  INT i,j;
	  
	  for (i=0; i<n; i++)
	    for (j=0; j<n-i-1; j++)
	      if (t[j] > t[j+1])
		{
		  INT tempt = t[j+1];
		  FLOAT tempv = v[j+1];
		  
		  t[j+1] = t[j];
		  v[j+1] = v[j];
		  t[j] = tempt;
		  v[j] = tempv;
		}
	}
      memFree_null(trowtab);
    }
#ifdef CSC_LOG
  fprintf(stdout, "<- CscOrder \n");
#endif
}

/*== Distribution/Remplissage ==*/
/******************************************************************************/
/* void CscDistrib(SymbolMatrix *symbmtx, CscMatrix *cscmtx,                  */
/*                 CscMatrix *cscdist, const FLOAT *transcsc,                 */
/*		   FLOAT **trandcsc)                                          */
/*                                                                            */
/* Distribution de la csc                                                     */
/*                                                                            */
/* symbmtx : symbol matrix locale                                             */
/* cscmtx : csc globale                                                       */
/* cscdist : csc locale                                                       */
/* transcsc : csr globale                                                     */
/* trandcsc : csr locale                                                      */
/******************************************************************************/
/* !!! FONCTION INUTILISEE !!! */
void CscDistrib(const SymbolMatrix *symbmtx, const CscMatrix *cscmtx,
		CscMatrix * cscdist, const FLOAT *transcsc, FLOAT **trandcsc)
{
  INT iterckcd=0; /* iter for cblk in csc(d) */
  INT itercol; /* iter for coltab */
  INT strdcol=0; /* stride for col */
  INT itervald=0; /* iter for valtab rowtab */
  INT itervalg=0;

#ifdef CSC_LOG
  fprintf(stdout, "-> CscDistrib \n");
#endif
  
  CSC_FNBR(cscdist) = symbmtx->cblknbr;
  CSC_FTAB(cscdist) = (CscFormat *) memAlloc(CSC_FNBR(cscdist)*sizeof(CscFormat));
  if (CSC_FTAB(cscdist) == NULL)
    errorPrint( "CscDistrib : Not enough memory for CSC_FTAB\n");
  
  
  for (iterckcd=0;
       iterckcd < symbmtx->cblknbr;
       iterckcd++)
    {
      CSC_COLNBR(cscdist,iterckcd) =
	symbmtx->cblktab[iterckcd].lcolnum-symbmtx->cblktab[iterckcd].fcolnum+1;
      CSC_COLTAB(cscdist,iterckcd) =
	(INT *) memAlloc((CSC_COLNBR(cscdist,iterckcd)+1)*sizeof(INT));
      if (CSC_COLTAB(cscdist,iterckcd) == NULL)
	errorPrint( "CscDistrib : Not enough memory for CSC_COLTAB\n");

      for (itercol = 0;
	   itercol < (CSC_COLNBR(cscdist,iterckcd)+1);
	   itercol++)
	{
	  CSC_COL(cscdist,iterckcd,itercol) =
	    CSC_COL(cscmtx,0,itercol+symbmtx->cblktab[iterckcd].fcolnum) -
	    CSC_COL(cscmtx,0,symbmtx->cblktab[iterckcd].fcolnum) +
	    strdcol;
	}
      strdcol = CSC_COL(cscdist,iterckcd,CSC_COLNBR(cscdist,iterckcd));
    }

  /* Remplissage des valtabs et rowtabs */
  CSC_ROWTAB(cscdist) = (INT *) memAlloc(strdcol*sizeof(INT));
  if (CSC_ROWTAB(cscdist) == NULL)
    errorPrint( "CscDistrib : Not enough memory for CSC_ROWTAB\n");
  CSC_VALTAB(cscdist) = (FLOAT *) memAlloc(strdcol*sizeof(FLOAT));
  if (CSC_VALTAB(cscdist) == NULL)
    errorPrint( "CscDistrib : Not enough memory for CSC_VALTAB\n");

  if (transcsc != NULL)
    {
      /* Rua */
      (*trandcsc) = (FLOAT *) memAlloc(strdcol*sizeof(FLOAT));
      if ((*trandcsc) == NULL)
	errorPrint( "CscDistrib : Not enough memory for trandcsc\n");
    }
  
  for (iterckcd = 0;
       iterckcd < CSC_FNBR(cscdist);
       iterckcd++)
    {
      for (itercol = 0;
	   itercol < CSC_COLNBR(cscdist,iterckcd);
	   itercol++)
	{
	  itervalg =
	    CSC_COL(cscmtx,0,itercol+symbmtx->cblktab[iterckcd].fcolnum);
	  
	  for (itervald = CSC_COL(cscdist,iterckcd,itercol);
	       itervald < CSC_COL(cscdist,iterckcd,itercol+1);
	       itervald++)
	    {
	      CSC_VAL(cscdist,itervald) = CSC_VAL(cscmtx,itervalg);

	      if (transcsc != NULL)
		(*trandcsc)[itervald] = transcsc[itervalg];
	      
	      CSC_ROW(cscdist,itervald) = CSC_ROW(cscmtx,itervalg);
	      itervalg++;
	    }
	}
    }

#ifdef CSC_LOG
  fprintf(stdout, "<- CscDistrib \n");
#endif
}


/******************************************************************************/
/* void Csc2solv(CscMatrix *cscmtx, SolverMatrix *solvmtx, FLOAT *trandcsc)   */
/*                                                                            */
/* Remplit la solvermatrix locale a partir de la csc locale, pour la partie   */
/* on utiliser la transpose de la csc globale qui est lu a partir d'un fichier*/
/*                                                                            */
/* cscmtx : csc locale                                                        */
/* solvmtx : solver locale                                                    */
/* stream : fichier de la csc transposee ouvert en lecture                    */
/******************************************************************************/
/* !!! FONCTION INUTILISEE !!! */
void Csc2symb(const CscMatrix *cscmtx, SymbolMatrix *symbmtx)
{
  INT itercblk;
  INT itercoltab;
  INT iterbloc;
  INT iterval;

#ifdef CSC_LOG
  fprintf(stdout, "-> Csc2symb \n");
#endif
  
  for (iterbloc=0; iterbloc < symbmtx->bloknbr; iterbloc++)
    symbmtx->bloktab[iterbloc].levfval=0;

  for (itercblk=0; itercblk < CSC_FNBR(cscmtx); itercblk++)
    {
      for (itercoltab=0;
	   itercoltab < CSC_COLNBR(cscmtx,itercblk);
	   itercoltab++)
	{
	  for (iterval = CSC_COL(cscmtx,itercblk,itercoltab);
	       iterval < CSC_COL(cscmtx,itercblk,itercoltab+1);
	       iterval++)
	    {
	      if (CSC_ROW(cscmtx,iterval) >=
		  symbmtx->cblktab[itercblk].fcolnum)
		{
		  iterbloc = symbmtx->cblktab[itercblk].bloknum;

		  while ((( symbmtx->bloktab[iterbloc].lrownum <
			    CSC_ROW(cscmtx,iterval)) ||
			  ( symbmtx->bloktab[iterbloc].frownum >
			    CSC_ROW(cscmtx,iterval))) &&
			 ( iterbloc < symbmtx->cblktab[itercblk+1].bloknum))
		    {
		      iterbloc++;
		    }

		  if ( iterbloc <
		       symbmtx->cblktab[itercblk+1].bloknum)
		    {
		      symbmtx->bloktab[iterbloc].levfval=1;
		    }
		  else printf("ILU: csc2symb drop coeff from CSC c=%ld(%ld) l=%ld(%ld) cblk=%ld fcol=%ld lcol=%ld\n",
			      (long)symbmtx->cblktab[itercblk].fcolnum+
			      (long)itercoltab,(long)itercoltab,
			      (long)CSC_ROW(cscmtx,iterval),(long)iterval,
			      (long)itercblk,
			      (long)symbmtx->cblktab[itercblk].fcolnum,
			      (long)symbmtx->cblktab[itercblk].lcolnum);
		}
	    }
	}
    }

#ifdef CSC_LOG
  fprintf(stdout, "<- Csc2symb \n");
#endif
}


/*== Divers ==*/
/******************************************************************************/
/* void CscTrans(CscMatrix *cscmtx, CscMatrix *csctrp)                        */
/*                                                                            */
/* Transpose une csc                                                          */
/*                                                                            */
/* cscmtx : csc                                                               */
/* csctrp : csc transposee                                                    */
/******************************************************************************/
/* !!! FONCTION INUTILISEE !!! */
void CscTrans(const CscMatrix *cscmtx, CscMatrix *csctrp)
{
  INT itercscf;
  INT itercol;
  INT valnbr;
  INT colcur;
  INT iterval;
  INT rowcur;
  INT iterval2;

#ifdef CSC_LOG
  fprintf(stdout, "-> CscTrans \n");
#endif

  /* Copie du coltab */
  CSC_FNBR(csctrp) = CSC_FNBR(cscmtx);
  CSC_FTAB(csctrp) = (CscFormat *) memAlloc(CSC_FNBR(cscmtx)*sizeof(CscFormat));
  if (CSC_FTAB(csctrp) == NULL)
    errorPrint( "CscTrans : Not enough memory for CSC_FTAB\n");

  for (itercscf = 0;
       itercscf < CSC_FNBR(csctrp);
       itercscf++)
    {
      CSC_COLNBR(csctrp,itercscf) = CSC_COLNBR(cscmtx,itercscf);
      CSC_COLTAB(csctrp,itercscf) =
	(INT *) memAlloc((CSC_COLNBR(csctrp,itercscf)+1)*sizeof(INT));
      if (CSC_COLTAB(csctrp,itercscf) == NULL)
	errorPrint( "CscTrans : Not enough memory for CSC_COLTAB\n");
      
      for (itercol = 0;
	   itercol < (CSC_COLNBR(csctrp,itercscf)+1);
	   itercol++)
	{
	  CSC_COL(csctrp,itercscf,itercol) = CSC_COL(cscmtx,itercscf,itercol);
	}
    }

  /* Allocation du valtab et rowtab */
  valnbr = CSC_VALNBR(cscmtx);
  printf("valnbr = %ld\n", (long)valnbr);
  CSC_ROWTAB(csctrp) = (INT *) memAlloc(valnbr*sizeof(INT));
  if (CSC_ROWTAB(csctrp) == NULL)
    errorPrint( "CscTrans : Not enough memory for CSC_ROWTAB\n");
  CSC_VALTAB(csctrp) = (FLOAT *) memAlloc(valnbr*sizeof(FLOAT));
  if (CSC_VALTAB(csctrp) == NULL)
    errorPrint( "CscTrans : Not enough memory for CSC_VALTAB\n");

  /* Renseignement des coeff */
  colcur = 0;
  for (itercscf = 0;
       itercscf < CSC_FNBR(cscmtx);
       itercscf++)
    {
      for (itercol = 0;
	   itercol < CSC_COLNBR(cscmtx,itercscf);
	   itercol++)
	{
	  for (iterval = CSC_COL(cscmtx,itercscf,itercol);
	       iterval < CSC_COL(cscmtx,itercscf,itercol+1);
	       iterval++)
	    {
	      INT cont = 0;
	      INT itercscf2 = 0;

	      rowcur = CSC_ROW(cscmtx,iterval);

	      /* Recherche dur coltab correspondant a rowcur dans trp */
	      do
		{
		  if (rowcur < CSC_COLNBR(csctrp,itercscf2))
		    {
		      iterval2 = CSC_COL(csctrp,itercscf2,rowcur);

		      CSC_ROW(csctrp,iterval2) = colcur;

		      CSC_VAL(csctrp,iterval2) = CSC_VAL(cscmtx,iterval);

		      (CSC_COL(csctrp,itercscf2,rowcur))++;

		      cont = 1;
		    }
		  else
		    {
		      rowcur -= CSC_COLNBR(csctrp,itercscf2);
		      itercscf2++;
		    }
		}
	      while (cont==0);
	    }
	  colcur++;
	}
    }

  rowcur = 0;
  /* Remet le coltab */
  for (itercscf = 0;
       itercscf < CSC_FNBR(csctrp);
       itercscf++)
    {
      for (itercol = 0;
	   itercol < (CSC_COLNBR(csctrp,itercscf));
	   itercol++)
	{
	  INT rowidx = CSC_COL(csctrp,itercscf,itercol);
	  CSC_COL(csctrp,itercscf,itercol) = rowcur;
	  rowcur = rowidx;
	}
    }

#ifdef CSC_LOG
  fprintf(stdout, "<- CscTrans \n");
#endif
}

/* !!! FONCTION INUTILISEE !!! */
void CscScaling2(char *Type, INT Ncol, INT *col, INT *row, FLOAT *val, FLOAT *rhs, FLOAT *rhs2)
{
  INT itercol;
  INT iterval;
  FLOAT *rowscal;
  FLOAT *colscal;
  const INT flag = (Type[1] != 'S');

#ifdef CSC_LOG
  fprintf(stdout, "-> CscScaling2 \n");
#endif

  if (flag)
    rowscal = (FLOAT *) memAlloc(Ncol*sizeof(FLOAT));
  else
    rowscal = NULL;
  colscal = (FLOAT *) memAlloc(Ncol*sizeof(FLOAT));

  if (flag)
    if (rowscal == NULL)
      errorPrint( "CscScaling2 : Not enough memory for rowscal\n");
  if (colscal == NULL)
    errorPrint( "CscScalign2 : Nto enough memory for colscal\n");
  
  for (itercol=0; itercol<Ncol; itercol++)
    {
      if (flag)
	rowscal[itercol] = 0;
      colscal[itercol] = 0;
    }

  for (itercol=0; itercol<Ncol; itercol++)
    {
      for (iterval=col[itercol]; iterval<col[itercol+1]; iterval++)
	{
	  colscal[itercol] += ABS_FLOAT(val[iterval-1]);
	  if (flag)
	    rowscal[row[iterval-1]-1] += ABS_FLOAT(val[iterval-1]);
	}
    }

  for (itercol=0; itercol<Ncol; itercol++)
    {
      for (iterval=col[itercol]; iterval<col[itercol+1]; iterval++)
	{
	  val[iterval-1] /= colscal[itercol];
	  if (flag)
	    val[iterval-1] /= rowscal[row[iterval-1]-1];
	  else
	    val[iterval-1] /= colscal[itercol];
	}
      if (rhs != NULL)
	{
	  if (flag)
	    rhs[itercol] /= rowscal[itercol];
	  else
	    rhs[itercol] /= colscal[itercol];
	}
      if (rhs2 != NULL)
	rhs2[itercol] *= colscal[itercol];
    }

  if (flag)
    memFree_null(rowscal);
  memFree_null(colscal);

#ifdef CSC_LOG
  fprintf(stdout, "<- CscScaling2 \n");
#endif
}

/******************************************************************************/
/* void CscScaling(CscMatrix *cscmtx, FLOAT *transcsc,                        */
/*                 FLOAT *rhs, FLOAT *rhs2)                                   */
/*                                                                            */
/* Scaling                                                                    */
/*                                                                            */
/* cscmtx : csc                                                               */
/* transcsc : transpose csc (NULL si on est en RSA)                           */
/* rhs : second membre                                                        */
/* rhs2 : solution                                                            */
/******************************************************************************/
/* !!! FONCTION INUTILISEE !!! */
void CscScaling(CscMatrix *cscmtx, FLOAT *transcsc, FLOAT *rhs, FLOAT *rhs2)
{
  const INT itercscf=0;
  INT itercol;
  INT iterval;
  FLOAT *rowscal;
  FLOAT *colscal;

#ifdef CSC_LOG
  fprintf(stdout, "-> CscScaling \n");
#endif
  
  rowscal = (FLOAT *) memAlloc(CSC_COLNBR(cscmtx,0)*sizeof(FLOAT));
  colscal = (FLOAT *) memAlloc(CSC_COLNBR(cscmtx,0)*sizeof(FLOAT));
  
  for (itercol=0; itercol<CSC_COLNBR(cscmtx,itercscf); itercol++)
    {
      rowscal[itercol] = 0;
      colscal[itercol] = 0;
    }

  /* Compute Scaling  */
  for (itercol=0; itercol<CSC_COLNBR(cscmtx,itercscf); itercol++)
    {
      for (iterval=CSC_COL(cscmtx,itercscf,itercol);
	   iterval<CSC_COL(cscmtx,itercscf,itercol+1);
	   iterval++)
	{
	  colscal[itercol] += ABS_FLOAT(CSC_VAL(cscmtx,iterval));
	  rowscal[CSC_ROW(cscmtx,iterval)] += ABS_FLOAT(CSC_VAL(cscmtx,iterval));
	}
    }

  /* Apply Scaling */
  for (itercol=0; itercol<CSC_COLNBR(cscmtx,itercscf); itercol++)
    {
      for (iterval=CSC_COL(cscmtx,itercscf,itercol);
	   iterval<CSC_COL(cscmtx,itercscf,itercol+1);
	   iterval++)
	{
	  CSC_VAL(cscmtx,iterval) /= colscal[itercol];
	  CSC_VAL(cscmtx,iterval) /= rowscal[CSC_ROW(cscmtx,iterval)];

	  if (transcsc != NULL)
	    {
	      transcsc[iterval] /= rowscal[itercol];
	      transcsc[iterval] /= colscal[CSC_ROW(cscmtx,iterval)];
	    }
	}
      if (rhs != NULL) rhs[itercol] /= rowscal[itercol];
      if (rhs2 != NULL) rhs2[itercol] *= colscal[itercol];
    }

  memFree_null(rowscal);
  memFree_null(colscal);

#ifdef CSC_LOG
  fprintf(stdout, "<- CscScaling \n");
#endif
}

/* !!! FONCTION INUTILISEE !!! */
void CscVerifUpdown(const UpDownVector *updovct, const SymbolMatrix *symbmtx,
		    const FLOAT *rhs2)
{
  INT itercblk;
  INT itercol;
  INT itersm2x;

#ifdef CSC_LOG
  fprintf(stdout, "-> CscVerifUpdown \n");
#endif

  for (itercblk=0; itercblk<symbmtx->cblknbr; itercblk++)
    {
      itersm2x = updovct->cblktab[itercblk].sm2xind;

      for (itercol=symbmtx->cblktab[itercblk].fcolnum;
	   itercol<symbmtx->cblktab[itercblk].lcolnum+1;
	   itercol++)
	{
#ifdef CPLX
	  errorPrint( "%ld (%10e,%10e) (%10e,%10e)\n",
		  (long)itercol, creal(rhs2[itercol]), cimag(rhs2[itercol]), creal(updovct->sm2xtab[itersm2x]), cimag(updovct->sm2xtab[itersm2x]));
#else
	  errorPrint( "%ld %10e %10e\n",
		  (long)itercol, rhs2[itercol], updovct->sm2xtab[itersm2x]);
#endif
	  itersm2x++;
	}
    }

#ifdef CSC_LOG
  fprintf(stdout, "<- CscVerifUpdown \n");
#endif
}

/* !!! FONCTION INUTILISEE !!! */
void CscUpdown(UpDownVector *updovct, /*const*/ SymbolMatrix *symbmtx,
	       const FLOAT *rhs)
{
  INT itercblk;
  INT itercol;
  INT itersm2x;

#ifdef CSC_LOG
  fprintf(stdout, "-> CscUpdown \n");
#endif
  
  for (itercblk=0; itercblk<symbmtx->cblknbr; itercblk++)
    {
      itersm2x = updovct->cblktab[itercblk].sm2xind;
      for (itercol=symbmtx->cblktab[itercblk].fcolnum;
	   itercol<symbmtx->cblktab[itercblk].lcolnum+1;
	   itercol++)
	{
	  updovct->sm2xtab[itersm2x] = rhs[itercol];
	  itersm2x++;
	}
    }

#ifdef CSC_LOG
  fprintf(stdout, "<- CscUpdown \n");
#endif
}

/* !!! FONCTION INUTILISEE !!! */
void CscUpdown2(UpDownVector *updovct, /*const*/ SymbolMatrix *symbmtx,
		const FLOAT *rhs)
{
  INT itercblk;
  INT itercol;
  INT itersm2x;

#ifdef CSC_LOG
  fprintf(stdout, "-> CscUpdown2 \n");
#endif

  for (itercblk=0; itercblk<symbmtx->cblknbr; itercblk++)
    {
      itersm2x = updovct->cblktab[itercblk].sm2xind;
      for (itercol=symbmtx->cblktab[itercblk].fcolnum;
	   itercol<symbmtx->cblktab[itercblk].lcolnum+1;
	   itercol++)
	{
	  updovct->sm2xtab[itersm2x] = 1.0;
	  itersm2x++;
	}
    }

#ifdef CSC_LOG
  fprintf(stdout, "<- CscUpdown2 \n");
#endif
}
/* !!! INUTILISEE, non testee !!! */
void Csc2updown_new(Sopalin_Data_t * sopalin_data, int me, const CscMatrix *cscmtx, UpDownVector *updovct,
		    /*const*/ SymbolMatrix *symbmtx, int n, MPI_Comm comm)
{
  SolverMatrix * datacode;
  FLOAT *tempy;
  INT    i;

#ifdef CSC_LOG
  fprintf(stdout, "-> Csc2updown_new \n");
#endif

  datacode = sopalin_data->datacode;

  MONOTHREAD_BEGIN;

  if (!(tempy = (FLOAT *)memAlloc(updovct->gnodenbr*sizeof(FLOAT)))) MALLOC_ERROR("tempy");
  sopalin_data->ptr_raff[0] = (void *)tempy;

  for (i=0;i<updovct->gnodenbr;i++)
    tempy[i] = (FLOAT)i;
  MONOTHREAD_END;

  SYNCHRO_THREAD;

  tempy = (FLOAT *)sopalin_data->ptr_raff[0];

  CscAx(sopalin_data, me, cscmtx, tempy, updovct->sm2xtab, symbmtx, updovct, comm);

  MONOTHREAD_BEGIN;
  memFree_null(tempy);
  MONOTHREAD_END;

#ifdef CSC_LOG
  fprintf(stdout, "<- Csc2updown_new \n");
#endif
}

/******************************************************************************/
/* void CscDiagDom(CscMatrix *cscmtx)                                         */
/*                                                                            */
/* Transforme la csc en csc a diagonale dominante                             */
/*                                                                            */
/* cscmtx : Csc                                                               */
/******************************************************************************/
/* !!! INUTILISEE, non testee !!! */
void CscDiagDom(CscMatrix *cscmtx)
{
  INT itercblk;
  INT itercoltab;
  INT iterval;
  INT itercol=0;
  INT diag=0;
  double sum41;

#ifdef CSC_LOG
  fprintf(stdout, "-> CscDiagDom \n");
#endif

  for (itercblk=0; itercblk<CSC_FNBR(cscmtx); itercblk++)
    {
      for (itercoltab=0;
	   itercoltab<CSC_COLNBR(cscmtx,itercblk);
	   itercoltab++)
	{
	  sum41 = 0;
	  for (iterval=CSC_COL(cscmtx,itercblk,itercoltab);
	       iterval<CSC_COL(cscmtx,itercblk,itercoltab+1);
	       iterval++)
	    {
	      if (CSC_ROW(cscmtx,iterval) != itercol)
		sum41 += ABS_FLOAT(CSC_VAL(cscmtx,iterval));
	      else
		diag = iterval;
	    }
	  CSC_VAL(cscmtx,diag) = sum41+1;
	  
	  itercol++;
	}
    }

#ifdef CSC_LOG
  fprintf(stdout, "<- CscDiagDom \n");
#endif
}

/* !!! INUTILISEE, non testee !!! */
INT CscStrucSym(CscMatrix *cscmtx)
{
  /* csc_fnbr = 1 */
  INT itercblk=0;
  INT itercoltab;
  INT iterval;
  INT iterval2;
  INT result=1;

#ifdef CSC_LOG
  fprintf(stdout, "-> CscStrucSym \n");
#endif

  for (itercoltab=0;
       itercoltab<CSC_COLNBR(cscmtx,itercblk);
       itercoltab++)
    {
      for (iterval=CSC_COL(cscmtx,itercblk,itercoltab);
	   iterval<CSC_COL(cscmtx,itercblk,itercoltab+1);
	   iterval++)
	{
	  if (CSC_ROW(cscmtx,iterval) != itercoltab)
	    {
	      INT rowidx=CSC_ROW(cscmtx, iterval);
	      INT flag=0;
	      
	      for(iterval2=CSC_COL(cscmtx,itercblk,rowidx);
		  iterval2<CSC_COL(cscmtx,itercblk,rowidx+1);
		  iterval2++)
		{
		  if (CSC_ROW(cscmtx,iterval2) == itercoltab)
		    {
		      flag=1;
		      break;
		    }
		}
	      result = result && flag;
	    }
	}
    }

#ifdef CSC_LOG
  fprintf(stdout, "<- CscStrucSym \n");
#endif
  
  return result;
}
