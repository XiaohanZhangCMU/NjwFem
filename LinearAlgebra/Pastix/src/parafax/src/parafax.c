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
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <pthread.h>

#ifdef FORCE_NOMPI
#include "nompi.h"
#else
#include <mpi.h>
#endif /* FORCE_NOMPI */

#include "ptscotch.h"
#include "common_pastix.h"
#include "cscd.h"
#include "parafax.h"

/*
  Function: pfax_vect_build

  Parameters:
    com_ind     -
    compo       -
    cscd_dgraph -
    collect     -
  
*/
int 
pfax_vect_build(int com_ind, int * compo, CSCD * cscd_dgraph, INT ** collect) 
{
  INT * tmp; 
  INT * vect;
  int   sizetmp = CSCD_IA(cscd_dgraph)[compo[com_ind]] - CSCD_IA(cscd_dgraph)[compo[com_ind+1]];
  INT   i;
  INT   j       = 2;
  INT   total   = 1;
  /* le cas critique etant un ensemble d'arete disjointe, ce qui donne des segments de longueur 1 */
  
  tmp  = (INT*)memAlloc(sizeof(INT) * sizetmp);	      
  vect = (INT*)memAlloc(sizeof(INT) * (2 * sizetmp +1) ); 

  /* copie du segment d'aretes concernées */
  memCpy(tmp,&(CSCD_JA(cscd_dgraph)[CSCD_IA(cscd_dgraph)[compo[com_ind]]-1]),sizetmp);

  /* on reordonne ce segment */
  cscd_quicksort(tmp, NULL, 0, sizetmp);
  
  vect[1] = tmp[0];
  vect[2] = tmp[0];

  for (i=0; i<sizetmp; i++) 
    {
      if (tmp[i] == vect[j]) 
	i++;
      else 
	{
	  if (tmp[i] == vect[j]+1) 
	    {
	      vect[j] = tmp[i];
	      i++;
	    }
	  else 
	    {
	      j++; vect[j] = tmp[i];
	      j++; vect[j] = tmp[i];
	      i++; total++;
	    }
	}
    }
  vect[0] = total;
  collect[com_ind] = vect;

  return 0;

}


/** Merge de deux vecteurs de facto
 ** Le resultat est envoye dans v1
 **/
int
pfax_vect_merge(INT * v1, INT * v2)
{
  INT *tmp;
  INT  size1 = v1[0];
  INT  size2 = v2[0];
  int  i = 0, j = 0, k = 0;
  
  tmp = (INT*)memAlloc(sizeof(INT) * (size1 + size2 +1));

  if ( v1[2*i+1] <= v2[2*j+1] ) 
    {
      tmp[2*k+1] = v1[2*i+1];
      tmp[2*k+2] = v1[2*i+2];
      i++;
    }
  else 
    {
      tmp[2*k+1] = v2[2*j+1];
      tmp[2*k+2] = v2[2*j+2];
      j++;
    }
      
  for (; i<size1 && j<size2;) 
    {
      if ( v1[2*i+1] <= v2[2*j+1] ) 
	{
	  if ( v1[2*i+1] <= tmp[2*k+2] +1 ) 
	    {
	      tmp[2*k+2] = v1[2*i+2]; i++;
	    }
	  else 
	    {
	      k++; 
	      tmp[2*k+1] = v1[2*i+1]; 
	      tmp[2*k+2] = v1[2*i+2];
	      i++;
	    }
	}
      else 
	{
	  if ( v2[2*j+1] <= tmp[2*k+2]+1 ) 
	    {
	      tmp[2*k+2] = v2[2*j+2];
	      j++;
	    }
	  else 
	    {
	      k++;
	      tmp[2*k+1] = v2[2*j+1]; 
	      tmp[2*k+2] = v2[2*j+2];
	      j++;
	    }
	}
    }

  if (i == size1) memCpy(tmp,&(v2[2*j+1]),(size2-j)*2);
  if (j == size2) memCpy(tmp,&(v1[2*i+1]),(size1-i)*2);
  
  return 0;
}



/** Fonction construisant un vecteur de booleen
 ** Ce vecteur dit si la composante est dans l'arbre des separateurs
 **/

INT
build_inTree(INT cblk, INT * Tree, INT * inTree) 
{
  int i=0;
  inTree[i] = 0;
  for (i=1; i<cblk; i++) 
    { 
      inTree[i] = 1;
      inTree[Tree[i]] = 0;
    }
  return 0;
}


/** Construit deux vecteurs donnant les limites de chaque composantes
 **
 **/

INT
build_compo(INT cblk, INT * beg_compo, INT * end_compo, INT * supnode, INT * inTree) 
{
  int i, count = 0;
  for (i=0; i<cblk; i++) 
    {
      if (inTree[i] == 1) 
	{
	  beg_compo[i] = count + 1;
	  count       += supnode[i];
	  end_compo[i] = count;
	}
      else 
	{
	  beg_compo[i] = -1;
	  end_compo[i] = -1; 
	}
    }
  return 0;
}


INT
compo(INT cblk, INT ind_vert, INT * beg_compo, INT * end_compo) 
{
  int i;
  for (i=0; i<cblk; i++)
    if (ind_vert >= beg_compo[i] && ind_vert <= end_compo[i])
      return i;
  
  return EXIT_FAILURE;
}

/* XL : commenté car compile pas et non utilisée */
/* INT */
/* a_une_contrib(INT *** contrib, INT ind_contributeur, INT ind_composante)  */
/* { */
  
/*   INT * v1 = contrib[ind_contributeur][ind_composante]; */
/*   INT i; */
  
/*   for(i=0; i<v1[0]; i++) */
/*     if (compo(v1[2*i+1]) == ind_composante) /\* XL : j'ai changé 2*&+1 en 2*i+1 au hazard ;p *\/ */
/*       return 1; */
/*   return 0; */
/* } */


INT
build_treetab(INT cblk, INT * Treetab, INT * Tree) 
{
  INT i, j;
  Treetab = memAlloc(sizeof(INT) * cblk);
  

  for (i=cblk-1; i>=0; i--) 
    {
      if (Treetab[i] < 0) 
	{
	  for (j=i-1; j>=0; j--)
	    if (Tree[j] == Tree[i])
	      Treetab[j] = i;
	  Treetab[i] = Tree[Tree[i]];
	}
    }
  return 0;
}

int
pfax_receive_merge(INT * tree, INT ** collect, int cblk, int my_comp, MPI_Comm mpi_comm) 
{
  INT *v1 = NULL;
  INT *v2 = NULL;
  int  i  = 0;
  
  while (v1 == NULL) 
    {
      if (tree[i] == my_comp) /* le premier fils rencontre est forcement local */
	v1 = collect[i];
      i++;
    }

  for (; i<cblk; i++) 
    {
      if (tree[i] == my_comp) 
	{
	  if (collect[i] == NULL) 
	    {
	      pfax_receive_vect(v2, i, mpi_comm);
	      pfax_vect_merge(v1, v2);
	    }
	}
    }
  
  return 0;
}


int 
pfax_receive_vect(INT * vect, INT comp_ind, MPI_Comm mpi_comm) 
{
  MPI_Status status;  
  INT        size_vect;
  int        target;

  MPI_Recv(&size_vect, 1, MPI_INT, MPI_ANY_SOURCE, comp_ind, mpi_comm, &status);

  vect   = memAlloc(sizeof(INT) * size_vect);
  target = status.MPI_SOURCE;

  MPI_Recv(vect, size_vect, MPI_INT, target, comp_ind, mpi_comm, &status);

  return 0;
}


int
pfax_vect_send(INT * vect, int dest, int comp_ind, MPI_Comm mpi_comm)
{
  int size = 2*vect[0]+1;

  MPI_Send(&size, 1,   MPI_INT, dest, comp_ind, mpi_comm);
  MPI_Send(vect, size, MPI_INT, dest, comp_ind, mpi_comm);

  return 0;
}


int
pfax_isLeaf_build(INT *isLeaf, INT *tree, INT cblk)
{
  int i;

  isLeaf = (INT*)memAlloc(sizeof(INT) * cblk);
  
  for (i=0; i<cblk; i++) 
    {
      isLeaf[i] = 1;
      isLeaf[tree[i]] = 0;
    }

  return 0;
}


int
pfax_vect_compo_build(INT *beg_compo, INT *end_compo, INT *isLeaf, INT *supnode, INT cblk)
{
  
  int i, count = 0;

  beg_compo = (INT*)memAlloc(sizeof(INT) * cblk);
  end_compo = (INT*)memAlloc(sizeof(INT) * cblk);

  for (i=0; i<cblk; i++) 
    {
      if (isLeaf[i] == 1) 
	{
	  beg_compo[i] = count + 1;
	  count       += supnode[i];
	  end_compo[i] = count;
	}
      else 
	{
	  beg_compo[i] = -1;
	  end_compo[i] = -1;
	}
    }

  return 0;
}

int
parafax(CSCD * cscd_dgraph, int myrank, INT * owner, INT cblk, INT * tree) 
{
  INT **collec = NULL;
  int  *compo  = NULL;
  INT  *child  = NULL;
  int   i;
  int   j;

  collec = (INT**)memAlloc(sizeof(INT*)*cblk);
  child  = (INT*) memAlloc(sizeof(INT) *cblk);

  for (i=1; i<cblk; i++)
    child[tree[i]]++;
  
  for(i=cblk-1; i>0; i++) 
    {
      if (owner[i] == myrank) 
	{
	  j = i;
	  if (child[i] == 0)
	    pfax_vect_build(i, compo, cscd_dgraph, collec);
	  if (child[i] >0)
	    pfax_receive_merge(tree, collec, cblk, j, cscd_dgraph->mpi_comm);
	}
    }
  
  /* tant que je suis le premier fils */
  while (tree[j] == j-1) 
    {
      j--;
      pfax_receive_merge(tree, collec, cblk, j, cscd_dgraph->mpi_comm);
      owner[j] = myrank;
    }

  if (tree[j] != -1)
    pfax_vect_send(collec[j], owner[tree[j]], j, cscd_dgraph->mpi_comm);
  else 
    {
      /* construire la symbol matrix */
    }
  
  return 0;
}


