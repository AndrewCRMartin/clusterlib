/*************************************************************************

   Program:    strucclus
   File:       strucclus.c
   
   Version:    V0.1
   Date:       12.10.12
   Function:   
   
   Copyright:  (c) UCL / Dr. Andrew C. R. Martin 2012
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Institute of Structural & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   EMail:      andrew.martin@ucl.ac.uk
               andrew@bioinf.org.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! 

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bioplib/macros.h"
#include "bioplib/general.h"
#include "bioplib/pdb.h"
#include "bioplib/array.h"

#include "cluster.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 256

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
REAL **ReadPDBData(FILE *fp, CLUSTERLABELS *labels, int nres, int *numVec, int *vecSize);
void Usage(void);
BOOL ParseCmdLine(int argc, char **argv, char *resfile, char *infile, 
                  char *outfile, int *type, BOOL *printMatrix, 
                  int *nClus, BOOL *verbose, BOOL *doFix);
CLUSTERLABELS *ReadLabels(FILE *fp);
BOOL WriteFixedResList(FILE *in, FILE *out, CLUSTERLABELS *labels);

/************************************************************************/
int main(int argc, char **argv)
{
   REAL **dataMatrix = NULL;
   int numVec;
   int vecSize;
   int **clusters = NULL;
   REAL *distances;
   FILE *in  = stdin,
        *res = NULL,
        *out = stdout;
   char infile[MAXBUFF],
        outfile[MAXBUFF],
        resfile[MAXBUFF];
   int  type = CLUSTER_SINGLE,
        nClus = 0;
   BOOL printMatrix = TRUE,
        verbose = FALSE,
        doFix = FALSE;
   CLUSTERLABELS *labels = NULL;

   if(ParseCmdLine(argc, argv, resfile, infile, outfile, &type, &printMatrix, &nClus, &verbose, &doFix))
   {
      if((res=fopen(resfile, "r")))
      {
         if((labels = ReadLabels(res)))
         {
            if(blOpenStdFiles(infile, outfile, &in, &out))
            {
               CLUSTERLABELS *r;
               int nres;

               if(doFix)
               {
                  WriteFixedResList(in, out, labels);
                  return(0);
               }
               
               /* Count the number of residues in our labels - that is maximum for 
                  numVec 
               */
               for(r=labels, nres=0; r!=NULL; NEXT(r))
               {
                  nres++;
               }

               if((dataMatrix = ReadPDBData(in, labels, nres, &numVec, &vecSize))==NULL)
               {
                  return(1);
               }
               
               clusters=cluster(dataMatrix, numVec, vecSize, type, 
                                &distances, verbose);
               if(printMatrix)
               {
                  PrintClusterMatrix(stdout, clusters, numVec, distances);
               }
               if(nClus)
               {
                  PrintClusters(stdout, clusters, nClus, dataMatrix, 
                                numVec, vecSize, labels);
               }
            }
         }
         else
         {
            fprintf(stderr,"Error (strucclus): Unable to read residue \
list from file (%s)\n",
                    resfile);
            return(1);
         }
      }
      else
      {
         fprintf(stderr,"Error (strucclus): Unable to open residue \
list file (%s)\n",
                 resfile);
         return(1);
      }
   }
   else
   {
      Usage();
   }
   
   return(0);
}


/************************************************************************/
/*>
   ---------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *infile      Input file (or blank string)
            char   *outfile     Output file (or blank string)
   Returns: BOOL                Success?

   Parse the command line
   
   08.11.96 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *resfile, char *infile, 
                  char *outfile, int *type, BOOL *printMatrix, 
                  int *nClus, BOOL *verbose, BOOL *doFix)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = resfile[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'm':
            *printMatrix = FALSE;
            break;
         case 'v':
            *verbose = TRUE;
            break;
         case 't':
            argc--;
            argv++;
            if(argc)
            {
               if(!strcmp(argv[0], "single"))
               {
                  *type = CLUSTER_SINGLE;
               }
               else if(!strcmp(argv[0], "complete"))
               {
                  *type = CLUSTER_COMPLETE;
               }
               else if(!strcmp(argv[0], "average"))
               {
                  *type = CLUSTER_AVERAGE;
               }
               else if(!strcmp(argv[0], "centroid"))
               {
                  *type = CLUSTER_CENTROID;
               }
            }
            else
            {
               return(FALSE);
            }
            break;
         case 's':
            *type = CLUSTER_SINGLE;
            break;
         case 'a':
            *type = CLUSTER_AVERAGE;
            break;
         case 'c':
            *type = CLUSTER_COMPLETE;
            break;
         case 'd':
            *type = CLUSTER_CENTROID;
            break;
         case 'f':
            *doFix = TRUE;
            break;
         case 'n':
            argc--;
            argv++;
            if(argc)
            {
               if(!sscanf(argv[0], "%d", nClus))
               {
                  return(FALSE);
               }
            }
            else
            {
               return(FALSE);
            }
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are only 1, 2, or 3 arguments left         */
         if(argc > 3)
            return(FALSE);

         /* First argument is the resfile                               */
         strcpy(resfile, argv[0]);

         /* If there's another, copy it to infile                       */
         argc--;
         argv++;
         if(argc)
         {
            strcpy(infile, argv[0]);
         
            /* If there's another, copy it to outfile                   */
            argc--;
            argv++;
            if(argc)
               strcpy(outfile, argv[0]);
         }
            
         return(TRUE);
      }
      argc--;
      argv++;
   }

   /* Check the resfile was specified                                   */
   if(!(resfile[0]))
      return(FALSE);
   
   return(TRUE);
}


void Usage(void)
{
   fprintf(stderr, "\nstrucclus V1.0 (c) 2012, UCL, Dr. Andrew C.R. \
Martin\n");
   fprintf(stderr, "\nUsage: strucclus [-n num] [-m] [-t type] \
[-s | -c | -a | -d] [-f] reslist.dat [input.pdb [output.txt]]\n");
   fprintf(stderr, "       -n     Specify a number of clusters and \
display the vectors associated\n");
   fprintf(stderr, "              with those clusters\n");
   fprintf(stderr, "       -m     Do NOT display the clustering \
matrix\n");
   fprintf(stderr, "       -t     Specify the clustering type:\n");
   fprintf(stderr, "                 single, complete, average, centroid\n");
   fprintf(stderr, "       -s     Do single linkage (same as -t \
single)\n");
   fprintf(stderr, "       -c     Do complete linkage (same as -t \
complete)\n");
   fprintf(stderr, "       -a     Do average linkage (same as -t \
average)\n");
   fprintf(stderr, "       -d     Do centroid linkage (same as -t \
centroid)\n");
   fprintf(stderr, "       -f     Fix the reslist file - generates a \
copy of the file with\n");
   fprintf(stderr, "              residues missing from the PDB file \
removed\n");
   fprintf(stderr, "\n");

   fprintf(stderr, "Performs hierarchical cluster analysis according to \
the specified method\n");
   fprintf(stderr, "to identify clusters of residues in a PDB file. \n");
   fprintf(stderr, "reslist.dat contains a list of residue identifiers \
(in the \n");
   fprintf(stderr, "form [c]nnn[i]) which are to be clustered.\n");
   fprintf(stderr, "\n");
   fprintf(stderr, "Output is a clustering matrix with the vector \
numbers arranged horizontally and\n");
   fprintf(stderr, "following rows showing the cluster numbers to which \
the vectors are assigned\n");
   fprintf(stderr, "for each number of clusters. This is suppressed \
with -m\n\n");
}

CLUSTERLABELS *ReadLabels(FILE *fp)
{
   char buffer[MAXBUFF];
   CLUSTERLABELS *labels = NULL,
           *r;
   
   while(fgets(buffer, MAXBUFF, fp))
   {
      char *b;
      
      TERMINATE(buffer);
      KILLLEADSPACES(b, buffer);
      KILLTRAILSPACES(b);
      
      if(strlen(b) && (*b != '#'))
      {
         if(labels == NULL)
         {
            INIT(labels, CLUSTERLABELS);
            r = labels;
         }
         else
         {
            ALLOCNEXT(r, CLUSTERLABELS);
         }
         
         if(r==NULL)
         {
            return(NULL);
         }
         
         strncpy(r->label, b, MAXBUFF);
      }
   }

   return(labels);
}

REAL **ReadPDBData(FILE *fp, CLUSTERLABELS *labels, int nres, int *numVec, int *vecSize)
{
   PDB *pdb = NULL;
   int natoms;
   REAL **data = NULL;
   BOOL missing;
   
   /* vecSize is 3 for the three dimensions                             */
   *vecSize = 3;
   
   /* Allocate the data array                                           */
   if((data = (REAL **)blArray2D(sizeof(REAL), nres, *vecSize))==NULL)
   {
      return(NULL);
   }
   
   /* open the PDB file                                                 */
   if((pdb=blReadPDB(fp, &natoms))!=NULL)
   {
      PDB *p;
      CLUSTERLABELS *r;
      BOOL found;
      char chain[8], insert[8];
      int resnum;

      nres = 0;
      missing = FALSE;
      
      /* We MUST search through the labels first so that the data items
         end up in the correct order for the labels to be valid
      */
      for(r=labels; r!=NULL; NEXT(r))
      {
         found = FALSE;
         for(p=pdb; p!=NULL; NEXT(p))
         {
            if(!strncmp(p->atnam, "CA  ", 4))
            {
               if(blParseResSpec(r->label, chain, &resnum, insert))
               {
                  if((chain[0]  == p->chain[0]) &&
                     (resnum    == p->resnum) &&
                     (insert[0] == p->insert[0]))
                  {
                     data[nres][0] = p->x;
                     data[nres][1] = p->y;
                     data[nres][2] = p->z;
                     nres++;
                     found = TRUE;
                     break;
                  }
               }
            }
         }
         if(!found)
         {
            fprintf(stderr,"Info (strucclus) Residue %c%d%c not found in PDB file\n", 
                    chain[0], resnum, insert[0]);
            missing = TRUE;
         }
      }
      *numVec = nres;
   }
   else
   {
      return(NULL);
   }

   if(missing)
   {
      fprintf(stderr,"\nError (strucclus) Remove missing residues from the residue file\n");
      return(NULL);
   }
   
   return(data);
}

BOOL WriteFixedResList(FILE *in, FILE *out, CLUSTERLABELS *labels)
{
   PDB *pdb = NULL;
   int natoms;
   
   /* read the PDB file                                                 */
   if((pdb=blReadPDB(in, &natoms))!=NULL)
   {
      PDB *p;
      CLUSTERLABELS *r;
      char chain[8], insert[8];
      int resnum;

      for(r=labels; r!=NULL; NEXT(r))
      {
         for(p=pdb; p!=NULL; NEXT(p))
         {
            if(!strncmp(p->atnam, "CA  ", 4))
            {
               if(blParseResSpec(r->label, chain, &resnum, insert))
               {
                  if((chain[0]  == p->chain[0]) &&
                     (resnum    == p->resnum) &&
                     (insert[0] == p->insert[0]))
                  {
                     fprintf(out, "%s\n", r->label);
                     break;
                  }
               }
            }
         }
      }
   }
   else
   {
      return(FALSE);
   }

   return(TRUE);
}

