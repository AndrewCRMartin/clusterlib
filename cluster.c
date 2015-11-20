#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bioplib/macros.h"
#include "bioplib/general.h"

#include "cluster.h"

#define MAXBUFF 256


/* Prototypes */
int main(int argc, char **argv);
REAL **ReadData(FILE *fp, int *numVec, int *vecSize);
void Usage(void);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  int *type, BOOL *printMatrix, int *nClus, BOOL *verbose);

int main(int argc, char **argv)
{
   REAL **dataMatrix = NULL;
   int numVec;
   int vecSize;
   int **clusters = NULL;
   REAL *distances;
   FILE *in  = stdin,
        *out = stdout;
   char infile[MAXBUFF],
        outfile[MAXBUFF];
   int  type = CLUSTER_SINGLE,
        nClus = 0;
   BOOL printMatrix = TRUE,
        verbose = FALSE;

   if(ParseCmdLine(argc, argv, infile, outfile, &type, &printMatrix, &nClus, &verbose))
   {
      if(OpenStdFiles(infile, outfile, &in, &out))
      {
         dataMatrix = ReadData(in, &numVec, &vecSize);
         clusters=cluster(dataMatrix, numVec, vecSize, type, &distances, verbose);
         if(printMatrix)
         {
            PrintClusterMatrix(stdout, clusters, numVec, distances);
         }
         if(nClus)
         {
            PrintClusters(stdout, clusters, nClus, dataMatrix, numVec, vecSize, NULL);
         }
      }
   }
   else
   {
      Usage();
   }
   
   return(0);
}


REAL **ReadData(FILE *fp, int *numVec, int *vecSize)
{
   char *line, *ptr, *pLine;
   char buffer[MAXBUFF];
   int  nFields = 0,
        nVec = 0,
        nf = 0;
   REAL **data = NULL;

   while((line=fgetsany(fp))!=NULL)
   {
      TERMINATE(line);
      KILLLEADSPACES(pLine, line);
      if(pLine[0] == '#')
      {
         continue;
      }
      
      if(strlen(pLine))
      {
         nVec++;
         
         if(nFields==0)
         {
            ptr = pLine;
            do
            {
               ptr=GetWord(ptr, buffer, MAXBUFF);
               nFields++;
            }  while(ptr!=NULL);
         }

         if((data = (REAL **)realloc(data, nVec*sizeof(REAL *)))==NULL)
         {
            return(NULL);
         }
         if((data[nVec-1] = (REAL *)malloc(nFields * sizeof(REAL)))==NULL)
         {
            return(NULL);
         }
         
         ptr = pLine;
         nf = 0;
         do
         {
            ptr=GetWord(ptr, buffer, MAXBUFF);
            if(nf >= nFields)
            {
               break;
            }
            if(!sscanf(buffer, "%lf", &(data[nVec-1][nf])))
            {
               return(NULL);
            }
            nf++;
         }  while(ptr!=NULL);
      }
      
      free(line);
   }

   *numVec = nVec;
   *vecSize = nFields;
   
   return(data);
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
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  int *type, BOOL *printMatrix, int *nClus, BOOL *verbose)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = '\0';
   
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
         /* Check that there are only 1 or 2 arguments left             */
         if(argc > 2)
            return(FALSE);
         
         /* Copy the first to infile                                    */
         strcpy(infile, argv[0]);
         
         /* If there's another, copy it to outfile                      */
         argc--;
         argv++;
         if(argc)
            strcpy(outfile, argv[0]);
            
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}


void Usage(void)
{
   fprintf(stderr, "\ncluster V2.0 (c) 2012, UCL, Dr. Andrew C.R. Martin\n");
   fprintf(stderr, "\nUsage: cluster [-n num] [-m] [-t type] [-s | -c | -a | -d] [input.dat [output.txt]]\n");
   fprintf(stderr, "       -n     Specify a number of clusters and display the vectors associated\n");
   fprintf(stderr, "              with those clusters\n");
   fprintf(stderr, "       -m     Do NOT display the clustering matrix\n");
   fprintf(stderr, "       -t     Specify the clustering type:\n");
   fprintf(stderr, "                 single, complete, average, centroid\n");
   fprintf(stderr, "       -s     Do single linkage (same as -t single)\n");
   fprintf(stderr, "       -c     Do complete linkage (same as -t complete)\n");
   fprintf(stderr, "       -a     Do average linkage (same as -t average)\n");
   fprintf(stderr, "       -d     Do centroid linkage (same as -t centroid)\n");
   fprintf(stderr, "\n");
   fprintf(stderr, "Performs hierarchical cluster analysis according to the specified method.\n");
   fprintf(stderr, "input.dat is a simple file containing vectors to be clustered. Blank lines and\n");
   fprintf(stderr, "comments introduced with a # character are ignored. e.g.\n");
   fprintf(stderr, "\n");
   fprintf(stderr, "   # Sample data file\n");
   fprintf(stderr, "   1.3 5.7 2.8\n");
   fprintf(stderr, "   1.9 5.2 3.4\n");
   fprintf(stderr, "   3.7 2.5 8.7\n");
   fprintf(stderr, "   3.2 2.6 9.2\n");
   fprintf(stderr, "   5.3 1.0 0.1\n");
   fprintf(stderr, "\n");
   fprintf(stderr, "Output is a clustering matrix with the vector numbers arranged horizontally and\n");
   fprintf(stderr, "following rows showing the cluster numbers to which the vectors are assigned\n");
   fprintf(stderr, "for each number of clusters. This is suppressed with -m\n\n");
}

