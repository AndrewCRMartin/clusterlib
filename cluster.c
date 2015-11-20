#include <stdio.h>
#include <stdlib.h>

#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"
#include "bioplib/macros.h"

typedef struct _cluster
{
   int level;
   /* These are used by the terminal cluster */
   REAL *data;
   int  id;
   /* These are used by parent clusters      */
   struct _cluster **children;
}  CLUSTER;
   typedef struct _node
   {
      CLUSTER *leaf;
      struct _node *next;
   }  NODE;
   


/* Prototypes */
int main(int argc, char **argv);
REAL **FillMatrix(int *numVec, int *vecSize);
REAL CalcDistance(CLUSTER *member0, CLUSTER *member1, int vecSize);
CLUSTER *cluster(REAL **dataMatrix, int numVec, int vecSize);
CLUSTER **DoClustering(CLUSTER **members, int level, int nClusters, int vecSize);
CLUSTER **PopulateLeaves(REAL **dataMatrix, int numVec, int vecSize);
CLUSTER *BuildTree(CLUSTER ***leaves, int numVec, int vecSize);
NODE *BuildNodeList(CLUSTER *member);
REAL VecDistance(NODE *set0, NODE *set1, int vecSize);



int main(int argc, char **argv)
{
   REAL **dataMatrix = NULL;
   int numVec;
   int vecSize;
   CLUSTER *tree = NULL;
   
   /* Dummy matrix population */
   dataMatrix = FillMatrix(&numVec, &vecSize);
   
   tree=cluster(dataMatrix, numVec, vecSize);

   return(0);
   
}

REAL VecDistance(NODE *set0, NODE *set1, int vecSize)
{
   return(0.0);
}

CLUSTER *BuildTree(CLUSTER ***leaves, int numVec, int vecSize)
{
   CLUSTER *cluster = NULL;
   return(cluster);
}

NODE *BuildNodeList(CLUSTER *member)
{
   NODE *node = NULL;
   return(node);
}


REAL **FillMatrix(int *numVec, int *vecSize)
{
   REAL **dataMatrix;
   int i;
   *numVec = 12;
   *vecSize = 10;
   
   dataMatrix = (REAL **)malloc(*numVec * sizeof(REAL *));
   for(i=0; i<*numVec; i++)
   {
      dataMatrix[i] = (REAL *)malloc(*vecSize * sizeof(REAL));
   }
   dataMatrix[0][0] = 1;
   dataMatrix[0][1] = 5;
   dataMatrix[0][2] = 6;
   dataMatrix[0][3] = 3;
   dataMatrix[0][4] = 5;
   dataMatrix[0][5] = 3;
   dataMatrix[0][6] = 6;
   dataMatrix[0][7] = 8;
   dataMatrix[0][8] = 9;
   dataMatrix[0][9] = 5;

   dataMatrix[1][0] = 5;
   dataMatrix[1][1] = 6;
   dataMatrix[1][2] = 7;
   dataMatrix[1][3] = 5;
   dataMatrix[1][4] = 3;
   dataMatrix[1][5] = 4;
   dataMatrix[1][6] = 5;
   dataMatrix[1][7] = 6;
   dataMatrix[1][8] = 8;
   dataMatrix[1][9] = 3;

   dataMatrix[2][0] = 5;
   dataMatrix[2][1] = 6;
   dataMatrix[2][2] = 7;
   dataMatrix[2][3] = 5;
   dataMatrix[2][4] = 4;
   dataMatrix[2][5] = 5;
   dataMatrix[2][6] = 3;
   dataMatrix[2][7] = 5;
   dataMatrix[2][8] = 7;
   dataMatrix[2][9] = 1;

   dataMatrix[3][0] = 3;
   dataMatrix[3][1] = 4;
   dataMatrix[3][2] = 5;
   dataMatrix[3][3] = 6;
   dataMatrix[3][4] = 7;
   dataMatrix[3][5] = 8;
   dataMatrix[3][6] = 6;
   dataMatrix[3][7] = 5;
   dataMatrix[3][8] = 6;
   dataMatrix[3][9] = 4;

   dataMatrix[4][0] = 5;
   dataMatrix[4][1] = 6;
   dataMatrix[4][2] = 5;
   dataMatrix[4][3] = 6;
   dataMatrix[4][4] = 7;
   dataMatrix[4][5] = 4;
   dataMatrix[4][6] = 5;
   dataMatrix[4][7] = 6;
   dataMatrix[4][8] = 4;
   dataMatrix[4][9] = 3;

   dataMatrix[5][0] = 4;
   dataMatrix[5][1] = 5;
   dataMatrix[5][2] = 6;
   dataMatrix[5][3] = 7;
   dataMatrix[5][4] = 4;
   dataMatrix[5][5] = 5;
   dataMatrix[5][6] = 6;
   dataMatrix[5][7] = 4;
   dataMatrix[5][8] = 3;
   dataMatrix[5][9] = 8;

   dataMatrix[6][0] = 7;
   dataMatrix[6][1] = 8;
   dataMatrix[6][2] = 7;
   dataMatrix[6][3] = 6;
   dataMatrix[6][4] = 7;
   dataMatrix[6][5] = 8;
   dataMatrix[6][6] = 6;
   dataMatrix[6][7] = 4;
   dataMatrix[6][8] = 5;
   dataMatrix[6][9] = 5;

   dataMatrix[7][0] = 6;
   dataMatrix[7][1] = 4;
   dataMatrix[7][2] = 3;
   dataMatrix[7][3] = 4;
   dataMatrix[7][4] = 5;
   dataMatrix[7][5] = 7;
   dataMatrix[7][6] = 8;
   dataMatrix[7][7] = 9;
   dataMatrix[7][8] = 2;
   dataMatrix[7][9] = 3;

   dataMatrix[8][0] = 4;
   dataMatrix[8][1] = 3;
   dataMatrix[8][2] = 5;
   dataMatrix[8][3] = 4;
   dataMatrix[8][4] = 5;
   dataMatrix[8][5] = 6;
   dataMatrix[8][6] = 7;
   dataMatrix[8][7] = 6;
   dataMatrix[8][8] = 5;
   dataMatrix[8][9] = 3;

   dataMatrix[9][0] = 5;
   dataMatrix[9][1] = 7;
   dataMatrix[9][2] = 4;
   dataMatrix[9][3] = 6;
   dataMatrix[9][4] = 7;
   dataMatrix[9][5] = 4;
   dataMatrix[9][6] = 6;
   dataMatrix[9][7] = 3;
   dataMatrix[9][8] = 2;
   dataMatrix[9][9] = 3;

   dataMatrix[10][0] = 4;
   dataMatrix[10][1] = 6;
   dataMatrix[10][2] = 7;
   dataMatrix[10][3] = 8;
   dataMatrix[10][4] = 5;
   dataMatrix[10][5] = 6;
   dataMatrix[10][6] = 8;
   dataMatrix[10][7] = 5;
   dataMatrix[10][8] = 4;
   dataMatrix[10][9] = 3;

   dataMatrix[11][0] = 4;
   dataMatrix[11][1] = 2;
   dataMatrix[11][2] = 3;
   dataMatrix[11][3] = 2;
   dataMatrix[11][4] = 4;
   dataMatrix[11][5] = 5;
   dataMatrix[11][6] = 6;
   dataMatrix[11][7] = 7;
   dataMatrix[11][8] = 8;
   dataMatrix[11][9] = 3;
   
   return(dataMatrix);
   
}


REAL CalcDistance(CLUSTER *member0, CLUSTER *member1, int vecSize)
{
   /* For the moment just does single linkage */
   NODE *set0 = NULL,
        *set1 = NULL,
        *n, *m;
   REAL dist, bestDist;

   if((set0 = BuildNodeList(member0))==NULL)
   {
      return(-1);
   }
   if((set1 = BuildNodeList(member1))==NULL)
   {
      return(-1);
   }

   bestDist = VecDistance(set0, set1, vecSize);

   for(m=set0; m!=NULL; NEXT(m))
   {
      for(n=set1; n!=NULL; NEXT(n))
      {
         dist = VecDistance(m, n, vecSize);
         if(dist < bestDist)
         {
            bestDist = dist;
         }
      }
   }
   
   FREELIST(set0, NODE);
   FREELIST(set1, NODE);
   
   return(bestDist);
}




CLUSTER *cluster(REAL **dataMatrix, int numVec, int vecSize)
{
   int i;
   
   CLUSTER ***leaves = NULL;

   /* Allocate array of leaf node pointers */
   if((leaves[0] = PopulateLeaves(dataMatrix, numVec, vecSize))==NULL)
   {
      return(NULL);
   }
   for(i=1; i<numVec; i++)
   {
      int nClusters = numVec - (i-1);
      leaves[i] = DoClustering(leaves[i-1], i, nClusters, vecSize);
   }
   return(BuildTree(leaves, numVec, vecSize));
}

CLUSTER **DoClustering(CLUSTER **members, int level, int nClusters, int vecSize)
{
   REAL minDist, dist;
   int i, j, iBest, jBest;
   CLUSTER **newClusters = NULL,
      *mergedCluster = NULL;
   
   /* Find the closest pair */
   minDist = CalcDistance(members[0], members[1], vecSize);
   iBest = 0;
   jBest = 1;
   
   for(i=0; i<nClusters; i++)
   {
      for(j=0; j<nClusters; j++)
      {
         if(i!=j)
         {
            dist = CalcDistance(members[i], members[j], vecSize);
            if(dist < minDist)
            {
               minDist = dist;
               iBest = i;
               jBest = j;
            }
         }
      }
   }

   /* Create the new set of clusters */
   if((newClusters = (CLUSTER **)malloc((nClusters-1)*sizeof(CLUSTER *)))==NULL)
   {
      return(NULL);
   }
   
   j=0;
   for(i=0;i<nClusters; i++)
   {
      if((i!=iBest)&&(i!=jBest))
      {
         newClusters[j++] = members[i];
      }
   }
   if((mergedCluster = (CLUSTER *)malloc(sizeof(CLUSTER)))==NULL)
   {
      return(NULL);
   }
   mergedCluster->id = (-1);
   mergedCluster->level = level;
   mergedCluster->data = NULL;
   if((mergedCluster->children = (CLUSTER **)malloc(2*sizeof(CLUSTER *)))==NULL)
   {
      return(NULL);
   }
   mergedCluster->children[0] = members[iBest];
   mergedCluster->children[1] = members[jBest];
   
   newClusters[j] = mergedCluster;
   
   return(newClusters);
}



CLUSTER **PopulateLeaves(REAL **dataMatrix, int numVec, int vecSize)
{
   CLUSTER **leaves = NULL;
   int i, j;
   
   if((leaves = (CLUSTER **)malloc(numVec * sizeof(CLUSTER *)))==NULL)
   {
      return(NULL);
   }
   
   /* Create those nodes and populate with data */
   for(i=0; i<numVec; i++)
   {
      if((leaves[i] = (CLUSTER *)malloc(sizeof(CLUSTER)))==NULL)
      {
         return(NULL);
      }
      leaves[i]->id = i;
      leaves[i]->level = numVec;
      leaves[i]->children = NULL;
      if((leaves[i]->data = (REAL *)malloc(vecSize * sizeof(REAL)))==NULL)
      {
         return(NULL);
      }
      
      for(j=0;j<vecSize;j++)
      {
         leaves[i]->data[j] = dataMatrix[i][j];
      }
   }
   
   return(leaves);
}
