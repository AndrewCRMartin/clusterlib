#include <stdio.h>
#include <stdlib.h>

#include "bioplib/TypeDefs.h"
#include "bioplib/MathType.h"

typedef struct _clusterdata
{
   APTR *data;
}  CLUSTERDATA;

typedef struct _cluster
{
   CLUSTERDATA *data;
   CLUSTER *clusters;
   REAL *distances;
   int level;
   int id;
}  CLUSTER;

#define MAXDATA 10


int main(int argc, char **argv)
{
   REAL **dMatrix = {
      {1,5,6,3,5,3,6,8,9,5},
      {5,5,6,7,5,3,4,5,6,8},
      {3,5,6,7,5,4,5,3,5,7},
      {1,3,4,5,6,7,8,6,5,6},
      {4,5,6,5,6,7,4,5,6,4},
      {3,4,5,6,7,4,5,6,4,3},
      {8,7,8,7,6,7,8,6,4,5},
      {5,6,4,3,4,5,7,8,9,2},
      {3,4,3,5,4,5,6,7,6,5},
      {3,4,2,3,2,4,5,6,7,8}
   }
   int matSize = 10;
   
   CLUSTER *tree = NULL;
   
   tree=cluster(dMatrix, matSize);
}

/* This should take the data not the distance matrix! */
CLUSTER *cluster(REAL **dMatrix, int matSize)
{
   int i;
   
   CLUSTER **leaves = NULL;

   /* Allocate array of leaf node pointers */
   if((leaves = (CLUSTER **)malloc(matSize * sizeof(CLUSTER *)))==NULL)
   {
      return(NULL);
   }
   
   /* Create those nodes and populate with distances to the other nodes */
   for(i=0; i<matSize; i++)
   {
      CLUSTER *leaf;
      if((leaf = (CLUSTER *)malloc(sizeof(CLUSTER)))==NULL)
      {
         return(NULL);
      }
      leaf->id = i;
      leaf->level = 0;
      leaf->distances = (REAL *)malloc(matSize * sizeof(REAL));
      for(j=0;j<matSize;j++)
      {
         /* Really we should be calculating this from the data not assuming a pre-created distance matrix */
         leaf->distances[j] = dMatrix[i][j];
      }
   }
   
}
