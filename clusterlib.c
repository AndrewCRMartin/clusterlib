#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bioplib/array.h"
#include "bioplib/macros.h"

#include "cluster.h"

static REAL CalcDistance(int *members0, int nMembers0, int *members1, int nMembers1, REAL **dataMatrix, int vecSize, int type, REAL **dataDistances);
static REAL doClustering(int *newClusters, int *oldClusters, REAL **dataMatrix, int level, int numVec, int vecSize, int type, REAL **dataDistances);
/* static REAL VecDistance(int item0, int item1, REAL **dataMatrix, int vecSize); */
static REAL doSingleDistance(int *members0, int nMembers0, int *members1, int nMembers1, REAL **dataDistances);
static REAL doCompleteDistance(int *members0, int nMembers0, int *members1, int nMembers1, REAL **dataDistances);
static REAL doCentroidDistance(int *members0, int nMembers0, int *members1, int nMembers1, REAL **dataMatrix, int vecSize);
static REAL doAverageDistance(int *members0, int nMembers0, int *members1, int nMembers1, REAL **dataDistances);
static REAL Distance(REAL *vector0, REAL *vector1, int vecSize);
static REAL **FillDataDistanceMatrix(REAL **dataMatrix, int numVec, int vecSize, BOOL verbose);



/*
labels is an optional linked list of labels to be used when one asks
to print a certain number of clusters and their members. If not used
it is set to NULL.
 */
void PrintClusters(FILE *fp, int **clusters, int nClus, REAL **data, int numVec, int vecSize, CLUSTERLABELS *labels)
{
   int row, i, j, k, 
       clusNum, 
       *members = NULL;
   CLUSTERLABELS **labelIndex = NULL;

   row = numVec-nClus;

   if(labels!=NULL)
   {
      int nLabels = 0;
      CLUSTERLABELS *l;
      for(l=labels; l!=NULL; NEXT(l))
      {
         nLabels++;
      }
      if((labelIndex = (CLUSTERLABELS **)malloc(nLabels * sizeof(CLUSTERLABELS *)))==NULL)
      {
         labels = NULL;
      }
      
      for(l=labels, nLabels=0; l!=NULL; NEXT(l))
      {
         labelIndex[nLabels++] = l;
      }
   }
   

   if((members=(int *)malloc(numVec * sizeof(int)))!=NULL)
   {
      /* Initialize members array */
      for(i=0; i<numVec; i++)
      {
         members[i] = 0;
      }
      /* Count the number of items in each cluster */
      for(i=0; i<numVec; i++)
      {
         members[clusters[row][i]]++;
      }
      
      clusNum = 1;
      for(i=0; i<numVec; i++)
      {
         if(members[i])
         {
            fprintf(fp, "Cluster %d (%d members)\n", clusNum++, members[i]);
            for(j=0; j<numVec; j++)
            {
               if(clusters[row][j]==i)
               {
                  for(k=0; k<vecSize; k++)
                  {
                     fprintf(fp, "%g ", data[j][k]);
                  }
                  if(labels != NULL)
                  {
                     fprintf(fp, "# %s", labelIndex[j]->label);
                  }
                  
                  fprintf(fp, "\n");
               }
            }
         }
      }
      free(members);
   }
}

void PrintClusterMatrix(FILE *fp, int **clusters, int numVec, REAL *distances)
{
   int i, j;
   
   fprintf(fp,"Vector   :");
   for(i=0; i<numVec; i++)
   {
      fprintf(fp, "%4d", i);
   }
   fprintf(fp, "\n");
   
   for(i=0; i<numVec; i++)
   {
      fprintf(fp,"Clust%-4d:",numVec-i);
      for(j=0; j<numVec; j++)
      {
         fprintf(fp, "%4d", clusters[i][j]);
      }
      if((distances != NULL) && (distances[i] >= 0.0))
      {
         fprintf(fp, " (%.3f)\n", distances[i]);
      }
      else
      {
         fprintf(fp, "\n");
      }
   }
}

int **cluster(REAL **dataMatrix, int numVec, int vecSize, int type, REAL **distances, BOOL verbose)
{
   int i;
   int **clusters;
   int printStep;
   REAL **dataDistances = NULL;
   
   if(verbose)
   {
      printStep = ((int)log10(numVec))-2;
      printStep = pow(10.0, printStep);
      
      if(printStep < 1)
      {
         printStep = 1;
      }
      if(printStep > 100)
      {
         printStep = 100;
      }
   }
   

   /* Allocate space for storing the clustering array */
   if(verbose)
   {
      fprintf(stderr,"Initializing clustering matrix...");
   }
   
   if((clusters = (int **)malloc(numVec * sizeof(int *)))==NULL)
   {
      return(NULL);
   }
   for(i=0; i<numVec; i++)
   {
      int j;

      if((clusters[i] = (int *)malloc(numVec * sizeof(int)))==NULL)
      {
         return(NULL);
      }
      for(j=0; j<numVec; j++)
      {
         clusters[i][j] = 0;
      }
   }

   /* Initially each vector is in its own cluster */
   for(i=0; i<numVec; i++)
   {
      clusters[0][i] = i;
   }
   if(verbose)
   {
      fprintf(stderr,"done\n");
   }

   /* Array to store precalculated distances between data points
      used by single, complete and average linkage
   */
   if((type == CLUSTER_SINGLE) ||
      (type == CLUSTER_COMPLETE) ||
      (type == CLUSTER_AVERAGE))
   {
      if((dataDistances = FillDataDistanceMatrix(dataMatrix, numVec, vecSize, verbose))==NULL)
      {
         return(NULL);
      }
   }

   /* Array to store the distances at each clustering step */
   if(verbose)
   {
      fprintf(stderr,"Initializing distance matrix...");
   }
   if(distances != NULL)
   {
      if((*distances = (REAL *)malloc(numVec * sizeof(REAL)))==NULL)
      {
         return(NULL);
      }
      (*distances)[0] = -1.0;
   }
   if(verbose)
   {
      fprintf(stderr,"done\n");
   }
   
   /* Now do the clustering */
   if(verbose)
   {
      fprintf(stderr, "Clustering step 1/%d\n", numVec);
   }
   
   for(i=1; i<numVec; i++)
   {
      REAL dist;

      if(verbose)
      {
         if(!(i%printStep))
         {
            fprintf(stderr, "Clustering step %d/%d\n", i, numVec);
         }
      }
      
      dist = doClustering(clusters[i], clusters[i-1], dataMatrix, i, numVec, vecSize, type, dataDistances);
      if(distances != NULL)
      {
         (*distances)[i] = dist;
      }
   }

   if(dataDistances != NULL)
   {
      FreeArray2D((char **)dataDistances, numVec, numVec);
      dataDistances = NULL;
   }

   /* Return the cluster tree */
   return(clusters);
}

static REAL **FillDataDistanceMatrix(REAL **dataMatrix, int numVec, int vecSize, BOOL verbose)
{
   REAL **dataDistances = NULL;
   int i, j;
   
   if(verbose)
   {
      fprintf(stderr,"Precalculating data distance matrix...");
   }
   if((dataDistances = (REAL **)Array2D(sizeof(REAL), numVec, numVec))==NULL)
   {
      return(NULL);
   }
   for(i=0; i<numVec; i++)
   {
      for(j=i; j<numVec; j++)
      {
         if(i==j)
         {
            dataDistances[i][j] = (REAL)0.0;
         }
         else
         {
            dataDistances[i][j] = dataDistances[j][i] = Distance(dataMatrix[i], dataMatrix[j], vecSize);
         }
      }
   }
   
   if(verbose)
   {
      fprintf(stderr,"done\n");
   }
   return(dataDistances);
}



CLUSTER *BuildTree(int **clusters, int numVec, int vecSize)
{
   CLUSTER *cluster = NULL;
   return(cluster);
}


static REAL doClustering(int *newClusters, int *oldClusters, REAL **dataMatrix, int level, int numVec, int vecSize, int type, REAL **dataDistances)
{
   REAL minDist, dist;
   int i, j, iBest, jBest;
   int **members, *nMembers;
   
   /* Allocate this number of arrays to hold the members of each cluster */
   if((members = (int **)malloc(numVec * sizeof(int *)))==NULL)
   {
      return(-1);
   }
   for(i=0; i<numVec; i++)
   {
      if((members[i] = (int *)malloc(numVec * sizeof(int)))==NULL)
      {
         return(-1);
      }
   }
   if((nMembers = (int *)malloc(numVec * sizeof(int)))==NULL)
   {
      return(-1);
   }
   for(i=0; i<numVec; i++)
   {
      nMembers[i] = 0;
   }

   /* Populate the members arrays */
   for(i=0; i<numVec; i++)
   {
      int n;
      
      n = nMembers[oldClusters[i]];
      members[oldClusters[i]][n] = i;
      nMembers[oldClusters[i]]++;
   }
   
   /* For each pair of clusters see which is nearest */
   minDist = (-1);
   for(i=0; i<numVec-1; i++)
   {
      for(j=i+1; j<numVec; j++)
      {
         if(nMembers[i] && nMembers[j])
         {
            dist = CalcDistance(members[i], nMembers[i], members[j], nMembers[j], dataMatrix, vecSize, type, dataDistances);
            if((dist < minDist)||(minDist < 0))
            {
               minDist = dist;
               iBest = i;
               jBest = j;
            }
            
         }
      }
   }

   /* Free memory */
   for(i=0; i<numVec; i++)
   {
      free(members[i]);
   }
   free(members);
   free(nMembers);

   /* Now merge those two clusters */
   for(i=0; i<numVec; i++)
   {
      if(oldClusters[i] == jBest)
      {
         newClusters[i] = iBest;
      }
      else
      {
         newClusters[i] = oldClusters[i];
      }
   }
   
   return(minDist);
}



static REAL CalcDistance(int *members0, int nMembers0, int *members1, int nMembers1, REAL **dataMatrix, int vecSize, int type, REAL **dataDistances)
{
   /* For the moment just does single linkage */
   REAL dist;

   /* Single and complete linkage */
   switch(type)
   {
   case CLUSTER_SINGLE:
      dist = doSingleDistance(members0, nMembers0, members1, nMembers1, dataDistances);
      break;
   case CLUSTER_COMPLETE:
      dist = doCompleteDistance(members0, nMembers0, members1, nMembers1, dataDistances);
      break;
   case CLUSTER_CENTROID:
      dist = doCentroidDistance(members0, nMembers0, members1, nMembers1, dataMatrix, vecSize);
      break;
   case CLUSTER_AVERAGE:
      dist = doAverageDistance(members0, nMembers0, members1, nMembers1, dataDistances);
      break;
   default:
      CLUSTERDIE("(clusterlib) Error: Illegal clustering type: %d\n", type);
      break;
   }
   return(dist);
}

static REAL doSingleDistance(int *members0, int nMembers0, int *members1, int nMembers1, REAL **dataDistances)
{
   REAL bestDist, dist;
   int  i, j;

   bestDist = dataDistances[members0[0]][members1[0]];
  
   for(i=0; i<nMembers0; i++)
   {
      for(j=0; j<nMembers1; j++)
      {
         dist = dataDistances[members0[i]][members1[j]];
         if(dist < bestDist)
         {
            bestDist = dist;
         }
      }
   }
   
   return(bestDist);
}

static REAL doCompleteDistance(int *members0, int nMembers0, int *members1, int nMembers1, REAL **dataDistances)
{
   REAL bestDist, dist;
   int  i, j;

   bestDist = dataDistances[members0[0]][members1[0]];
   for(i=0; i<nMembers0; i++)
   {
      for(j=0; j<nMembers1; j++)
      {
         dist = dataDistances[members0[i]][members1[j]];
         if(dist > bestDist)
         {
            bestDist = dist;
         }
      }
   }
   
   return(bestDist);
}

static REAL doCentroidDistance(int *members0, int nMembers0, int *members1, int nMembers1, REAL **dataMatrix, int vecSize)
{
   REAL dist;
   int  dim, vec;

   REAL *centroid0 = NULL,
        *centroid1 = NULL;

   if((centroid0 = (REAL *)malloc(vecSize * sizeof(REAL)))==NULL)
   {
      return(-1);
   }
   if((centroid1 = (REAL *)malloc(vecSize * sizeof(REAL)))==NULL)
   {
      return(-1);
   }

   for(dim=0; dim<vecSize; dim++)
   {
      centroid0[dim] = (REAL)0.0;
      for(vec=0; vec<nMembers0; vec++)
      {
         centroid0[dim] += dataMatrix[members0[vec]][dim];
      }
      centroid0[dim] /= nMembers0;
   }

   for(dim=0; dim<vecSize; dim++)
   {
      centroid1[dim] = (REAL)0.0;
      for(vec=0; vec<nMembers1; vec++)
      {
         centroid1[dim] += dataMatrix[members1[vec]][dim];
      }
      centroid1[dim] /= nMembers1;
   }

   dist = Distance(centroid0, centroid1, vecSize);
   
   free(centroid0);
   free(centroid1);

   return(dist);
}


static REAL doAverageDistance(int *members0, int nMembers0, int *members1, int nMembers1, REAL **dataDistances)
{
   REAL dist = 0;
   int  divisor = 0,
      i, j;

   for(i=0; i<nMembers0; i++)
   {
      for(j=0; j<nMembers1; j++)
      {
         divisor++;
         dist += dataDistances[members0[i]][members1[j]];
      }
   }
   dist /= divisor;
   return(dist);
}



/*
 * static REAL VecDistance(int item0, int item1, REAL **dataMatrix, int vecSize)
 * {
 *    return(Distance(dataMatrix[item0], dataMatrix[item1], vecSize));
 * }
 */

static REAL Distance(REAL *vector0, REAL *vector1, int vecSize)
{
   int i;
   REAL distSq = (REAL)0.0;
   
   for(i=0; i<vecSize; i++)
   {
      distSq += (vector0[i] - vector1[i]) * (vector0[i] - vector1[i]);
   }
   return(sqrt(distSq));
}

/* Ward's:

D = ||X_K - X_L||^2 / (1/n_k + 1/n_l)

Sum of square deviations from the centroid

For the kth cluster, define the Error Sum of Squares as 
ESS_k = sum of squared deviations from the cluster centroid
If there are C clusters, define the Total ESS as
ESS_{tot} = \sum_{k=1}^C ESS_k

Combine the 2 cluster that minimize the increase in ESS_{tot}

The division converts this to variances

Using a Lance-Williams recursive algorith, (as in the 'cluster' program)
At each stage, the distances between the clusters are calculated by
testing a new cluster k composed of two sub-clusters, i and j as
{[D_ik * (N_i + N_k)] + [D_jk * (N_j + N_k)] - [D_ij * N_k]} / (N_i + N_j + N_k)

*/

