#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"

#define MAXCLUSLABEL 256

typedef struct _cluster
{
   int level;
   /* These are used by the terminal cluster */
   REAL *data;
   int  id;
   /* These are used by parent clusters      */
   struct _cluster **children;
}  CLUSTER;

typedef struct _clusterlabels
{
   char                  label[MAXCLUSLABEL];
   struct _clusterlabels *next;
}  CLUSTERLABELS;



int **cluster(REAL **dataMatrix, int numVec, int vecSize, int type, REAL **distances, BOOL verbose);
void PrintClusterMatrix(FILE *fp, int **clusters, int numVec, REAL *distances);
CLUSTER *BuildTree(int **clusters, int numVec, int vecSize);
void PrintClusters(FILE *fp, int **clusters, int nClus, REAL **data, int numVec, int vecSize, CLUSTERLABELS *labels);

#define CLUSTER_SINGLE   0
#define CLUSTER_COMPLETE 1
#define CLUSTER_AVERAGE  2
#define CLUSTER_CENTROID 3

#define CLUSTERDIE(x, y) fprintf(stderr, x, y); exit(1)
