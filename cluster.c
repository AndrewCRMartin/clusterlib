typedef (void *) DATA;
typedef struct _cluster
{
   CLUSTER children[2];
   DATA    data;
   REAL    *dists;
   int     level;
} CLUSTER;
