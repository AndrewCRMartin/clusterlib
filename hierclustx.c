// I N C L U D E S ///////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

// D E F I N E S /////////////////////////////////////////////////////////////

#define MAXVAL (double)1.0e+38f

#define MIN(a, b) ( a < b ? a : b )
#define MAX(a, b) ( a > b ? a : b )

// T Y P E S /////////////////////////////////////////////////////////////////

// S T R U C T U R E S ///////////////////////////////////////////////////////

typedef struct
{
	double **m;
	int   num_rows;
	int   num_cols;
} fmatrix_t;

typedef struct
{
	int   **m;
	int   num_rows;
	int   num_cols;
} imatrix_t;

// G L O B A L S //////////////////////////////////////////////////////////////


// F U N C T I O N S //////////////////////////////////////////////////////////

void FMatrixCreate(fmatrix_t *mat, int num_rows, int num_cols)
{
	int row;

	mat->num_rows = num_rows;
	mat->num_cols = num_cols;

	mat->m = (double**)malloc(mat->num_rows * sizeof(double*));

	for (row = 0; row < mat->num_rows; row++)
	{
		mat->m[row] = (double*)malloc(mat->num_cols * sizeof(double));
		memset(mat->m[row], 0, mat->num_cols * sizeof(double));
	}
} // end FMatrixCreate

///////////////////////////////////////////////////////////////////////////////

void FMatrixDestroy(fmatrix_t *mat)
{
	int row ;

	for (row = 0; row < mat->num_rows; row++)
	{
		free(mat->m[row]);
	}
	free(mat->m);
} // end FMatrixDestroy

///////////////////////////////////////////////////////////////////////////////

void IMatrixCreate(imatrix_t *mat, int num_rows, int num_cols)
{
	int row;

	mat->num_rows = num_rows;
	mat->num_cols = num_cols;

	mat->m = (int**)malloc(mat->num_rows * sizeof(int*));

	for (row = 0; row < mat->num_rows; row++)
	{
		mat->m[row] = (int*)malloc(mat->num_cols * sizeof(int));
		memset(mat->m[row], 0, mat->num_cols * sizeof(int));
	}
} // end IMatrixCreate

///////////////////////////////////////////////////////////////////////////////

void IMatrixDestroy(imatrix_t *mat)
{
	int row ;

	for (row = 0; row < mat->num_rows; row++)
	{
		free(mat->m[row]);
	}
	free(mat->m);
} // end IMatrixDestroy

///////////////////////////////////////////////////////////////////////////////

void Dissim(fmatrix_t *diss, const double *data, const double *wts, int num_rows, int num_cols)
{
	int i1, i2, j;
	double tmp, diff_a;

	// create a nxn matrix
	FMatrixCreate(diss, num_rows, num_rows);

	for (i1 = 0; i1 < diss->num_rows; i1++)
	{	
		for (i2 = 0; i2 < i1; i2++)
		{
			for (j = 0; j < num_cols; j++)
			{
				// We use the squared Euclidean distance, weighted.
				diff_a = data[i1 * num_cols + j] - data[i2 * num_cols + j];

				tmp = (wts[i1] * wts[i2]) / (wts[i1] + wts[i2]) * (double)(diff_a * diff_a);

				diss->m[i1][i2] += tmp;
			} // end for j
			diss->m[i2][i1] = diss->m[i1][i2];
		} // end for i2
	} // end for i1

} // end Dissim

///////////////////////////////////////////////////////////////////////////////

void GetNNs(int *nn, double *nndiss, const fmatrix_t *diss, const int *flag)
{
	int minobs, i1, i2;
	double mindist;

	for (i1 = 0; i1 < diss->num_rows; i1++)
	{
		if (flag[i1] == 1)
		{
			minobs = -1;
			mindist = MAXVAL;
			for (i2 = 0; i2 < diss->num_cols; i2++)
			{
				if ((diss->m[i1][i2] < mindist) && (i1 != i2))
				{
					mindist = diss->m[i1][i2];
					minobs  = i2;
				} // end if
			} // end for i2
			nn[i1] = minobs;
			nndiss[i1] = mindist;
		} // end if flag == 1
	} // end for i1

} // // end GetNNs

///////////////////////////////////////////////////////////////////////////////

void HierClust(int *order, int *ia, int *ib, double *levels, 
			   const double *data, const double *wts, 
			   const int *nrows_ptr, const int *ncols_ptr)
{
	int num_rows, num_cols;

	fmatrix_t diss;

	imatrix_t clusmat;

	double *mass = NULL,
		  *nndiss = NULL;

	int *flag = NULL,
	    *a = NULL,
		*b = NULL,
	    *card = NULL,
		*tmp_buffer = NULL,
		*nn = NULL;

	int i, i2, ncl, 
		lastind, 
		order_list_size, 
		minobs, 
		clus1, clus2, 
		left, right, 
		tobeexp,
		tmp;

	double mindis, tmp1, tmp2, tmp3;
	
	num_rows = (*nrows_ptr);
	num_cols = (*ncols_ptr);

	memset(ia, 0, sizeof(int) * (num_rows-1));
	memset(ib, 0, sizeof(int) * (num_rows-1));
	memset(levels, 0, sizeof(double) * (num_rows-1));
	memset(order, 0, sizeof(int) * num_rows);

	mass = (double*)malloc(sizeof(double) * num_rows);
	memcpy(mass, wts, sizeof(double)* num_rows);


	Dissim(&diss, data, mass, num_rows, num_cols);
    

	a = (int*)malloc(sizeof(int) * (num_rows-1));
	b = (int*)malloc(sizeof(int) * (num_rows-1));
	memset(a, 0, sizeof(int) * (num_rows-1));
	memset(b, 0, sizeof(int) * (num_rows-1));

	card  = (int*)malloc(sizeof(int) * num_rows);
	flag  = (int*)malloc(sizeof(int) * num_rows);
	
	for (i = 0; i < num_rows; i++)
	{
		flag[i] = 1;
		card[i] = 1;
	}

	nn     = (int*)malloc(sizeof(int) * num_rows);
	nndiss = (double*)malloc(sizeof(double) * num_rows);
	memset(nn, 0, sizeof(int) * num_rows);
	memset(nndiss, 0, sizeof(double) * num_rows);

	GetNNs(nn, nndiss, &diss, flag);


	// create a nxn matrix
	IMatrixCreate(&clusmat, num_rows, num_rows);
	for (i = 0; i < num_rows; i++) 
	{ 
		clusmat.m[i][num_rows-1] = i; 
	}

	
	for (ncl = (num_rows - 2); ncl >= 0; ncl--)
	{
		// check for  agglomerable pair
		minobs = -1;
		mindis = MAXVAL;
		for (i = 0; i < num_rows; i++)
		{
			if (flag[i] == 1)
			{
				if (nndiss[i] < mindis)
				{
					mindis = nndiss[i];
					minobs = i;
				} // end if
			} // end if flag == 1
		} // end for i

		// find aglomerands clus1 and clus2 with former < latter
		if (minobs < nn[minobs])
		{
			clus1 = minobs;
			clus2 = nn[minobs];
		} // end if
		else
		if (minobs > nn[minobs])
		{
			clus2 = minobs;
			clus1 = nn[minobs];
		} // end if

		// So, agglomeration of pair clus1 < clus2 defines cluster ncl

		//------------------------------------ Block for subnode labels
		a[ncl] = clus1; // aine, or left child
		b[ncl] = clus2; // benjamin, or right


		// Now build up ia, ib as version of a, b which is R-complement
		if (card[clus1] == 1) ia[ncl] = (-clus1); // singleton
		if (card[clus2] == 1) ib[ncl] = (-clus2); // singleton

		if (card[clus1] > 1 )
		{
			// left child is non singleton
			lastind = 0;

			for (i2 = (num_rows-2); i2 >= (ncl+1); i2--)
			{
				if (a[i2] == clus1) lastind = i2; // Only concerns a[i2]
			}
			ia[ncl] = num_rows - lastind; // label of non-singleton
		}

		if (card[clus2] > 1)
		{
			// right child is non singleton
			lastind  = 0;

			for (i2 = (num_rows-2); i2 >= (ncl+1); i2--)
			{
				if (a[i2] == clus2) lastind = i2; // Only concerns a[i2]
			} // end for i2
			ib[ncl] = num_rows - lastind; // label of non-singleton
		}	

		if (ia[ncl] > 0 || ib[ncl] > 0)
		{
			// Check that left < right
			left  = MIN(ia[ncl], ib[ncl]);
			right = MAX(ia[ncl], ib[ncl]);

			ia[ncl] = left;
			ib[ncl] = right;
		} // end if


	    //----------------------------------------------------------------

		levels[ncl] = mindis;

		for (i = 0; i < num_rows; i++)
		{
			clusmat.m[i][ncl] = clusmat.m[i][ncl+1];

			if (clusmat.m[i][ncl] == clus2) clusmat.m[i][ncl] = clus1;
		} // end for i

		// Next we need to update diss array

		for (i = 0; i < num_rows; i++)
		{
			if ( (i != clus1) && (i != clus2) && (flag[i] == 1) )
			{
				tmp1 = ((mass[clus1] + mass[i]) / (mass[clus1] + mass[clus2] + mass[i])) * diss.m[clus1][i];
				tmp2 = ((mass[clus2] + mass[i]) / (mass[clus1] + mass[clus2] + mass[i])) * diss.m[clus2][i]; 
				tmp3 = (mass[i] / (mass[clus1] + mass[clus2] + mass[i])) * diss.m[clus1][clus2];

				diss.m[clus1][i] = tmp1 + tmp2 - tmp3;
				diss.m[i][clus1] = diss.m[clus1][i];
				
			} // end if
		} // end for i

		mass[clus1] += mass[clus2];  // Update mass of new cluster
		card[clus1] += card[clus2];  // Update card of new cluster
		

		// Cluster label clus2 is knocked out; following not nec. but no harm

		flag[clus2]   = 0;
		nndiss[clus2] = MAXVAL;
		mass[clus2]   = 0.0f;

		for (i=0; i<num_rows; i++)
		{
			diss.m[clus2][i] = MAXVAL;
			diss.m[i][clus2] = diss.m[clus2][i];
		} // end for i

		// finaly update nn and nndiss
		// i.e. nearest neighbors and the nearest neigh. dissomilarity
		GetNNs(nn, nndiss, &diss, flag);
		
	} // end for ncl


	for(i = 0; i < (num_rows-1); i++)
	{
		ia[i] -= 1;
		ib[i] -= 1;
	}

	// reverse ia, ib
	for(i = 0; i < (num_rows-1) / 2; i++)
	{

		tmp = ia[num_rows-2-i];
		ia[num_rows-2-i] = ia[i];
		ia[i] = tmp;

		tmp = ib[num_rows-2-i];
		ib[num_rows-2-i] = ib[i];
		ib[i] = tmp;

	}

	// merge is R-compliant; later suppress merge2

	//------------------------------------- Build R-compatiblr order from ia, ib

	tmp_buffer = (int*)malloc(sizeof(int) * num_rows);

	order[0] = ia[num_rows - 2];
	order[1] = ib[num_rows - 2];

	order_list_size = 2;

	for (i = 0; i < (num_rows - 2); i++) // for precisely n-2 further node expansions
	{
		for (i2 = 0; i2 < order_list_size; i2++)
		{
			if (order[i2] > 0)
			{
				tobeexp = (order[i2]) - 1;
				
				if (i2 == 0)
				{
					tmp_buffer[0] = ia[tobeexp];
					tmp_buffer[1] = ib[tobeexp];

					memcpy(&tmp_buffer[2], &order[1], sizeof(int)*(order_list_size-1));
	
				} // end if
				else
				if (i2 == (order_list_size-1))
				{
					memcpy(tmp_buffer, order, sizeof(int)*(order_list_size - 1));

					tmp_buffer[order_list_size - 1] = ia[tobeexp];
					tmp_buffer[order_list_size]     = ib[tobeexp];
					
				} // end  if
				else
				if (i2 > 0 && i2 < (order_list_size-1))
				{
					memcpy(tmp_buffer, order, sizeof(int)*(i2-1));

					tmp_buffer[i2]   = ia[tobeexp];
					tmp_buffer[i2+1] = ib[tobeexp];

					memcpy(&tmp_buffer[i2+2], &order[i2+1], sizeof(int)*(order_list_size - (i2 + 1)));
						
				} // end if
		
				order_list_size++;
				memcpy(order, tmp_buffer, sizeof(int)*order_list_size);

			} // end if order[i2] > 0
		} // end for i2
	} // end for i


	// invert order sign
	for (i=0; i<num_rows; i++)
	{
		order[i] = -order[i];
	}


	// reverse levels
	for (i=0; i<(num_rows - 1)/2; i++)
	{
		tmp1 = levels[num_rows-2-i];
		levels[num_rows-2-i] = levels[i];
		levels[i] = tmp1;
	}

	free(tmp_buffer);
	free(mass);
	free(flag);
	free(a);
	free(b);
	free(card);
	free(nn);
	free(nndiss);
	
	IMatrixDestroy(&clusmat);
	FMatrixDestroy(&diss);

} // end HierClust
