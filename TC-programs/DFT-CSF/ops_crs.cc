#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

using namespace std;

long int  find_nonzeros(int dim, int umo, int omo, int llim);
void rows_cols_crs(int cis_size, long int* col_ind, long int* row_ptr);
void sparse_matrix(int nx, int ny, double* matrix, long int nz, double* nonzero, long int* col_ind, long int* row_ptr);
double get_matrix_element(int row, int column, double* nonzero, long int* col_ind, long int* row_ptr);

long int find_nonzeros(int dim, int umo, int omo, int llim){
   long int x = 0;
   int i1 = 0;
   int i2 = 0;
   int f1 = 0;
   int f2 = 0;
   for(int i = 0; i < dim; i++){
       if(i != 0){
	   i1 = (i-1)/umo+llim;
	   f1 = (i-1)%umo+omo+llim;
       }else{
	   i1 = 0;
	   f1 = 0;
       }
       for(int j = 0; j < dim; j++){
	   if(j != 0){
	       i2 = (j-1)/umo+llim;
	       f2 = (j-1)%umo+omo+llim;
	   }else{
	       i2 = 0;
	       f2 = 0;
	   }
	   if(i == 0 && j == 0){
	       x++;
	   }else if(i == 0 && j != 0){
	       x++;
	   }else if(i != 0 && j == 0){
	       x++;
	   }else if(i == j){
	       x++;
	   }else if(i1 == i2 && f1 != f2){
	       x++;
	   }else if(i1 != i2 && f1 == f2){
	       x++;
	   }
       }
   }
   return(x);
}

void rows_cols_crs(int cis_size, int umo, int omo, int llim, long int* col_ind, long int* row_ptr){
    int i1 = 0;
    int i2 = 0;
    int f1 = 0;
    int f2 = 0;
    long int z = 0;
    int x = 0;
    for(int i = 0; i < cis_size; i++){
	row_ptr[x] = z;
	x++;
	if(i != 0){
	i1 = (i-1)/umo+llim;
	f1 = (i-1)%umo+omo+llim;
	}
	for(int j = 0; j < cis_size; j++){
	    if(j != 0){
	    i2 = (j-1)/umo+llim;
	    f2 = (j-1)%umo+omo+llim;
	    }
	    if(i == 0 && j == 0){
		col_ind[z] = (long int) j;
		z++;
	    }else if(i == 0 && j != 0){
		col_ind[z] = (long int) j;
		z++;
	    }else if(i != 0 && j == 0){
		col_ind[z] = (long int) j;
		z++;
	    }else if(i == j){
		col_ind[z] = (long int) j;
		z++;
	    }else if(i1 == i2 && f1 != f2){
		col_ind[z] = (long int) j;
		z++;
	    }else if(i1 != i2 && f1 == f2){
		col_ind[z] = (long int) j;
		z++;
	    }
	}
    }
    row_ptr[x] = z;
}


void sparse_matrix(int nx, int ny, double* matrix, long int nz, double* nonzero, long int* col_ind, long int* row_ptr){
  nz = 0;
  for(int i = 0; i < nx; i++){
     for(int j = 0; j < ny; j++){
      if(matrix[i*ny+j] != 0.) nz++;
     }
  }
  long int z = 0;
  int x = 0;
  for(int i = 0; i < nx; i++){
     row_ptr[x] = z;
     x++;
     for(int j = 0; j < ny; j++){
       if(matrix[i*ny+j] != 0.){
	   nonzero[z] = matrix[i*ny+j];
	   col_ind[z] = (long int) j;
	   z++;
       }
     }  
  }
  row_ptr[x] = z;
}

double get_matrix_element(int row, int column, double* nonzero, long int* col_ind, long int* row_ptr){
    int curr_row = (int) row_ptr[row];
    int next_row = (int) row_ptr[row+1];
    double value = 0.;
    for(int i = curr_row; i < next_row; i++){
	if(col_ind[i] == (long int) column){
	    value = nonzero[i];
	    break;
	}
    }
    return(value);
}

