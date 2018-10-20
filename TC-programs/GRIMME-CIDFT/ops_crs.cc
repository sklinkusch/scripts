#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

using namespace std;

uint32_t  find_nonzeros(uint32_t dim, double* matrix);
void sparse_matrix(uint32_t nx, uint32_t ny, double* matrix, uint32_t nz, double* nonzero, uint32_t* col_ind, uint32_t* row_ptr);
double get_matrix_element(uint32_t row, uint32_t column, double* nonzero, uint32_t* col_ind, uint32_t* row_ptr);

uint32_t find_nonzeros(uint32_t dim, double* matrix){
   uint32_t x = 0;
   for(uint32_t i = 0; i < dim; i++){
       if(matrix[i] != 0) x++;
   }
   return(x);
}


void sparse_matrix(uint32_t nx, uint32_t ny, double* matrix, uint32_t nz, double* nonzero, uint32_t* col_ind, uint32_t* row_ptr){
  nz = 0;
  for(uint32_t i = 0; i < nx; i++){
     for(uint32_t j = 0; j < ny; j++){
      if(matrix[i*ny+j] != 0.) nz++;
     }
  }
  uint32_t z = 0;
  uint32_t x = 0;
  for(uint32_t i = 0; i < nx; i++){
     row_ptr[x] = z;
     x++;
     for(uint32_t j = 0; j < ny; j++){
       if(matrix[i*ny+j] != 0.){
	   nonzero[z] = matrix[i*ny+j];
	   col_ind[z] = j;
	   z++;
       }
     }  
  }
  row_ptr[x] = z;
}

double get_matrix_element(uint32_t row, uint32_t column, double* nonzero, uint32_t* col_ind, uint32_t* row_ptr){
    uint32_t curr_row = row_ptr[row];
    uint32_t next_row = row_ptr[row+1];
    double value = 0.;
    for(uint32_t i = curr_row; i < next_row; i++){
	if(col_ind[i] == column){
	    value = nonzero[i];
	    break;
	}
    }
    return(value);
}

