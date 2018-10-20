# include <iostream>
# include <math.h>
# include <stdio.h>
# include <stdlib.h>

using namespace std;

extern long long int find_nonzeros(long long int dim, double* matrix);
extern void sparse_matrix(long long int nx, long long int ny, double* matrix, long long int nz, double* nonzero, long long int* col_ind, long long int* row_ptr);
extern double get_matrix_element(long long int row, long long int column, double* vals, long long int* col_ind, long long int* row_ptr); 

int main(void){
    long long int nx = 4;
    long long int ny = 4;
    double* A = new double[nx*ny];
    A[0] = 1; A[1] = 0; A[2] = 0; A[3] = 4;
    A[4] = 0; A[5] = 3; A[6] = 5; A[7] = 6;
    A[8] = 0; A[9] = 1; A[10] = 5; A[11] = 0;
    A[12] = 2; A[13] = 3; A[14] = 0; A[15] = 7;

    long long int nz = find_nonzeros(nx*ny, A);
    double* val = new double[nz];
    long long int* col_ind = new long long int[nz];
    long long int* row_ptr = new long long int[nx+1];

    sparse_matrix(nx, ny, A, nz, val, col_ind, row_ptr);
    cout << "Number of rows: " << nx << "\n";
    cout << "Number of columns: " << ny << "\n";
    cout << "Number of non-zero elements: " << nz << "\n";
    cout << "Matrix in dense representation: \n";
    for(long long int i = 0; i < nx; i++){
	for(long long int j = 0; j < ny; j++){
	    cout << A[i*ny+j] << " ";
	}
	cout << "\n";
    }
    cout << "Non-zero elements array: ";
    for(long long int i = 0; i < nz; i++) cout << val[i] << " ";
    cout << "\n";
    cout << "Column indices: ";
    for(long long int i = 0; i < nz; i++) cout << col_ind[i] << " ";
    cout << "\n";
    cout << "Row pointer: ";
    for(long long int i = 0; i <= nx; i++) cout << row_ptr[i] << " ";
    cout << "\n";

    long long int elem;
    cout << "Which element do you want?\n";
    cin >> elem;
    if(elem < 0 || elem > nx*ny){
       	cerr << "Element not in the relephant range\n";
    }else{
	long long int row = elem/ny;
	long long int column = elem%ny;
	double value = get_matrix_element(row, column, val, col_ind, row_ptr);
	cout << value << "\n";
    }
}

