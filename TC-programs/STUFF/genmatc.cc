#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

void   diag_mat(int nroao, double* mat, double* vals);
extern "C" void dsyev_(char* jobz, char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* info);

int main(void){
    int nbasis = 100;
    int nstates = 100;
    double* configurations = new double[nbasis];
    double* energy = new double[nbasis];
    configurations[0] = 0.;
    for(int i = 1; i < nbasis; i++){
	configurations[i] = 0.1+((double) i * 0.001);
    }
    double* ionization = new double[nbasis];
    for(int i = 0; i < 4; i++) ionization[i] = 0.;
    for(int i = 4; i < nbasis; i++) ionization[i] = 1.e-5;

    double* eigenvectors = new double[nbasis*nbasis];
    for(int i = 0; i < nbasis; i++){
	for(int j = 0; j < nbasis; j++){
	    eigenvectors[i*nbasis+j] = 0.;
	}
    }
    eigenvectors[0] = configurations[0];
    eigenvectors[1*nbasis+1] = configurations[1];
    int j;
    for(int i = 2; i < nbasis; i++){
	j = i-1;
	eigenvectors[i*nbasis+j] = 0.01;
	eigenvectors[i*nbasis+i] = configurations[i];
	eigenvectors[j*nbasis+i] = 0.01;
    }
    diag_mat(nbasis,eigenvectors,energy);
//    dsyev_('V','U',nbasis,eigenvectors,nbasis,energy,work,lwork,info);
    clog << "Hier\n";

    double* dx = new double[nbasis*nbasis];
    double* dy = new double[nbasis*nbasis];
    double* dz = new double[nbasis*nbasis];
    for(int i = 0; i < nbasis; i++){
	for(int j = 0; j < nbasis; j++){
	    dx[i*nbasis+j] = 0.;
	    dy[i*nbasis+j] = 0.;
	    dz[i*nbasis+j] = 0.;
	}
    }
    for(int i = 0; i < nbasis; i++){
	int ipo = i+1;
	int imo = i-1;
	if(ipo < nbasis){
	    dx[i*nbasis+ipo] = 0.1;
	    dy[i*nbasis+ipo] = 0.1;
	    dz[i*nbasis+ipo] = 0.1;
	}
	if(imo >= 0){
	    dx[i*nbasis+imo] = 0.1;
	    dy[i*nbasis+imo] = 0.1;
	    dz[i*nbasis+imo] = 0.1;
	}
    }
/*	    if(fabs(i-j) == 1){
		dx[i*nbasis+j] = 0.1;
		dy[i*nbasis+j] = 0.1;
		dz[i*nbasis+j] = 0.1;
	    }else{
	        dx[i*nbasis+j] = 0.;
	        dy[i*nbasis+j] = 0.;
	        dz[i*nbasis+j] = 0.;
	    }
	}
    }*/
    double dum = 1.e-6*(pow(0.01,2.));
    double* rates = new double[nbasis*nbasis];
    double* dephasing = new double[nbasis*nbasis];
    for(int j = 0; j < nbasis; j++){
	for(int i = 0; i < nbasis; i++){
	    if(i != j){
		rates[i*nbasis+j] = dum/(pow((configurations[i]-configurations[j]),2.));
	    }else{
		rates[i*nbasis+j] = 0.;
	    }
	    dephasing[i*nbasis+j] = 0.;
	}
    }

    ofstream outf;
    outf.open("benchmark.sys");
    outf.write((char *) &nbasis, sizeof(int));
    outf.write((char *) &nstates, sizeof(int));
    outf.write((char *) configurations, sizeof(double)*nbasis);
    outf.write((char *) ionization, sizeof(double)*nbasis);
    outf.write((char *) dx, sizeof(double)*nbasis*nbasis);
    outf.write((char *) dy, sizeof(double)*nbasis*nbasis);
    outf.write((char *) dz, sizeof(double)*nbasis*nbasis);
    outf.write((char *) rates, sizeof(double)*nbasis*nbasis);
    outf.write((char *) dephasing, sizeof(double)*nbasis*nbasis);
    outf.write((char *) energy, sizeof(double)*nstates);
    outf.write((char *) eigenvectors, sizeof(double)*nbasis*nbasis);
    outf.close();
}

void   diag_mat(int nroao, double* mat, double* vals){

 static int N = -1 ;
  static double* WORK;
  
  int LDA = nroao, LWORK=3*nroao-1, INFO;
  
  if(N != nroao){
    if(WORK != NULL){
      delete [] WORK; 
     }
    
    WORK = new double[LWORK];
    N = nroao;
  }

  char JOBZ, UPLO;
  double* A    = mat; //    = vecs;
  double* W    = vals;
  JOBZ = 'V';
  UPLO = 'U';  

  for(int x = 0; x <  nroao*nroao; x++)
    A[x] = mat[x];
  dsyev_( &JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK, &INFO );
  for(int x = 0; x < nroao*nroao; x++) mat[x] = A[x];
}

/*
program genmat

    implicit none
    integer,parameter :: Nbasis=100,Nstates=100
    real(kind=8) :: energy(Nbasis),ionization(Nbasis),configurations(Nbasis),eigenvectors(Nbasis,Nbasis)
    real(kind=8) :: work(4*Nbasis),dipole(Nbasis,Nbasis,3),dephasing(Nbasis,Nbasis),rates(Nbasis,Nbasis)
    real(kind=8) :: dum
    integer :: i,j,lwork=4*Nbasis

! create fictive system
    configurations(1) = 0d0
    do i = 2,Nbasis
       configurations(i) = 0.1d0+(i-1)*0.001d0
    enddo

    ionization(:) = 0d0
    ionization(5:Nbasis) = 1d-5

    eigenvectors(:,:) = 0d0
    eigenvectors(1,1) = configurations(1)
    eigenvectors(2,2) = configurations(2)
    do i = 3,Nbasis
       eigenvectors(i,i-1) = 0.01d0
       eigenvectors(i  ,i) = configurations(i)
       eigenvectors(i-1,i) = 0.01d0
    enddo
    call dsyev('V','U',Nbasis,eigenvectors,Nbasis,energy,work,lwork,i)

! generate fictive dipole coupling
    dipole(:,:,:) = 0d0
    do i = 2,Nbasis
       dipole(i,i-1,1:3) = 0.1d0
       dipole(i-1,i,1:3) = 0.1d0
    enddo

! generate dubious rates
    dum = 1d-6*(0.01d0**2)
    do j = 1,Nbasis
    do i = 1,Nbasis
      if(i/=j) then
         rates(i,j) = dum/(configurations(i)-configurations(j))**2
      else
         rates(i,j) = 0d0
      endif
      dephasing(i,j) = 0d0
    enddo
    enddo

! output
    open(2,file='benchmark.sys',form='unformatted',access="stream")
    write(2) Nbasis,Nstates
    write(2) configurations
    write(2) ionization
    write(2) dipole(:,:,1)
    write(2) dipole(:,:,2)
    write(2) dipole(:,:,3)
    write(2) rates
    write(2) dephasing
    write(2) (energy(i),i=1,Nstates)
    do i=1,Nbasis
      write(2) eigenvectors(:,i)
    enddo
    close(2)

end program genmat
*/
