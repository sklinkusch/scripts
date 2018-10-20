#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

//void   diag_mat(int nroao, double* mat, double* vals);
//extern "C" void dsyev_(char* jobz, char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* info);
void read_wav_HF(char* wavfile, int nroao, double* MOens, double* MOs);

int main(int argc, char* argv[]){
    if(argc != 6){
	cerr << "Usage ./gensys <eigenenergy/vector file> <configuration file> <dipole file> <MO file> <output prefix>\n";
	exit(1);
    }
    double gsenergy;
    int omo, umo, llim, nroao, nbasis, nstates;
    ifstream ensf;
    ensf.open(argv[1]);
    ensf.read((char *) &gsenergy, sizeof(double));
    ensf.read((char *) &omo, sizeof(int));
    ensf.read((char *) &umo, sizeof(int));
    ensf.read((char *) &llim, sizeof(int));
    ensf.read((char *) &nroao, sizeof(int));
    ensf.read((char *) &nstates, sizeof(int));
    ensf.read((char *) &nbasis, sizeof(int));
    double* configurations = new double[nbasis];
    double* energy = new double[nstates];
    double* eigenvectors = new double[nbasis*nstates];
    double* downvectors = new double[nbasis*nstates];
    ensf.read((char *) energy, sizeof(double)*nstates);
    ensf.read((char *) eigenvectors, sizeof(double)*nstates*nbasis);
    ensf.read((char *) downvectors, sizeof(double)*nstates*nbasis);
    ensf.close();
    double* dx = new double[nbasis*nbasis];
    double* dy = new double[nbasis*nbasis];
    double* dz = new double[nbasis*nbasis];
    ifstream conf;
    ifstream dipf;
    int idum, jdum, kdum, ldum;
    conf.open(argv[2]);
    dipf.open(argv[3]);
    for(int i = 0; i < (nbasis*nbasis); i++){
	dx[i] = 0.;
	dy[i] = 0.;
	dz[i] = 0.;
    }
    configurations[0] = 0.;
    for(int i = 1; i < nbasis; i++) conf >> idum >> jdum >> configurations[i];
    for(int i = 0; i < nbasis; i++){
	for(int j = 0; j < nbasis; j++){
	    if((j != 0 && i != 0) || i == 0){
		dipf >> idum >> jdum >> kdum >> ldum >> dx[i*nbasis+j] >> dy[i*nbasis+j] >> dz[i*nbasis+j];
	    }else{
		dx[i*nbasis+j] = dx[j*nbasis+i];
		dy[i*nbasis+j] = dy[j*nbasis+i];
		dz[i*nbasis+j] = dz[j*nbasis+i];
	    }
	}
    }
    conf.close();
    dipf.close();
    double* MOens = new double[nroao];
    double* MOs = new double[nroao*nroao];
    char mosfile[512];
    sprintf(mosfile, "%s", argv[4]);
    read_wav_HF(mosfile, nroao, MOens, MOs);
    delete [] MOs;
    int MO_r = 0;
    double damp = 1.;
    double moenerg = 0.;
    double* ionization = new double[nbasis];
    for(int i = 0; i < nbasis; i++) ionization[i] = 0.;
    for(int i = 0; i < nbasis; i++){
	MO_r = (i-1)%umo+omo+llim;
	moenerg = MOens[MO_r];
	if(moenerg > 0.) ionization[i] = damp*sqrt(moenerg);
	cout << i << " " << moenerg << " " << ionization[i] << "\n";
    }
    double* rates = new double[nbasis*nbasis];
    double* dephasing = new double[nbasis*nbasis];
    double* mu_tot = new double[nbasis*nbasis];
    double* omega = new double[nbasis*nbasis];
    for(int i = 0; i < nbasis; i++){
	for(int j = 0; j < nbasis; j++){
	    mu_tot[i*nbasis+j] = sqrt(pow(dx[i*nbasis+j],2.)+pow(dy[i*nbasis+j],2.)+pow(dz[i*nbasis+j],2.));
	    omega[i*nbasis+j] = configurations[i] - configurations[j]; 
	}
    }
    for(int i = 0; i < nbasis; i++){
	for(int j = 0; j < nbasis; j++){
	    rates[i*nbasis+j] = (4.*(double) j*pow(fabs(mu_tot[i*nbasis+j]),2.)*pow(fabs(omega[i*nbasis+j]),3.))/(3.*pow(137.036,3.));
	}
    }
    double sum = 0.;
    double dum = (10./1.7776791e+09)*(1/pow(0.2254,2.));
    for(int i = 0; i << nbasis; i++){
	for(int j = 0; j < nbasis; j++){
	    sum = 0.;
	    for(int k = 0; k < nbasis; k++){
		sum += (rates[i*nbasis+k]+rates[j*nbasis+k]);
	    }
	    dephasing[i*nbasis+j] = sum/2. + dum*omega[i*nbasis+j];
	}
    }
    double* eigenvectors_t = new double[nbasis*nstates];
    double* rates_t = new double[nbasis*nbasis];
    double* dephasing_t = new double[nbasis*nbasis];
    int count = 0;
    for(int i = 0; i < nbasis; i++){
	for(int j = 0; j < nstates; j++){
	    eigenvectors_t[count] = eigenvectors[j*nbasis+i];
	    count++;
	}
    }
    for(int i = 0; i < nbasis; i++){
	for(int j = 0; j < nbasis; j++){
	    rates_t[i*nbasis+j] = rates[j*nbasis+i];
	    dephasing_t[i*nbasis+j] = dephasing[j*nbasis+i];
	}
    }

    ofstream outf;
    char dumc[512];
    sprintf(dumc, "%s.sys", argv[5]);
    outf.open(dumc);
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
    outf.write((char *) eigenvectors, sizeof(double)*nbasis*nstates);
    outf.close();

    ofstream logf;
    sprintf(dumc, "%s.log", argv[5]);
    logf.open(dumc);
    logf << "Nr of configuration state functions: " << nbasis << "\n";
    logf << "Nr of electronic states: " << nstates << "\n";
    logf << "Energies of configurations: \n";
    for(int i = 0; i < nbasis; i++){
	logf << configurations[i] << "  ";
	if(i%10 == 9) logf << "\n";
    }
    logf << "\nIonization rates: \n";
    for(int i = 0; i < nbasis; i++){
	logf << ionization[i] << "  ";
	if(i%10 == 9) logf << "\n";
    }
    logf << "\nDipole moments (x): \n";
    for(int i = 0; i < nbasis; i++){
	for(int j = 0; j < nbasis; j++){
	    logf << dx[i*nbasis+j] << " ";
	    if(j%10 == 9) logf << "\n";
	}
	logf << "\n";
    }
    logf << "Dipole moments (y): \n";
    for(int i = 0; i < nbasis; i++){
	for(int j = 0; j < nbasis; j++){
	    logf << dy[i*nbasis+j] << " ";
	    if(j%10 == 9) logf << "\n";
	}
	logf << "\n";
    }
    logf << "Dipole moments (z): \n";
    for(int i = 0; i < nbasis; i++){
	for(int j = 0; j < nbasis; j++){
	    logf << dz[i*nbasis+j] << " ";
	    if(j%10 == 9) logf << "\n";
	}
	logf << "\n";
    }
    logf << "Relaxation rates: \n";
    for(int i = 0; i < nbasis; i++){
	for(int j = 0; j < nbasis; j++){
	    logf << rates[i*nbasis+j] << " ";
	    if(j%10 == 9) logf << "\n";
	}
	logf << "\n";
    }
    logf << "Dephasing rates: \n";
    for(int i = 0; i < nbasis; i++){
	for(int j = 0; j < nbasis; j++){
	    logf << dephasing[i*nbasis+j] << " ";
	    if(j%10 == 9) logf << "\n";
	}
	logf << "\n";
    }
    logf << "Eigenenergies: \n";
    for(int i = 0; i < nstates; i++){
	logf << energy[i] << " ";
	if(i%10 == 9) logf << "\n";
    }
    logf << "\nEigenvectors: \n";
    for(int i = 0; i < nstates; i++){
	for(int j = 0; j < nbasis; j++){
	    logf << eigenvectors[i*nbasis+j] << " ";
	    if(j%10 == 9) logf << "\n";
	}
	logf << "\n";
    }
    logf.close();
}
/*
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
}*/

void read_wav_HF(char* wavfile, int nroao, double* MOens, double* MOs){
  ifstream inf(wavfile);
  int real_nroao = 0;
  inf.read((char *) &real_nroao, sizeof(int));
  if(real_nroao != nroao){
    cerr << "Wrong HF wavefunction size in read_wav_HF!\n";
    exit(3);
  }
  inf.read((char *) MOens,  sizeof(double)*nroao);
  inf.read((char *) MOs,    sizeof(double)*nroao*nroao);

  inf.close();
}

