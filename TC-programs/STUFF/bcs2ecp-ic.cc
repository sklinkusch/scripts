#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sstream>

using namespace std;

//Functions
int rem_com(char* filename, char* streamstring, int string_length);
void zero(int length, double* mat);
void ecp_read(int nros, double* inens, double* indipx, double* indipy, double* indipz, double* invr, double* inion, double* inkin, ifstream* inf);
void bcs_size(ifstream* inf, int &nroao, int &nroe, int &llim, int &ulim, int &omo, int &umo, int &nstates, int &cis_size, int switch_nstates);
void bcs_read(int nstates, int cis_size, double* cisens, double* muvalsx, double* muvexx, double* muvalsy, double* muvexy, double* muvalsz, double* muvexz, double* corr, ifstream* inf);
int calc_cissel(int cis_size, double* cisens, double* corr, double offset, double maxens);
void transform_energies(int cis_size, int cissel, int* enarray, double offset, double limit, double* ens, double* corr, double* newens);
inline bool test_array(int nroel, int testnr, int* array);
void transform_dipoles(int cis_size, int cissel, int* enarray, double* vals, double* vecs, double* mat);
void join_energies(int nros, int cissel, double* in, double* bcs, double* join);
void join_matrix(int nros, int cissel, int joined, double* ina, double* inb, double* join);
void join_dipoles(int nros, int cissel, int joined, double* inax, double* inay, double* inaz, double* inbx, double* inby, double* inbz, double* joinx, double* joiny, double* joinz);
void read_ions(int cis_size, int cissel, int* enarray, double* ionrat, ifstream* inf);
int count_nonion(int joined, double* joinion);
void presort(int joined, double* joinens, double* joinion, int* enssort, int cnt_nonion, int cnt_ion);
void sort_energies(int joined, int* enssort, double* joinens, double* newens);
void sort_dipoles(int joined, int* enssort, double* joindips, double* newdips);
void calc_vrs(int cissel, double* vr, double* ens, double matsubara);
double deltafkt(double x, double width);
void calc_ics(int nros, int cissel, double* ic, double* ens, double matsubara, double theta);
double calc_olapnorm(double* ens, int nros, int cissel, int i, double matsubara, double theta, int mode);
double calc_olap(double* ens, int i, int k, double matsubara, double theta, double norm);
void calc_kin(int cissel, double* ens, double* ion, double ip, double* kin);
void calc_kins(int joined, int nros, double* ens, double* ion, double ip, double* inkin, double* kin);
void write_ecp(ofstream* outf, int joined, double* ens, double* dipx, double* dipy, double* dipz, double* vr, double* ion, double* kin);

int main(int argc, char* argv[]){
    if(argc != 3){
	cerr << "Usage: ./bcs2ecp-ic <input file> <output prefix>\n";
	exit(1);
    }
    char infile[256];
    char outfile[256];
    sprintf(infile,"%s",argv[1]);
    sprintf(outfile,"%s.ecp",argv[2]);
    int buff_length = 65536;
    char* file_buff = new char[buff_length];
    rem_com(infile,file_buff,buff_length);
    istringstream ist(file_buff);
    int mode;
    char finecp[256];
    char finbcs[256];
    int switch_nstates;
    double offset;
    double maxen;
    char finirx[256];
    double ip;
    double matsubara;
    double theta;
    ist >> mode;
    if(mode == 1) ist >> finecp;
    ist >> finbcs >> switch_nstates >> offset >> maxen >> finirx >> ip >> matsubara;
    if(mode == 1) ist >> theta;
    ofstream outf;
    if(mode == 1){
	ifstream datf;                                                                                                         // define input filestream datf
	datf.open(finecp);                                                                                                     // connect filestream datf with filename finecp and open file
	int nros;                                                                                                              // number of electronic states
	datf.read((char *) &nros, sizeof(int));                                                                                // read (binary) nros
	cout << "Nr of states read\n";
	double* inens = new double[nros];                                                                                      // array of energies
	double* indipx = new double[nros*nros];                                                                                // array for "x" dipole matrix
	double* indipy = new double[nros*nros];                                                                                // array for "y" dipole matrix
	double* indipz = new double[nros*nros];                                                                                // array for "z" dipole matrix
	double* invr = new double[nros*nros];                                                                                  // vibrational relaxation rates
	double* inion = new double[nros];                                                                                      // array for ionization rates
	double* inkin = new double[nros];                                                                                      // array for kinetic energies
	ecp_read(nros, inens, indipx, indipy, indipz, invr, inion, inkin, &datf);                                              // function to read the previously defined arrays from the file finecp
	cout << "ECP file read\n";
	datf.close();                                                                                                          // close filestream datf
	int inbound = count_nonion(nros,inion);
	cout << "ECP-File non-ionizing: " << inbound << "\n";
	ifstream bcsf;                                                                                                         // define input filestream bcsf
	bcsf.open(finbcs);                                                                                                     // connect filestream bcsf with filename finbcs and open file
	int nroao, nroe, llim, ulim, omo, umo, nstates, cis_size;                                                              // define integers (number of atomic orbitals, number of electrons, lower/upper limit, occ./unocc. orb., number of states)
	bcs_size(&bcsf, nroao, nroe, llim, ulim, omo, umo, nstates, cis_size, switch_nstates);                                 // read size parameters nroao, nroe, llim, and ulim from bcs file and calculate omo, umo, and cis_size
	double* cisens = new double[nstates];                                                                                  // CIS energy array
	double* muvalsx = new double[nstates];                                                                                 // Array for eigenvalues of the "x" dipole matrix
	double* muvexx = new double[nstates*nstates];                                                                          // Array for eigenvectors of the "x" dipole matrix
	double* muvalsy = new double[nstates];                                                                                 // Array for eigenvalues of the "y" dipole matrix
	double* muvexy = new double[nstates*nstates];                                                                          // Array for eigenvectors of the "y" dipole matrix
	double* muvalsz = new double[nstates];                                                                                 // Array for the eigenvalues of the "z" dipole matrix
	double* muvexz = new double[nstates*nstates];                                                                          // Array for the eigenvectors of the "z" dipole matrix
	double* corr = new double[nstates];                                                                                    // Array for CIS(D) corrections to CIS energies
	bcs_read(nstates, cis_size, cisens, muvalsx, muvexx, muvalsy, muvexy, muvalsz, muvexz, corr, &bcsf);                   // read the previously defined arrays from the bcs file
	bcsf.close();                                                                                                          // close bcs file
	int cissel = calc_cissel(nstates, cisens, corr, offset, maxen);                                                        // calculate number of selected CIS(D) states
	int* enarray = new int[cissel];                                                                                        // energy array containing the numbers of selected states
	double* bcsens = new double[cissel];                                                                                   // array with the reduced number of CIS(D) energies
	transform_energies(nstates, cissel, enarray, offset, maxen, cisens, corr, bcsens);                                     // determine the elements of enarray and bcsens
	delete [] cisens; delete [] corr;                                                                                      // delete the arrays cisens and corr (not needed any more)
	double* bcsdipx = new double[cissel*cissel];                                                                           // array for the "x" dipole matrix
	transform_dipoles(nstates, cissel, enarray, muvalsx, muvexx, bcsdipx);                                                 // calculate "x" dipole matrix for all cis_size states and reduce it to the cissel*cissel matrix
	delete [] muvalsx; delete [] muvexx;                                                                                   // delete the arrays muvalsx and muvexx (not needed any more)
	double* bcsdipy = new double[cissel*cissel];                                                                           // repeat the prcedure for "y"
	transform_dipoles(nstates, cissel, enarray, muvalsy, muvexy, bcsdipy); 
	delete [] muvalsy; delete [] muvexy;
	double* bcsdipz = new double[cissel*cissel];                                                                           // repeat the procedure for "z"
	transform_dipoles(nstates, cissel, enarray, muvalsz, muvexz, bcsdipz);
	delete [] muvalsz; delete [] muvexz;
	double* bcsvr = new double[cissel*cissel];
	int joined = nros + cissel;                                                                                            // number of joined states
	double* joinens = new double[joined];                                                                                  // common energy array
	join_energies(nros, cissel, inens, bcsens, joinens);                                                                   // join energies to the common array joinens
	double* joindipx = new double[joined*joined];                                                                          // joined "x" dipole matrix
	double* joindipy = new double[joined*joined];                                                                          // joined "y" dipole matrix
	double* joindipz = new double[joined*joined];                                                                          // joined "z" dipole matrix
	join_matrix(nros, cissel, joined, indipx, bcsdipx, joindipx);
	join_matrix(nros, cissel, joined, indipy, bcsdipy, joindipy);
	join_matrix(nros, cissel, joined, indipz, bcsdipz, joindipz);
//	join_dipoles(nros, cissel, joined, indipx, indipy, indipz, bcsdipx, bcsdipy, bcsdipz, joindipx, joindipy, joindipz);   // join dipole matrices to the common arrays joindipx, joindipy, joindipz
	double* ionrat = new double[cissel];                                                                                   // array to store the selected ionization rates from the "new" states
	ifstream irxf;                                                                                                         // filestream to read ionization rates
	irxf.open(finirx);                                                                                                     // connect filestream irxf with filenume finirxf and open file
	read_ions(nstates, cissel, enarray, ionrat, &irxf);                                                                    // read the ionization rates for all states and store the ones for the selected states in ionrat
	irxf.close();                                                                                                          // close the filei
	int bcsbound = count_nonion(cissel,ionrat);
	cout << "BCS-File non-ionizing: " << bcsbound << "\n";
	calc_vrs(cissel, bcsvr, bcsens, matsubara);
	double* joinvr = new double[joined*joined];
	zero(joined*joined, joinvr);
	join_matrix(nros,cissel,joined,invr,bcsvr,joinvr);
	calc_ics(nros,cissel,joinvr,joinens,matsubara, theta);
	double* joinion = new double[joined];                                                                                  // common ionization rate array
	join_energies(nros, cissel, inion, ionrat, joinion);                                                                   // join ionization rates to the common array joinion
	delete [] inion; delete [] ionrat;                                                                                     // delete arrays that are not needed any more
	int cnt_nonion = count_nonion(joined, joinion);                                                                        // calculate the number of nonionizing states
	cout << "Total non-ionizing: " << cnt_nonion << "\n";
	int cnt_ion = joined - cnt_nonion;                                                                                     // calculate the number of ionizing states
	int* enssort = new int[joined];                                                                                        // integer array to show the sequence of the "original" states
	presort(joined, joinens, joinion, enssort, cnt_nonion, cnt_ion);                                                       // subroutine that sorts the "original" states in the integer array
	ip += offset;
	double* okin = new double[joined];
	calc_kins(joined, nros, joinens, joinion, ip, inkin, okin);
	double* newens = new double[joined];                                                                                   // "sorted" energy array
	sort_energies(joined, enssort, joinens, newens);                                                                       // sort the energies according to the order in the enssort array
	delete [] joinens;                                                                                                     // delete the old energy array
	double* newdipx = new double[joined*joined];                                                                           // "sorted" "x" dipole array
	sort_dipoles(joined, enssort, joindipx, newdipx);                                                                      // sort the dipole moments
	delete [] joindipx;                                                                                                    // delete the old dipole array
	double* newdipy = new double[joined*joined];                                                                           // repeat the procedure for y
	sort_dipoles(joined, enssort, joindipy, newdipy);
	delete [] joindipy;
	double* newdipz = new double[joined*joined];                                                                           // repeat the procedure for z
	sort_dipoles(joined, enssort, joindipz, newdipz);
	delete [] joindipz;
	double* newvr = new double[joined*joined];
	sort_dipoles(joined, enssort, joinvr, newvr);
	delete [] joinvr;
	double* newion = new double[joined];                                                                                   // repeat the procedure for the ionization rates
	sort_energies(joined, enssort, joinion, newion);
	delete [] joinion;
	double* kin = new double[joined];
	sort_energies(joined, enssort, okin, kin);
	delete [] okin;
	ofstream outf;                                                                                                         // define output filestream
	outf.open(outfile);                                                                                                    // connect filestream outf with filename foutecp and open the file
	write_ecp(&outf, joined, newens, newdipx, newdipy, newdipz, newvr, newion, kin);                                       // subroutine to write out the generated data
	outf.close();                                                                                                          // close the file
    }else{
        ifstream bcsf;                                                                                                         // define input filestream bcsf
	bcsf.open(finbcs);                                                                                                     // connect filestream bcsf with filename finbcs and open file
	int nroao, nroe, llim, ulim, omo, umo, nstates, cis_size;                                                              // define integers (number of atomic orbitals, number of electrons, lower/upper limit, occ./unocc. orb., number of states)
	bcs_size(&bcsf, nroao, nroe, llim, ulim, omo, umo, nstates, cis_size, switch_nstates);                                 // read size parameters nroao, nroe, llim, and ulim from bcs file and calculate omo, umo, and cis_size
	cout << "Reading bcs data from: " << finbcs << "\n";
	cout << "Number of atomic orbitals: " << nroao << "\n";
	cout << "Number of electrons: " << nroe << "\n";
	cout << "Lower limit: " << llim << "\n";
	cout << "Upper limit: " << ulim << "\n";
	cout << "Occupied MOs: " << omo << "\n";
	cout << "Unoccupied MOs: " << umo << "\n";
	cout << "Number of states: " << nstates << "\n";
	cout << "Number of configurations: " << cis_size << "\n";
	double* cisens = new double[nstates];                                                                                  // CIS energy array
	double* muvalsx = new double[nstates];                                                                                 // Array for eigenvalues of the "x" dipole matrix
	double* muvexx = new double[nstates*nstates];                                                                          // Array for eigenvectors of the "x" dipole matrix
	double* muvalsy = new double[nstates];                                                                                 // Array for eigenvalues of the "y" dipole matrix
	double* muvexy = new double[nstates*nstates];                                                                          // Array for eigenvectors of the "y" dipole matrix
	double* muvalsz = new double[nstates];                                                                                 // Array for the eigenvalues of the "z" dipole matrix
	double* muvexz = new double[nstates*nstates];                                                                          // Array for the eigenvectors of the "z" dipole matrix
	double* corr = new double[nstates];                                                                                    // Array for CIS(D) corrections to CIS energies
	bcs_read(nstates, cis_size, cisens, muvalsx, muvexx, muvalsy, muvexy, muvalsz, muvexz, corr, &bcsf);                   // read the previously defined arrays from the bcs file
	cout << "BCS data read\n";
	bcsf.close();                                                                                                          // close bcs file
	int cissel = calc_cissel(nstates, cisens, corr, offset, maxen);                                                        // calculate number of selected CIS(D) states
	cout << "Number of selected CIS states: " << cissel << "\n";
	int* enarray = new int[cissel];                                                                                        // energy array containing the numbers of selected states
	double* bcsens = new double[cissel];                                                                                   // array with the reduced number of CIS(D) energies
	transform_energies(nstates, cissel, enarray, offset, maxen, cisens, corr, bcsens);                                     // determine the elements of enarray and bcsens
	cout << "Energies transformed\n";
	delete [] cisens; delete [] corr;                                                                                      // delete the arrays cisens and corr (not needed any more)
	double* bcsdipx = new double[cissel*cissel];                                                                           // array for the "x" dipole matrix
	transform_dipoles(nstates, cissel, enarray, muvalsx, muvexx, bcsdipx);                                                 // calculate "x" dipole matrix for all cis_size states and reduce it to the cissel*cissel matrix
	delete [] muvalsx; delete [] muvexx;                                                                                   // delete the arrays muvalsx and muvexx (not needed any more)
	double* bcsdipy = new double[cissel*cissel];                                                                           // repeat the prcedure for "y"
	transform_dipoles(nstates, cissel, enarray, muvalsy, muvexy, bcsdipy);
	delete [] muvalsy; delete [] muvexy;
	double* bcsdipz = new double[cissel*cissel];                                                                           // repeat the procedure for "z"
	transform_dipoles(nstates, cissel, enarray, muvalsz, muvexz, bcsdipz);
	cout << "Dipole moments transformed\n";
	delete [] muvalsz; delete [] muvexz;
	double* ionrat = new double[cissel];                                                                                   // array to store the selected ionization rates from the "new" states
	ifstream irxf;                                                                                                         // filestream to read ionization rates
	irxf.open(finirx);
	read_ions(nstates, cissel, enarray, ionrat, &irxf);                                                                   // read the ionization rates for all states and store the ones for the selected states in ionrat
	cout << "Ionization rates read\n";
	irxf.close();
	int cnt_nonion = count_nonion(cissel, ionrat);                                                                         // calculate the number of nonionizing states
	cout << "Number of nonionizing states: " << cnt_nonion << "\n";
	int cnt_ion = cissel - cnt_nonion;                                                                                     // calculate the number of ionizing states
	cout << "Number of ionizing states: " << cnt_ion << "\n";
	int* enssort = new int[cissel];                                                                                        // integer array to show the sequence of the "original" states
	presort(cissel, bcsens, ionrat, enssort, cnt_nonion, cnt_ion);                                                         // subroutine that sorts the "original" states in the integer array
	double* newens = new double[cissel];                                                                                   // "sorted" energy array
	sort_energies(cissel, enssort, bcsens, newens);                                                                        // sort the energies according to the order in the enssort array
	cout << "Energies sorted\n";
	delete [] bcsens;                                                                                                      // delete the old energy array
	double* newdipx = new double[cissel*cissel];                                                                           // "sorted" "x" dipole array
	sort_dipoles(cissel, enssort, bcsdipx, newdipx);                                                                       // sort the dipole moments
	delete [] bcsdipx;                                                                                                     // delete the old dipole array
	double* newdipy = new double[cissel*cissel];                                                                           // repeat the procedure for y
	sort_dipoles(cissel, enssort, bcsdipy, newdipy);
	delete [] bcsdipy;
	double* newdipz = new double[cissel*cissel];                                                                           // repeat the procedure for z
	sort_dipoles(cissel, enssort, bcsdipz, newdipz);
	cout << "Dipole moments sorted\n";
	delete [] bcsdipz;
	double* newion = new double[cissel];                                                                                   // repeat the procedure for the ionization rates
	sort_energies(cissel, enssort, ionrat, newion);
	cout << "Ionization rates sorted\n";
	delete [] ionrat;
	ip += offset;
	double* vr  = new double[cissel*cissel];
	calc_vrs(cissel,vr,newens,matsubara);
	double* kin = new double[cissel];
	calc_kin(cissel, newens, newion, ip, kin);
	cout << "Kinetic energies calculated\n";
	outf.open(outfile);                                                                                                    // connect filestream outf with filename foutecp and open the file
	write_ecp(&outf, cissel, newens, newdipx, newdipy, newdipz, vr, newion, kin);                                          // subroutine to write out the generated data
	cout << "ECP file written\n";
	outf.close();                                                                                                          // close the file
    }
}

int rem_com(char* filename, char* streamstring, int string_length){
  const char com_B = '#';
  const char com_E = '\n';
  
  int pos = 0;
  char cc;

  ifstream inf(filename);
  
  while(inf.get(cc)&& pos < string_length-1){
    if(cc != com_B) 
      streamstring[pos++] = cc;
    else{
      while(cc != com_E && inf.get(cc));
      streamstring[pos++] = com_E;
    }
  }
  streamstring[pos] = 0;
  if(pos == string_length-1){
    cerr << "Buffer size exceeded !\n"; exit(0);
  }
  return(strlen(streamstring));
}

void zero(int length, double* mat){
#pragma omp parallel for
    for(int i = 0; i < length; i++) mat[i] = 0.;
}

void ecp_read(int nros, double* inens, double* indipx, double* indipy, double* indipz, double* invr, double* inion, double* inkin, ifstream* inf){
    inf->read((char *) inens, sizeof(double)*nros);
    inf->read((char *) indipx, sizeof(double)*nros*nros);
    inf->read((char *) indipy, sizeof(double)*nros*nros);
    inf->read((char *) indipz, sizeof(double)*nros*nros);
    inf->read((char *) invr, sizeof(double)*nros*nros);
    inf->read((char *) inion, sizeof(double)*nros);
    inf->read((char *) inkin, sizeof(double)*nros);
}

void bcs_size(ifstream* inf, int &nroao, int &nroe, int &llim, int &ulim, int &omo, int &umo, int &nstates, int &cis_size, int switch_nstates){
    inf->read((char *) &nroao, sizeof(int));
    inf->read((char *) &nroe, sizeof(int));
    inf->read((char *) &llim, sizeof(int));
    inf->read((char *) &ulim, sizeof(int));
    if(switch_nstates == 1) inf->read((char *) &nstates, sizeof(int));
    omo = nroe/2 - llim;
    umo = ulim + 1 - omo - llim;
    cis_size = omo * umo + 1;
    if(switch_nstates != 1) nstates = cis_size;
}

void bcs_read(int nstates, int cis_size, double* cisens, double* muvalsx, double* muvexx, double* muvalsy, double* muvexy, double* muvalsz, double* muvexz, double* corr, ifstream* inf){
    inf->read((char *) cisens, sizeof(double)*nstates);
    double* cisvex = new double[nstates*cis_size];
    inf->read((char *) cisvex, sizeof(double)*nstates*cis_size);
    delete [] cisvex;
    inf->read((char *) muvalsx, sizeof(double)*nstates);
    inf->read((char *) muvexx, sizeof(double)*nstates*nstates);
    inf->read((char *) muvalsy, sizeof(double)*nstates);
    inf->read((char *) muvexy, sizeof(double)*nstates*nstates);
    inf->read((char *) muvalsz, sizeof(double)*nstates);
    inf->read((char *) muvexz, sizeof(double)*nstates*nstates);
    inf->read((char *) corr, sizeof(double)*nstates);
}

int calc_cissel(int cis_size, double* cisens, double* corr, double offset, double maxens){
    double dumdum;
    int cissel = 0;
    for(int i = 0; i < cis_size; i++){
	dumdum = offset + cisens[i] + corr[i];
	if(dumdum <= maxens) cissel += 1;
    }
    return(cissel);
}

void transform_energies(int cis_size, int cissel, int* enarray, double offset, double limit, double* ens, double* corr, double* newens){
    double dumdum;
    int k = 0;
    for(int i = 0; i < cis_size; i++){
	dumdum = offset + ens[i] + corr[i];
	if(dumdum <= limit){
  	    newens[k] = dumdum;
	//    cout << i << " " << k << " " << newens[k] << "\n";
	    enarray[k] = i;
	    k += 1;
	}else{
	    continue;
	}
    }
}

inline bool test_array(int nroel, int testnr, int* array){
    bool dumbool = false;
    for(int i = 0; i < nroel; i++){
	if(array[i] == testnr){
	    dumbool = true;
	    break;
	}else{
	    continue;
	}
    }
    return(dumbool);
}

void transform_dipoles(int cis_size, int cissel, int* enarray, double* vals, double* vecs, double* mat){
    double* premat = new double[cis_size*cis_size];
    double dumval = 0.;
#pragma omp parallel for reduction(+:dumval)
    for(int i = 0; i < cis_size; i++){
	dumval = 0.;
	for(int j = 0; j < cis_size; j++){
	    dumval += vecs[j*cis_size+i]*vals[j]*vecs[j*cis_size+i];
	}
	premat[i*cis_size+i] = dumval;
    }
#pragma omp parallel for reduction(+:dumval)
    for(int i = 0; i < cis_size; i++){
	for(int j = i+1; j < cis_size; j++){
	    dumval = 0.;
	    for(int k = 0; k < cis_size; k++){
		dumval += vecs[k*cis_size+i]*vals[k]*vecs[k*cis_size+j];
	    }
	    premat[i*cis_size+j] = dumval;
	    premat[j*cis_size+i] = dumval;
	}
    }
    bool testi, testj;
    int k = 0;
    for(int i = 0; i < cis_size; i++){
	testi = test_array(cissel, i, enarray);
	if(testi == false) continue;
	for(int j = 0; j < cis_size; j++){
	    testj = test_array(cissel, j, enarray);
	    if(testj == false) continue;
	    mat[k] = premat[i*cis_size+j];
	    k++;
	}
    }
}

void join_energies(int nros, int cissel, double* in, double* bcs, double* join){
    for(int i = 0; i < nros; i++){
	join[i] = in[i];
    }
    for(int i = 0; i < cissel; i++){
	join[nros+i] = bcs[i];
    }
}

void join_matrix(int nros, int cissel, int joined, double* ina, double* inb, double* join){
#pragma omp parallel for
    for(int i = 0; i < joined; i++){
	for(int j = 0; j < joined; j++){
	    if(i < nros && j < nros){
		join[i*joined+j] = ina[i*nros+j];
	    }else if(i < nros || j < nros){
		join[i*joined+j] = 0.;
	    }else{
		join[i*joined+j] = inb[(i-nros)*cissel+(j-nros)];
	    }
	}
    }
}

void join_dipoles(int nros, int cissel, int joined, double* inax, double* inay, double* inaz, double* inbx, double* inby, double* inbz, double* joinx, double* joiny, double* joinz){
#pragma omp parallel for
    for(int i = 0; i < joined; i++){
	for(int j = 0; j < joined; j++){
	    if(i < nros && j < nros){
		joinx[i*joined+j] = inax[i*nros+j];
		joiny[i*joined+j] = inay[i*nros+j];
		joinz[i*joined+j] = inaz[i*nros+j];
	    }else if(i < nros || j < nros){
		joinx[i*joined+j] = 0.;
		joiny[i*joined+j] = 0.;
		joinz[i*joined+j] = 0.;
	    }else{
		joinx[i*joined+j] = inbx[(i-nros)*cissel+(j-nros)];
		joiny[i*joined+j] = inby[(i-nros)*cissel+(j-nros)];
		joinz[i*joined+j] = inbz[(i-nros)*cissel+(j-nros)];
	    }
	}
    }
}

void read_ions(int cis_size, int cissel, int* enarray, double* ionrat, ifstream* inf){
    int current;
    double* premat = new double[cis_size];
    for(int i = 0; i < cis_size; i++){
	*inf >> current >> premat[i];
    }
    bool testi;
    int k = 0;
    for(int i = 0; i < cis_size; i++){
	testi = test_array(cissel, i, enarray);
	if(testi == true){
	    ionrat[k] = premat[i];
	    k++;
	}
    }
}

int count_nonion(int joined, double* joinion){
    int nonion = 0;
    for(int i = 0; i < joined; i++){
	if(fabs(joinion[i]) < 1.e-12){
	   nonion++;
	}
    }
    return(nonion);
}

void presort(int joined, double* joinens, double* joinion, int* enssort, int cnt_nonion, int cnt_ion){
    // Bring the non-ionizing states to front
    int k = 0;
    for(int i = 0; i < joined; i++){
	if(joinion[i] == 0){
	    enssort[k] = i;
	    k++;
	}
    }
    // then the ionizing states
    for(int i = 0; i < joined; i++){
	if(joinion[i] != 0){
	    enssort[k] = i;
	    k++;
	}
    }
    int tmp;
    for(int z = 0; z < (cnt_nonion-2); z++){
	if(z%2 == 0){
//#pragma omp parallel for
	    for(int i = 0; i < (cnt_nonion/2); i++){
		if(joinens[enssort[2*i+1]] < joinens[enssort[2*i]]){
		    tmp = enssort[2*i];
		    enssort[2*i] = enssort[2*i+1];
		    enssort[2*i+1] = tmp;
		}
	    }
	}else{
//#pragma omp parallel for
	    for(int i = 0; i < (cnt_nonion/2-1); i++){
		if(joinens[enssort[2*i+2]] < joinens[enssort[2*i+1]]){
		    tmp = enssort[2*i+1];
		    enssort[2*i+1] = enssort[2*i+2];
		    enssort[2*i+2] = tmp;
		}
	    }
	}
    }
    for(int z = 0; z < (cnt_ion-2); z++){
	if(z%2 == 0){
//#pragma omp parallel for
	    for(int i = (cnt_nonion/2); i < (joined/2); i++){
		if(joinens[enssort[2*i+1]] < joinens[enssort[2*i]]){
		    tmp = enssort[2*i];
		    enssort[2*i] = enssort[2*i+1];
		    enssort[2*i+1] = tmp;
		}
	    }
	}else{
//#pragma omp parallel for
	    for(int i = (cnt_nonion/2); i < (joined/2-1); i++){
		if(joinens[enssort[2*i+2]] < joinens[enssort[2*i+1]]){
		    tmp = enssort[2*i+1];
		    enssort[2*i+1] = enssort[2*i+2];
		    enssort[2*i+2] = tmp;
		}
	    }
	}
    }
}

void sort_energies(int joined, int* enssort, double* joinens, double* newens){
#pragma omp parallel for
    for(int i = 0; i < joined; i++){
	newens[i] = joinens[enssort[i]];
    }
}

void sort_dipoles(int joined, int* enssort, double* joindips, double* newdips){
#pragma omp parallel for
    for(int i = 0; i < joined; i++){
	for(int j = 0; j < joined; j++){
	    newdips[i*joined+j] = joindips[enssort[i]*joined+enssort[j]];
	}
    }
}

void calc_vrs(int cissel, double* vr, double* ens, double matsubara){
    int i,j;
    double diffen, delta_abs;
    zero(cissel*cissel, vr);
#pragma omp parallel for default(shared) private(j)
    for(i = 0; i < cissel; i++){
	for(j = (i+1); j < cissel; j++){
	    diffen = fabs(ens[i] - ens[j]);
	    delta_abs = deltafkt(diffen,matsubara);
	    vr[i*cissel+j] = delta_abs;
	    vr[j*cissel+i] = delta_abs;
	}
    }
}

double deltafkt(double x, double width){
    double delta = 0.;
    delta = (width/(pow(x,2.) + pow(width,2.)));
    return(delta);
}

void calc_ics(int nros, int cissel, double* ic, double* ens, double matsubara, double theta){
    int joined = nros + cissel;
    double norm, lortz, deltaE, summand, sum, olap;
//    printf("%4s %4s %4s %18s %8s %20s %20s %20s\n","i","j","k","overlap","deltaE","lorentz","product","summe");
//    cout << "----------------------------------------------------------------------------------------------------------------------------------------------\n";
    for(int i = 0; i < nros; i++){
	for(int j = nros; j < joined; j++){
	    if(ens[i] > ens[j]){
		norm = calc_olapnorm(ens,nros,cissel,i,matsubara, theta,2);
		sum = 0.;
#pragma omp parallel for reduction(+:sum)
		for(int k = nros; k < joined; k++){
		    if(k != j){
			olap = calc_olap(ens,i,k,matsubara, theta,norm);
			deltaE = fabs(ens[j] - ens[k]);
			lortz = deltafkt(deltaE,matsubara);
			summand = lortz*pow(olap,2.);
			sum += summand;
//			printf("A%3d B%3d B%3d %18.16f %8.6f %20.16f %20.16f %20.16f\n",i,j-nros,k-nros,olap,deltaE,lortz,summand,sum);
		    }
		}
		ic[i*joined+j] = sum;
		ic[j*joined+i] = sum;
	    }else{
		norm = calc_olapnorm(ens,nros,cissel,j,matsubara, theta,1);
		sum = 0.;
#pragma omp parallel for reduction(+:sum)
		for(int k = 0; k < nros; k++){
		    if(k != i){
			olap = calc_olap(ens,j,k,matsubara, theta,norm);
			deltaE = fabs(ens[i] - ens[k]);
			lortz = deltafkt(deltaE,matsubara);
			summand = lortz*pow(olap,2.);
			sum += summand;
//			printf("B%3d A%3d A%3d %18.16f %8.6f %20.16f %20.16f %20.16f\n",j-nros,i,k,olap,deltaE,lortz,summand,sum);
		    }
		}
		ic[i*joined+j] = sum;
		ic[j*joined+i] = sum;
	    }
	}
    }
}

double calc_olapnorm(double* ens, int nros, int cissel, int i, double matsubara, double theta, int mode){
    int joined = nros + cissel;
    double sum = 0.;
    if(mode == 1){
	for(int x = 0; x < nros; x++) sum += calc_olap(ens,i,x,matsubara, theta,1.);
    }else{
	for(int x = nros; x < joined; x++) sum += calc_olap(ens,i,x,matsubara, theta,1.);
    }
    return(sum);
}

double calc_olap(double* ens, int i, int k, double matsubara, double theta, double norm){
    double sigma = theta * matsubara;
    double olap = (exp(-0.5*pow(ens[i]-ens[k],2.)/pow(sigma,2.)))/norm;
    return(olap);
}

void calc_kin(int cissel, double* ens, double* ion, double ip, double* kin){
    double eps = 1.e-10;
#pragma omp parallel for
    for(int i = 0; i < cissel; i++){
	if(ion[i] < eps){
	    kin[i] = 0.;
	}else{
	    kin[i] = ens[i] - ip;
	}
    }
}

void calc_kins(int joined, int nros, double* ens, double* ion, double ip, double* inkin, double* kin){
    double eps = 1.e-10;
#pragma omp parallel for
    for(int i = 0; i < nros; i++){
	kin[i] = inkin[i];
    }
#pragma omp parallel for
    for(int i = nros; i < joined; i++){
	if(ion[i] < eps){
	    kin[i] = 0.;
	}else{
	    kin[i] = ens[i] - ip;
	}
    }
}

void write_ecp(ofstream* outf, int joined, double* ens, double* dipx, double* dipy, double* dipz, double* vr, double* ion, double* kin){
    outf->write((char *) &joined, sizeof(int));
    outf->write((char *) ens, sizeof(double)*joined);
    outf->write((char *) dipx, sizeof(double)*joined*joined);
    outf->write((char *) dipy, sizeof(double)*joined*joined);
    outf->write((char *) dipz, sizeof(double)*joined*joined);
    outf->write((char *) vr, sizeof(double)*joined*joined);
    outf->write((char *) ion, sizeof(double)*joined);
    outf->write((char *) kin, sizeof(double)*joined);
    outf->flush();
}

