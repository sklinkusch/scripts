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
void ecp_read(int nros, double* inens, double* indipx, double* indipy, double* indipz, double* inion, double* inkin, ifstream* inf);
void bcs_size(ifstream* inf, int &nroao, int &nroe, int &llim, int &ulim, int &omo, int &umo, int &cis_size);
void bcs_read(int cis_size, double* cisens, double* muvalsx, double* muvexx, double* muvalsy, double* muvexy, double* muvalsz, double* muvexz, double* corr, ifstream* inf);
int calc_cissel(int cis_size, double* cisens, double* corr, double offset, double maxens);
void transform_energies(int cis_size, int cissel, int* enarray, double offset, double limit, double* ens, double* corr, double* newens);
inline bool test_array(int nroel, int testnr, int* array);
void transform_dipoles(int cis_size, int cissel, int* enarray, double* vals, double* vecs, double* mat);
void join_energies(int nros, int cissel, double* in, double* bcs, double* join);
void join_dipoles(int nros, int cissel, int joined, double* inax, double* inay, double* inaz, double* inbx, double* inby, double* inbz, double* joinx, double* joiny, double* joinz);
void read_ions(int cis_size, int cissel, int* enarray, double* ionrat, ifstream* inf);
int count_nonion(int joined, double* joinion);
void presort(int joined, double* joinens, double* joinion, int* enssort, int cnt_nonion, int cnt_ion);
void sort_energies(int joined, int* enssort, double* joinens, double* newens);
void sort_dipoles(int joined, int* enssort, double* joindips, double* newdips);
void calc_kin(int cissel, double* ens, double* ion, double ip, double* kin);
void calc_kins(int joined, int nros, double* ens, double* ion, double ip, double* inkin, double* kin);
void write_ecp(ofstream* outf, int joined, double* ens, double* dipx, double* dipy, double* dipz, double* ion, double* kin);

int main(int argc, char* argv[]){
    if(argc != 7 && argc != 8){
	cerr << "Usage: ./bcs2ecp <input ecp file> <input bcs file> <offset energy> <maximal energy> <input irx file> <IP> <output prefix>\n";
	cerr << "or: ./bcs2ecp <input bcs file> <offset energy> <maximal energy> <input irx file> <IP> <output prefix>\n";
	exit(1);
    }
    if(argc == 9){
	char finecp[128];                                                                                                      // name of input ecp file
	sprintf(finecp, "%s", argv[1]);                                                                                        // store commandline argument 2 in finecp
	ifstream datf;                                                                                                         // define input filestream datf
	datf.open(finecp);                                                                                                     // connect filestream datf with filename finecp and open file
	int nros;                                                                                                              // number of electronic states
	datf.read((char *) &nros, sizeof(int));                                                                                // read (binary) nros
	double* inens = new double[nros];                                                                                      // array of energies
	double* indipx = new double[nros*nros];                                                                                // array for "x" dipole matrix
	double* indipy = new double[nros*nros];                                                                                // array for "y" dipole matrix
	double* indipz = new double[nros*nros];                                                                                // array for "z" dipole matrix
	double* inion = new double[nros];                                                                                      // array for ionization rates
	double* inkin = new double[nros];                                                                                      // array for kinetic energies
	ecp_read(nros, inens, indipx, indipy, indipz, inion, inkin, &datf);                                                    // function to read the previously defined arrays from the file finecp
	datf.close();                                                                                                          // close filestream datf
	cout << "ECP-File non-ionizing: " << count_nonion(nros, inion) << "\n";
	ifstream bcsf;                                                                                                         // define input filestream bcsf
	char finbcs[128];                                                                                                      // name of input bcs file
	sprintf(finbcs, "%s", argv[2]);                                                                                        // store commandline argument 3 in finbcs
	bcsf.open(finbcs);                                                                                                     // connect filestream bcsf with filename finbcs and open file
	int nroao, nroe, llim, ulim, omo, umo, cis_size;                                                                       // define integers (number of atomic orbitals, number of electrons, lower/upper limit, occ./unocc. orb., number of states)
	bcs_size(&bcsf, nroao, nroe, llim, ulim, omo, umo, cis_size);                                                          // read size parameters nroao, nroe, llim, and ulim from bcs file and calculate omo, umo, and cis_size
	double* cisens = new double[cis_size];                                                                                 // CIS energy array
	double* muvalsx = new double[cis_size];                                                                                // Array for eigenvalues of the "x" dipole matrix
	double* muvexx = new double[cis_size*cis_size];                                                                        // Array for eigenvectors of the "x" dipole matrix
	double* muvalsy = new double[cis_size];                                                                                // Array for eigenvalues of the "y" dipole matrix
	double* muvexy = new double[cis_size*cis_size];                                                                        // Array for eigenvectors of the "y" dipole matrix
	double* muvalsz = new double[cis_size];                                                                                // Array for the eigenvalues of the "z" dipole matrix
	double* muvexz = new double[cis_size*cis_size];                                                                        // Array for the eigenvectors of the "z" dipole matrix
	double* corr = new double[cis_size];                                                                                   // Array for CIS(D) corrections to CIS energies
	bcs_read(cis_size, cisens, muvalsx, muvexx, muvalsy, muvexy, muvalsz, muvexz, corr, &bcsf);                            // read the previously defined arrays from the bcs file
	bcsf.close();                                                                                                          // close bcs file
	double offset = strtod(argv[3], NULL);                                                                                 // energy offset (added to each electronic state) in Hartree
	double maxen = strtod(argv[4], NULL);                                                                                  // maximal energy (higher states are cut off) in Hartree
	int cissel = calc_cissel(cis_size, cisens, corr, offset, maxen);                                                       // calculate number of selected CIS(D) states
	int* enarray = new int[cissel];                                                                                        // energy array containing the numbers of selected states
	double* bcsens = new double[cissel];                                                                                   // array with the reduced number of CIS(D) energies
	transform_energies(cis_size, cissel, enarray, offset, maxen, cisens, corr, bcsens);                                    // determine the elements of enarray and bcsens
	delete [] cisens; delete [] corr;                                                                                      // delete the arrays cisens and corr (not needed any more)
	double* bcsdipx = new double[cissel*cissel];                                                                           // array for the "x" dipole matrix
	transform_dipoles(cis_size, cissel, enarray, muvalsx, muvexx, bcsdipx);                                                // calculate "x" dipole matrix for all cis_size states and reduce it to the cissel*cissel matrix
	delete [] muvalsx; delete [] muvexx;                                                                                   // delete the arrays muvalsx and muvexx (not needed any more)
	double* bcsdipy = new double[cissel*cissel];                                                                           // repeat the prcedure for "y"
	transform_dipoles(cis_size, cissel, enarray, muvalsy, muvexy, bcsdipy);
	delete [] muvalsy; delete [] muvexy;
	double* bcsdipz = new double[cis_size*cis_size];                                                                       // repeat the procedure for "z"
	transform_dipoles(cis_size, cissel, enarray, muvalsz, muvexz, bcsdipz);
	delete [] muvalsz; delete [] muvexz;
	int joined = nros + cissel;                                                                                            // number of joined states
	double* joinens = new double[joined];                                                                                  // common energy array
	join_energies(nros, cissel, inens, bcsens, joinens);                                                                   // join energies to the common array joinens
	double* joindipx = new double[joined*joined];                                                                          // joined "x" dipole matrix
	double* joindipy = new double[joined*joined];                                                                          // joined "y" dipole matrix
	double* joindipz = new double[joined*joined];                                                                          // joined "z" dipole matrix
	join_dipoles(nros, cissel, joined, indipx, indipy, indipz, bcsdipx, bcsdipy, bcsdipz, joindipx, joindipy, joindipz);   // join dipole matrices to the common arrays joindipx, joindipy, joindipz
	double* ionrat = new double[cissel];                                                                                   // array to store the selected ionization rates from the "new" states
	ifstream irxf;                                                                                                         // filestream to read ionization rates
	char finirx[128];                                                                                                      // filename with ionization rates
	sprintf(finirx, "%s", argv[5]);                                                                                        // store commandline argument 6 in finirx
	irxf.open(finirx);                                                                                                     // connect filestream irxf with filenume finirxf and open file
	read_ions(cis_size, cissel, enarray, ionrat, &irxf);                                                                   // read the ionization rates for all states and store the ones for the selected states in ionrat
	irxf.close();                                                                                                          // close the filei
	cout << "BCS-File non-ionizing: " << count_nonion(cissel, ionrat) << "\n";
	double* joinion = new double[joined];                                                                                  // common ionization rate array
	join_energies(nros, cissel, inion, ionrat, joinion);                                                                   // join ionization rates to the common array joinion
	delete [] inion; delete [] ionrat;                                                                                     // delete arrays that are not needed any more
	int cnt_nonion = count_nonion(joined, joinion);                                                                        // calculate the number of nonionizing states
	cout << "Total non-ionizing: " << cnt_nonion << "\n";
	int cnt_ion = joined - cnt_nonion;                                                                                     // calculate the number of ionizing states
	int* enssort = new int[joined];                                                                                        // integer array to show the sequence of the "original" states
	presort(joined, joinens, joinion, enssort, cnt_nonion, cnt_ion);                                                       // subroutine that sorts the "original" states in the integer array
	double ip = strtod(argv[6], NULL);
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
	double* newion = new double[joined];                                                                                   // repeat the procedure for the ionization rates
	sort_energies(joined, enssort, joinion, newion);
	delete [] joinion;
	double* kin = new double[joined];
	sort_energies(joined, enssort, okin, kin);
	delete [] okin;
	char foutecp[128];                                                                                                     // filename of output ecp file
	sprintf(foutecp, "%s.ecp", argv[7]);                                                                                   // store commandline argument 9 in foutecp
	ofstream outf;                                                                                                         // define output filestream
	outf.open(foutecp);                                                                                                    // connect filestream outf with filename foutecp and open the file
	write_ecp(&outf, joined, newens, newdipx, newdipy, newdipz, newion, kin);                                              // subroutine to write out the generated data
	outf.close();                                                                                                          // close the file
    }else{
        ifstream bcsf;                                                                                                         // define input filestream bcsf
	char finbcs[128];                                                                                                      // name of input bcs file
	sprintf(finbcs, "%s", argv[1]);                                                                                        // store commandline argument 3 in finbcs
	bcsf.open(finbcs);                                                                                                     // connect filestream bcsf with filename finbcs and open file
	int nroao, nroe, llim, ulim, omo, umo, cis_size;                                                                       // define integers (number of atomic orbitals, number of electrons, lower/upper limit, occ./unocc. orb., number of states)
	bcs_size(&bcsf, nroao, nroe, llim, ulim, omo, umo, cis_size);                                                          // read size parameters nroao, nroe, llim, and ulim from bcs file and calculate omo, umo, and cis_size
	double* cisens = new double[cis_size];                                                                                 // CIS energy array
	double* muvalsx = new double[cis_size];                                                                                // Array for eigenvalues of the "x" dipole matrix
	double* muvexx = new double[cis_size*cis_size];                                                                        // Array for eigenvectors of the "x" dipole matrix
	double* muvalsy = new double[cis_size];                                                                                // Array for eigenvalues of the "y" dipole matrix
	double* muvexy = new double[cis_size*cis_size];                                                                        // Array for eigenvectors of the "y" dipole matrix
	double* muvalsz = new double[cis_size];                                                                                // Array for the eigenvalues of the "z" dipole matrix
	double* muvexz = new double[cis_size*cis_size];                                                                        // Array for the eigenvectors of the "z" dipole matrix
	double* corr = new double[cis_size];                                                                                   // Array for CIS(D) corrections to CIS energies
	bcs_read(cis_size, cisens, muvalsx, muvexx, muvalsy, muvexy, muvalsz, muvexz, corr, &bcsf);                            // read the previously defined arrays from the bcs file
	bcsf.close();                                                                                                          // close bcs file
	double offset = strtod(argv[2], NULL);                                                                                 // energy offset (added to each electronic state) in Hartree
	double maxen = strtod(argv[3], NULL);                                                                                  // maximal energy (higher states are cut off) in Hartree
	int cissel = calc_cissel(cis_size, cisens, corr, offset, maxen);                                                       // calculate number of selected CIS(D) states
	int* enarray = new int[cissel];                                                                                        // energy array containing the numbers of selected states
	double* bcsens = new double[cissel];                                                                                   // array with the reduced number of CIS(D) energies
	transform_energies(cis_size, cissel, enarray, offset, maxen, cisens, corr, bcsens);                                    // determine the elements of enarray and bcsens
	delete [] cisens; delete [] corr;                                                                                      // delete the arrays cisens and corr (not needed any more)
	double* bcsdipx = new double[cissel*cissel];                                                                           // array for the "x" dipole matrix
	transform_dipoles(cis_size, cissel, enarray, muvalsx, muvexx, bcsdipx);                                                // calculate "x" dipole matrix for all cis_size states and reduce it to the cissel*cissel matrix
	delete [] muvalsx; delete [] muvexx;                                                                                   // delete the arrays muvalsx and muvexx (not needed any more)
	double* bcsdipy = new double[cissel*cissel];                                                                           // repeat the prcedure for "y"
	transform_dipoles(cis_size, cissel, enarray, muvalsy, muvexy, bcsdipy);
	delete [] muvalsy; delete [] muvexy;
	double* bcsdipz = new double[cis_size*cis_size];                                                                       // repeat the procedure for "z"
	transform_dipoles(cis_size, cissel, enarray, muvalsz, muvexz, bcsdipz);
	delete [] muvalsz; delete [] muvexz;
	double* ionrat = new double[cissel];                                                                                   // array to store the selected ionization rates from the "new" states
	ifstream irxf;                                                                                                         // filestream to read ionization rates
	char finirx[128];                                                                                                      // filename with ionization rates
	sprintf(finirx, "%s", argv[4]);                                                                                        // store commandline argument 6 in finirx
	irxf.open(finirx);                                                                                                     // connect filestream irxf with filenume finirxf and open file
	read_ions(cis_size, cissel, enarray, ionrat, &irxf);                                                                   // read the ionization rates for all states and store the ones for the selected states in ionrat
	irxf.close();
	int cnt_nonion = count_nonion(cissel, ionrat);                                                                         // calculate the number of nonionizing states
	int cnt_ion = cissel - cnt_nonion;                                                                                     // calculate the number of ionizing states
	int* enssort = new int[cissel];                                                                                        // integer array to show the sequence of the "original" states
	presort(cissel, bcsens, ionrat, enssort, cnt_nonion, cnt_ion);                                                         // subroutine that sorts the "original" states in the integer array
	double* newens = new double[cissel];                                                                                   // "sorted" energy array
	sort_energies(cissel, enssort, bcsens, newens);                                                                        // sort the energies according to the order in the enssort array
	delete [] bcsens;                                                                                                      // delete the old energy array
	double* newdipx = new double[cissel*cissel];                                                                           // "sorted" "x" dipole array
	sort_dipoles(cissel, enssort, bcsdipx, newdipx);                                                                       // sort the dipole moments
	delete [] bcsdipx;                                                                                                     // delete the old dipole array
	double* newdipy = new double[cissel*cissel];                                                                           // repeat the procedure for y
	sort_dipoles(cissel, enssort, bcsdipy, newdipy);
	delete [] bcsdipy;
	double* newdipz = new double[cissel*cissel];                                                                           // repeat the procedure for z
	sort_dipoles(cissel, enssort, bcsdipz, newdipz);
	delete [] bcsdipz;
	double* newion = new double[cissel];                                                                                   // repeat the procedure for the ionization rates
	sort_energies(cissel, enssort, ionrat, newion);
	delete [] ionrat;
	char foutecp[128];                                                                                                     // filename of output ecp file
	double ip = strtod(argv[5], NULL);
	double* kin = new double[cissel];
	calc_kin(cissel, newens, newion, ip, kin);
	sprintf(foutecp, "%s.ecp", argv[6]);                                                                                   // store commandline argument 9 in foutecp
	ofstream outf;                                                                                                         // define output filestream
	outf.open(foutecp);                                                                                                    // connect filestream outf with filename foutecp and open the file
	write_ecp(&outf, cissel, newens, newdipx, newdipy, newdipz, newion, kin);                                              // subroutine to write out the generated data
	outf.close();                                                                                                          // close the file
    }
}

void ecp_read(int nros, double* inens, double* indipx, double* indipy, double* indipz, double* inion, double* inkin, ifstream* inf){
    inf->read((char *) inens, sizeof(double)*nros);
    inf->read((char *) indipx, sizeof(double)*nros*nros);
    inf->read((char *) indipy, sizeof(double)*nros*nros);
    inf->read((char *) indipz, sizeof(double)*nros*nros);
    inf->read((char *) inion, sizeof(double)*nros);
    inf->read((char *) inkin, sizeof(double)*nros);
}

void bcs_size(ifstream* inf, int &nroao, int &nroe, int &llim, int &ulim, int &omo, int &umo, int &cis_size){
    inf->read((char *) &nroao, sizeof(int));
    inf->read((char *) &nroe, sizeof(int));
    inf->read((char *) &llim, sizeof(int));
    inf->read((char *) &ulim, sizeof(int));
    omo = nroe/2 - llim;
    umo = ulim + 1 - omo - llim;
    cis_size = omo * umo + 1;
}

void bcs_read(int cis_size, double* cisens, double* muvalsx, double* muvexx, double* muvalsy, double* muvexy, double* muvalsz, double* muvexz, double* corr, ifstream* inf){
    inf->read((char *) cisens, sizeof(double)*cis_size);
    double* cisvex = new double[cis_size*cis_size];
    inf->read((char *) cisvex, sizeof(double)*cis_size*cis_size);
    delete [] cisvex;
    inf->read((char *) muvalsx, sizeof(double)*cis_size);
    inf->read((char *) muvexx, sizeof(double)*cis_size*cis_size);
    inf->read((char *) muvalsy, sizeof(double)*cis_size);
    inf->read((char *) muvexy, sizeof(double)*cis_size*cis_size);
    inf->read((char *) muvalsz, sizeof(double)*cis_size);
    inf->read((char *) muvexz, sizeof(double)*cis_size*cis_size);
    inf->read((char *) corr, sizeof(double)*cis_size);
}

int calc_cissel(int cis_size, double* cisens, double* corr, double offset, double maxens){
    double dumdum;
    int cissel = 0;
#pragma omp parallel for reduction(+:cissel)
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
#pragma omp parallel for
	    for(int i = 0; i < (cnt_nonion/2); i++){
		if(joinens[enssort[2*i+1]] < joinens[enssort[2*i]]){
		    tmp = enssort[2*i];
		    enssort[2*i] = enssort[2*i+1];
		    enssort[2*i+1] = tmp;
		}
	    }
	}else{
#pragma omp parallel for
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
#pragma omp parallel for
	    for(int i = (cnt_nonion/2); i < (joined/2); i++){
		if(joinens[enssort[2*i+1]] < joinens[enssort[2*i]]){
		    tmp = enssort[2*i];
		    enssort[2*i] = enssort[2*i+1];
		    enssort[2*i+1] = tmp;
		}
	    }
	}else{
#pragma omp parallel for
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

void write_ecp(ofstream* outf, int joined, double* ens, double* dipx, double* dipy, double* dipz, double* ion, double* kin){
    outf->write((char *) &joined, sizeof(int));
    outf->write((char *) ens, sizeof(double)*joined);
    outf->write((char *) dipx, sizeof(double)*joined*joined);
    outf->write((char *) dipy, sizeof(double)*joined*joined);
    outf->write((char *) dipz, sizeof(double)*joined*joined);
    outf->write((char *) ion, sizeof(double)*joined);
    outf->write((char *) kin, sizeof(double)*joined);
    outf->flush();
}

