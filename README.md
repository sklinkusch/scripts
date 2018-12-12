# scripts

This is a small description of programs within this folder

## data

data files for the scripts in the fahrinfo and finfo folders

- **fahrinfo-elinks.dat**: numbers of stations and bus stops (station names contain Umlaute)
- **fahrinfo-elinks2.dat**: same as _fahrinfo-elinks.dat_, but Umlaute are transformed (ä → ae, ö → oe, ü → ue, ß → ss)

## fahrinfo

programs to read the data from fahrinfo.vbb.de, manipulate and display them in a perl-Tk window

- **fahrinfo-dual.pl**: read arrivals and departure of a station/bus stop
- **fahrinfo-duo-sort.pl**: read departures or arrivals of two stations/bus stops
- **fahrinfo-duo-sort-filter.pl**: read departures or arrivals of two stations/bus stops and filter them
- **fahrinfo-ffilter-ddual.pl**: read departures or arrivals of a station/bus stop and filter them according to two different schemes
- **fahrinfo-ffilter-dual.pl**: read departures and arrivals of a station/bus stop and filter them according to two different schemes
- **fahrinfo-filter-dual.pl**: read departures and arrivals of a station/bus stop and filter them
- **fahrinfo-filter.pl**: read departures or arrivals of a station/bus stop and filter them
- **fahrinfo-fstrict-dual.pl**: read departures and arrivals of a station/bus stop (w/o neighboring stops)
- **fahrinfo-fstrict-duo.pl**: read departures or arrivals of two different stations/bus stops (w/o neighboring stops)
- **fahrinfo-fstrict-filter.pl**: read departures or arrivals of a station/bus stop (w/o neighboring stops) and filter them
- **fahrinfo-fstrict.pl**: read departures or arrivals of a station/bus stop (w/o neighboring stops)
- **fahrinfo.pl**: read departures or arrivals of a station/bus stop
- **fahrinfo-strict.pl**: read departures or arrivals of a station/bus stop (with information on neighboring stops/track)

further programs

- **fahrinfo-dlinks.sh**: print timetable of fahrinfo.vbb.de into the command line (usually in less)
- **fahrinfo-elinks.sh**: show timetable of fahrinfo.vbb.de in elinks
- **fahrinfo-elinksd.sh**: print timetable of fahrinfo.vbb.de into the command line (not in less)
- **fahrinfo-w3m.sh**: show timetable of fahrinfo.vbb.de in w3m

data file

- **Fahrinfo.param**: including the update interval in seconds

routines file

- **Fahrinfo.pm**: routines for the fahrinfo.pl program suite

## finfo

programs to read the data from fahrinfo.vbb, manipulate and display them in the command line

- **finfo-dual.pl**: read arrivals and departure of a station/bus stop
- **finfo-duo-sort.pl**: read departures or arrivals of two stations/bus stops
- **finfo-duo-sort-filter.pl**: read departures or arrivals of two stations/bus stops and filter them
- **finfo-ffilter-ddual.pl**: read departures or arrivals of a station/bus stop and filter them according to two different schemes
- **finfo-ffilter-dual.pl**: read departures and arrivals of a station/bus stop and filter them according to two different schemes
- **finfo-filter-dual.pl**: read departures and arrivals of a station/bus stop and filter them
- **finfo-filter.pl**: read departures or arrivals of a station/bus stop and filter them
- **finfo-fstrict-dual.pl**: read departures and arrivals of a station/bus stop (w/o neighboring stops)
- **finfo-fstrict-duo.pl**: read departures or arrivals of two different stations/bus stops (w/o neighboring stops)
- **finfo-fstrict-filter.pl**: read departures or arrivals of a station/bus stop (w/o neighboring stops) and filter them
- **finfo-fstrict.pl**: read departures or arrivals of a station/bus stop (w/o neighboring stops)
- **finfo.pl**: read departures or arrivals of a station/bus stop
- **finfo-strict.pl**: read departures or arrivals of a station/bus stop (with information on neighboring stops/track)

shell scripts to include the previously mentioned Perl scripts into a _watch_ command

- **finfo-dual.sh**
- **finfo-duo-sort.sh**
- **finfo-duo-sort-filter.sh**
- **finfo-ffilter-ddual.sh**
- **finfo-ffilter-dual.sh**
- **finfo-filter-dual.sh**
- **finfo-filter.sh**
- **finfo-fstrict-dual.sh**
- **finfo-fstrict-duo.sh**
- **finfo-fstrict-filter.sh**
- **finfo-fstrict.sh**
- **finfo.sh**
- **finfo-strict.sh**

## further-programs

some unrelated, but useful Perl scripts

- **average.pl**: calculate the average of at least one value (more is better)
- **bs-filter.pl**: take the _browser-sync_ output and print only the information about IP addresses, informations about refreshing the page are swallowed
- **margin.pl**: calculate the percentage of margin if only the two items with the highest percentage are regarded

some unrelated, but useful C++ programs (source files are included - .cc)

- **cardano**: calculate the roots of cubic equations
- **dewpoint**: calculate the dewpoint according to temperature and humidity of air
- **lunar**: calculate the current lunar phase
- **ostern**: calculate the (German) holidays for a certain year (weekdays and dates)
- **sofortrente**: calculate the remaining balance for a bank-hosted pension funds
- **sun**: calculate sunrise, sunset, and twilight for a certain day
- **sun2**: calculate sunrise, sunset, and twilight for a certain day (better algorithm)
- **sun3**: calculate the angle (azimuth and elevation) the sunbeam hits the earth

some shell scripts

- **bsf.sh**: execute _browser-sync_ and pipe the output to _bs-filter.pl_

### calendar

writes a calendar tex file with some dates marked

- **red**: Sundays and public holidays
- **grey**: Saturdays, Christmas Eve and New Year's Eve
- **other colors**: specified in certain files  
  Keep in mind, that no date can be marked with more than one color.

### Dok

includes a pdf file with the source code of the _sofortrente_ program

### Entsorgung

calculates regular dates for the disposal of different things (paper, plastics, remaining waste, Christmas trees, glass, ...); the
dates can be used within the _calendar_ program to make a calendar

### Sudoku

can solve Sudokus of arbitrary difficulty in a fast trial-and-error algorithm

### TV-Masse

a Perl script than can calculate the width and height of a television device according to its diagonal in inch or cm

# TC-programs

programs intended for the use in theoretical chemistry

## CIS-CSF

prepatory program for stochasticTDCI based on CIS calculations

- **istochcis**: reads a binary sys file and a binary dat.hfw or a binary hwf file ⇒ binary sys file
- **read_bcs**: reads the binary bcs file and writes data out to the shell
- **read_wav**: reads the binary wav file and writes MO energies and coefficients to an output file

## CIS-D

calculation of the doubles correction on top of configuration interaction singles (CIS) calculations

- **calc_norm**: calculates the norm of a certain electronic system
- **cis-d**: calculates the doubles correction within CIS(D) calculations
- **cis-d-orig**: older version of _cis-d_
- **motwoelint**: calculate the molecular two-electron integrals ⇒ binary moint file
- **readbin**: read the binary bin-file containing ground state energy and excitation energies ⇒ arbitrary out file
- **readint**: read the non-zero integrals from file ⇒ out file
- **readlog**: reads number of states, RHF energy, and CIS excitation energies from GAMESS log file ⇒ bin.dat file
- **readmat**: reads CI eigenenergies and CI matrix from a binary GAMESS matrix file
- **readsys**: reads number of atomic orbitals, of atoms, of integrals, further coordinates, charges, masses, Hamiltonian matrix, kinetic energy matrix, overlap matrix, dipole moments along _x_, _y_, and _z_, two-electron integral indices and values
- **readwav**: reads MO energies and vectors from a binary dat.hwf file ⇒ log file
- **renorm**: calculates the norm of a certain electronic system and renormalizes it to one ⇒ trnorm file
- **rgw**: reads MO energies and vectors from a GAMESS punch file ⇒ binary dat.hwf file
- **writemat**: writes a binary CI matrix file from a human readable ASCII file

## DFT-CSF

prepatory programs for stochasticTDCI based on DFT calculations

- **iirxdft2**: calculates ionization rates for each electronic state ⇒ irx file
- **irhodft**: older version of _irhodft_
- **irhodft2**: reads three binary files: sys, dat.hwf and vec ⇒ binary bcs file and human readable irx file
- **istochdft**: reads three binary files: sys, dat.hwf and vec ⇒ binary sys file and binary bcs file
- **read_bcs**: reads the binary bcs file and writes data out to the shell
- **read_wav**: reads the binary wav file and writes MO energies and coefficients to an output file

## DMP

programs that perform a ρ-TDCI calculation

- **rhoTDCI.bin**: reads a binary ecp file and (maybe) a binary rst file ⇒ new binary rst file and several human readable files: log, pop, fld, dip, enpop
- **rhoTDCI.x**: older version of _rhoTDCI.bin_ that reads data from ASCII files instead of the binary ecp file

## GRIMME-CIDFT

performs calculation according to an algorithm by Stefan Grimme

- **grimme**
- **grimme-20**

## IONRATS

a program calculating the ionzation rates for CIS or CIS(D) states

- **ionrats**: reading data from human readable ens and vex file as well as from the binary dat.hwf file ⇒ human readable irx file

## OLD_CIS3D

the original CIS3D program suite by Tillmann Klamroth and my extensions

- **fft** and **fft_n**: perform a Fourier transformation of a signal given in an ASCII file
- **mk_ff**: reading electric field data from a binary file and performs a Fourier transformation
- **norm**: reading populations during a laser pulse excitation and calculates the norm as a function of time
- **quadfkt**: reads (experimental) values and performs two types of regression: linear and quadratical
- **readwrite_wav**: reads time-dependent wave function and writes it out for a certain timestep
- **sep**: separates specified degenerate states from a population file

### CID

- **cid**: performs a configuration interaction doubles (CID) calculation

### CIS

- **read_bcs**: reading energies and dipole moments from a binary bcs file and writes data to stdout
- **read_bcs-new**: same as read_bcs, but bcs file contains only a subset of states
- **read_corr**: reads the correlation energy from a bcs file and prints it to stdout
- **read_wav**: reads MO energies and coefficients from a binary file and writes them to a log file

### CISD-OCT

programs performing laser pulse optimizations according to the _Optimal Control Theory_, based on a CIS or CISD calculation:

- **cis_oct**: OCT calculation, based on CIS
- **mk_field**: make an electric field for the OCT calculation to start with
- **mk_init_wav**: make an initial CIS wave function
- **mk_init_wav-cisd**: make an initial CISD wave function

### GAMMA-CSFS

programs to perform CIS calculations with an ionization scheme based on configuration state functions

- **cis**: make a CIS calculation
- **mk_field**: make an initial field
- **mk_in**: read integrals from a GAMESS calculation
- **mk_init_wav**: make initial wave function
- **rgw**: read the molecular orbitals' data from an ASCII file and write it to a binary file
- **rhf**: perform a restricted Hartree-Fock calculation to get a binary file containing the molecular orbitals' data
- **tdcis**: propagate the CIS wave function in time within a possible laser field
- **tdpop**: time-dependent populations
- **tdrhf**: propagate the HF wave function in time

### GAMMA-OCT

programs to perform OCT calculations, based on a CIS calculation including the photoionization algorithm according to S. Klinkusch, P. Saalfrank, and T. Klamroth, _The Journal of Chemical Physics_ **131**, 114304 (2009)

- **cis_oct**: main program
- **mk_ef**: read binary electronic field file
- **mk_field**: generate binary file containing the initial field
- **mk_init_wav**: generate binary file containing the initial wave function

## stochasticTDCI

## STUFF

## VOOCIS

## WRITEBCS
