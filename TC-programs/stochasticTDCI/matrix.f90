!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
module matrix

  implicit none
  private
  integer,public :: Nstates,Nbasis,Nbm1
  integer :: Nnz,mapdip(3)
  integer,allocatable :: ind(:)
  integer(kind=8),allocatable :: row_ptr(:),col_ind(:)
  real(kind=4),parameter :: zerof = 0.0
  real(kind=8),parameter :: zero = 0d0,half = 0.5d0,one = 1d0,two = 2d0
  real(kind=8),allocatable :: dipole(:,:)                                     ! the dipole matrices in CRS
  real(kind=8),allocatable :: rates(:,:),dephasing(:,:),chi(:)                ! the transition rates matrix
  real(kind=8),public,allocatable :: energy(:),configurations(:)              ! energies: states,configs
  real(kind=8),public,allocatable :: eigenvectors(:,:)                        ! eigenvectors for transformation
  real(kind=8),public,allocatable :: ionization(:)                            ! ionization rates
  real(kind=8),public,allocatable :: Et(:)                                    ! local time field
  integer,public :: Npol
  complex(kind=8),parameter :: zeroc = dcmplx(0d0,0d0),icmplx=dcmplx(0d0,1d0),onec=dcmplx(1d0,0d0)

  real(kind=8),public :: Tfinal
  public matvecI,normalize,initialize,writePOP,cdipole

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  subroutine initialize(restart,polarization,N_traj,psiINIT,temperature)

    implicit none
    logical :: restart
    character(len=3),intent(in) :: polarization
    integer,intent(in) :: N_traj
    real(kind=8),intent(in) :: Temperature
    complex(kind=8),intent(out) :: psiINIT(0:Nbasis,N_traj)
    integer :: i,j,idum,NstatesIN,ipol,io
    real(kind=8),parameter :: tinyFLOAT = 1d-10
    real(kind=8),parameter :: Knst = 3.1577465d5                ! K/hartree
    real(kind=8),parameter :: fs2au = 41.341373337d0
    real(kind=8),parameter :: eps = 1d-10
    real(kind=8) :: dum,bose_einstein,deltaE

! convert time to atomic units
    Tfinal = Tfinal*fs2au

! consistency check
    read(2) idum,NstatesIN
    !sk 2015.07.09
    !sk end 2015.07.09
    if(idum /= Nbasis) stop 'CSF basis has the wrong size'
    if(NstatesIN > Nbasis .or. Nstates > NstatesIN) stop 'Eigenstate basis too big'

! read the energies of the configurations, as well as their ionization rates
    allocate(energy(Nstates),eigenvectors(0:Nbm1,Nstates),configurations(0:Nbm1),ionization(0:nbm1))
    read(2) configurations
    read(2) ionization

! read the dipole matrix elements
    select case (polarization)
      case('xyz') 
        Npol = 3
        mapdip(1) = 1
        mapdip(2) = 2
        mapdip(3) = 3
      case('xy0') 
        Npol = 2
        mapdip(1) = 1
        mapdip(2) = 2
      case('x0z') 
        Npol = 2
        mapdip(1) = 1
        mapdip(2) = 3
      case('0yz') 
        Npol = 2
        mapdip(1) = 2
        mapdip(2) = 3
      case('x00') 
        Npol = 1
        mapdip(1) = 1
      case('0y0') 
        Npol = 1
        mapdip(1) = 2
      case('00z') 
        Npol = 1
        mapdip(1) = 3
      case('000') 
        Npol = 0
      case default
        stop 'Polarization keyword is soooo WRONG'
    end select

! allocate and initialize
    allocate(rates(0:Nbm1,0:Nbm1),dephasing(0:Nbm1,0:Nbm1),chi(0:nbm1),Et(npol))
    psiINIT(:,:) = zero
    
! store the components of dipole: using Compressed Row Storage
!jct 2015.09.02
    allocate(row_ptr(0:Nbasis))
! read mapping functions
    do i = 0,Nbasis; read(2) row_ptr(i);  enddo
    Nnz = row_ptr(Nbasis)-1
    allocate(dipole(0:Nnz,3),col_ind(0:Nnz))
    do i=0,Nnz;      read(2) col_ind(i);  enddo
! read x component
    do i=0,Nnz;      read(2) dipole(i,1); enddo
! read y component
    do i=0,Nnz;      read(2) dipole(i,2); enddo
! read z component
    do i=0,Nnz;      read(2) dipole(i,3); enddo

! read relaxation and dephasing rates
    rates(:,:) = zero
    read(2) rates
    dephasing(:,:) = zero
    read(2) dephasing

! thermalize rates
    dum = Knst/max(temperature,tinyfloat)
    do i = 0,Nbm1
    do j = i+1,Nbm1
    !sk 2015.07.01
      deltaE = abs(configurations(i)-configurations(j))
      if(temperature < eps .AND. deltaE >= zero) then
          bose_einstein = zero
      elseif(temperature < eps .AND. deltaE < zero) then
          bose_einstein = -one
      else
          bose_einstein = one/(exp(dum*deltaE)-one)
      endif
      !sk end 2015.07.01
      rates(i,j) = rates(i,j)*(bose_einstein)
      rates(j,i) = rates(j,i)*(one+bose_einstein)
    enddo
    enddo

! auxiliary vector for accelerating MV: summing rates
    do i = 0,nbm1
      dum = zero
      do j = 0,nbm1
        dum = dum + rates(i,j)
      enddo
      chi(i) = half*(dum+ionization(i))
    enddo

! read eigenvectors and eigenvalues defining the interaction picture
    read(2) (energy(i),i=1,Nstates)
    read(2) (eigenvectors(i,1),i=Nstates+1,NstatesIN) ! place reader to correct position
    do j=1,Nstates
      read(2) eigenvectors(:,j)
    enddo

! shift energies
    dum = minval(energy)
    energy(:) = energy(:) - dum

! generate starting density matrix
    if(restart) then
! restart file
      open(8,file='psi.rst',status='old',form='unformatted',access='stream',iostat=io)
      read(8,iostat=io) idum,i
      if(io == 0 .and. idum == nbasis .and. i == N_traj) then
        do i = 1,N_traj
          read(8,iostat=io) psiINIT(:,i)
        enddo
        if(io /= 0) stop 'problem reading the psi.rst file'
      else
        write(1,*) 'WARNING: wrong restart file. initializing density to the ground state'
        do i = 1,N_traj
          psiINIT(0,i) = onec
          psiINIT(1:Nbasis,i) = zeroc
        enddo
      endif
      close(8)
    else
! initialize in ground state
      do i = 1,N_traj
!!! debug begin !
!         psiINIT(0:Nbasis,i) = zeroc
!         psiINIT(5,i) = onec
!!! debug end !!!
        psiINIT(0,i) = onec
        psiINIT(1:Nbasis,i) = zeroc
      enddo
    endif

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  subroutine matvecI(T,psiIN,psiOUT,N_traj)

    use field, only: ElectricField
    implicit none
    integer,intent(in) :: N_traj
    real(kind=8),intent(in) :: T
    complex(kind=8),intent(in) :: psiIN(0:Nbasis,N_traj)
    complex(kind=8),intent(out) :: psiOUT(0:Nbasis,N_traj)

    integer :: j,m,mp,ipol,i
    complex(kind=8) :: dumc,dumd,psiS(0:nbasis),tmp(0:nbasis)

! field intensity at time T
    call ElectricField(T,Et,npol)

! for each trajectory 
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,ipol,m,mp,dumc,psiS,tmp)
    do j = 1,N_traj

!! conversion to the Schroedinger picture
    psiOUT(:,j) = zeroc
    call convertIntSchroe(psiIN(0,j),psiS,-T)

! coupling to the electric field (CRS)
      tmp(:) = zeroc
      do ipol = 1,Npol
        do m = 0,nbm1
          dumc = zeroc
          do mp = row_ptr(m),row_ptr(m+1)-1
            dumc = dumc + dipole(mp,mapdip(ipol))*psiS(col_ind(mp))
          enddo
          tmp(m) = tmp(m) + icmplx*Et(ipol)*dumc
        enddo
      enddo
      tmp(Nbasis) = psiS(Nbasis)

! conversion to the Interaction picture
      call convertIntSchroe(tmp,psiOUT(0,j),T)

! matrix-vector product w/o laser field: dissipation
      do m = 0,Nbm1
        psiOUT(m,j) = psiOUT(m,j) - chi(m)*psiIN(m,j)
      enddo
      psiOUT(Nbasis,j) = zeroc
    enddo
!$OMP END PARALLEL DO

  end subroutine matvecI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  subroutine convertIntSchroe(vecIN,vecOUT,time)
    implicit none
    real(kind=8),intent(in) :: time
    complex(kind=8),intent(in) :: vecIN(0:Nbasis)
    complex(kind=8),intent(out) :: vecOUT(0:Nbasis)
    integer :: m,jvec
    complex(kind=8) :: dumc,tmp(0:Nstates)

    do jvec = 1,Nstates
      dumc = zeroc
      do m = 0,Nbm1
        dumc = dumc + eigenvectors(m,jvec)*vecIN(m)
      enddo
      tmp(jvec) = exp(icmplx*time*energy(jvec))*dumc
    enddo

    vecOUT(:) = zeroc
    do jvec = 1,Nstates
      do m = 0,Nbm1
        vecOUT(m) = vecOUT(m) + eigenvectors(m,jvec)*tmp(jvec)
      enddo
    enddo

! continuum state does not evolve
    vecOUT(Nbasis) = vecIN(Nbasis)

  end subroutine convertIntSchroe

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  subroutine normalize(N_traj,psiIN,psiOUT)

    implicit none
    integer,intent(in) :: N_traj
    complex(kind=8),intent(in) :: psiIN(0:Nbasis,N_traj)
    complex(kind=8),intent(out) :: psiOUT(0:Nbasis,N_traj)

    integer :: j,m,mp,nseed
    integer :: seed(2*N_traj)
    real(kind=8) :: epsilon0(2,N_traj),prob(0:Nbasis),cumul
    real(kind=8) :: norme,loss

! generate uniformly distributed random numbers
    nseed = 2*N_traj
    do j = 1,nseed
      call system_clock(count=m)
      seed(j) = modulo(min(m,m*19937528),214748)+j
!      seed(j) = modulo(m*(j+nseed)+10983651,1111111111)
    enddo
    call random_seed(size = nseed)
    call random_seed(put = seed(1:nseed))
    call random_number(epsilon0)

    do j = 1,N_traj

! compute the loss of norm
      Norme = zero
      do m = 0,Nbasis
        Norme = Norme + dble(psiIN(m,j))**2+dimag(psiIN(m,j))**2
      enddo
      loss = one - Norme

      if(loss > epsilon0(1,j)) then
! Stochastic jump: collapse to a zeroth-order state
! 1) compute probabilities
        norme = zero
        do m = 0,Nbm1
          prob(m) = zero
          do mp = 0,Nbm1
            prob(m) = prob(m) + rates(mp,m)*conjg(psiIN(mp,j))*psiIN(mp,j)
          enddo
          norme = norme + prob(m)
        enddo
! ionization probability 
        prob(Nbasis) = zero
        do mp = 0,Nbm1
            prob(Nbasis) = prob(Nbasis) + ionization(mp)*conjg(psiIN(mp,j))*psiIN(mp,j)
        enddo
! normalize probability
        norme = norme + prob(Nbasis)
! renormalize probability
        prob(:) = prob(:)/norme
! 2) transfer the WP to the most probable state and normalize
        mp = 0
        cumul = prob(0)
        findmax : do
          if(mp == Nbasis) exit findmax
          if(cumul >= epsilon0(2,j)) exit findmax
          mp = mp + 1
          cumul = cumul + prob(mp)
        enddo findmax
        psiOUT(:,j) = zeroc
        psiOUT(mp,j) = one

      else
! no jump occured: renormalize the WF
        norme = sqrt(Norme)
        do m = 0,Nbasis
          psiOUT(m,j) = psiIN(m,j)/norme
        enddo
      endif

    enddo

  end subroutine normalize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  subroutine writePOP(time,psi,err,N_traj,Nwrite)

    implicit none
    integer,intent(in) :: N_traj,Nwrite 
    real(kind=8),intent(in) :: time,err(N_traj)
    complex(kind=8),intent(in) :: psi(0:Nbasis,N_traj)
    integer :: i,j,m,jvec
    real(kind=4) :: population(0:Nbasis,N_traj),average(0:Nbasis),rdum
    real(kind=4) :: statepop(0:Nstates)
    complex(kind=8) :: dumc

! compute state population for individual trajectories
    rdum = zero
    do j = 1,N_traj
      do i = 0,Nbasis
        population(i,j) = dble(psi(i,j))**2+dimag(psi(i,j))**2
        rdum = rdum + population(i,j)
      enddo
    enddo

    statepop(:) = zerof
    do jvec = 1,Nstates
      do j = 1, N_traj
        dumc = zeroc
        do m = 0,Nbm1
          dumc = dumc + eigenvectors(m,jvec)*psi(m,j)
        enddo
        statepop(jvec) =  statepop(jvec)+real(dumc*conjg(dumc))
      enddo
    enddo
    statepop(:) = statepop(:)/real(N_traj)
    do jvec = 1,Nstates
      statepop(jvec-1) = statepop(jvec)
    enddo



! norm of stochastic density matrix
    write(1,'(f15.2,f10.6,e15.4)') time,rdum/real(N_traj),maxval(err)
    call flush(1)

! average population
    average(:) = zero
    do j = 1,N_traj
      do i = 0,Nbasis
        average(i) = average(i) + population(i,j)
      enddo
    enddo
    average(:) = average(:)/real(N_traj)
    statepop(Nstates) = average(Nbasis)

!    write(8) time,psi ; call flush(8)       ! output state vector
    write(9,*) time,population(:,1:Nwrite)   ! output individual populations
    write(10,*) time,average                 ! output averaged population
    write(12,*) time,statepop

  end subroutine writePOP

  function cdipole(i,j,ipol)
  
  implicit none
  integer,intent(in) :: i, j, ipol
  integer :: curr_row, next_row, x
  real(kind=8) :: cdipole

  curr_row = row_ptr(i)
  next_row = row_ptr(i+1)
  cdipole = zero
  do x=curr_row, next_row-1
    if(j .EQ. col_ind(x)) then
      cdipole = dipole(x, ipol)
      exit
    endif
  enddo
  return
  end function cdipole
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
end module matrix
