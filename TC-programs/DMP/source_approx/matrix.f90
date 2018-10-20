!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
module matrix

  implicit none
  private
  integer :: NfieldIN,Npol
  integer,allocatable :: map_dip(:)
  logical,allocatable :: deltaK(:,:)
  logical,public :: readfield
  real(kind=8):: f0,carrier,width,tpeak
  real(kind=8),parameter :: zero = 0d0,half = 0.5d0,one = 1d0, two = 2d0
  real(kind=8),parameter :: mtwo = -2d0, mone = -1d0
  real(kind=8),allocatable :: dipole(:,:,:),rates(:,:),dephasing(:,:)              ! the dipole and the transition rate matrices
  real(kind=8),allocatable :: Tin(:),field(:,:),coeff_lin(:,:),coeff_quad(:,:),coeff_quart(:,:),fieldQUASI(:,:)
  complex(kind=8),allocatable :: energies(:), trans_energy(:,:)     ! the transition energies
  complex(kind=8),public,parameter :: zeroc = dcmplx(0d0,0d0),icmplx=dcmplx(0d0,1d0), halfc=dcmplx(0.5d0,0d0)

  integer,public :: Nbasis,Nsys,Nbm1,Nsm1
  public matvec,matvecI,writePOP,initialize,trpesrk,trpesa,writekin

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  subroutine initialize(restart,Frst,polarization,rho,temperature,scaleION,Ffield,Tinitial,Tfinal,freq,amplitude)

    use scalarOP
    implicit none
    logical :: restart
    character(len=3),intent(in) :: polarization
    character(len=40),intent(in) :: Ffield, Frst
    real(kind=8),intent(in) :: Temperature,Tinitial,Tfinal,amplitude,freq,scaleION
    complex(kind=8),intent(out) :: rho(0:Nbm1,0:Nsys)
    integer :: i,j,idum,jdum,ipol,io
    real(kind=8),parameter :: Knst = 3.1577465d5
    real(kind=8) :: dum(3),dum1,dum2,dum3,dum4,dum5
    real(kind=8) :: tmp(0:Nbm1),pop(0:Nbm1)
    complex(kind=8) :: dumc

! read the energies and evaluate the exponential transition energies
    allocate(fieldQUASI(0:nbm1,0:Nbm1),energies(0:Nbm1),deltaK(0:nbm1,0:Nbm1))
    fieldQUASI(:,:) = zeroc
    deltaK(:,:) = .false.
    read(2,*)
    do i = 0,Nbm1
      read(2,*,iostat=io) idum,dum1,dum2,dum3
      if(io /= 0) stop 'problem while reading energies of the subsystem'
      energies(idum) = dcmplx(dum1,-scaleION*dum2)
      pop(idum) = dum3
    enddo
    allocate(trans_energy(0:Nbm1,0:Nsys))
    do j = 0,Nsm1
      do i = j,Nbm1
        trans_energy(i,j) = energies(i)-conjg(energies(j))
        idum = nint(dble(trans_energy(i,j))/freq)
        fieldQUASI(i,j) = dble(trans_energy(i,j))-idum*freq
        fieldQUASI(j,i) = fieldQUASI(i,j)
        if(abs(idum) == 1) then
          deltaK(i,j) = .true.
          deltaK(i,j) = .true.
        endif
      enddo
    enddo
    do i = Nsys,Nbm1
      trans_energy(i,Nsys) = energies(i)-conjg(energies(i))
    enddo

! read the input density and the dipole matrix elements
    
    select case (polarization)
      case('xyz') 
        Npol = 3
      case('xy0','x0z','0yz') 
        Npol = 2
      case('x00','0y0','00z') 
        Npol = 1
      case('000') 
        Npol = 0
      case default
        stop 'Polarization keyword WRONG'
    end select

! initialize matrices
    allocate(map_dip(Npol),dipole(3,0:Nbm1,0:Nbm1),rates(0:Nbm1,0:Nbm1),dephasing(0:Nbm1,0:Nbm1))
    do j = 0,Nsm1
      do i = j,Nbm1
        rho(i,j) = zeroc
      enddo
    enddo
    do j = Nsys,Nbm1
      rho(j,Nsys) =  zeroc
    enddo
    do j = 0,Nbm1
      do i = 0,Nbm1
        do idum = 1,3
          dipole(idum,i,j) = zero
        enddo
         rates(i,j) = zero
      enddo
    enddo

! create maping for dipole moment
    map_dip(:) = 0
    ipol = 0
    if(polarization(1:1) == 'x') then
      ipol = ipol + 1
      map_dip(ipol) = 1
    endif
    if(polarization(2:2) == 'y') then
      ipol = ipol + 1
      map_dip(ipol) = 2
    endif
    if(polarization(3:3) == 'z') then
      ipol = ipol + 1
      map_dip(ipol) = 3
    endif

    inmu: do
      read(3,*,iostat=io) jdum,idum,dum1,dum2,dum3,dum4,dum5
      if(io /=0 .or. idum > nbm1) exit inmu
! store the components of dipole
      dipole(1,idum,jdum) = dum1; dipole(1,jdum,idum) = dum1
      dipole(2,idum,jdum) = dum2; dipole(2,jdum,idum) = dum2
      dipole(3,idum,jdum) = dum3; dipole(3,jdum,idum) = dum3
      rates(idum,jdum) = dum4
      dephasing(idum,jdum) = dum5; dephasing(jdum,idum) = dum5
      if(temperature > 1e-10)  rates(jdum,idum) = dum4*exp(-Knst*dble(energies(idum)-conjg(energies(jdum)))/temperature)
    enddo inmu

! compute rates
    do j = 0,Nbm1
      do i = 0,Nbm1
        dum1 = zero
        do idum = 0,Nbm1
          dum1 = dum1 + rates(i,idum) + rates(j,idum)
        enddo
        dephasing(i,j) = dephasing(i,j) + half*dum1
      enddo
    enddo

! generate starting density matrix
    if(restart) then
! restart file
      open(11,file=Frst,status='old',form='unformatted')
      read(11,iostat=io) rho  !((rho(i,j), i=0,Nbm1), j=0,Nsys)  ! (0:Nbm1,0:Nsys)
      if(io /= 0) stop 'problem reading the rst-file '
      close(11)
    else
! selected density matrix
      if(temperature < 1e-10) then
        do i = 0,Nbm1
          tmp(i) = sqrt(pop(i))
        enddo
      else
! thermalize density
        do i = 0,Nbm1
          tmp(i) = sqrt(exp(-Knst*dble(trans_energy(i,0))/temperature))
        enddo
      endif
      do j = 0,Nsm1
        do i = j,Nbm1
          rho(i,j) = tmp(i)*tmp(j)
        enddo
      enddo
      do j = Nsys,Nbm1
        rho(j,Nsys) = tmp(j)*tmp(j)
      enddo
! normalize density
      dumc = trace(rho,Nbasis,Nsys)
      call normalize(dumc,rho,Nbasis,Nsys)
    endif

    if(readfield) then
      open(4,file=Ffield,status='old',iostat=io)
! read the field
      read(4,*) NfieldIN
      allocate(Tin(NfieldIN),field(NfieldIN,Npol),coeff_lin(NfieldIN,Npol),coeff_quad(NfieldIN,Npol),coeff_quart(NfieldIN,Npol))
      do i = 1,NfieldIN
        read(4,*) dum1,(dum(j),j=1,ipol)
        Tin(i) = dum1
! store the components of the electric field
          ipol = 0
          if(polarization(1:1) == 'x') then
            ipol = ipol + 1
            field(i,ipol) = dum(ipol)
          endif
          if(polarization(2:2) == 'y') then
            ipol = ipol + 1
            field(i,ipol) = dum(ipol)
          endif
          if(polarization(3:3) == 'z') then
            ipol = ipol + 1
            field(i,ipol) = dum(ipol)
          endif
      enddo
      if(abs(Tfinal-Tin(NfieldIN))/tfinal > 1d-7) stop 'final time and field mismatch'
      write(1,*) 'Fitting the field components using cubic splines'
      call flush(1)
      do ipol = 1,Npol
        call spline(NfieldIN,Tin,field(1,ipol),coeff_lin(1,ipol),coeff_quad(1,ipol),coeff_quart(1,ipol))
      enddo
      close(4)
! to prevent numerical explosion in the polarizability calculation
      carrier = 1d300
    else
      carrier = freq
      f0 = amplitude
      width = Tfinal-Tinitial
      tpeak = 0.5d0*(Tfinal-Tinitial)
    endif
! evaluate memory requirement
   write(1,*) 'Memory requirement:',(5*Nbasis**2+22*Nbasis*(Nsys+1)+(4*npol+1)*NfieldIN)*8./(1024.**3), 'Gb'

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  subroutine writePOP(interaction,time,rhoIN)

    implicit none
    logical,intent(in) :: interaction
    real(kind=8),intent(in) :: time
    complex(kind=8),intent(in) :: rhoIN(0:Nbm1,0:Nsys)
    integer :: m,mp,k,l,j,ipol
    real(kind=4) :: population(0:Nbasis)
    real(kind=8) :: time_tau,norm
    real(kind=8) :: expectDIP(3),IM_pol(3)
    complex(kind=8) :: rhoS(0:Nbm1,0:Nsys),tmp(0:Nbm1,0:Nsys)

    if(interaction) then
! convert density to Schroedinger representation
!$OMP PARALLEL DO
      do mp = 0,Nsm1
      do m = mp,Nbm1
        rhoS(m,mp) = exp(-icmplx*dble(trans_energy(m,mp))*time)*rhoIN(m,mp)
      enddo
      enddo
      !$OMP END PARALLEL DO
    else
        !$OMP PARALLEL DO
      do mp = 0,Nsm1
      do m = mp,Nbm1
        rhoS(m,mp) = rhoIN(m,mp)
      enddo
      enddo
      !$OMP END PARALLEL DO
    endif
    !$OMP PARALLEL DO
    do mp = Nsys,Nbm1
      rhoS(mp,Nsys) = rhoIN(mp,Nsys)
    enddo
    !$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute population in states below ionization threshold
    norm = zero
    !$OMP PARALLEL DO REDUCTION(+:norm)
    do m = 0,Nsm1
      population(m) = rhoS(m,m)
      norm = norm + population(m)
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO REDUCTION(+:norm)
    do m = Nsys,Nbm1
      population(m) = rhoS(m,Nsys)
      norm = norm + population(m)
    enddo
    !$OMP END PARALLEL DO
    population(Nbasis) = 1e0 - norm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! COMPUTE dipole expectation value

    expectDIP(:) = zero
    do mp = 0,Nsm1
      do ipol = 1,3
        expectDIP(ipol) = expectDIP(ipol) + dipole(ipol,mp,mp)*dble(rhoS(mp,mp))
      enddo
      do m = mp+1,Nbm1
        do ipol = 1,3
          expectDIP(ipol) = expectDIP(ipol) + two*dipole(ipol,m,mp)*dble(rhoS(m,mp))
        enddo
      enddo
    enddo
    do mp = Nsys,Nbm1
      do ipol = 1,3
        expectDIP(ipol) = expectDIP(ipol) + dipole(ipol,mp,mp)*dble(tmp(mp,Nsys))
      enddo
    enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute the imaginary part of the polarizability
    time_tau = time + one/carrier
    call matvec(time_tau,rhoS,tmp)
    IM_pol(:) = zero
    do mp = 0,Nsm1
      do ipol = 1,3
        IM_pol(ipol) = IM_pol(ipol) + dipole(ipol,mp,mp)*dble(tmp(mp,mp))
      enddo
      do m = mp+1,Nbm1
        do ipol = 1,3
          IM_pol(ipol) = IM_pol(ipol) + two*dipole(ipol,m,mp)*dble(tmp(m,mp))
        enddo
      enddo
    enddo
    do mp = Nsys,Nbm1
      do ipol = 1,3
        IM_pol(ipol) = IM_pol(ipol) + dipole(ipol,mp,mp)*dble(tmp(mp,Nsys))
      enddo
    enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! output
    write(7,*) real(time),(evalfield(time,ipol),ipol=1,npol)   ! field
    write(8,*) real(time),norm,expectDIP,real(time_tau),IM_pol ! dipole components
    write(9,*) real(time),population(:)                        ! populations

  end subroutine writePOP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function evalfield(T,ipol)

    implicit none
    integer,intent(in) :: ipol
    real(kind=8) :: evalfield
    real(kind=8),intent(in) :: T
    real(kind=8),parameter :: pi = 3.1415926535897931d0
    real(kind=8),external :: seval

    if(readfield) then
      evalfield = seval(NfieldIN,T,Tin,field(1,ipol),coeff_lin(1,ipol),coeff_quad(1,ipol),coeff_quart(1,ipol))
    else
      evalfield = f0*sin(pi*T/width)**2*cos(carrier*(T-tpeak))
    endif
    return

  end function evalfield
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

  subroutine matvec(T,rhoIN,rhoOUT)

    implicit none
    real(kind=8),intent(in) :: T
    complex(kind=8),intent(in) :: rhoIN(0:Nbm1,0:Nsys)
    complex(kind=8),intent(out) :: rhoOUT(0:Nbm1,0:Nsys)

    integer :: i,j,m,mp,ipol
    real(kind=8) :: Et(Npol)
    complex(kind=8) :: dumc

! field intensity at time T
    !$OMP PARALLEL DO
    do ipol = 1,Npol
      Et(ipol) = evalfield(T,ipol)
    enddo
    !$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! diagonal elements
!!!!!!!!!!!!!!!!!!!
! bound states
!$OMP PARALLEL DO REDUCTION(+:dumc)
    do mp = 0,Nsm1
      dumc = zeroc
      do i = 0,mp-1
        do ipol = 1,Npol
          dumc = dumc + two*Et(ipol)*dipole(map_dip(ipol),i,mp)*dimag(rhoIN(mp,i))
        enddo
      enddo
      do i = mp,Nbm1
        do ipol = 1,Npol
          dumc = dumc + mtwo*Et(ipol)*dipole(map_dip(ipol),i,mp)*dimag(rhoIN(i,mp))
        enddo
      enddo
      do i = 0,Nsm1
        dumc = dumc + (rates(i,mp)*rhoIN(i,i) - rates(mp,i)*rhoIN(mp,mp))
      enddo
      do i = Nsys,Nbm1
        dumc = dumc + (rates(i,mp)*rhoIN(i,Nsys) - rates(mp,i)*rhoIN(mp,mp))
      enddo
      rhoOUT(mp,mp) = dumc + dimag(trans_energy(mp,mp))*rhoIN(mp,mp)
    enddo
    !$OMP END PARALLEL DO
!!!!!!!!!!!!!!!!!!
! ionizing states
    !$OMP PARALLEL DO REDUCTION(+:dumc)
    do mp = Nsys,Nbm1
      dumc = zeroc
      do i = 0,Nsm1
        do ipol = 1,Npol
          dumc = dumc + two*Et(ipol)*dipole(map_dip(ipol),i,mp)*dimag(rhoIN(mp,i))
        enddo
      enddo
      do ipol = 1,Npol
        dumc = dumc + mtwo*Et(ipol)*dipole(map_dip(ipol),mp,mp)*dimag(rhoIN(mp,Nsys))
      enddo
      do i = 0,Nsm1
        dumc = dumc + (rates(i,mp)*rhoIN(i,i) - rates(mp,i)*rhoIN(mp,Nsys))
      enddo
      do i = Nsys,Nbm1
        dumc = dumc + (rates(i,mp)*rhoIN(i,Nsys) - rates(mp,i)*rhoIN(mp,Nsys))
      enddo
      rhoOUT(mp,Nsys) = dumc + dimag(trans_energy(mp,Nsys))*rhoIN(mp,Nsys)
    enddo
    !$OMP END PARALLEL DO
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! off-diagonal elements
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! coupling between bound states
    !$OMP PARALLEL DO REDUCTION(+:dumc)
    do mp = 0,Nsm1
      do m = mp+1,Nsm1
        dumc = zeroc
        do i = 0,mp-1
          do ipol = 1,Npol
            dumc = dumc + Et(ipol)*dipole(map_dip(ipol),i,m)*conjg(rhoIN(mp,i))
          enddo
        enddo
        do i = mp,Nbm1
          do ipol = 1,Npol
            dumc = dumc + Et(ipol)*dipole(map_dip(ipol),i,m)*rhoIN(i,mp)
          enddo
        enddo
        do i = 0,m-1
          do ipol = 1,Npol
            dumc = dumc + mone*Et(ipol)*dipole(map_dip(ipol),i,mp)*rhoIN(m,i)
          enddo
        enddo
        do i = m,Nbm1
          do ipol = 1,Npol
            dumc = dumc + mone*Et(ipol)*dipole(map_dip(ipol),i,mp)*conjg(rhoIN(i,m))
          enddo
        enddo
        rhoOUT(m,mp) = icmplx*(dumc - trans_energy(m,mp)*rhoIN(m,mp)) - dephasing(m,mp)*rhoIN(m,mp)
      enddo
    enddo
    !$OMP END PARALLEL DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! coupling between bound and ionizing states
    !$OMP PARALLEL DO REDUCTION(+:dumc)
    do mp = 0,Nsm1
      do m = Nsys,Nbm1
        dumc = zeroc
        do i = 0,mp-1
          do ipol = 1,Npol
            dumc = dumc + Et(ipol)*dipole(map_dip(ipol),i,m)*conjg(rhoIN(mp,i))
          enddo
        enddo
        do i = mp,Nsm1
          do ipol = 1,Npol
            dumc = dumc + Et(ipol)*dipole(map_dip(ipol),i,m)*rhoIN(i,mp)
          enddo
        enddo
        do i = 0,Nsm1
          do ipol = 1,Npol
            dumc = dumc + mone*Et(ipol)*dipole(map_dip(ipol),i,mp)*rhoIN(m,i)
          enddo
        enddo
        do ipol = 1,Npol
          dumc = dumc + mone*Et(ipol)*dipole(map_dip(ipol),m,mp)*rhoIN(m,Nsys)
        enddo
        rhoOUT(m,mp) = icmplx*(dumc - trans_energy(m,mp)*rhoIN(m,mp)) - dephasing(m,mp)*rhoIN(m,mp)
      enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine matvec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

  subroutine matvecI(T,rhoIN,rhoOUT)

    use scalarOP
    implicit none
    real(kind=8),intent(in) :: T
    complex(kind=8),intent(in) :: rhoIN(0:Nbm1,0:Nsys)
    complex(kind=8),intent(out) :: rhoOUT(0:Nbm1,0:Nsys)

    integer :: i,j,m,mp,ipol
    real(kind=8) :: field,Et(Npol)
    complex(kind=8) :: dumc,rhoS(0:Nbm1,0:Nbm1)

! field intensity at time T
!$OMP PARALLEL DO
    do ipol = 1,Npol
      Et(ipol) = evalfield(T,ipol)
    enddo
    !$OMP END PARALLEL DO
! convert density to Schroedinger representation
!$OMP PARALLEL DO
    do mp = 0,Nsm1
    do m = mp,Nbm1
      rhoS(m,mp) = exp(-icmplx*dble(trans_energy(m,mp))*T)*rhoIN(m,mp)
    enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO
    do mp = Nsys,Nbm1
      rhoS(mp,Nsys) = rhoIN(mp,Nsys)
    enddo
    !$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! diagonal elements
!!!!!!!!!!!!!!!
! bound states
!$OMP PARALLEL DO REDUCTION(+:dumc)
    do mp = 0,Nsm1
      dumc = zeroc
      do i = 0,mp-1
        do ipol = 1,Npol
          dumc = dumc + two*Et(ipol)*dipole(map_dip(ipol),i,mp)*dimag(rhoS(mp,i))
        enddo
      enddo
      do i = mp,Nbm1
        do ipol = 1,Npol
          dumc = dumc + mtwo*Et(ipol)*dipole(map_dip(ipol),i,mp)*dimag(rhoS(i,mp))
        enddo
      enddo
      do i = 0,Nsm1
        dumc = dumc + (rates(i,mp)*rhoS(i,i) - rates(mp,i)*rhoS(mp,mp))
      enddo
      do i = Nsys,Nbm1
        dumc = dumc + (rates(i,mp)*rhoS(i,Nsys) - rates(mp,i)*rhoS(mp,mp))
      enddo
      rhoOUT(mp,mp) = dumc + dimag(trans_energy(mp,mp))*rhoIN(mp,mp)
    enddo
    !$OMP END PARALLEL DO
!!!!!!!!!!!!!!!!!!
! ionizing states
!$OMP PARALLEL DO REDUCTION (+:dumc)
    do mp = Nsys,Nbm1
      dumc = zeroc
      do i = 0,Nsm1
        do ipol = 1,Npol
          dumc = dumc + two*Et(ipol)*dipole(map_dip(ipol),i,mp)*dimag(rhoS(mp,i))
        enddo
      enddo
      do ipol = 1,Npol
        dumc = dumc + two*Et(ipol)*dipole(map_dip(ipol),mp,mp)*dimag(rhoS(mp,Nsys))
      enddo
      do i = 0,Nsm1
        dumc = dumc + (rates(i,mp)*rhoS(i,i) - rates(mp,i)*rhoS(mp,Nsys))
      enddo
      do i = Nsys,Nbm1
        dumc = dumc + (rates(i,mp)*rhoS(i,Nsys) - rates(mp,i)*rhoS(mp,Nsys))
      enddo
      rhoOUT(mp,Nsys) = dumc + dimag(trans_energy(mp,Nsys))*rhoIN(mp,Nsys)
    enddo
    !$OMP END PARALLEL DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! off-diagonal elements
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! coupling between bound states
!$OMP PARALLEL DO REDUCTION(+:dumc)
    do mp = 0,Nsm1
      do m = mp+1,Nsm1
        dumc = zeroc
        do i = 0,mp-1
          do ipol = 1,Npol
            dumc = dumc + Et(ipol)*dipole(map_dip(ipol),i,m)*conjg(rhoS(mp,i))
          enddo
        enddo
        do i = mp,Nbm1
          do ipol = 1,Npol
            dumc = dumc + Et(ipol)*dipole(map_dip(ipol),i,m)*rhoS(i,mp)
          enddo
        enddo
        do i = 0,m-1
          do ipol = 1,Npol
            dumc = dumc + mone* Et(ipol)*dipole(map_dip(ipol),i,mp)*rhoS(m,i)
          enddo
        enddo
        do i = m,Nbm1
          do ipol = 1,Npol
            dumc = dumc + mone*Et(ipol)*dipole(map_dip(ipol),i,mp)*conjg(rhoS(i,m))
          enddo
        enddo
        rhoOUT(m,mp) = exp(icmplx*dble(trans_energy(m,mp))*T)*(icmplx*dumc - dephasing(m,mp)*rhoS(m,mp)) &
                  + dimag(trans_energy(m,mp))*rhoIN(m,mp)
      enddo
    enddo
    !$OMP END PARALLEL DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! coupling between bound and ionizing states
!$OMP PARALLEL DO REDUCTION(+:dumc)
    do mp = 0,Nsm1
      do m = Nsys,Nbm1
        dumc = zeroc
        do i = 0,mp-1
          do ipol = 1,Npol
            dumc = dumc + Et(ipol)*dipole(map_dip(ipol),i,m)*conjg(rhoS(mp,i))
          enddo
        enddo
        do i = mp,Nsm1
          do ipol = 1,Npol
            dumc = dumc + Et(ipol)*dipole(map_dip(ipol),i,m)*rhoS(i,mp)
          enddo
        enddo
        do i = 0,Nsm1
          do ipol = 1,Npol
            dumc = dumc + mone*Et(ipol)*dipole(map_dip(ipol),i,mp)*rhoS(m,i)
          enddo
        enddo
        do ipol = 1,Npol
          dumc = dumc + mone*Et(ipol)*dipole(map_dip(ipol),m,mp)*rhoS(m,Nsys)
        enddo

        rhoOUT(m,mp) = exp(icmplx*dble(trans_energy(m,mp))*T)*(icmplx*dumc - dephasing(m,mp)*rhoS(m,mp)) &
                  + dimag(trans_energy(m,mp))*rhoIN(m,mp)
      enddo
    enddo
    !$OMP END PARALLEL DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine matvecI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!
!  subroutine matvecRW(T,rhoIN,rhoOUT)
!
!    implicit none
!    real(kind=8),intent(in) :: T
!    complex(kind=8),intent(in) :: rhoIN(0:Nbm1,0:Nbm1)
!    complex(kind=8),intent(out) :: rhoOUT(0:Nbm1,0:Nbm1)
!
!    integer :: i,j,m,mp,ipol
!    real(kind=8) :: Et(Npol)
!    complex(kind=8) :: dumc,RWfield(0:Nbm1,0:Nbm1)
!
!! field intensity at time T
!    do ipol = 1,Npol
!      Et(ipol) = evalfield(T,ipol)
!    enddo
!
!! evaluate the iquasiresonant field 
!    do mp = 0,Nbm1
!    do m = 0,Nbm1
!      RWfield(m,mp) = exp(-icmplx*fieldQUASI(m,mp)*T)
!    enddo
!    enddo
!      
!    do mp = 0,Nbm1
!! diagonal elements
!      dumc = zeroc
!      do i = 0,Nbm1
!        if(deltaK(i,mp)) then
!          do ipol = 1,Npol
!            dumc = dumc + Et(ipol)*dimag(RWfield(i,mp)*dipole(map_dip(ipol),i,mp)*rhoIN(i,mp))
!          enddo
!        endif
!      enddo
!      dumc = -dumc
!      do i = 0,Nbm1
!        dumc = dumc + (rates(i,mp)*rhoIN(i,i) - rates(mp,i)*rhoIN(mp,mp)) 
!      enddo
!      rhoOUT(mp,mp) = dumc + dimag(trans_energy(mp,mp))*rhoIN(mp,mp)
!! off-diagonal elements
!      do m = mp+1,Nbm1
!        dumc = zeroc
!        do i = 0,Nbm1
!          if(deltaK(m,i)) then
!            do ipol = 1,Npol
!              dumc = dumc + Et(ipol)*RWfield(m,i)*dipole(map_dip(ipol),i,m)*rhoIN(i,mp)
!            enddo
!          endif
!          if(deltaK(i,mp)) then
!            do ipol = 1,Npol
!              dumc = dumc - Et(ipol)*RWfield(i,mp)*dipole(map_dip(ipol),i,mp)*conjg(rhoIN(i,m))
!            enddo
!          endif
!        enddo
!        dumc = half*icmplx*dumc - dephasing(m,mp)*rhoIN(m,mp) + dimag(trans_energy(m,mp))*rhoIN(m,mp)
!        rhoOUT(m,mp) = dumc
!        rhoOUT(mp,m) = conjg(dumc)
!      enddo
!    enddo
!
!  end subroutine matvecRW

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

subroutine trpesrk(dT, rho, rhokin)
    integer :: m
    real(kind=8), intent(in) :: dT
    complex(kind=8), intent(in) :: rho(0:Nbm1,0:Nbm1)
    real(kind=8), intent(inout) :: rhokin(Nsys:Nbm1)
    real(kind=8) :: dum
    dum = half*mone*dT
    do m = Nsys, Nbm1
    rhokin(m) = rhokin(m) + dum*dimag(energies(m))*dble(rho(m,Nsys))
    enddo

end subroutine trpesrk

subroutine trpesa(dTdid, rhoold, rho, rhokin)
    integer :: m
    real(kind=8), intent(in) :: dTdid
    complex(kind=8), intent(in) :: rho(0:Nbm1,0:Nbm1)
    real(kind=8), intent(in) :: rhoold(Nsys:Nbm1)
    real(kind=8), intent(inout) :: rhokin(Nsys:Nbm1)
    real(kind=8) :: dum
    dum =  half*mone*dTdid
    do m = Nsys, Nbm1
    rhokin(m) = rhokin(m) + dum*dimag(energies(m))*(rhoold(m)+dble(rho(m,Nsys)))
    enddo
end subroutine trpesa

subroutine writekin(Fkin, rhokin)
    integer :: m
    character(len=400), intent(in) :: Fkin
    real(kind=8), intent(in) :: rhokin(Nsys:Nbm1)
    open(17,file=Fkin,status='replace')
    do m = Nsys, Nbm1
     write(17,*) dble(energies(m)), dble(rhokin(m))
     call flush(17)
    enddo
    close(17)
end subroutine writekin

end module matrix
