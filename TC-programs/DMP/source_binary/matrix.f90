!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
module matrix

  implicit none
  private
  integer :: NfieldIN,Npol
  integer,allocatable :: map_dip(:)
  logical,allocatable :: deltaK(:,:)
  logical,public :: readfield
  real(kind=8):: f0,carrier,width,tpeak,dumreal,esig,bosefkt
  real(kind=8),parameter :: zero = 0d0,half = 0.5d0,one = 1d0, two = 2d0
  real(kind=8),parameter :: mtwo = -2d0, mone = -1d0, highest = 0.001d0
  real(kind=8),parameter :: eps =  0.00000001d0
  real(kind=8),allocatable :: dipole(:,:,:),rates(:,:),dephasing(:,:)              ! the dipole and the transition rate matrices
  real(kind=8),allocatable :: Tin(:),field(:,:),coeff_lin(:,:),coeff_quad(:,:),coeff_quart(:,:),fieldQUASI(:,:)
  real(kind=8) :: permdip(3)
  real(kind=8), allocatable :: ens(:), dx(:,:), dy(:,:), dz(:,:), ion(:), vr(:,:)
  complex(kind=8),allocatable :: energies(:), trans_energy(:,:)     ! the transition energies
  complex(kind=8),public,parameter :: zeroc = dcmplx(0d0,0d0),icmplx=dcmplx(0d0,1d0), halfc=dcmplx(0.5d0,0d0)

  integer,public :: Nbasis,Nsys,Nbm1,Nsm1,Ntotal,Ntm1
  public matvec,matvec_d,matvecI,matvecI_d,writePOP,writePOP_d,initialize,trpesrk,trpesa,writekin,calc_cpl 
  public writeenpop

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  subroutine initialize(restart,Frst,Fecp,Fpop,polarization,rho,temperature,scaleION,&
          dep_pref,kopplung,koppl_pref,Ffield,Tinitial,Tfinal,freq,amplitude,kinens,readkin,&
          readpop,ionstates,matsubara,read_coupl,tempweight,zerogs,printrates)

    use scalarOP
    implicit none
    logical :: restart,kopplung,readkin,readpop,ionstates,read_coupl,zerogs
    logical :: printrates
    character(len=3),intent(in) :: polarization
    character(len=40),intent(in) :: Ffield, Frst, Fecp, Fpop
    character(len=400),intent(in) :: tempweight
    real(kind=8),intent(in) :: Temperature,Tinitial,Tfinal,amplitude,freq,scaleION
    real(kind=8),intent(in) :: dep_pref,koppl_pref,matsubara
    real(kind=8),intent(out) :: kinens(0:Nbm1)
    complex(kind=8),intent(out) :: rho(0:Nbm1,0:Nsys)
    integer :: i,j,idum,ipol,io,cpos
    real(kind=8),parameter :: Knst = 3.1577465d5
    real(kind=8),parameter :: minrate = 3.778036d-5
    real(kind=8) :: dum(3),dum1,deltaE
    real(kind=8) :: tmp(0:Nbm1),pop(0:Nbm1)
    complex(kind=8) :: dumc

! read the energies and evaluate the exponential transition energies
    allocate(energies(0:Nbm1))
    open(71,file=Fecp,access='stream',status='old')
    cpos = 1
! read number of states in the ecp file
    read(71,pos=cpos) Ntotal
    Ntm1 = Ntotal - 1
    if(Ntotal -  Nbasis < 0) stop 'ecp file has a wrong system size'
    allocate(ens(0:Nbm1), dx(0:Nbm1,0:Nbm1), dy(0:Nbm1,0:Nbm1), dz(0:Nbm1,0:Nbm1), ion(0:Nbm1), &
        vr(0:Nbm1,0:Nbm1))
    cpos = cpos + 4
! read energies
    do i=0, Nbm1
      read(71,pos=cpos) ens(i)
      cpos = cpos + 8
    enddo
    do i=Nbasis,Ntm1
      cpos = cpos + 8
    enddo
!    write(*,*) 'Energies:'
!    do i=0, Nbm1
!      write(*,*) i, ens(i)
!    enddo
! read dipole moment (x)
    do i=0, Nbm1
      do j=0, Nbm1
        read(71,pos=cpos) dx(i,j)
        cpos = cpos + 8
      enddo
      do j =Nbasis,Ntm1
        cpos = cpos + 8
      enddo
    enddo
    do i=Nbasis,Ntm1
      do j=0,Ntm1
        cpos = cpos + 8
      enddo
    enddo
! read dipole moment (y)
    do i=0, Nbm1
      do j=0, Nbm1
        read(71,pos=cpos) dy(i,j)
        cpos = cpos + 8
      enddo
      do j =Nbasis,Ntm1
        cpos = cpos + 8
      enddo
    enddo
    do i=Nbasis,Ntm1
      do j=0,Ntm1
        cpos = cpos + 8
      enddo
    enddo
! read dipole moment (z)
    do i=0, Nbm1
      do j=0, Nbm1
        read(71,pos=cpos) dz(i,j)
        cpos = cpos + 8
      enddo
      do j =Nbasis,Ntm1
        cpos = cpos + 8
      enddo
   enddo
    do i=Nbasis,Ntm1
      do j=0,Ntm1
        cpos = cpos + 8
      enddo
    enddo
!    write(*,*) 'Dipole moments (xyz): '
!    do i=0, Nbm1
!      do j=0, Nbm1
!        write(*,*) i, '->', j, dx(i,j), dy(i,j), dz(i,j)
!      enddo
!    enddo
if(read_coupl) then
    !read vibrational relaxation rates
    do i=0, Nbm1
      do j=0, Nbm1
        read(71,pos=cpos) vr(i,j)
        cpos = cpos + 8
      enddo
      do j=Nbasis,Ntm1
        cpos = cpos + 8
      enddo
    enddo
    do i=Nbasis, Ntm1
      do j=0, Ntm1
        cpos = cpos + 8
      enddo
    enddo
    !read internal conversion rates
!    do i=0, Nbm1
!      do j=0, Nbm1
!        read(71,pos=cpos) ic(i,j)
!        cpos = cpos + 8
!      enddo
!      do j=Nbasis,Ntm1
!        cpos = cpos + 8
!      enddo
!    enddo
!    do i=Nbasis, Ntm1
!      do j=0, Ntm1
!        cpos = cpos + 8
!      enddo
!    enddo
endif
! read relaxation rates
!    do i=0, Nbm1
!      do j=0, Nbm1
!        read(71,pos=cpos) cpl(i,j)
!        cpos = cpos + 8
!      enddo
!      do j=Nbasis,Ntm1
!        cpos = cpos + 8
!      enddo
!    enddo
!    do i=Nbasis, Ntm1
!      do j=0,Ntm1
!        cpos = cpos + 8
!      enddo
!    enddo
!endif
!    write(*,*) 'Relaxation rates: '
!    do i=0, Nbm1
!      do j=i, Nbm1
!        write(*,*) i, '<-', j, cpl(i,j)
!      enddo
!    enddo
! read kinetic energies
! read ionization rates
    do i=0,Nbm1
      read(71,pos=cpos) ion(i)
      cpos = cpos + 8
    enddo
    do i=Nbasis,Ntm1
      cpos = cpos + 8
    enddo
!    write(*,*) 'Ionization rates: '
!    do i=0, Nbm1
!      write(*,*) i, ion(i)
!    enddo
    if(readkin) then
      do i=0, Nbm1
        read(71,pos=cpos) kinens(i)
        cpos = cpos + 8
      enddo
      do i=Nbasis,Ntm1
        cpos = cpos + 8
      enddo
!      write(*,*) 'Kinetic energies: '
!      do i=0, Nbm1
!        write(*,*) i, kinens(i)
!      enddo
    endif
    close(71)
    !$OMP PARALLEL DO
    do i=0, Nbm1
     energies(i) = dcmplx(ens(i),mone*scaleION*ion(i))
    enddo
    !$OMP END PARALLEL DO
! read initial populations
    if(readpop) then
    open(73, file=Fpop, status='old', iostat=io)
    do i=0, Nbm1
    read(73,*) pop(i)
    enddo
    close(73)
    endif
!    write(*,*) 'Initial pops: '
!    do i=0, Nbm1
!      write(*,*) i, pop(i)
!    enddo
!    stop

    if(ionstates) then
        allocate(trans_energy(0:Nbm1,0:Nsys))
        do j = 0,Nsm1
          do i = 0,Nbm1
            trans_energy(i,j) = energies(i)-conjg(energies(j))
          enddo
        enddo
        do i = Nsys,Nbm1
          trans_energy(i,Nsys) = energies(i)-conjg(energies(i))
        enddo
    else
        allocate(trans_energy(0:Nbm1,0:Nbm1))
        do j=0,Nbm1
          do i=0,Nbm1
            trans_energy(i,j) = energies(i)-conjg(energies(j))
          enddo
        enddo
    endif

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
    write(1,*) "Number of laser polarizations: ", Npol
    call flush(1)

! initialize matrices
    allocate(map_dip(Npol),dipole(3,0:Nbm1,0:Nbm1),rates(0:Nbm1,0:Nbm1),dephasing(0:Nbm1,0:Nbm1))
    do j = 0,Nsm1
      do i = j,Nbm1
        rho(i,j) = zeroc
      enddo
    enddo
    if(ionstates) then
      do j = Nsys,Nbm1
        rho(j,Nsys) =  zeroc
      enddo
    endif
    do j = 0,Nbm1
      do i = 0,Nbm1
        do idum = 1,3
          dipole(idum,i,j) = zero
        enddo
         rates(i,j) = zero
      enddo
    enddo

! create mapping for dipole moment
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

    do i=0, Nbm1
      do j=0, Nbm1
        dipole(1,i,j) = dx(i,j)
        dipole(2,i,j) = dy(i,j)
        dipole(3,i,j) = dz(i,j)
      enddo
    enddo
!    do i=1,3
!      permdip(i) = dipole(i,0,0)
!      do j=0, Nbm1
!        dipole(i,j,j) = dipole(i,j,j) - permdip(i)
!      enddo
!    enddo

    if(ionstates) then
    do i=0, Nbm1
        do j=0, Nsys
           dephasing(i,j) = dep_pref * (dble(trans_energy(i,j)))**2
           dephasing(j,i) = dephasing(i,j)
        enddo
    enddo
    else
    do i=0, Nbm1
      do j=0, Nbm1
        dephasing(i,j) = dep_pref * (dble(trans_energy(i,j)))**2
        dephasing(j,i) = dephasing(i,j)
        enddo
      enddo
    endif

    if(kopplung) then
      if ( .NOT. read_coupl ) then 
          call calc_cpl(matsubara)
      endif
          do i=0, Nbm1
            do j=(i+1), Nbm1
               deltaE = abs(dble(energies(i)) - dble(energies(j)))
!               if(deltaE < eps) then
!                   rates(i,j) = zero
!                   rates(j,i) = zero
                   !rates(i,j) = koppl_pref*minrate
                   !rates(j,i) = koppl_pref*minrate
!               else
                  if(trim(tempweight) == 'bose') then
                    call bose(bosefkt, deltaE, Temperature)
                    rates(i,j) = koppl_pref*vr(i,j)*bosefkt
                    rates(j,i) = vr(j,i)*(one+bosefkt)
                  elseif(trim(tempweight) == 'boltz') then
                    call boltz(bosefkt, deltaE, Temperature)
                    rates(i,j) = koppl_pref*vr(i,j)*(bosefkt)
                    rates(j,i) = koppl_pref*vr(j,i)
                  else
                    stop 'Choose between bose and boltz'
                  endif
!                endif
        enddo
      enddo
      if(printrates) then
          do i=0,Nbm1
             do j=i,Nbm1
                 write(*,*) i, '->', j, ': ', abs(dble(energies(i) - energies(j))), rates(i,j), rates(j,i)
             enddo
          enddo
          stop
      endif
!      do i=0,Nbm1
!        do j=0, Nbm1
!          write(*,*) i, '->', j, ': ', rates(i,j), rates(j,i), cpl(i,j),&
!              cpl(j,i)
!        enddo
!      enddo
!write(*,'(6f20.10)') rates
!do i=0,Nbm1
!write(*,*) i, ens(i), kinens(i)
!enddo
!stop
    endif
! generate starting density matrix
    if(restart) then
! restart file
      open(11,file=Frst,status='old',form='unformatted')
      read(11,iostat=io) rho  !((rho(i,j), i=0,Nbm1), j=0,Nsys)  ! (0:Nbm1,0:Nsys)
      if(io /= 0) stop 'problem reading the rst-file '
      close(11)
      if(zerogs) then
          tmp(0) = sqrt(zero)
          do i = 0, Nsm1
             rho(0,i) = zero
             rho(i,0) = zero
          enddo
      endif
    else
! selected density matrix
      if(temperature < 1e-10 .AND. readpop) then
          ! generate density matrix from read populations at zero temperature
        do i = 0,Nbm1
          tmp(i) = sqrt(pop(i))
        enddo
        do i = 0,Nsm1
          do j = i,Nsm1
            rho(i,j) = tmp(i)*tmp(j)
            rho(j,i) = rho(i,j)
          enddo
        enddo
        if(ionstates) then
          do j = Nsys,Nbm1
            rho(j,Nsys) = tmp(j)**2
          enddo
        endif
      else if(.NOT. readpop .AND. temperature < 1e-10) then
          ! populate ground state at zero temperature (others are not populated)
        tmp(0) = one
        do i=1, Nbm1
          tmp(i) = zero
        enddo
        do i=0, Nsm1
          do j=i, Nsm1
            rho(i,j) = tmp(i)*tmp(j)
            rho(j,i) = rho(i,j)
          enddo
        enddo
        if(ionstates) then
            do j = Nsys, Nbm1
              rho(j, Nsys) = tmp(j)**2
            enddo
        endif
      else if(temperature .GE. 1e-10 .AND. readpop) then
          ! generate density matrix from read populations at arbitrary
          ! temperature
        do i=0,Nbm1
          tmp(i) = sqrt(pop(i))
        enddo
        do i=0,Nsm1
          do j=0,Nsm1
            rho(i,j) = tmp(i)*tmp(j)
            rho(j,i) = rho(i,j)
          enddo
        enddo
      else
! thermalize density
        do i = 0,Nbm1
          tmp(i) = exp(-Knst*dble(trans_energy(i,0))/temperature)
        enddo
        do i = 0,Nsm1
          rho(i,i) = tmp(i)
        enddo
        if(ionstates) then
          do j = Nsys,Nbm1
            rho(j,Nsys) = tmp(j)
          enddo
        endif
      endif
! normalize density
      if(ionstates) then
          dumc = trace(rho,Nbm1,Nsys)
          call normalize(dumc,rho,Nbasis,Nsys)
      else
          dumc = trace_d(rho,Nbm1)
          call normalize_d(dumc,rho,Nbasis)
      endif
    endif
!    do j = 0, Nsm1
!    write(*,*) j, real(rho(j,j))
!    enddo
!    do j = Nsys, Nbm1
!    write(*,*) j, real(rho(j,Nsys))
!    enddo
!    stop

    if(readfield) then
      open(4,file=Ffield,status='old',iostat=io)
! read the field
      write(*,*) Tfinal
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
      write(*,*) Tin(NfieldIN)
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
      width = (Tfinal-Tinitial)
      tpeak = 0.5d0*(Tfinal-Tinitial)
    endif
! evaluate memory requirement
!   write(1,*) 'Memory requirement:',(5*Nbasis**2+22*Nbasis*(Nsys+1)+(4*npol+1)*NfieldIN)*8./(1024.**3), 'Gb'
   call flush(1)

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  subroutine writePOP(interaction,time,rhoIN)

    implicit none
    logical,intent(in) :: interaction
    real(kind=8),intent(in) :: time
    complex(kind=8),intent(in) :: rhoIN(0:Nbm1,0:Nsys)
    integer :: m,mp,ipol
    real(kind=8) :: population(0:Nbasis)
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
      population(m) = dble(rhoS(m,m))
      norm = norm + population(m)
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO REDUCTION(+:norm)
    do m = Nsys,Nbm1
      population(m) = dble(rhoS(m,Nsys))
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
    call flush(7)
    call flush(8)
    call flush(9)

  end subroutine writePOP

  subroutine writePOP_d(interaction, time, rhoIN)

    implicit none
    logical,intent(in) :: interaction
    real(kind=8),intent(in) :: time
    complex(kind=8),intent(in) :: rhoIN(0:Nbm1,0:Nbm1)
    integer :: m,mp,ipol
    real(kind=8) :: population(0:Nbasis)
    real(kind=8) :: time_tau,norm
    real(kind=8) :: expectDIP(3),IM_pol(3)
    complex(kind=8) :: rhoS(0:Nbm1,0:Nbm1),tmp(0:Nbm1,0:Nbm1)

    if(interaction) then
! convert density to Schroedinger representation
!$OMP PARALLEL DO
      do mp = 0,Nbm1
      do m = mp,Nbm1
        rhoS(m,mp) = exp(-icmplx*dble(trans_energy(m,mp))*time)*rhoIN(m,mp)
      enddo
      enddo
      !$OMP END PARALLEL DO
    else
        !$OMP PARALLEL DO
      do mp = 0,Nbm1
      do m = mp,Nbm1
        rhoS(m,mp) = rhoIN(m,mp)
      enddo
      enddo
      !$OMP END PARALLEL DO
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute population in states below ionization threshold
    norm = zero
    !$OMP PARALLEL DO REDUCTION(+:norm)
    do m = 0,Nbm1
      population(m) = dble(rhoS(m,m))
      norm = norm + population(m)
    enddo
    !$OMP END PARALLEL DO
    population(Nbasis) = 1e0 - norm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! COMPUTE dipole expectation value

    expectDIP(:) = zero
    do mp = 0,Nbm1
      do ipol = 1,3
        expectDIP(ipol) = expectDIP(ipol) + dipole(ipol,mp,mp)*dble(rhoS(mp,mp))
      enddo
      do m = mp+1,Nbm1
        do ipol = 1,3
          expectDIP(ipol) = expectDIP(ipol) + two*dipole(ipol,m,mp)*dble(rhoS(m,mp))
        enddo
      enddo
    enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute the imaginary part of the polarizability
    time_tau = time + one/carrier
    call matvec_d(time_tau,rhoS,tmp)
    IM_pol(:) = zero
    do mp = 0,Nbm1
      do ipol = 1,3
        IM_pol(ipol) = IM_pol(ipol) + dble(tmp(mp,mp))
      enddo
      do m = mp+1,Nbm1
        do ipol = 1,3
          IM_pol(ipol) = IM_pol(ipol) + two*dble(tmp(m,mp))
        enddo
      enddo
    enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! output
    write(7,*) real(time),(evalfield(time,ipol),ipol=1,npol)   ! field
    write(8,*) real(time),norm,expectDIP,real(time_tau),IM_pol ! dipole components
    write(9,*) real(time),population(:)                        ! populations
    call flush(7)
    call flush(8)
    call flush(9)

  end subroutine writePOP_d
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

    integer :: i,m,mp,ipol
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
    dumc = zeroc
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

  subroutine matvec_d(T,rhoIN,rhoOUT)

    implicit none
    real(kind=8),intent(in) :: T
    complex(kind=8),intent(in) :: rhoIN(0:Nbm1,0:Nbm1)
    complex(kind=8),intent(out) :: rhoOUT(0:Nbm1,0:Nbm1)

    integer :: i,m,mp,ipol
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
    dumc = zeroc
!$OMP PARALLEL DO REDUCTION(+:dumc)
    do mp = 0,Nbm1
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
      do i = 0,Nbm1
        dumc = dumc + (rates(i,mp)*rhoIN(i,i) - rates(mp,i)*rhoIN(mp,mp))
      enddo
      rhoOUT(mp,mp) = dumc + dimag(trans_energy(mp,mp))*rhoIN(mp,mp)
    enddo
    !$OMP END PARALLEL DO
!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! off-diagonal elements
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! coupling between bound states
    !$OMP PARALLEL DO REDUCTION(+:dumc)
    do mp = 0,Nbm1
      do m = mp+1,Nbm1
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
  end subroutine matvec_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

  subroutine matvecI(T,rhoIN,rhoOUT)

    use scalarOP
    implicit none
    real(kind=8),intent(in) :: T
    complex(kind=8),intent(in) :: rhoIN(0:Nbm1,0:Nsys)
    complex(kind=8),intent(out) :: rhoOUT(0:Nbm1,0:Nsys)

    integer :: i,m,mp,ipol
    real(kind=8) :: Et(Npol)
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
    dumc = zeroc
    !$OMP PARALLEL DO REDUCTION (+:dumc)
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

  subroutine matvecI_d(T,rhoIN,rhoOUT)

    use scalarOP
    implicit none
    real(kind=8),intent(in) :: T
    complex(kind=8),intent(in) :: rhoIN(0:Nbm1,0:Nbm1)
    complex(kind=8),intent(out) :: rhoOUT(0:Nbm1,0:Nbm1)

    integer :: i,m,mp,ipol
    real(kind=8) :: Et(Npol)
    complex(kind=8) :: dumc,rhoS(0:Nbm1,0:Nbm1)

! field intensity at time T
!$OMP PARALLEL DO
    do ipol = 1,Npol
      Et(ipol) = evalfield(T,ipol)
    enddo
    !$OMP END PARALLEL DO
! convert density to Schroedinger representation

!$OMP PARALLEL DO
    do mp = 0,Nbm1
    do m = mp,Nbm1
      rhoS(m,mp) = exp(-icmplx*dble(trans_energy(m,mp))*T)*rhoIN(m,mp)
    enddo
    enddo
    !$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! diagonal elements
!!!!!!!!!!!!!!!
! bound states
    dumc = zeroc
    !$OMP PARALLEL DO REDUCTION (+:dumc)
    do mp = 0,Nbm1
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
      rhoOUT(mp,mp) = dumc + dimag(trans_energy(mp,mp))*rhoIN(mp,mp)
    enddo
    !$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!
! off-diagonal elements
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! coupling between bound states
!$OMP PARALLEL DO REDUCTION(+:dumc)
    do mp = 0,Nbm1
      do m = mp+1,Nbm1
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

    end subroutine matvecI_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

subroutine writekin(Fkin, rhokin,kinens,readkin)
    logical :: readkin
    integer :: m
    character(len=400), intent(in) :: Fkin
    real(kind=8), intent(in) :: rhokin(Nsys:Nbm1),kinens(0:Nbm1)

    open(17,file=Fkin,status='replace')
    do m = Nsys, Nbm1
     if(readkin) then
         write(17,*) dble(kinens(m)), dble(rhokin(m))
     else
         write(17,*) dble(energies(m)),dble(rhokin(m))
     endif
     call flush(17)
    enddo
    close(17)
end subroutine writekin

subroutine calc_cpl(matsubara)
    integer :: i,j
    real(kind=8), intent(in) :: matsubara

    do i=0, Nbm1
      do j=0, Nbm1
        vr(i,j) = zero
      enddo
    enddo

    do i=0,Nsm1
      do j=i+1, Nsm1
        vr(i,j) = (matsubara / (matsubara**2 + (dble(energies(i)-conjg(energies(j))))**2))
!        write(*,*) i,j,dble(energies(i)),dble(energies(j)),cpl(i,j)
        vr(j,i) = vr(i,j)
      enddo
    enddo
!    stop

    end subroutine calc_cpl

    subroutine writeenpop(interaction, time, rhoIN)
        logical :: interaction
        integer :: m,mp
        complex(kind=8) :: rhoS(0:Nbm1,0:Nsys) 
        real(kind=8) :: pop(0:Nbm1)
        real(kind=8), intent(in) :: time
        complex(kind=8), intent(in) :: rhoIN(0:Nbm1,0:Nsys)

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
    !$OMP PARALLEL DO
    do m = 0,Nsm1
      pop(m) = dble(rhoS(m,m))
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO
    do m = Nsys,Nbm1
      pop(m) = dble(rhoS(m,Nsys))
    enddo
    !$OMP END PARALLEL DO
    do m=0,Nbm1
    write(415,*) m, dble(energies(m)), pop(m)
    enddo



    end subroutine writeenpop

end module matrix
