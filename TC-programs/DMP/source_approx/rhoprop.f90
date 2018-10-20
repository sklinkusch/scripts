!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
program RKprop

  use scalarOP
  use matrix
  use rungekutta
  implicit none
  logical :: interaction,restart
  integer :: i,j,cnt,io,Nsteps,Nskip,m,mdl
  character(len=3) :: polarization
  character(len=400) :: Fin,Fmu,Fenergies,Ffield,fout,foutpop,foutlog,foutdip,foutfld,method
  character(len=400) :: Fkin, Frst
  real(kind=8) :: temperature,freq,amplitude,scaleION
  real(kind=8) :: Tinitial,Tfinal,time,dT,dTdid,dTnext,eps,err
  complex(kind=8),allocatable :: rho(:,:),drho(:,:),tmp(:,:)            ! the density matrix and its first derivative
  real(kind=8),allocatable :: rhokin(:),rhoold(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! get input file name

  call getarg(1,Fin)

! read input parameters
  open(11,file=Fin,status='old',iostat=io)
  if(io /= 0) stop 'Problem opening input file'
  read(11,*)
  read(11,*) restart
  read(11,*) Frst
  read(11,*) Fenergies
  read(11,*) Fmu
  read(11,*) Nbasis, Nsys
  read(11,*) scaleION
  read(11,*) Temperature
  read(11,*)
  read(11,*)
  read(11,*) Tinitial,Tfinal; Tinitial = Tinitial*41.341373337d0; Tfinal = Tfinal*41.341373337d0
  read(11,*) Nsteps
  read(11,*) method
  read(11,*) interaction
  read(11,*) eps
  read(11,*)
  read(11,*)
  read(11,*) readfield
  read(11,*) Ffield
  read(11,*) polarization
  read(11,*) freq
  read(11,*) amplitude
  read(11,*)
  read(11,*)
  read(11,*) fout
  read(11,*) Nskip
  close(11)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! define problem dimensions
  Nbm1 = Nbasis - 1
  Nsm1 = Nsys - 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! I/O
  foutfld = trim(fout)//'.fld'
  foutlog = trim(fout)//'.log'
  foutpop = trim(fout)//'.pop'
  foutdip = trim(fout)//'.dip'
  Fkin    = trim(fout)//'.kin'
  open(1,file=foutlog,status='replace',iostat=io)
  if(io /= 0) stop 'problem opening the log file'
  write(1,*) '----------------------------------------------------------------------'
  write(1,*) 'Number of eigenstates included in the basis: ',Nbasis
  write(1,*) 'Number of eigenstates below the ionization barrier: ',Nsys
  write(1,*) 'Temperature of the system: ',real(Temperature)
  write(1,*) 'Energies of the isolated subsystem read from file: ',Fenergies
  write(1,*) 'Dipole matrix element and initial density read from file: ',Fmu
  write(1,*) 'Field read from file: ',Ffield
  write(1,*) 'Polarization of the electric field :: ',polarization
  call flush(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! initializing the matrices for the propagation
  open(2,file=Fenergies,status='old',iostat=io)
  if(io /= 0) stop 'problem opening the energy file'
  open(3,file=Fmu,status='old',iostat=io)
  if(io /= 0) stop 'problem opening the dipole file'
  allocate(rho(0:Nbm1,0:Nsys),tmp(0:Nbm1,0:Nsys),drho(0:Nbm1,0:Nsys))
  allocate(rhokin(Nsys:Nbm1),rhoold(Nsys:Nbm1))

  call initialize(restart,Frst,polarization,rho,temperature,scaleION,Ffield,Tinitial,Tfinal,freq,amplitude)
  dt = (Tfinal-Tinitial)/(Nsteps+1)
  close(2)
  close(3)

! summerizing propagation parameters
  write(1,*) 'Propagation time: ',real((Tfinal-Tinitial)/41.341375583725707),' fs'
  write(1,*) 'using ', Nsteps, ' timesteps of length ', dt, ' fs'
  write(1,*) 'writing out every ', Nskip, ' timesteps'
  if(trim(method) == 'adaptive') then
    write(1,*) 'Adaptive time-step Runge-Kutta Cash-Karp scheme'
  elseif(trim(method) == 'rungekutta') then
    write(1,*) 'Simple Runge-Kutta Cash-Karp scheme'
  elseif(trim(method) == 'quasires') then
    write(1,*) 'Rotating wave approximation + Adaptive time-step Runge-Kutta Cash-Karp scheme'
  else
    stop 'No such method. Choose between rungekutta, adaptive and quasires.'
  endif
  if(interaction) then
    write(1,*) 'Propagating in the interaction representation'
  else
    write(1,*) 'Propagating in the Schroedinger representation'
  endif
  write(1,*) 'Population output in file: ',foutpop
  write(1,*) 'Dipole and field output in file: ',foutdip
  write(1,*) '----------------------------------------------------------------------'
  write(1,*) ''
  call flush(1)

  close(2)
  close(3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! main propagation loop
  open(7,file=foutfld,status='replace')
  open(8,file=foutdip,status='replace')
  open(9,file=foutpop,status='replace')
! write initial population
  time = Tinitial
  do m = Nsys, Nbm1
   rhokin(m) = 0d0
  enddo
  call writePOP(interaction,time,rho)
  write(1,*) real(time), real(trace(rho,Nbasis,Nsys))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! adaptive step-size Runge-Kutta Cash-Karp propagation
  if(trim(method) == 'adaptive') then
    cnt = 1
    timeloop: do 
! propagate for one step
      if(time+dT > Tfinal) exit timeloop
      do m = Nsys, Nbm1
        rhoold(m) = dble(rho(m,Nsys))
      enddo
      if(interaction) then
        call matvecI(time,rho,drho)
        call adapt(Nbasis,Nsys,rho,drho,time,dT,eps,err,dTdid,dTnext,matvecI)
      else
        call matvec(time,rho,drho)
        call adapt(Nbasis,Nsys,rho,drho,time,dT,eps,err,dTdid,dTnext,matvec)
      endif
      dT = dTnext
      write(1,*) real(time), dT, real(trace(rho,Nbasis, Nsys)), err
      call trpesa(dTdid, rhoold, rho, rhokin)
      mdl = modulo(cnt,Nskip)
      if(mdl == 0) call writePOP(interaction,time,rho)
      cnt = cnt + 1
    enddo timeloop
! propagate last step
    dT = Tfinal - time
    do m = Nsys, Nbm1
     rhoold(m) = dble(rho(m,m))
    enddo
    if(interaction) then
      call matvecI(time,rho,drho)
      call CashKarp_step(Nbasis,Nsys,rho,drho,time,dT,tmp,err,matvecI)
    else
      call matvec(time,rho,drho)
      call CashKarp_step(Nbasis,Nsys,rho,drho,time,dT,tmp,err,matvec)
    endif
    call trpesa(dT, rhoold, rho, rhokin)
    call writePOP(interaction,Tfinal,rho)
    call writekin(Fkin,rhokin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! simple Runge-Kutta Cash-Karp propagation
  elseif(trim(method) == 'rungekutta') then
    cnt = 1
    ckloop : do
! propagate for one step
      time = time + dT
      if(time > Tfinal+eps) then
          call writekin(Fkin, rhokin)
          exit ckloop
      endif
      call trpesrk(dT,rho,rhokin)
      if(interaction) then
        call matvecI(time,rho,drho)
        call CashKarp_step(Nbasis,Nsys,rho,drho,time,dT,tmp,err,matvecI)
      else
        call matvec(time,rho,drho)
        call CashKarp_step(Nbasis,Nsys,rho,drho,time,dT,tmp,err,matvec)
      endif
      rho(0:Nbm1,0:Nsys) = tmp(0:Nbm1,0:Nsys)
      call trpesrk(dT,rho,rhokin)
      write(1,*) real(time),real(trace(tmp,Nbasis,Nsys)),err ! tmp = rho
      if(modulo(cnt,Nskip)==0) call writePOP(interaction,time,rho)
      cnt = cnt + 1
    enddo ckloop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! propagation within the rotating wave approximation
  else
    if(Nsys /= Nbasis) stop 'Quasiresonant approximation not implemented'
!    cnt = 1
!    rwloop: do 
!! propagate for one step
!      if(time+dT > Tfinal) exit rwloop
!      call matvecRW(time,rho,drho)
!      call adapt(Nbasis,Nsys,rho,drho,time,dT,eps,err,dTdid,dTnext,matvecRW)
!      dT = dTnext
!      if(modulo(cnt,Nskip)==0) call writePOP(interaction,time,rho)
!      cnt = cnt + 1
!    enddo rwloop
!! propagate last step
!    dT = Tfinal - time
!    call matvecRW(time,rho,drho)
!    call adapt(Nbasis,Nsys,rho,drho,time,dT,eps,err,dTdid,dTnext,matvecRW)
!    call writePOP(interaction,Tfinal,rho)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  endif
  close(7)
  close(8)
  close(9)
  close(1)

! save density matrix for restarting
  foutlog = trim(fout)//'.rst'
  open(11,file=foutlog,status='replace',form='unformatted')
  write(11) rho(0:Nbm1,0:Nsys)
  close(11)

end program RKprop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
