!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
program RKprop

  use scalarOP
  use matrix
  use rungekutta
  implicit none
  logical :: interaction,restart,kopplung,readkin,outkin,readpop,ionstates
  logical :: zerogs
  logical :: read_coupl,printrates
  integer :: cnt,io,Nsteps,Nskip,m,mdl
  integer :: today(3), now(3)
  character(len=3) :: polarization
  character(len=400) :: Fin,Ffield,fout,foutpop,foutlog,foutdip,foutfld,method
  character(len=400) :: Fkin, Frst, Fecp, Fpop,tempweight,Fenpop
  real(kind=8) :: temperature,freq,amplitude,scaleION, dep_pref
  real(kind=8) :: Tinitial,Tfinal,time,dT,dTdid,dTnext,eps,err
  real(kind=8) :: koppl_pref,matsubara
  real(kind=8),parameter :: zero=0.0d0
  complex(kind=8),allocatable :: rho(:,:),drho(:,:),tmp(:,:)            ! the density matrix and its first derivative
  real(kind=8),allocatable :: kinens(:),rhokin(:),rhoold(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! get input file name

  call getarg(1,Fin)
  kopplung = .true.

! read input parameters
  open(11,file=Fin,status='old',iostat=io)
  if(io /= 0) stop 'Problem opening input file'
  read(11,*)
  read(11,*) restart
  read(11,*) Frst
  read(11,*) Fecp
  read(11,*) readpop
  read(11,*) Fpop
  read(11,*) Nbasis, Nsys
  read(11,*) scaleION
  read(11,*) Temperature
  if(Temperature < zero) Temperature = zero
  read(11,*) tempweight
  read(11,*) readkin
  read(11,*) outkin
  read(11,*) zerogs
  read(11,*)
  read(11,*)
  read(11,*) Tinitial,Tfinal; Tinitial = Tinitial*41.341373337d0; Tfinal = Tfinal*41.341373337d0
  read(11,*) Nsteps
  read(11,*) method
  read(11,*) interaction
  read(11,*) eps
  read(11,*) dep_pref
  read(11,*) koppl_pref
  if(koppl_pref == zero) kopplung = .false.
  read(11,*) read_coupl
  read(11,*) matsubara
  read(11,*) printrates
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

  if(Nbasis > Nsys) then
      ionstates = .true.
  else if(Nbasis .EQ. Nsys) then
      ionstates = .false.
      scaleION = zero
  else
      stop "Nr of nonionizing states cannot be larger than nr of total states"
  endif
  
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
  Fenpop  = trim(fout)//'.enpop'
  open(1,file=foutlog,status='replace',iostat=io)
  if(io /= 0) stop 'problem opening the log file'
  call idate(today) ! today(1) = day, today(2) = month, today(3) = year
  call itime(now) ! now(1) = hour, now(2) = minute, now(3) = second
1055 format ('Current date and time: ',i2.2,'.',i2.2,'.',i4.4,' ', i2.2,':',&
    i2.2,':',i2.2)
  write(1,1055) today, now
  write(1,*) '----------------------------------------------------------------------'
  if(restart) then
      write(1,*) 'Reading initial density matrix from ', Frst
  else
      write(1,*) 'Constructing initial density matrix from ', Fpop
  endif
  write(1,*) 'Reading electronic states and dipoles from ', Fecp
  write(1,*) 'Reading initial populations from ', Fpop
  write(1,*) 'Number of eigenstates included in the basis: ',Nbasis
  write(1,*) 'Number of eigenstates below the ionization barrier: ',Nsys
  write(1,*) 'damping factor for ionization: ', scaleION
  write(1,*) 'Temperature of the system: ',real(Temperature)
  write(1,*) 'Initial and final propagation time (a.u.): ', Tinitial, Tfinal
  write(1,*) 'Number of timesteps for propagation: ', Nsteps
  write(1,*) 'Method used for propagation: ', method
  if(interaction) then
      write(1,*) 'Propagating in the interaction picture'
  else
      write(1,*) 'Propagating in the Schroedinger picture'
  endif
  write(1,*) 'Relative tolerance for adaptive step size: ', eps
  if(readfield) then
      write(1,*) 'Field read from file: ',Ffield
  else
      write(1,*) 'Polarization of the electric field: ',polarization
      write(1,*) 'Frequency of the electric field: ', freq
      write(1,*) 'Amplitude of the electric field: ', amplitude
  endif
  write(1,*) 'Writing populations to ', foutpop
  write(1,*) 'Writing dipoles to ', foutdip
  write(1,*) 'Writing electric field to ', foutfld
  write(1,*) 'Writing TRPES information to ', Fkin
  call flush(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! initializing the matrices for the propagation
  if(ionstates) then
      allocate(rho(0:Nbm1,0:Nsys),tmp(0:Nbm1,0:Nsys),drho(0:Nbm1,0:Nsys))
      allocate(rhoold(Nsys:Nbm1))
      if(outkin) then
          allocate(rhokin(Nsys:Nbm1))
      endif
  else
      allocate(rho(0:Nbm1,0:Nbm1),tmp(0:Nbm1,0:Nbm1),drho(0:Nbm1,0:Nbm1))
  endif
  allocate(kinens(0:Nbm1))
!  if(outkin) then
!      allocate(rhokin(Nsys:Nbm1))
!  endif
!  allocate(kinens(0:Nbm1))


  call initialize(restart,Frst,Fecp,Fpop,polarization,rho,temperature,&
      scaleION,dep_pref,kopplung,koppl_pref,Ffield,Tinitial,Tfinal,&
      freq,amplitude, kinens, readkin,readpop,ionstates,matsubara,&
      read_coupl,tempweight,zerogs,printrates)
  dt = (Tfinal-Tinitial)/(Nsteps+1)
  write(1,*) '----------------------------------------------------------------------'
  write(1,*) ''
  call flush(1)
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
  else
    stop 'No such method. Choose between rungekutta and  adaptive.'
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! main propagation loop
  open(7,file=foutfld,status='replace')
  open(8,file=foutdip,status='replace')
  open(9,file=foutpop,status='replace')
  open(415,file=Fenpop,status='replace')
! write initial population
  time = Tinitial
  if(ionstates .AND. outkin) then
    do m = Nsys, Nbm1
      rhokin(m) = zero
    enddo
  endif
  if(ionstates) then
      call writePOP(interaction,time,rho)
      write(1,*) real(time), real(trace(rho,Nbm1,Nsys))
  else
      call writePOP_d(interaction,time,rho)
      write(1,*) real(time), real(trace_d(rho,Nbm1))
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! adaptive step-size Runge-Kutta Cash-Karp propagation
  if(trim(method) == 'adaptive') then
    cnt = 1
    timeloop: do 
! propagate for one step
      if(time+dT > Tfinal) exit timeloop
      if(ionstates) then
          do m = Nsys, Nbm1
             rhoold(m) = dble(rho(m,Nsys))
          enddo
      endif
      if(interaction) then
          if(ionstates) then
              call matvecI(time,rho,drho)
              call adapt(Nbasis,Nsys,rho,drho,time,dT,eps,err,dTdid,dTnext,matvecI)
          else
              call matvecI_d(time,rho,drho)
              call adapt_d(Nbasis,rho,drho,time,dT,eps,err,dTdid,dTnext,&
                  matvecI_d)
          endif
      else
          if(ionstates) then
              call matvec(time,rho,drho)
              call adapt(Nbasis,Nsys,rho,drho,time,dT,eps,err,dTdid,dTnext,matvec)
          else
              call matvec_d(time,rho,drho)
              call adapt_d(Nbasis,rho,drho,time,dT,eps,err,dTdid,dTnext,&
                  matvec_d)
          endif
      endif
      dT = dTnext
      if(ionstates) then
          write(1,*) real(time), dT, real(trace(rho,Nbm1, Nsys)), err
      else
          write(1,*) real(time), dT, real(trace_d(rho,Nbm1)), err
      endif
      if(ionstates .AND. outkin) then
      call trpesa(dTdid, rhoold, rho, rhokin)
      endif
      mdl = modulo(cnt,Nskip)
      if(mdl == 0) then 
          if(ionstates) then
              call writePOP(interaction,time,rho)
          else
              call writePOP_d(interaction,time,rho)
          endif
      endif
      cnt = cnt + 1
    enddo timeloop
! propagate last step
    dT = Tfinal - time
    if(ionstates) then
        do m = Nsys, Nbm1
           rhoold(m) = dble(rho(m,Nsys))
        enddo
    endif
    if(interaction) then
        if(ionstates) then
            call matvecI(time,rho,drho)
            call CashKarp_step(Nbasis,Nsys,rho,drho,time,dT,tmp,err,matvecI)
        else
            call matvecI_d(time,rho,drho)
            call CashKarp_step_d(Nbasis,rho,drho,time,dT,tmp,err,matvecI_d)
        endif
    else
        if(ionstates) then
            call matvec(time,rho,drho)
            call CashKarp_step(Nbasis,Nsys,rho,drho,time,dT,tmp,err,matvec)
        else
            call matvec_d(time,rho,drho)
            call CashKarp_step_d(Nbasis,rho,drho,time,dT,tmp,err,matvec_d)
        endif
    endif
    if(ionstates .AND. outkin) then
       call trpesa(dT, rhoold, rho, rhokin)
    endif
    if(ionstates) then
        call writePOP(interaction,Tfinal,rho)
    else
        call writePOP_d(interaction,Tfinal,rho)
    endif
    call writeenpop(interaction, time, rho)
    if(ionstates .AND. outkin) then
        call writekin(Fkin,rhokin,kinens,readkin)
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! simple Runge-Kutta Cash-Karp propagation
  elseif(trim(method) == 'rungekutta') then
    cnt = 1
    ckloop : do
! propagate for one step
      time = time + dT
      if(time > Tfinal+eps) then
          if(ionstates .AND. outkin) then
              call writekin(Fkin, rhokin, kinens,readkin)
          endif
          call writeenpop(interaction, time, rho)
          exit ckloop
      endif
      if(ionstates .AND. outkin) then
         call trpesrk(dT,rho,rhokin)
      endif
      if(interaction) then
          if(ionstates) then
              call matvecI(time,rho,drho)
              call CashKarp_step(Nbasis,Nsys,rho,drho,time,dT,tmp,err,matvecI)
          else
              call matvecI_d(time,rho,drho)
              call CashKarp_step_d(Nbasis,rho,drho,time,dT,tmp,err,&
                  matvecI_d)
          endif
      else
          if(ionstates) then
              call matvec(time,rho,drho)
              call CashKarp_step(Nbasis,Nsys,rho,drho,time,dT,tmp,err,matvec)
          else
              call matvec_d(time,rho,drho)
              call CashKarp_step_d(Nbasis,rho,drho,time,dT,tmp,err,&
                  matvec_d)
          endif
      endif
      if(ionstates) then
          rho(0:Nbm1,0:Nsys) = tmp(0:Nbm1,0:Nsys)
      else
          rho(0:Nbm1,0:Nbm1) = tmp(0:Nbm1,0:Nbm1)
      endif
      if(ionstates .AND. outkin) then
          call trpesrk(dT,rho,rhokin)
      endif
      if(ionstates) then
          write(1,*) real(time),real(trace(tmp,Nbm1,Nsys)),err ! tmp = rho
      else
          write(1,*) real(time),real(trace_d(tmp,Nbm1)), err
      endif
      if(modulo(cnt,Nskip)==0) then
          if(ionstates) then
              call writePOP(interaction,time,rho)
          else
              call writePOP_d(interaction,time,rho)
          endif
      endif
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
  close(415)
  call idate(today) ! today(1) = day, today(2) = month, today(3) = year
  call itime(now) ! now(1) = hour, now(2) = minute, now(3) = second
  write(1, 1055) today, now
  close(1)

! save density matrix for restarting
  foutlog = trim(fout)//'.rst'
  open(11,file=foutlog,status='replace',form='unformatted')
  if(ionstates) then
      write(11) rho(0:Nbm1,0:Nsys)
  else
      write(11) rho(0:Nbm1,0:Nbm1)
  endif
  close(11)

end program RKprop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
