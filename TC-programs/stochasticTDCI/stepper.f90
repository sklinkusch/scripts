!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
subroutine stepper(psi,N_traj,method,Nsteps,Nskip,Nwrite,eps,OPv)

  use matrix
  use rungekutta
  implicit none
  character(len=30),intent(in) :: method
  integer,intent(in) :: N_traj,Nskip,Nwrite,Nsteps
  real(kind=8),intent(in) :: eps
  complex(kind=8),intent(inout) :: psi(0:Nbasis,N_traj)
  external OPv

  integer :: cnt,Nbp1,i,j
  real(kind=8) :: time,dT,dTdid,dTnext
  real(kind=8) :: err(N_traj)
  real(kind=8),parameter :: zero = 0d0
  complex(kind=8) :: dpsi(0:Nbasis,N_traj),psiNEW(0:Nbasis,N_traj)   ! state vectors
  Nbp1 = Nbasis+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! write initial population
  dT = Tfinal/(Nsteps+1)
  time = zero
  err(:) = zero
  call writePOP(time,psi,err,N_traj,Nwrite)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! adaptive step-size Runge-Kutta Cash-Karp propagation
  if(trim(method) == 'adaptive') then

    cnt = 0
    timeloop: do 
      if(time+dT > Tfinal) exit timeloop
! compute gradient
      call OPv(time,psi,dpsi,N_traj)
! non unitary evolution
      call adapt(Nbp1,psi,dpsi,psiNEW,N_traj,time,dT,eps,err,dTdid,dTnext,OPv)
! normalization
      call normalize(N_traj,psiNEW,psi)
! prepar for next step
      dT = dTnext
      if(modulo(cnt,Nskip)==0) call writePOP(time,psi,err,N_traj,Nwrite)
      cnt = cnt + 1
    enddo timeloop

! compute gradient and optimize field
    call OPv(time,psi,dpsi,N_traj)
! propagate last step
    dT = Tfinal - time
    call CashKarp_step(Nbp1,psi,dpsi,N_traj,time,dT,psiNEW,err,OPv)
! normalization
    call normalize(N_traj,psiNEW,psi)
    Time = Tfinal
    call writePOP(Tfinal,psi,err,N_traj,Nwrite)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! simple Runge-Kutta Cash-Karp propagation
  elseif(trim(method) == 'rungekutta') then

    cnt = 0
    ckloop : do
      if(time + dT > Tfinal+eps) exit ckloop
! compute gradient
      call OPv(time,psi,dpsi,N_traj)
! non unitary evolution
      call CashKarp_step(Nbp1,psi,dpsi,N_traj,time,dT,psiNEW,err,OPv)
! normalization
      call normalize(N_traj,psiNEW,psi)
! prepare for next iteration
      time = time + dT
      if(modulo(cnt,Nskip)==0) call writePOP(time,psi,err,N_traj,Nwrite)
      cnt = cnt + 1
    enddo ckloop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  else
    write(1,*) 'Oups! No such propagator available.'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  endif

end subroutine stepper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
