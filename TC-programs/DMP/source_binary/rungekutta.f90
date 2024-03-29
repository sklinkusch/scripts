!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! contains the subroutines for propagating a set of N density matrices
! over a time interval dT using a Runge-Kutta Cash-Karp integrator 
! with adaptive time stpe
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
MODULE RungeKutta

  implicit none
  private
  public adapt,adapt_d,cashKarp_step,cashKarp_step_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
CONTAINS
  SUBROUTINE adapt(Nbasis,Nsys,rho,drho,time,dTtry,eps,err,dTdid,dTnext,derivatives)

    implicit none
    integer,intent(in) :: Nbasis,Nsys
    real(kind=8),intent(in) :: dTtry,eps
    real(kind=8),intent(out) :: dTdid,dTnext
    real(kind=8),intent(inout) :: time,err
    complex(kind=8),intent(inout) :: rho(Nbasis,0:Nsys),drho(Nbasis,0:Nsys)
    external derivatives

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! local variables

    real(kind=8) :: dT,Tnew,errmax
    complex(kind=8) :: tmp(Nbasis,0:Nsys)
    real(kind=8),parameter :: SAFETY = 0.9d0
    real(kind=8),parameter :: PGROW = 0.2d0
    real(kind=8),parameter :: PSHRNK = 0.25d0
    real(kind=8),parameter :: ERRCON = 1.89d-4
    real(kind=8),parameter :: minDT = 1d-5        ! minimal time step in a.u.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! adaptive step loop
    dT = dTtry
    step: do 
! try a first step
      call CashKarp_step(Nbasis,Nsys,rho,drho,time,dT,tmp,err,derivatives)
      errmax = err/eps
      dTdid = dT
      if(errmax < 1.) exit step
! if the error is too big, retake the step with a smaller time interval
      dT = min(SAFETY*dTdid*(errmax**PSHRNK),0.1d0*dTdid)
      Tnew = time + dT
      if(dT < minDT) stop 'stepsize underflow'
      write(1,*) 'Rejected step...'
    enddo step

! guess the time step for the next iteration
    if(errmax > ERRCON)then
      dTnext = SAFETY*dT*(errmax**PGROW)
    else
      dTnext = 5d0*dT
    endif
    time = time + dT

! save output density
    rho(1:Nbasis,0:Nsys) = tmp(1:Nbasis,0:Nsys)

  END SUBROUTINE adapt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  SUBROUTINE adapt_d(Nbasis,rho,drho,time,dTtry,eps,err,dTdid,dTnext,derivatives)

    implicit none
    integer,intent(in) :: Nbasis
    real(kind=8),intent(in) :: dTtry,eps
    real(kind=8),intent(out) :: dTdid,dTnext
    real(kind=8),intent(inout) :: time,err
    complex(kind=8),intent(inout) :: rho(Nbasis,Nbasis),drho(Nbasis,Nbasis)
    external derivatives

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! local variables

    real(kind=8) :: dT,Tnew,errmax
    complex(kind=8) :: tmp(Nbasis,Nbasis)
    real(kind=8),parameter :: SAFETY = 0.9d0
    real(kind=8),parameter :: PGROW = 0.2d0
    real(kind=8),parameter :: PSHRNK = 0.25d0
    real(kind=8),parameter :: ERRCON = 1.89d-4
    real(kind=8),parameter :: minDT = 1d-5        ! minimal time step in a.u.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! adaptive step loop
    dT = dTtry
    step: do 
! try a first step
      call CashKarp_step_d(Nbasis,rho,drho,time,dT,tmp,err,derivatives)
      errmax = err/eps
      dTdid = dT
      if(errmax < 1.) exit step
! if the error is too big, retake the step with a smaller time interval
      dT = min(SAFETY*dTdid*(errmax**PSHRNK),0.1d0*dTdid)
      Tnew = time + dT
      if(dT < minDT) stop 'stepsize underflow'
      write(1,*) 'Rejected step...'
    enddo step

! guess the time step for the next iteration
    if(errmax > ERRCON)then
      dTnext = SAFETY*dT*(errmax**PGROW)
    else
      dTnext = 5d0*dT
    endif
    time = time + dT

! save output density
    rho(1:Nbasis,1:Nbasis) = tmp(1:Nbasis,1:Nbasis)

  END SUBROUTINE adapt_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  SUBROUTINE CashKarp_step(Nbasis,Nsys,rho,drho,time,dT,rho_out,err,derivatives)

    use scalarOP, only : trace
    implicit none
    integer,intent(in) :: Nbasis,Nsys
    real(kind=8),intent(in) :: dT,time
    real(kind=8),intent(out) :: err
    complex(kind=8),intent(in) :: rho(Nbasis,0:Nsys),drho(Nbasis,0:Nsys)
    complex(kind=8),intent(out) :: rho_out(Nbasis,0:Nsys)
    external derivatives

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! parameters and variables for the embedded Cash-Karp Runge-Kutta method
    complex(kind=8) :: ak2(Nbasis,0:Nsys),ak3(Nbasis,0:Nsys),ak4(Nbasis,0:Nsys),ak5(Nbasis,0:Nsys),&
                       ak6(Nbasis,0:Nsys),tmp(Nbasis,0:Nsys)
    real(kind=8) :: A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53,&
         B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
    parameter (A2=.2d0,A3=.3d0,A4=.6d0,A5=1.d0,A6=.875d0,B21=.2d0,B31=3.d0/40.d0,B32=9.d0/40.d0,&
         B41=.3d0,B42=-.9d0,B43=1.2d0,B51=-11.d0/54.d0,B52=2.5d0,B53=-70.d0/27.d0,B54=35.d0/27.d0,&
         B61=1631.d0/55296.d0,B62=175.d0/512.d0,B63=575.d0/13824.d0,B64=44275.d0/110592.d0,B65=253.d0/4096.d0,&
         C1=37.d0/378.d0,C3=250.d0/621.d0,C4=125.d0/594.d0,C6=512.d0/1771.d0,&
         DC1=C1-2825.d0/27648.d0,DC3=C3-18575.d0/48384.d0,DC4=C4-13525.d0/55296.d0,DC5=-277.d0/14336.d0,DC6=C6-.25d0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! generating approximations of the N functions on the interval
    tmp(1:Nbasis,0:Nsys) = rho(1:Nbasis,0:Nsys) + dT*(B21*drho(1:Nbasis,0:Nsys))
    call derivatives(time+A2*dT,tmp,ak2)
    tmp(1:Nbasis,0:Nsys) = rho(1:Nbasis,0:Nsys) + dT*(B31*drho(1:Nbasis,0:Nsys)+B32*ak2(1:Nbasis,0:Nsys))
    call derivatives(time+A3*dT,tmp,ak3)
    tmp(1:Nbasis,0:Nsys) = rho(1:Nbasis,0:Nsys) + dT*(B41*drho(1:Nbasis,0:Nsys)+&
      B42*ak2(1:Nbasis,0:Nsys)+B43*ak3(1:Nbasis,0:Nsys))
    call derivatives(time+A4*dT,tmp,ak4)
    tmp(1:Nbasis,0:Nsys) = rho(1:Nbasis,0:Nsys) + dT*(B51*drho(1:Nbasis,0:Nsys)+&
      B52*ak2(1:Nbasis,0:Nsys)+B53*ak3(1:Nbasis,0:Nsys)+B54*ak4(1:Nbasis,0:Nsys))
    call derivatives(time+A5*dT,tmp,ak5)
    tmp(1:Nbasis,0:Nsys) = rho(1:Nbasis,0:Nsys) + dT*(B61*drho(1:Nbasis,0:Nsys)+&
      B62*ak2(1:Nbasis,0:Nsys)+B63*ak3(1:Nbasis,0:Nsys)+B64*ak4(1:Nbasis,0:Nsys)+B65*ak5(1:Nbasis,0:Nsys))
    call derivatives(time+A6*dT,tmp,ak6)

! combining the approximations to integrate the N functions at the final time
    rho_out(1:Nbasis,0:Nsys) = rho(1:Nbasis,0:Nsys) + dT*(C1*drho(1:Nbasis,0:Nsys)+&
      C3*ak3(1:Nbasis,0:Nsys)+C4*ak4(1:Nbasis,0:Nsys)+C6*ak6(1:Nbasis,0:Nsys))
!       normc = trace(rho_out(1,1),Nbasis,Nsys)
!       rho_out(1:Nbasis,0:Nsys) = rho_out(1:Nbasis,0:Nsys)/normc

! evaluating the errors
    tmp(1:Nbasis,0:Nsys) = dT*(DC1*drho(1:Nbasis,0:Nsys)+DC3*ak3(1:Nbasis,0:Nsys)+&
      DC4*ak4(1:Nbasis,0:Nsys)+DC5*ak5(1:Nbasis,0:Nsys)+DC6*ak6(1:Nbasis,0:Nsys))
      err = maxval(abs(tmp(1:Nbasis,0:Nsys)))/maxval(abs(rho_out(1:Nbasis,0:Nsys)))
!     err = 0d0
!     do i = 1,Nbasis
!     do j = 0,Nsys
!       err = max(abs(tmp(i,j))/max(abs(dT*drho(i,j)),1d-7),err)
!     enddo
!     enddo
!    err = maxval(abs(tmp(1:Nbasis,0:Nsys)))/maxval(abs(rho_out(1:Nbasis,0:Nsys)))

  END SUBROUTINE CashKarp_step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  SUBROUTINE CashKarp_step_d(Nbasis,rho,drho,time,dT,rho_out,err,derivatives)

    use scalarOP, only : trace
    implicit none
    integer,intent(in) :: Nbasis
    real(kind=8),intent(in) :: dT,time
    real(kind=8),intent(out) :: err
    complex(kind=8),intent(in) :: rho(Nbasis,Nbasis),drho(Nbasis,Nbasis)
    complex(kind=8),intent(out) :: rho_out(Nbasis,Nbasis)
    external derivatives

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! parameters and variables for the embedded Cash-Karp Runge-Kutta method
    complex(kind=8) :: ak2(Nbasis,Nbasis),ak3(Nbasis,Nbasis),ak4(Nbasis,Nbasis),ak5(Nbasis,Nbasis),&
                       ak6(Nbasis,Nbasis),tmp(Nbasis,Nbasis)
    real(kind=8) :: A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53,&
         B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
    parameter (A2=.2d0,A3=.3d0,A4=.6d0,A5=1.d0,A6=.875d0,B21=.2d0,B31=3.d0/40.d0,B32=9.d0/40.d0,&
         B41=.3d0,B42=-.9d0,B43=1.2d0,B51=-11.d0/54.d0,B52=2.5d0,B53=-70.d0/27.d0,B54=35.d0/27.d0,&
         B61=1631.d0/55296.d0,B62=175.d0/512.d0,B63=575.d0/13824.d0,B64=44275.d0/110592.d0,B65=253.d0/4096.d0,&
         C1=37.d0/378.d0,C3=250.d0/621.d0,C4=125.d0/594.d0,C6=512.d0/1771.d0,&
         DC1=C1-2825.d0/27648.d0,DC3=C3-18575.d0/48384.d0,DC4=C4-13525.d0/55296.d0,DC5=-277.d0/14336.d0,DC6=C6-.25d0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! generating approximations of the N functions on the interval
    tmp(1:Nbasis,1:Nbasis) = rho(1:Nbasis,1:Nbasis) + dT*(B21*drho(1:Nbasis,1:Nbasis))
    call derivatives(time+A2*dT,tmp,ak2)
    tmp(1:Nbasis,1:Nbasis) = rho(1:Nbasis,1:Nbasis) + dT*(B31*drho(1:Nbasis,1:Nbasis)+B32*ak2(1:Nbasis,1:Nbasis))
    call derivatives(time+A3*dT,tmp,ak3)
    tmp(1:Nbasis,1:Nbasis) = rho(1:Nbasis,1:Nbasis) + dT*(B41*drho(1:Nbasis,1:Nbasis)+&
      B42*ak2(1:Nbasis,1:Nbasis)+B43*ak3(1:Nbasis,1:Nbasis))
    call derivatives(time+A4*dT,tmp,ak4)
    tmp(1:Nbasis,1:Nbasis) = rho(1:Nbasis,1:Nbasis) + dT*(B51*drho(1:Nbasis,1:Nbasis)+&
      B52*ak2(1:Nbasis,1:Nbasis)+B53*ak3(1:Nbasis,1:Nbasis)+B54*ak4(1:Nbasis,1:Nbasis))
    call derivatives(time+A5*dT,tmp,ak5)
    tmp(1:Nbasis,1:Nbasis) = rho(1:Nbasis,1:Nbasis) + dT*(B61*drho(1:Nbasis,1:Nbasis)+&
      B62*ak2(1:Nbasis,1:Nbasis)+B63*ak3(1:Nbasis,1:Nbasis)+B64*ak4(1:Nbasis,1:Nbasis)+B65*ak5(1:Nbasis,1:Nbasis))
    call derivatives(time+A6*dT,tmp,ak6)

! combining the approximations to integrate the N functions at the final time
    rho_out(1:Nbasis,1:Nbasis) = rho(1:Nbasis,1:Nbasis) + dT*(C1*drho(1:Nbasis,1:Nbasis)+&
      C3*ak3(1:Nbasis,1:Nbasis)+C4*ak4(1:Nbasis,1:Nbasis)+C6*ak6(1:Nbasis,1:Nbasis))
!       normc = trace(rho_out(1,1),Nbasis,Nsys)
!       rho_out(1:Nbasis,0:Nsys) = rho_out(1:Nbasis,0:Nsys)/normc

! evaluating the errors
    tmp(1:Nbasis,1:Nbasis) = dT*(DC1*drho(1:Nbasis,1:Nbasis)+DC3*ak3(1:Nbasis,1:Nbasis)+&
      DC4*ak4(1:Nbasis,1:Nbasis)+DC5*ak5(1:Nbasis,1:Nbasis)+DC6*ak6(1:Nbasis,1:Nbasis))
      err = maxval(abs(tmp(1:Nbasis,1:Nbasis)))/maxval(abs(rho_out(1:Nbasis,1:Nbasis)))
!     err = 0d0
!     do i = 1,Nbasis
!     do j = 0,Nsys
!       err = max(abs(tmp(i,j))/max(abs(dT*drho(i,j)),1d-7),err)
!     enddo
!     enddo
!    err = maxval(abs(tmp(1:Nbasis,0:Nsys)))/maxval(abs(rho_out(1:Nbasis,0:Nsys)))

  END SUBROUTINE CashKarp_step_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
END MODULE RungeKutta
