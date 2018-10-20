!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! contains the subroutines for propagating a set of N state vectors
! over a time interval dT using a Runge-Kutta Cash-Karp integrator 
! with adaptive time stpe
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
MODULE rungekutta

  implicit none
  private
  public adapt,cashKarp_step

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
CONTAINS
  SUBROUTINE adapt(BasisSize,psi,dpsi,psi_out,Nequations,time,dTtry,eps,err,dTdid,dTnext,derivatives)

    implicit none
    integer,intent(in) :: BasisSize,Nequations
    real(kind=8),intent(in) :: dTtry,eps
    real(kind=8),intent(out) :: dTdid,dTnext
    real(kind=8),intent(inout) :: time,err(Nequations)
    complex(kind=8),intent(in) :: psi(BasisSize,Nequations),dpsi(BasisSize,Nequations)
    complex(kind=8),intent(out) :: psi_out(BasisSize,Nequations)
    external derivatives

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! local variables

    real(kind=8) :: dT,Tnew,errmax
    real(kind=8),parameter :: SAFETY = 0.9d0
    real(kind=8),parameter :: PGROW = 0.2d0
    real(kind=8),parameter :: PSHRNK = 0.25d0
    real(kind=8),parameter :: ERRCON = 1.89d-4
    real(kind=8),parameter :: minDT = 1d-10
!    real(kind=8),parameter :: minDT = 1d-5        ! minimal time step in a.u.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! adaptive step loop
    dT = dTtry
    step: do 
! try a first step
      call CashKarp_step(BasisSize,psi,dpsi,Nequations,time,dT,psi_out,err,derivatives)
      errmax = maxval(err)/eps
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

  END SUBROUTINE adapt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  SUBROUTINE CashKarp_step(BasisSize,psi,dpsi,Nequations,time,dT,psi_out,err,derivatives)

    implicit none
    integer,intent(in) :: BasisSize,Nequations
    real(kind=8),intent(in) :: dT,time
    real(kind=8),intent(out) :: err(Nequations)
    complex(kind=8),intent(in) :: psi(BasisSize,Nequations),dpsi(BasisSize,Nequations)
    complex(kind=8),intent(out) :: psi_out(BasisSize,Nequations)
    external derivatives

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! parameters and variables for the embedded Cash-Karp Runge-Kutta method
    integer :: i
    complex(kind=8) :: ak2(BasisSize,Nequations),ak3(BasisSize,Nequations),ak4(BasisSize,Nequations),&
                       ak5(BasisSize,Nequations),ak6(BasisSize,Nequations),tmp(BasisSize,Nequations)
    real(kind=8) :: A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53,&
         B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
    parameter (A2=.2d0,A3=.3d0,A4=.6d0,A5=1.d0,A6=.875d0,B21=.2d0,B31=3.d0/40.d0,B32=9.d0/40.d0,&
         B41=.3d0,B42=-.9d0,B43=1.2d0,B51=-11.d0/54.d0,B52=2.5d0,B53=-70.d0/27.d0,B54=35.d0/27.d0,&
         B61=1631.d0/55296.d0,B62=175.d0/512.d0,B63=575.d0/13824.d0,B64=44275.d0/110592.d0,B65=253.d0/4096.d0,&
         C1=37.d0/378.d0,C3=250.d0/621.d0,C4=125.d0/594.d0,C6=512.d0/1771.d0,&
         DC1=C1-2825.d0/27648.d0,DC3=C3-18575.d0/48384.d0,DC4=C4-13525.d0/55296.d0,DC5=-277.d0/14336.d0,DC6=C6-.25d0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! generating approximations of the Nequations functions on the interval
     tmp(:,:) = psi(:,:) + dT*(B21*dpsi(:,:))
     call derivatives(time+A2*dT,tmp,ak2,Nequations)
     tmp(:,:) = psi(:,:) + dT*(B31*dpsi(:,:)+B32*ak2(:,:))
     call derivatives(time+A3*dT,tmp,ak3,Nequations)
     tmp(:,:) = psi(:,:) + dT*(B41*dpsi(:,:)+B42*ak2(:,:)+B43*ak3(:,:))
     call derivatives(time+A4*dT,tmp,ak4,Nequations)
     tmp(:,:) = psi(:,:) + dT*(B51*dpsi(:,:)+B52*ak2(:,:)+B53*ak3(:,:)+B54*ak4(:,:))
     call derivatives(time+A5*dT,tmp,ak5,Nequations)
     tmp(:,:) = psi(:,:) + dT*(B61*dpsi(:,:)+B62*ak2(:,:)+B63*ak3(:,:)+B64*ak4(:,:)+B65*ak5(:,:))
     call derivatives(time+A6*dT,tmp,ak6,Nequations)

! combining the approximations to integrate the N functions at the final time
     psi_out(:,:) = psi(:,:) + dT*(C1*dpsi(:,:)+C3*ak3(:,:)+C4*ak4(:,:)+C6*ak6(:,:))

! evaluating the errors
     tmp(:,:) = dT*(DC1*dpsi(:,:)+DC3*ak3(:,:)+DC4*ak4(:,:)+DC5*ak5(:,:)+DC6*ak6(:,:))
     do i = 1,Nequations
       err(i) = maxval(abs(tmp(:,i)))/maxval(abs(psi_out(:,i)))
     enddo

  END SUBROUTINE CashKarp_step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
END MODULE RungeKutta
