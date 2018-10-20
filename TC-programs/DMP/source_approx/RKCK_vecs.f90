!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! contains the subroutines for propagating a set of N vectors
! over a time interval dT using a Runge-Kutta Cash-Karp integrator 
! with adaptive time stpe
!
! original routines from Numerical Recipes, adapted and bundled by 
! Jean Christophe Tremblay (22.06.2009)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
MODULE RungeKutta

  implicit none
  private
  public adapt,cashKarp_step

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
CONTAINS
  SUBROUTINE adapt(Nbasis,vecs,dvecs,N,time,dTtry,eps,err,dTdid,dTnext,derivatives)
!
!  Nbasis = size of a vector;
!  vecs = array of dimension (Nbasis,N) containing the N vectors at the beginning of the interval;
!  dvecs = array of dimension (Nbasis,N) containing the N vector derivatives at the beginning of the interval;
!  N = number of vectors;
!  time = actual time spent since the begining of the propagation;
!  vecs_out = array of dimension (Nbasis,N) containing the N vectors at the end of the interval;
!  err = vector of dimension (N) containing the error estimate of each vector
!  derivatives = name of the subroutine for computing the derivatives
!  dTtry = trial time interval;
!  dTdid = actual time interval;
!  dTnext = trial time interval for the next step;
!  eps = relative error tolerance
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none
    integer,intent(in) :: Nbasis,N
    real(kind=8),intent(in) :: dTtry,eps
    real(kind=8),intent(out) :: dTdid,dTnext
    real(kind=8),intent(inout) :: time,err(N)
    complex(kind=8),intent(inout) :: vecs(Nbasis,N),dvecs(Nbasis,N)
    external derivatives

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! local variables

    integer :: i
    real(kind=8) :: dT,Tnew,errmax
    complex(kind=8) :: tmp(Nbasis,N)
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
      call CashKarp_step(Nbasis,vecs,dvecs,N,time,dT,tmp,err,derivatives)
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

! save output density
    do i = 1,N
      vecs(:,i) = tmp(:,i)
    enddo

  END SUBROUTINE adapt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  SUBROUTINE CashKarp_step(Nbasis,vecs,dvecs,N,time,dT,vecs_out,err,derivatives)
!
!  Nbasis = size of a vector;
!  vecs = array of dimension (Nbasis,N) containing the N vectors at the beginning of the interval;
!  dvecs = array of dimension (Nbasis,N) containing the N vector derivatives at the beginning of the interval;
!  N = number of vectors;
!  time = actual time spent since the begining of the propagation;
!  dT = time interval for the RK step;
!  vecs_out = array of dimension (Nbasis,N) containing the N vectors at the end of the interval;
!  err = vector of dimension (N) containing the error estimate of each vector
!  derivatives = name of the subroutine for computing the derivatives
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use scalarOP, only : trace
    implicit none
    integer,intent(in) :: Nbasis,N
    real(kind=8),intent(in) :: dT,time
    real(kind=8),intent(out) :: err(N)
    complex(kind=8),intent(in) :: vecs(Nbasis,N),dvecs(Nbasis,N)
    complex(kind=8),intent(out) :: vecs_out(Nbasis,N)
    external derivatives

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! parameters and variables for the embedded Cash-Karp Runge-Kutta method
    integer :: i
    complex(kind=8) :: ak2(Nbasis,N),ak3(Nbasis,N),ak4(Nbasis,N),ak5(Nbasis,N),&
                       ak6(Nbasis,N),tmp(Nbasis,N)
    real(kind=8) :: A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53,&
         B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
    parameter (A2=.2d0,A3=.3d0,A4=.6d0,A5=1.d0,A6=.875d0,B21=.2d0,B31=3.d0/40.d0,B32=9.d0/40.d0,&
         B41=.3d0,B42=-.9d0,B43=1.2d0,B51=-11.d0/54.d0,B52=2.5d0,B53=-70.d0/27.d0,B54=35.d0/27.d0,&
         B61=1631.d0/55296.d0,B62=175.d0/512.d0,B63=575.d0/13824.d0,B64=44275.d0/110592.d0,B65=253.d0/4096.d0,&
         C1=37.d0/378.d0,C3=250.d0/621.d0,C4=125.d0/594.d0,C6=512.d0/1771.d0,&
         DC1=C1-2825.d0/27648.d0,DC3=C3-18575.d0/48384.d0,DC4=C4-13525.d0/55296.d0,DC5=-277.d0/14336.d0,DC6=C6-.25d0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! generating approximations of the N functions on the interval
    do i = 1,N
       tmp(:,i) = vecs(:,i) + dT*(B21*dvecs(:,i))
    end do
    call derivatives(time+A2*dT,tmp,ak2,Nbasis,N)
    do i = 1,N
       tmp(:,i) = vecs(:,i) + dT*(B31*dvecs(:,i)+B32*ak2(:,i))
    end do
    call derivatives(time+A3*dT,tmp,ak3,Nbasis,N)
    do i = 1,N
       tmp(:,i) = vecs(:,i) + dT*(B41*dvecs(:,i)+B42*ak2(:,i)+B43*ak3(:,i))
    end do
    call derivatives(time+A4*dT,tmp,ak4,Nbasis,N)
    do i = 1,N
       tmp(:,i) = vecs(:,i) + dT*(B51*dvecs(:,i)+B52*ak2(:,i)+B53*ak3(:,i)+B54*ak4(:,i))
    end do
    call derivatives(time+A5*dT,tmp,ak5,Nbasis,N)
    do i = 1,N
       tmp(:,i) = vecs(:,i) + dT*(B61*dvecs(:,i)+B62*ak2(:,i)+B63*ak3(:,i)+B64*ak4(:,i)+B65*ak5(:,i))
    end do
    call derivatives(time+A6*dT,tmp,ak6,Nbasis,N)

! combining the approximations to integrate the N functions at the final time
    do i = 1,N
       vecs_out(:,i) = vecs(:,i) + dT*(C1*dvecs(:,i)+C3*ak3(:,i)+C4*ak4(:,i)+C6*ak6(:,i))
    end do

! evaluating the errors
    do i = 1,N
       tmp(:,i) = dT*(DC1*dvecs(:,i)+DC3*ak3(:,i)+DC4*ak4(:,i)+DC5*ak5(:,i)+DC6*ak6(:,i))
       err(i) = maxval(abs(tmp(:,i)))/maxval(abs(vecs_out(:,i)))
    end do

  END SUBROUTINE CashKarp_step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
END MODULE RungeKutta
