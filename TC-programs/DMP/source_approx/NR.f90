MODULE rkqsmod
  USE rkckmod
  implicit none
  REAL*8::SAFETY,PGROW,PSHRNK,ERRCON
  PARAMETER (SAFETY=0.9d0,PGROW=-.2d0,PSHRNK=-.25d0,ERRCON=1.89d-4)
CONTAINS
  SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
    INTEGER,intent(in)::n
    REAL*8,intent(in)::htry,eps,yscal(:)
    REAL*8,intent(inout)::dydx(:),y(:),x
    REAL*8,intent(out)::hdid,hnext
    INTERFACE
       SUBROUTINE derivs(x,y,dydx)
         REAL*8,intent(in)::x,y(:)
         REAL*8,intent(out)::dydx(:)
       END SUBROUTINE derivs
    END INTERFACE

    !CU    USES derivs,rkck
    INTEGER i
    REAL*8 errmax,h,xnew,yerr(n),ytemp(n)
    h=htry
1   call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)
    errmax=0.
    do 11 i=1,n
       errmax=max(errmax,abs(yerr(i)/yscal(i)))
11  continue
    errmax=errmax/eps
    hdid=h
    if(errmax.gt.1.)then
       h=SAFETY*hdid*(errmax**PSHRNK)
       if(h.lt.0.1*hdid)then
          h=.1*hdid
       endif
       xnew=x+h
       if(xnew.eq.x) then 
          print*,'x=',x,'h=',h,'htry=',htry,'errmax=',errmax
          print*,'yerr=',yerr(1:n)
          print*,'yscal=',yscal(1:n)
          pause 'stepsize underflow in rkqs'
       end if
       goto 1
    else
       if(errmax.gt.ERRCON)then
          hnext=SAFETY*h*(errmax**PGROW)
       else
          hnext=5.*h
       endif
       x=x+h
       do 12 i=1,n
          y(i)=ytemp(i)
12     continue
       return
    endif
  END SUBROUTINE rkqs
  !C     (C) Copr. 1986-92 Numerical Recipes Software 0"0)P.

  SUBROUTINE rkqsC(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
    INTEGER,intent(in)::n
    COMPLEX*16,intent(in)::htry
    real*8,intent(in)::yscal(:),eps
    COMPLEX*16,intent(inout)::dydx(:),y(:),x
    COMPLEX*16,intent(out)::hdid,hnext
    INTERFACE
       SUBROUTINE derivs(x,y,dydx)
         COMPLEX*16,intent(in)::x,y(:)
         COMPLEX*16,intent(out)::dydx(:)
       END SUBROUTINE derivs
    END INTERFACE

    INTEGER::i
    REAL*8::errmax
    COMPLEX*16:: h,xnew,yerr(n),ytemp(n)

    h=htry
1   call rkckC(y,dydx,n,x,h,ytemp,yerr,derivs)
    errmax=0.
    do 11 i=1,n
       errmax=max(errmax,abs(yerr(i)/yscal(i)))
11  continue
    errmax=errmax/eps
    hdid=h
    if(errmax.gt.1.)then
       h=SAFETY*hdid*(errmax**PSHRNK)
       if(abs(h)<0.1*abs(hdid)) then
          h=.1*hdid
       endif
       xnew=x+h
       if(xnew.eq.x)pause 'stepsize underflow in rkqs'
       goto 1
    else
       if(errmax.gt.ERRCON)then
          hnext=SAFETY*h*(errmax**PGROW)
       else
          hnext=5.*h
       endif
       x=x+h
       do 12 i=1,n
          y(i)=ytemp(i)
12     continue
       return
    endif
  END SUBROUTINE rkqsC
  !C     (C) Copr. 1986-92 Numerical Recipes Software 0"0)P.
END MODULE rkqsmod





MODULE rkckmod
CONTAINS
  SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr,derivs)
    INTEGER::n
    REAL*8::h,x,dydx(n),y(n),yerr(n),yout(n)
    INTERFACE
       SUBROUTINE derivs(x,y,dydx)
         REAL*8,intent(in)::x,y(:)
         REAL*8,intent(out)::dydx(:)
       END SUBROUTINE derivs
    END INTERFACE

    !CU    USES derivs
    INTEGER::i
    REAL*8::ak2(n),ak3(n),ak4(n),ak5(n),ak6(n),ytemp(n)
    REAL*8,PARAMETER::A2=.2d0,A3=.3d0,A4=.6d0,A5=1.d0,A6=.875d0,B21=.2d0,&
         B31=3./40.d0,&
         B32=9.d0/40.d0,B41=.3d0,B42=-.9d0,B43=1.2d0,B51=-11.d0/54.d0,&
         B52=2.5d0,&
         B53=-70.d0/27.d0,B54=35.d0/27.d0,B61=1631.d0/55296.d0,&
         B62=175.d0/512.d0,&
         B63=575.d0/13824.d0,B64=44275.d0/110592.d0,B65=253.d0/4096.d0&
         ,C1=37.d0/378.d0,&
         C3=250.d0/621.d0,C4=125.d0/594.d0,C6=512.d0/1771.d0,&
         DC1=C1-2825.d0/27648.d0,&
         DC3=C3-18575.d0/48384.d0,DC4=C4-13525.d0/55296.d0&
         ,DC5=-277.d0/14336.d0,&
         DC6=C6-.25d0
    do 11 i=1,n
       ytemp(i)=y(i)+B21*h*dydx(i)
11  continue
    call derivs(x+A2*h,ytemp,ak2)
    do 12 i=1,n
       ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
12  continue
    call derivs(x+A3*h,ytemp,ak3)
    do 13 i=1,n
       ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
13  continue
    call derivs(x+A4*h,ytemp,ak4)
    do 14 i=1,n
       ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
14  continue
    call derivs(x+A5*h,ytemp,ak5)
    do 15 i=1,n
       ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+&
            B65*ak5(i))
15  continue
    call derivs(x+A6*h,ytemp,ak6)
    do 16 i=1,n
       yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
16  continue
    do 17 i=1,n
       yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*&
            ak6(i))
17  continue
    return
  END SUBROUTINE rkck
!C  (C) Copr. 1986-92 Numerical Recipes Software 0"0)P.

  SUBROUTINE rkckC(y,dydx,n,x,h,yout,yerr,derivs)
    INTEGER,intent(in)::n
    complex*16,intent(inout)::h,x,dydx(n),y(n),yerr(n),yout(n)
    INTERFACE
       SUBROUTINE derivs(x,y,dydx)
         COMPLEX*16,intent(in)::x,y(:)
         COMPLEX*16,intent(out)::dydx(:)
       END SUBROUTINE derivs
    END INTERFACE

    !CU    USES derivs
    INTEGER::i
    COMPLEX*16::ak2(n),ak3(n),ak4(n),ak5(n),ak6(n),ytemp(n)
    REAL*8::A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53,&
         B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
    PARAMETER (A2=.2d0,A3=.3d0,A4=.6d0,A5=1.d0,A6=.875d0,B21=.2d0,&
         B31=3./40.d0,&
         B32=9.d0/40.d0,B41=.3d0,B42=-.9d0,B43=1.2d0,B51=-11.d0/54.d0,&
         B52=2.5d0,&
         B53=-70.d0/27.d0,B54=35.d0/27.d0,B61=1631.d0/55296.d0,&
         B62=175.d0/512.d0,&
         B63=575.d0/13824.d0,B64=44275.d0/110592.d0,B65=253.d0/4096.d0&
         ,C1=37.d0/378.d0,&
         C3=250.d0/621.d0,C4=125.d0/594.d0,C6=512.d0/1771.d0,&
         DC1=C1-2825.d0/27648.d0,&
         DC3=C3-18575.d0/48384.d0,DC4=C4-13525.d0/55296.d0&
         ,DC5=-277.d0/14336.d0,&
         DC6=C6-.25d0)

!    PARAMETER (A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40.,&
!         B32=9./40.,B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5,&
!         B53=-70./27.,B54=35./27.,B61=1631./55296.,B62=175./512.,&
!         B63=575./13824.,B64=44275./110592.,B65=253./4096.,C1=37./378.,&
!         C3=250./621.,C4=125./594.,C6=512./1771.,DC1=C1-2825./27648.,&
!         DC3=C3-18575./48384.,DC4=C4-13525./55296.,DC5=-277./14336.,&
!         DC6=C6-.25)
    do i=1,n
       ytemp(i)=y(i)+B21*h*dydx(i)
    end do
    call derivs(x+A2*h,ytemp,ak2)
    do i=1,n
       ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
    end do
    call derivs(x+A3*h,ytemp,ak3)
    do i=1,n
       ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
    end do
    call derivs(x+A4*h,ytemp,ak4)
    do i=1,n
       ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
    end do
    call derivs(x+A5*h,ytemp,ak5)
    do i=1,n
       ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+&
            B65*ak5(i))
    end do
    call derivs(x+A6*h,ytemp,ak6)
    do i=1,n
       yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
    end do
    do i=1,n
       yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*&
            ak6(i))
    end do
    return
  END SUBROUTINE rkckC
!C  (C) Copr. 1986-92 Numerical Recipes Software 0"0)P.
END MODULE rkckmod
