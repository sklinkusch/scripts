!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module scalarOP

  implicit none
  private
  integer :: i,j,k
  real(kind=8) :: zero = 0d0, eps=1d-10, one=1.0d0
  real(kind=8),parameter :: kBm1=3.1577465d5

  complex(kind=8) :: dumc

  public trace,trace_d,normalize,normalize_d,bose,boltz

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind=8) function trace(mat,Nb,Ns)

    implicit none
    integer,intent(in) :: Nb,Ns
    complex(kind=8),intent(in) :: mat(0:Nb,0:Ns)

    dumc = dcmplx(zero,zero)
    do i = 0,Ns-1
      dumc = dumc + mat(i,i)
    enddo
    do i = Ns,Nb
      dumc = dumc + mat(i,Ns)
    enddo
    trace = dumc

  end function trace

  complex(kind=8) function trace_d(mat,Nb)

      implicit none
      integer,intent(in) :: Nb
      complex(kind=8),intent(in) :: mat(0:Nb,0:Nb)

      dumc = dcmplx(zero,zero)
      do i=0,Nb
        dumc = dumc + mat(i,i)
      enddo
      trace_d = dumc

      end function trace_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine normalize(normc,mat,Nb,Ns)

    implicit none
    integer,intent(in) :: Nb,Ns
    complex(kind=8),intent(in) :: normc
    complex(kind=8),intent(inout) :: mat(0:Nb,0:Ns)

    dumc = 1d0/normc
    do j = 0,Ns-1
    do i = 0,Nb
      mat(i,j) = dumc*mat(i,j)
    enddo
    enddo
    do i = Ns,Nb
      mat(i,Ns) = dumc*mat(i,Ns)
    enddo

  end subroutine normalize

  subroutine normalize_d(normc,mat,Nb)

      implicit none
      integer,intent(in) :: Nb
      complex(kind=8),intent(in) :: normc
      complex(kind=8),intent(inout) :: mat(0:(Nb-1),0:(Nb-1))

      dumc = 1d0/normc
      do j = 0, Nb-1
        do i = 0, Nb-1
          mat(i,j) = dumc*mat(i,j)
        enddo
      enddo

    end subroutine normalize_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine bose(bosefkt, deltaE, T)
     implicit none
     real(kind=8), intent(in) :: deltaE, T
     real(kind=8), intent(out) :: bosefkt
     real(kind=8) :: dum

     dum = Kbm1 / T
     if(T < eps .and. deltaE >= zero) then
         bosefkt = zero
     elseif(T < eps .and. deltaE < zero) then
         bosefkt = -one
     else
         bosefkt = one/(exp(dum*deltaE)-one)
     endif

  end subroutine bose

  subroutine boltz(bosefkt, deltaE, T)
     implicit none
     real(kind=8), intent(in) :: deltaE, T
     real(kind=8), intent(out) :: bosefkt

     if(T < eps) then
         bosefkt = zero
     else
         bosefkt = exp(-deltaE*kBm1/T)
     endif

  end subroutine boltz


end module scalarOP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
