!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module scalarOP

  implicit none
  private
  integer :: i,j,k
  real(kind=8) :: zero = 0d0
  complex(kind=8) :: dumc

  public trace,normalize

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind=8) function trace(mat,Nb,Ns)

    implicit none
    integer,intent(in) :: Nb,Ns
    complex(kind=8),intent(in) :: mat(0:Nb,0:Ns)
    complex(kind=8) :: dum

    dumc = dcmplx(zero,zero)
    do i = 0,Ns-1
      dumc = dumc + mat(i,i)
    enddo
    do i = Ns,Nb
      dumc = dumc + mat(i,Ns)
    enddo
    trace = dumc

  end function trace

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine normalize(normc,mat,Nb,Ns)

    implicit none
    integer,intent(in) :: Nb,Ns
    complex(kind=8),intent(in) :: normc
    complex(kind=8),intent(inout) :: mat(0:Nb,0:Ns)
    complex(kind=8) :: dum

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module scalarOP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
