!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module field

  private
  integer :: NfieldIN
  real(kind=8),parameter :: zero = 0d0
  real(kind=8),public :: Emax
  real(kind=8),allocatable :: Tin(:),FIELDin(:,:),coeff_lin(:,:),coeff_quad(:,:),coeff_cub(:,:)
  public ElectricField,setfield

contains  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ElectricField(time,fieldOUT,npol)

    implicit none
    integer :: ipol
    integer,intent(in) :: npol
    real(kind=8),intent(in) :: time
    real(kind=8),intent(out) :: fieldOUT(npol)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! security check
    if(time > Tin(NfieldIN)) then
      fieldOUT(:) = zero
    else
      do ipol = 1,npol
        fieldOUT(ipol) = seval(NfieldIN,time,Tin,FIELDin(1,ipol),&
         coeff_lin(1,ipol),coeff_quad(1,ipol),coeff_cub(1,ipol))
      enddo
    endif

  end subroutine ElectricField
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine setfield(polarization,Npol)

    implicit none
    integer,intent(in) :: Npol
    character(len=3),intent(in) :: polarization
    integer :: i,ipol
    real(kind=8) :: Tdum,dum(npol)

! read the field
    read(4,*) NfieldIN,Tdum
    allocate(Tin(NfieldIN),FIELDin(NfieldIN,Npol))
    allocate(coeff_lin(NfieldIN,Npol),coeff_quad(NfieldIN,Npol),coeff_cub(NfieldIN,Npol))
    do i = 1,NfieldIN
      read(4,*) Tdum,(dum(ipol),ipol=1,Npol)
      Tin(i) = Tdum
! store the components of the electric field
      ipol = 0
      if(polarization(1:1) == 'x') then
        ipol = ipol + 1
        FIELDin(i,ipol) = dum(ipol)
      endif
      if(polarization(2:2) == 'y') then
        ipol = ipol + 1
        FIELDin(i,ipol) = dum(ipol)
      endif
      if(polarization(3:3) == 'z') then
        ipol = ipol + 1
        FIELDin(i,ipol) = dum(ipol)
      endif
    enddo

! fit the field
    write(1,*) 'Fitting the field components using cubic splines'
    call flush(1)
    do ipol = 1,Npol
      call spline(NfieldIN,Tin,FIELDin(1,ipol),coeff_lin(1,ipol),coeff_quad(1,ipol),coeff_cub(1,ipol))
    enddo
    write(1,*) 'Finished fit using cubic splines'
    call flush(1)

  end subroutine setfield
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine spline (n, x, y, b, c, d)
    integer n
    real(kind=8) ::x(n), y(n), b(n), c(n), d(n)
!
!  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
!  for a cubic interpolating spline
!
!    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!
!    for  x(i) .le. x .le. x(i+1)
!
!  input..
!
!    n = the number of data points or knots (n.ge.2)
!    x = the abscissas of the knots in strictly increasing order
!    y = the ordinates of the knots
!
!  output..
!
!    b, c, d  = arrays of spline coefficients as defined above.
!
!  using  p  to denote differentiation,
!
!    y(i) = s(x(i))
!    b(i) = sp(x(i))
!    c(i) = spp(x(i))/2
!    d(i) = sppp(x(i))/6  (derivative from the right)
!
!  the accompanying function subprogram  seval  can be used
!  to evaluate the spline.
!
!
    integer nm1, ib, i
    real(kind=8) ::t

    nm1 = n-1
    if (n < 2) then
      return
    elseif ( n == 2 ) then
      b(1) = (y(2)-y(1))/(x(2)-x(1))
      c(1) = 0d0
      d(1) = 0d0
      b(2) = b(1)
      c(2) = 0d0
      d(2) = 0d0
      return
    else
!
!  set up tridiagonal system
!  b = diagonal, d = offdiagonal, c = right hand side.
!
      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      do i = 2, nm1
        d(i) = x(i+1) - x(i)
        b(i) = 2d0*(d(i-1) + d(i))
        c(i+1) = (y(i+1) - y(i))/d(i)
        c(i) = c(i+1) - c(i)
      enddo
!
!  end conditions.  third derivatives at  x(1)  and  x(n)
!  obtained from divided differences
!
      b(1) = -d(1)
      b(n) = -d(nm1)
      c(1) = 0d0
      c(n) = 0d0
      if ( n /= 3 ) then
        c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
        c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
        c(1) = c(1)*d(1)**2/(x(4)-x(1))
        c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
      endif
!
!  forward elimination
!
      do i = 2, n
        t = d(i-1)/b(i-1)
        b(i) = b(i) - t*d(i-1)
        c(i) = c(i) - t*c(i-1)
      enddo
!
!  back substitution
!
      c(n) = c(n)/b(n)
      do ib = 1, nm1
        i = n-ib
        c(i) = (c(i) - d(i)*c(i+1))/b(i)
      enddo
!
!  c(i) is now the sigma(i) of the text
!
!  compute polynomial coefficients
!
      b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2d0*c(n))
      do i = 1, nm1
        b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2d0*c(i))
        d(i) = (c(i+1) - c(i))/d(i)
        c(i) = 3d0*c(i)
      enddo
      c(n) = 3d0*c(n)
      d(n) = d(nm1)
      return
    endif

  end subroutine spline
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function seval( n, u, x, y, b, c, d)

    implicit none
    integer n
    real(kind=8),intent(in) :: u, x(n), y(n), b(n), c(n), d(n)
    real(kind=8) :: seval
!
!  this subroutine evaluates the cubic spline function
!
!    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
!
!    where  x(i) .lt. u .lt. x(i+1), using horner's rule
!
!  if  u .lt. x(1) then  i = 1  is used.
!  if  u .ge. x(n) then  i = n  is used.
!
!  input..
!
!    n = the number of data points
!    u = the abscissa at which the spline is to be evaluated
!    x,y = the arrays of data abscissas and ordinates
!    b,c,d = arrays of spline coefficients computed by spline
!
!  if  u  is not in the same interval as the previous call, then a
!  binary search is performed to determine the proper interval.
!
    integer i, j, k
    real(kind=8) ::dx
    data i /1/
    if ( i >= n ) i = 1
!
!  binary search
!
    if( u < x(i) ) then
      j = i+1
      i = 1
    elseif( u >= x(i+1) ) then
      j = n + 1
    else
      j = i+1
    endif
    loop20 : do 
      if ( j == i+1 ) exit loop20
      k = (i+j)/2
      if ( u < x(k) ) j = k
      if ( u >= x(k) ) i = k
    enddo loop20
!
!  evaluate spline
!
    dx = u - x(i)
    seval = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
    return

  end function seval
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module field
