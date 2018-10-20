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
    c(1) = 0.
    d(1) = 0.
    b(2) = b(1)
    c(2) = 0.
    d(2) = 0.
    return
  else

!
!  set up tridiagonal system
!
!  b = diagonal, d = offdiagonal, c = right hand side.
!
    d(1) = x(2) - x(1)
    c(2) = (y(2) - y(1))/d(1)
    do i = 2, nm1
      d(i) = x(i+1) - x(i)
      b(i) = 2.*(d(i-1) + d(i))
      c(i+1) = (y(i+1) - y(i))/d(i)
      c(i) = c(i+1) - c(i)
    enddo
!
!  end conditions.  third derivatives at  x(1)  and  x(n)
!  obtained from divided differences
!
    b(1) = -d(1)
    b(n) = -d(n-1)
    c(1) = 0.
    c(n) = 0.
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
    b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.*c(n))
    do i = 1, nm1
      b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.*c(i))
      d(i) = (c(i+1) - c(i))/d(i)
      c(i) = 3.*c(i)
    enddo
    c(n) = 3.*c(n)
    d(n) = d(n-1)
    return
  endif

end subroutine spline

real(kind=8) function seval(n, u, x, y, b, c, d)

  implicit none
  integer n
  real(kind=8) :: u, x(n), y(n), b(n), c(n), d(n)
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
  data i/1/
  if ( i >= n ) i = 1
  if ( u < x(i) ) go to 10
  if ( u <= x(i+1) ) go to 30
!
!  binary search
!
 10 i = 1
  j = n+1
 20 k = (i+j)/2
  if ( u < x(k) ) j = k
  if ( u >= x(k) ) i = k
  if ( j > i+1 ) go to 20
!
!  evaluate spline
!
 30 dx = u - x(i)
  seval = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
  return

end function seval
