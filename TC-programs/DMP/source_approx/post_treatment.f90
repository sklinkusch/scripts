program post_treatment
! finds max dipole and max field and spits out the dynamic 
! polarizability with some error measure
  implicit none

  character(len=50) :: filename
  integer :: io,Nin,indF,indD,indM,indI,indx,indy,jdum
  real(kind=8) :: t2dum,tdum,fdum,ddum,cddum,iddum,u1(2),u2(2)
  real(kind=8) :: dphs,t_m,omega,period,dphs2
  real(kind=4) :: a1,a2,a_I
  real(kind=8),allocatable :: tin(:),ttau(:),fieldIN(:),dipoleIN(:),IdipoleIN(:)

  call getarg(1,filename)
  open(11,file=filename,status='old')
! count number of elements and allocate memory
  Nin = 0
  do
    read(11,*,iostat=io) jdum,tdum,fdum,ddum,t2dum,iddum
    if(io /= 0) exit
    Nin = Nin + 1
  enddo
  allocate(tIN(Nin),ttau(Nin),fieldIN(Nin),dipoleIN(Nin),IdipoleIN(Nin))
  rewind(11)

! store the field and vectors
  do io = 1,Nin
    read(11,*) jdum,tdum,fdum,ddum,t2dum,iddum
    tin(io) = tdum
    fieldIN(io) = fdum
    dipoleIN(io) = ddum
    ttau(io) = t2dum
    IdipoleIN(io) = iddum
  enddo
  


! find the max field, the max dipole, and the maximum of the ellipse
  indI = index(filename,".dip",.false.) - 1
  read(filename(1:indI),*) omega
  period = 1d0/omega
  indI = Nin/2
  indF = indI; indD = indI; indM = indI
  do io = max(indI-10,1),min(indI+10,Nin)
    if(abs(tin(io)-0.5d0*tin(Nin)) > period) cycle
    if(fieldIN(io) > fieldIN(indF)) indF = io
    if((dipoleIN(io)**2+fieldIN(io)**2) > (dipoleIN(indM)**2+fieldIN(indM)**2)) indM = io
    if(abs(dipoleIN(io)) > abs(dipoleIN(indD))) indD = io
  enddo

! polarizability at different critical points
  a1 = dipoleIN(indF)/fieldIN(indF)
  a2 = dipoleIN(indD)/fieldIN(indD)
  a_I = IdipoleIN(indF)/fieldIN(indF)
!  dphs = sign(acos(abs(dipoleIN(indF)/maxval(abs(dipoleIN)))),dble(a1))
  dphs = sign(acos(abs(dipoleIN(indF)/dipoleIN(indD))),dble(a1))
  dphs2 = (tin(indD)-0.5d0*tin(Nin))*omega
!  dphs = sign(acos(dipoleIN(indF)/dipoleIN(indD)),dble(a1))
!  x_dum = sign(acos(fieldIN(indD)/fieldIN(indF)),dble(a1))

! output
  write(*,'(a50,5e15.7)') filename,a1,a2,dphs,a_I,dphs2

  close(11)

end program post_treatment
