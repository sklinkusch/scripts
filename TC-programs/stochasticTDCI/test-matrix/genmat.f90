program genmat

    implicit none
    integer,parameter :: Nbasis=100,Nstates=100
    real(kind=8) :: energy(Nbasis),ionization(Nbasis),configurations(Nbasis),eigenvectors(Nbasis,Nbasis)
    real(kind=8) :: work(4*Nbasis),dipole(Nbasis,Nbasis,3),dephasing(Nbasis,Nbasis),rates(Nbasis,Nbasis)
    real(kind=8) :: dum
    integer :: i,j,lwork=4*Nbasis

! create fictive system
    configurations(1) = 0d0
    do i = 2,Nbasis
       configurations(i) = 0.1d0+(i-1)*0.001d0
    enddo

    ionization(:) = 0d0
    ionization(5:Nbasis) = 1d-5

    eigenvectors(:,:) = 0d0
    eigenvectors(1,1) = configurations(1)
    eigenvectors(2,2) = configurations(2)
    do i = 3,Nbasis
       eigenvectors(i,i-1) = 0.01d0
       eigenvectors(i  ,i) = configurations(i)
       eigenvectors(i-1,i) = 0.01d0
    enddo
    call dsyev('V','U',Nbasis,eigenvectors,Nbasis,energy,work,lwork,i)

! generate fictive dipole coupling
    dipole(:,:,:) = 0d0
    do i = 2,Nbasis
       dipole(i,i-1,1:3) = 0.1d0
       dipole(i-1,i,1:3) = 0.1d0
    enddo

! generate dubious rates
    dum = 1d-6*(0.01d0**2)
    do j = 1,Nbasis
    do i = 1,Nbasis
      if(i/=j) then
         rates(i,j) = dum/(configurations(i)-configurations(j))**2
      else
         rates(i,j) = 0d0
      endif
      dephasing(i,j) = 0d0
    enddo
    enddo

! output
    open(2,file='benchmark.sys',form='unformatted',access="stream")
    write(2) Nbasis,Nstates
    write(2) configurations
    write(2) ionization
    write(2) dipole(:,:,1)
    write(2) dipole(:,:,2)
    write(2) dipole(:,:,3)
    write(2) rates
    write(2) dephasing
    write(2) (energy(i),i=1,Nstates)
    do i=1,Nbasis
      write(2) eigenvectors(:,i)
    enddo
    close(2)

end program genmat
