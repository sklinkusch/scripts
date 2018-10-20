program fbcs

      integer :: nros,io,i,j
      real(kind=8),allocatable :: ens(:), vecs(:,:), dxvals(:)
      real(kind=8),allocatable :: dxvecs(:,:), dyvals(:)
      real(kind=8),allocatable :: dyvecs(:,:), dzvals(:)
      real(kind=8),allocatable :: dzvecs(:,:), corr(:)
      character(len=400) :: fin, fout

      nros = 3461

      allocate(ens(1:nros), vecs(1:nros,1:nros), dxvals(1:nros),&
          dxvecs(1:nros,1:nros),dyvals(1:nros),dyvecs(1:nros,1:nros),&
          dzvals(1:nros), dzvecs(1:nros,1:nros), corr(1:nros))

      call getarg(1,fin)
      call getarg(2,fout)

      open(31,file=fin,access='stream',status='old')
      open(32,file=fout,status='replace',iostat=io)
      read(31,pos=17) ens
      write(32,*) 'Energies: '
      do i=1, nros
        write(32,*) ens(i)
      enddo
      call flush(32)
      read(31) vecs
      write(32,*) 'Vectors: '
      do i=1, nros
       do j=1, nros
        write(32,*) vecs(j,i)
       enddo
      enddo
      call flush(32)
      read(31) dxvals
      write(32,*) 'x-Dipoles: '
      do i=1, nros
        write(32,*) dxvals(i)
      enddo
      call flush(32)
      read(31) dxvecs
      write(32,*) 'x-Dipole Vectors: '
      do i=1, nros
       do j=1, nros
        write(32,*) dxvecs(j,i)
       enddo
      enddo
      call flush(32)
      read(31) dyvals
      write(32,*) 'y-Dipoles: '
      do i=1, nros
        write(32,*) dyvals(i)
      enddo
      call flush(32)
      read(31) dyvecs
      write(32,*) 'y-Dipole Vectors: '
      do i=1, nros
       do j=1, nros
        write(32,*) dyvecs(j,i)
       enddo
      enddo
      call flush(32)
      read(31) dzvals
      write(32,*) 'z-Dipoles: '
      do i=1, nros
        write(32,*) dzvals(i)
      enddo
      call flush(32)
      read(31) dzvecs
      write(32,*) 'z-Dipole Vectors: '
      do i=1, nros
       do j=1, nros
        write(32,*) dzvecs(j,i)
       enddo
      enddo
      call flush(32)
      read(31) corr
      close(31)
      write(32,*) 'Energy Corrections: '
      do i=1, nros
        write(32,*) corr(i)
      enddo
      call flush(32)
      close(32)

      end program fbcs

