program fecp

      integer :: nros,nsm1, io,i,j
      real(kind=8),allocatable :: ens(:), dx(:,:), dy(:,:)
      real(kind=8),allocatable :: dz(:,:), ion(:), cpl(:,:)
      real(kind=8),allocatable :: kin(:)
      real(kind=8) :: trel
      real(kind=8),parameter :: limit=1.0d-30, zero=0.0d0, at2fs=2.4188843d-2
      character(len=400) :: fin, fout

      call getarg(1,fin)
      call getarg(2,fout)

      open(31,file=fin,access='stream',status='old')
      open(32,file=fout,status='replace',iostat=io)

      read(31) nros
  !    write(32,*) 'Number of states: ', nros
      nsm1 = nros - 1

      allocate(ens(0:nsm1), dx(0:nsm1,0:nsm1), dy(0:nsm1,0:nsm1),&
          dz(0:nsm1,0:nsm1),ion(0:nsm1),cpl(0:nsm1,0:nsm1),&
          kin(0:nsm1))

      read(31) ens
      write(32,*) 'Energies: '
 !     do i=0, nsm1
 !       write(32,*) ens(i)
 !     enddo
 !     call flush(32)
      read(31) ((dx(j,i),i=0,nsm1),j=0,nsm1)
      read(31) ((dy(j,i),i=0,nsm1),j=0,nsm1)
      read(31) ((dz(j,i),i=0,nsm1),j=0,nsm1)
!      write(32,*) 'Dipole (x,y,z): '
!      do i=0, nsm1
!        do j=i, nsm1
!          if(abs(dx(i,j)) >= limit .OR. abs(dy(i,j)) >= limit&
!             .OR. abs(dz(i,j)) >= limit) then
!                write(32,*) i, '->', j, dx(i,j), dy(i,j), dz(i,j)
!                call flush(32)
!            endif
!        enddo
!      enddo
      read(31) ion
!      write(32,*) 'Ionization rates: '
!      do i=0, nsm1
!        write(32,*) ion(i)
!        call flush(32)
!      enddo
      read(31) ((cpl(j,i),i=0,nsm1),j=0,nsm1)
      write(32,*) 'Non-Dipole Couplings: '
      do i=0, nsm1
       if(ion(i) .eq. zero) then
        do j=(i+1), nsm1
          if(ion(j) .eq. zero) then
              trel = at2fs/cpl(i,j)
              write(32,*) i, j, cpl(i,j), trel, dx(i,j), dy(i,j), dz(i,j)
              call flush(32)
          endif
        enddo
       endif
      enddo
      read(31) kin
!      write(32,*) 'Kinetic energies: '
!      do i=0, nsm1
!        write(32,*) kin(i)
!        call flush(32)
!      enddo
      close(31)
      close(32)

      end program fecp

