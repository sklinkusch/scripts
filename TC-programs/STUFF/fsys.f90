
      program fsys

          implicit none
          
          integer :: i,j,Nbasis,Nstates,Nbm1,Nsm1,nonzero, nnzm1
          integer(kind=8), allocatable :: row_ptr(:), col_ind(:)
          real(kind=8),allocatable :: configurations(:),ionization(:)
          real(kind=8),allocatable :: dx(:),dy(:),dz(:),rates(:,:)
          real(kind=8),allocatable :: dephasing(:,:),eigenvectors(:,:)
          real(kind=8),allocatable :: energy(:)
          character(len=50) :: Fenergies

          call getarg(1,Fenergies)

          open(2,File=Fenergies,status='old',access='stream',&
            form='unformatted')
          read(2) Nbasis
          Nbm1 = Nbasis - 1
          write(*,*) 'Nbasis: ', Nbasis
          read(2) Nstates
          Nsm1 = Nstates - 1
          write(*,*) 'Nstates: ', Nstates
          allocate(configurations(0:Nbm1), ionization(0:Nbm1),&
              rates(0:Nbm1,0:Nbm1),dephasing(0:Nbm1,0:Nbm1),&
              energy(0:Nsm1),eigenvectors(Nbasis,Nstates),&
              row_ptr(0:Nbasis))
          write(*,*) 'Configuration energies:'
          do i=0, Nbm1
            read(2) configurations(i)
            write(*,*) i, configurations(i)
          enddo
          write(*,*) 'Ionization rates: '
          do i=0, Nbm1
            read(2) ionization(i)
            write(*,*) i, ionization(i)
          enddo
          read(2) row_ptr
          write(*,*) 'Row pointer: '
          do i=0, Nbasis
           write(*,*) i, row_ptr(i)
          enddo
          nonzero = row_ptr(Nbasis)
          nnzm1 = nonzero - 1
          allocate(col_ind(0:nnzm1), dx(0:nnzm1), dy(0:nnzm1),&
              dz(0:nnzm1))
          read(2) col_ind
          write(*,*) 'Column indices: '
          do i=0, nnzm1
            write(*,*) i, col_ind(i)
          enddo
          read(2) dx
          read(2) dy
          read(2) dz
          write(*,*) 'Dipole(x, y, z): '
          do i=0, nnzm1
            do j=0, Nbasis
              if(i .LT. row_ptr(j)) then
                  write(*,*) j-1, col_ind(i), dx(i), dy(i), dz(i)
                  exit
              endif
            enddo
          enddo
          read(2) rates
          write(*,*) 'Relaxation rates: '
          do i=0, Nbm1
            do j=0, Nbm1
              write(*,*) i, j, rates(i,j)
            enddo
          enddo
          read(2) dephasing
          write(*,*) 'Dephasing: '
          do i=0, Nbm1
            do j=0, Nbm1
              write(*,*) i, j, dephasing(i,j)
            enddo
          enddo
          read(2) energy
          write(*,*) 'Energies: '
          do i=0, Nsm1
            write(*,*) i, energy(i)
          enddo
          do j=1,Nstates
            read(2) eigenvectors(:,j)
          enddo
          write(*,*) 'Eigenvectors: '
          do i=1, Nbasis
            do j=1, Nstates
              write(*,*) i, j, eigenvectors(i,j)
            enddo
          enddo
          close(2)
          end program fsys
