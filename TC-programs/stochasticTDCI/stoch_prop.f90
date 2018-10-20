!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
program stoch_prop

  use matrix
  use field
  implicit none
  logical :: restart
  integer :: i,io,N_traj,Nsteps,Nskip,Nwrite
  character(len=3) :: polarization
  character(len=30) :: Fin,Fenergies,Ffield,method
  character(len=30) :: fout,foutpop,foutavg,foutlog,foutrst,foutsts
  real(kind=8) :: temperature
  real(kind=8) :: eps
  real(kind=8),allocatable :: err(:)
  complex(kind=8),allocatable :: psi(:,:)            ! the state vectors

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! get input file name

  call getarg(1,Fin)

! read input parameters
  open(11,file=Fin,status='old',iostat=io)
  if(io /= 0) stop 'Problem opening input file'
  read(11,*)
  read(11,*) restart
  read(11,*) Fenergies
  read(11,*) Nbasis, Nstates
  Nbm1 = Nbasis - 1 ! auxiliary basis length
  read(11,*) N_traj, Nwrite
  read(11,*) Temperature
  read(11,*)
  read(11,*)
  read(11,*) Tfinal
  read(11,*) Nsteps
  read(11,*) method
  read(11,*) eps
  read(11,*)
  read(11,*)
  read(11,*) Ffield
  read(11,*) polarization
  read(11,*)
  read(11,*)
  read(11,*) fout
  read(11,*) Nskip
  close(11)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!! initializing matrices for the propagation
! name of the output files
  foutlog = trim(fout)//'.log'
  foutpop = trim(fout)//'.pop'
  foutavg = trim(fout)//'.avg'
  foutrst = trim(fout)//'.rst'
  foutsts = trim(fout)//'.sts'

! output to logfile for consistency check
  open(1,file=foutlog,status='replace',iostat=io)
  if(io /= 0) stop 'problem opening the log file'
  write(1,*) '----------------------------------------------------------------------'
  write(1,*) 'Number of configurations included in the basis: ',Nbasis
  write(1,*) 'Number of eigenstates included in the basis: ',Nstates
  write(1,*) 'Number of quantum trajectories performed: ',N_traj
  write(1,*) 'Temperature of the system: ',real(Temperature)

! initialize the matrix: configurations, eigenstates, transition dipoles, rates
  write(1,*) 'Information in the isolated subsystem read from file: ',Fenergies
  open(2,file=Fenergies,status='old',form='unformatted',access='stream')
  if(io /= 0) stop 'problem opening the system file'

  allocate(psi(0:Nbasis,N_traj),err(N_traj))
  call initialize(restart,polarization,N_traj,psi,temperature)
  close(2)
  close(3)

! read and interpolate electric field
  write(1,*) 'Field read from file: ',Ffield
  write(1,*) 'Polarization of the electric field :: ',polarization
  open(4,file=Ffield,status='old',iostat=io)
  if(io /= 0) stop 'problem opening the field file'
  call setfield(polarization,Npol)
  close(4)

! summarizing propagation parameters
  write(1,*) 'Propagation time: ',real(Tfinal/41.341373337d0),' fs'
  if(trim(method) == 'adaptive') then
    write(1,*) 'Adaptive time-step Runge-Kutta Cash-Karp scheme'
  elseif(trim(method) == 'rungekutta') then
    write(1,*) 'Simple Runge-Kutta Cash-Karp scheme'
  else
    stop 'No such method. Choose between rungekutta or adaptive.'
  endif
  write(1,*) 'Propagating in the interaction representation'
  write(1,*) 'Individual trajectory population output in file: ',foutpop
  write(1,*) 'Incoherently averaged population output in file: ',foutavg
  write(1,*) '----------------------------------------------------------------------'
  write(1,*) ''
  call flush(1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! main propagation loop
  open(9,file=foutpop,status='replace')
  open(10,file=foutavg,status='replace')
  open(12,file=foutsts,status='replace')

  call stepper(psi,N_traj,method,Nsteps,Nskip,Nwrite,eps,matvecI)
! end of propagation
  close(1)   ! close log file
  close(9)   ! close population file
  close(10)  ! close average file
  close(12)  ! close statepop file

! write restart file
  open(8,file=foutrst,status='replace',form='unformatted',access='stream')
  write(8) Nbasis,N_traj
  do i = 1,N_traj
    write(8) psi(:,i)
  enddo
  write(8) Tfinal
  close(8)

end program stoch_prop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
