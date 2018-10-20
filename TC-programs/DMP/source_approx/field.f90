!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program genfield

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  character(len=5) :: output,lineshape
  integer :: Nsteps
  real(kind=8) :: dt,time,field
  real(kind=8) :: t1on,t1off,tfin,tpulse1,mu1,freq1,E01

! parameter of the Lorentzian
  real(kind=8) :: nu0,nu,numin,numax,weight,dnu,area,E
  real(kind=8) :: fwhm,const,const2,norm

! constants
  real(kind=8),parameter :: pi = 3.141592653589793d0
  real(kind=8),parameter :: fs2au = 41.341373337d0
  real(kind=8),parameter :: light = 2.99792458d10
  real(kind=8),parameter :: zero = 0d0
  real(kind=8),parameter :: conv= 2.756923555732303547d-4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read input

  open(11,file='pulse.param',status='old')
  read(11,*) output
  read(11,*) lineshape
  read(11,*) t1on,t1off,tfin
  read(11,*) Nsteps
  read(11,*) mu1
  read(11,*) freq1
  read(11,*) E01
  read(11,*) dnu
  close(11)

! set pulse parameters
  tpulse1 = t1off-t1on           ! pulse duration
  if(E01 < 0.) then
     E01 = 2d0*pi/(tpulse1*fs2au*abs(mu1))
     write(*,*) 'Field fluence ::',real(E01/conv)**2
  else
     write(*,*) 'Field fluence ::',real(E01)
     E01 = conv*sqrt(E01)
  endif
  write(*,*) 'Field amplitude ::',real(E01)
  dt = tfin/Nsteps
  nu0 = 2d0*pi*freq1*light*1d-15 ! conversion WN to frequency

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! field envelope for a sin2 pulse

  if(lineshape == "shape") then
     open(11,file='shape.dat',status='replace')
     if(output == "mctdh") then
        write(11,22)
     else
        write(11,*) Nsteps+1,dt*fs2au
     endif
     time = zero
     field = zero
! compute the intensity of a PI pulse
     do
        if(time >= t1on .and. time <= t1off) then
           field = E01*sin(pi*(time-t1on)/tpulse1)**2
        endif
        if(output == "mctdh") then
           write(11,'(f20.10,f20.12)') time*fs2au,field
        else
           write(11,'(e20.7)') field
        endif
        time = time + dt
        if(time > tfin) exit
     enddo
     close(11)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! field for a delta pulse

  elseif(lineshape == "delta") then
     open(11,file='delta.dat',status='replace')
     if(output == "mctdh") then
        write(11,22)
     else
        write(11,*) Nsteps+1,dt*fs2au
     endif
     time = zero
     field = zero
! compute the intensity of a PI pulse
     do
        if(time >= t1on .and. time <= t1off) then
           field = E01*sin(pi*(time-t1on)/tpulse1)**2*cos(nu0*time)
        endif
        if(output == "mctdh") then
           write(11,'(f20.10,f20.12)') time*fs2au,field
        else
           write(11,'(e20.7)') field
        endif
        time = time + dt
        if(time > tfin) exit
     enddo
     close(11)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! field for a chirped pulse

  elseif(lineshape == "chirp") then
     dnu = 2d0*pi*dnu*light*1d-15 ! conversion WN to frequency
     const = 9d0*pi/(mu1*E01*tpulse1)
     const2 = dnu/(1d0 + const)

     open(11,file='chirp.dat',status='replace')
     if(output == "mctdh") then
        write(11,22)
     else
        write(11,*) Nsteps+1,dt*fs2au
     endif
     time = zero
     field = zero
! compute the field
     do
        if(time <= t1on) then
           field = E01*sin(0.5d0*pi*time/t1on)**2*cos(nu0*time)
        elseif(time > t1on .and. time <= t1off) then
           nu = nu0 + const2*(((time-t1on)/tpulse1)**2+const*(time-t1on)/tpulse1)
           field = E01*cos(nu*time)
        else
           field = E01*cos(0.5d0*pi*(time-t1off)/(tfin-t1off))**2*cos(nu*time)
        endif
        if(output == "mctdh") then
           write(11,'(f20.10,f20.12)') time*fs2au,field
        else
           write(11,'(e20.7)') field
        endif
        time = time + dt
        if(time > tfin) exit
     enddo
     close(11)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! field for a lorentzian-shaped line

  elseif(lineshape == "lortz") then

     open(11,file='lorentz.dat',status='replace')
     fwhm = 2d0*pi*(50d0)*light*1d-15 ! conversion WN to frequency
     numin = nu0 - 5d0*fwhm
     numax = nu0 + 5d0*fwhm
     dnu = (numax-numin)/1d3
     const = (fwhm*0.5d0)/pi
     const2 = (fwhm*0.5d0)**2

! compute the normalizing constant
     nu = numin
     norm = zero
     norme: do
        if(nu > numax) exit norme
        weight = const/((nu-nu0)**2+const2)
        norm = norm + weight*dnu
        nu = nu + dnu
     enddo norme

! set pulse area equal to the area of an equivalent delta pulse
     area = zero
     do
        if(time >= t1on .and. time <= t1off) then
           const = E01*sin(pi*(time-t1on)/tpulse1)**2*cos(nu0*time)
           area = area + abs(const)**2
        endif
        time = time + dt
        if(time > tfin) exit
     enddo
     area = area*dt
  
! compute the intensity of the desired pulse
     E = zero
     time = zero
     field = zero
     do
        if(time >= t1on .and. time <= t1off) then
           nu = numin
           field = zero
           aire: do
              if(nu > numax) exit aire
              weight = const/((nu-nu0)**2+const2)
              field = field + (sin(pi*(time-t1on)/tpulse1)**2*cos(nu*time)*weight/norm)*dnu
              nu = nu + dnu
           enddo aire
           E = E + abs(field)**2
        endif
        time = time + dt
        if(time > tfin) exit
     enddo
     E = sqrt(area/(E*dt))

! compute the field
     if(output == "mctdh") then
        write(11,22)
     else
        write(11,*) Nsteps+1,dt*fs2au
     endif
     time = zero
     do
        field = zero
        if(time >= t1on .and. time <= t1off) then
           nu = numin
           pulse: do
              if(nu > numax) exit pulse
              weight = const/((nu-nu0)**2+const2)
              field = field + E*sin(pi*(time-t1on)/tpulse1)**2*cos(nu*time)*weight*dnu/norm
              nu = nu + dnu
           enddo pulse
        endif
        if(output == "mctdh") then
           write(11,'(f20.10,f20.12)') time*fs2au,field
        else
           write(11,'(e20.7)') field
        endif
        time = time + dt
        if(time > tfin) exit
     enddo
     close(11)
  endif
 22 format('# Time [atomic units]   E(t) [e a_0]')

end program genfield
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
