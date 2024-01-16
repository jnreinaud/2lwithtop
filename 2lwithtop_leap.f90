! ============================================= !
! = Pseudo-spectral code for a QG 2layer flow = ! 
! =    -- with free-surface                   = !
! =    -- with topography                     = !
! =       JNR @ St Andrews 06/09/23           = !
! =           - St Andrews 18/10/23           = ! 
! ============================================= !

program two
  implicit none
  integer, parameter :: dp = kind(1.d0)
  real(dp), parameter :: pi = acos(-1.d0)
  integer :: nx, ny ! grid size in x and y
  integer :: nwx , nwy ! number of wave numbers (e.g. nwx = nx/2)
  real(dp), dimension(:,:), allocatable :: q1,q2,q2a,p1,p2,eta ! PV layer 1,2; Streamfunction layer 1,2; bathymetry
                                                               ! q2 includes topography, q2a is the anomaly 
  real(dp), dimension(:,:), allocatable :: q1i,q1f,q2i,q2f ! auxilary arrays for integration PV with RK4
  real(dp), dimension(:,:), allocatable :: u1,u2,v1,v2 ! horizontal velocity
  real(dp), dimension(:), allocatable :: rkx,rky,rkxf,rkyf ! wave numbers and filtered ones
  real(dp), dimension(:,:), allocatable :: rksi,rks ! Laplacian in spectral and square of wave numbers
  real(dp), dimension(:), allocatable :: trig,work ! used by FFT init 
  real(dp), dimension(:,:), allocatable :: emo,epo ! used for diffusion (integrating factor)
  real(dp), dimension(:,:), allocatable :: q1x,q2x,q1y,q2y ! spatial derivatives of PV 1,2 w.r.t. x,y
  real(dp) :: t,tc,dt,dt2,dt3,dt6,dfac ! time and time-steps
  real(dp) :: dtn ! temporary time step used to make sure data saved at the right time
  real(dp) :: tmax,dtsav ! maximum time and interval between 
  real(dp) :: nu ! viscosity
  real(dp) :: cfl,umax ! clf and used for clf
  real(dp) :: q1max,q2max,q2amax ! monitor max values of q1,q2,q2a
  real(dp) :: hx,hy,dg ! grid size (dg = min(hx,hy))
  integer :: ix,iy,kx,ky ! index in physical arrays and in spectral arrays
  real(dp) :: kd1,kd2,kd0,keta ! inverse deformation lengths (squared after reading) for interface in layer 1,2 (kd1,2)/free surface(kd0)/bathymetry(keta)
  ! kd1 = f_0^2/(g'H_1), g'=g(\rho_2-\rho_1)\rho_1
  ! kd2 = f_0^2/(g'H_2)
  ! kd0 = f_0^2/(gH_1)
  ! keta = f_0/H_2
  ! Vallis p213, eqn 5.85a,b
  real(dp) :: drc ! shift for merger
  
  logical :: idump ! logical for datasaves
  integer :: ite
  
  call readinit
  call allo
  call specinit
!  call dea
  call topo
  call PVinit
  
  t = 0.d0
  tc = 0.d0
  ite = 0
  call dump

  call stepping_euler_2
  ite = 1
  t = t+dt
  tc = tc+dt

  do while (t < tmax)
     write(*,*) 't = ',t,umax
!     call stepping
!     call stepping_euler
     !     call stepping_sem
     if ( mod(ite,50) == 0) then
        call steppint_euler_2
     else   
        call stepping_leap
     endif   
     t = t+dt
     tc = tc+dt
     if (tc > dtsav) then
        write(*,*) '   -- saved at t = ',t
        tc = tc - dtsav
        ite = ite+1
        call dump
     endif   
  enddo   
  
Contains
  
  ! ================================================= !
  ! = Read in parameter (grid size, cfl, viscosity) = !
  ! ================================================= !  
  subroutine readinit
    open(unit=10,file='param.dat',action='read',status='old')
    read(10,*) nx,ny
    read(10,*) nu,cfl
    read(10,*) kd0,kd1,kd2,keta
    read(10,*) tmax,dtsav
    read(10,*) drc
    kd0 = kd0**2  ! Square all inverse deformation lengths
    kd1 = kd1**2
    kd2 = kd2**2
    keta = keta**2
    hx = 2.d0*pi/dble(nx)
    hy = 2.d0*pi/dble(ny)
    dg = min(hx,hy)
    nwx = nx/2
    nwy = ny/2
    close(10)
  end subroutine readinit
  
  ! ======================================== !
  ! = Allocate memory for dynamical arrays = !
  ! ======================================== !
  subroutine allo

!    allocate(debug(nx,ny)) 
    
    allocate(q1(nx,ny))
    allocate(q2(nx,ny))
    allocate(q2a(nx,ny))
    allocate(p1(nx,ny))
    allocate(p2(nx,ny))
    allocate(eta(nx,ny))
    allocate(q1i(nx,ny))
    allocate(q1f(nx,ny))
    allocate(q2i(nx,ny))
    allocate(q2f(nx,ny))
    allocate(u1(nx,ny))
    allocate(u2(nx,ny))
    allocate(v1(nx,ny))
    allocate(v2(nx,ny))
    allocate(rksi(nx,ny))
    allocate(rks(nx,ny))
    allocate(emo(nx,ny))
    allocate(epo(nx,ny))
    allocate(q1x(nx,ny))
    allocate(q2x(nx,ny))
    allocate(q1y(nx,ny))
    allocate(q2y(nx,ny))
    allocate(rkx(nx))
    allocate(rkxf(nx))
    allocate(trig(2*max(nx,ny)))
    allocate(work(nx*ny))
    allocate(rky(ny))
    allocate(rkyf(ny))
  end subroutine allo

  ! ============================== !    
  ! = Define all spectral arrays = !
  ! ============================== !      
  subroutine specinit
    integer :: ifail
    real(dp) :: small = 1.d-13

!     Use array q1 to init FFTS
    do kx=1,nx
       do ky=1,ny
          q1(kx,ky) = 0.d0
       enddo
    enddo

    do kx=1,nwx
       rkx(kx) = dble(kx-1)
       rkx(nx+1-kx) = dble(kx)
    enddo
    do ky=1,nwy
       rky(ky) = dble(ky-1)
       rky(ny+1-ky) = dble(ky)
    enddo
    do kx=1,nx
       rkxf(kx) = rkx(kx)*exp(-36.d0*(rkx(kx)/dble(nwx))**36.d0)
    enddo   
    do ky=1,ny
       rkyf(ky) = rky(ky)*exp(-36.d0*(rky(ky)/dble(nwy))**36.d0)
    enddo

    do kx=1,nx
       do ky=1,ny
          rks(kx,ky) = rkx(kx)**2 + rky(ky)**2
          if (rks(kx,ky) .gt. small) then
             rksi(kx,ky) = 1.d0/rks(kx,ky)
          else
             rksi(kx,ky) = 0.d0
          endif
       enddo   
    enddo

    call c06fpf(nx,ny,q1,'i',trig,work,ifail)

  end subroutine specinit

  ! ==================================== !
  ! =  Dealiasing  using the 2/3 rule  = !
  ! ==================================== !

  subroutine dea
    integer :: kdx,kdy
    kdx = int(2.0*dble(nwx)/3.0)
    kdy = int(2.0*dble(nwy)/3.0)
    do kx=1,nx
       if (rkxf(kx) > kdx) rkxf(kx) = 0.d0
    enddo
    do ky=1,ny
       if (rkyf(ky) > kdy) rkyf(ky) = 0.d0
    enddo   
  end subroutine dea
    
  ! =========================================== !
  ! ==  FFT from physical to spectral space   = !
  ! =========================================== !
  subroutine ptosp(array)
    real(dp), dimension(:,:), allocatable :: array
    real(dp), dimension(:,:), allocatable :: trans
    integer :: ifail

    allocate(trans(ny,nx))

    call c06fpf(nx,ny,array,'s',trig,work,ifail)
    do ky=1,ny
       do kx=1,nx
          trans(ky,kx) = array(kx,ky)
       enddo
    enddo
    call c06fpf(ny,nx,trans,'s',trig,work,ifail)
    do ky=1,ny
       do kx=1,nx
          array(kx,ky) = trans(ky,kx)
       enddo
    enddo

    deallocate(trans)

  end subroutine ptosp
 
  ! =========================================== !
  ! ==  FFT from physical to spectral space   = !
  ! =========================================== !

  subroutine sptop(array)
      real(dp), dimension(:,:), allocatable :: array
      real(dp), dimension(:,:), allocatable :: trans
      integer :: ifail

      allocate(trans(ny,nx))

      call c06gqf(nx,ny,array,ifail)
      call c06fqf(nx,ny,array,'s',trig,work,ifail)
      do ix=1,nx
         do iy=1,ny
            trans(iy,ix) = array(ix,iy)
         enddo
      enddo
      call c06gqf(ny,nx,trans,ifail)
      call c06fqf(ny,nx,trans,'s',trig,work,ifail)
      do ix=1,nx
         do iy=1,ny
            array(ix,iy) = trans(iy,ix)
         enddo
      enddo
      deallocate(trans)

    end subroutine sptop

  ! ============================================ !
  ! = Defines topography in the physical space = !
  ! =    and moves it in spectral space        = !
  ! ============================================ !
  subroutine topo
    real(dp) :: x,y
    real(dp), parameter :: eta0 = 1.0d0, eps = 1.0d0
    do ix=1,nx
       x = -pi + dble(ix-1)*hx
       do iy=1,ny
          y = -pi + dble(iy-1)*hy
          eta(ix,iy) = eta0*exp(-(x/eps)**2-y**2)
          !eta(ix,iy) = 0.5d0*eta0*(sin(5*x)+sin(5*y))
       enddo
    enddo
    call ptosp(eta)
  end subroutine topo

  ! ====================================== !
  ! =  Defines initial PV anomaly field  = !
  ! =   and moves it in spectral space   = !
  ! ====================================== !
  subroutine PVinit
    real(dp) :: x,y,rloc
    real(dp), parameter :: q10= 1.d0,q20=-1.d0,r10=1.d0/3.d0,r20=r10
    real(dp), parameter :: fac = 8.d0
    drc = 0.5d0*drc*r10
    do ix=1,nx
       x = -pi + dble(ix-1)*hx
       do iy=1,ny
          y = -pi + dble(iy-1)*hy
! ---- CASE 1: Heton with mode 2 ---!          
!          rloc = (x/1.1d0)**2+y**2
!          q1(ix,iy)  = q10*exp(-fac*(rloc/r10**2)**fac)
!          q2a(ix,iy) = q20*exp(-fac*(rloc/r20**2)**fac)
! ---- CASE 2: Merger layer 1 --- !          
          q1(ix,iy)  = q10*exp(-fac*(((x-r10-drc)**2+y**2)/r10**2)**fac)
          q1(ix,iy)  = q1(ix,iy)+q10*exp(-fac*(((x+r10+drc)**2+y**2)/r10**2)**fac)
          q2a(ix,iy) = 0.d0          
! ---- CASE 3: Merger layer 2 --- !          
!          q2a(ix,iy)  = q20*exp(-fac*(((x-r20)**2+y**2)/r20**2)**fac)
!          q2a(ix,iy)  = q2a(ix,iy)+q20*exp(-fac*(((x+r20)**2+y**2)/r20**2)**fac)
!          q1(ix,iy) = 0.d0
! ---- CASE 2: dipole layer 1 --- !          
!          q1(ix,iy)  = q10*exp(-fac*(((x-r10)**2+(y-pi+r10)**2)/r10**2)**fac)
!          q1(ix,iy)  = q1(ix,iy)-q10*exp(-fac*(((x+r10)**2+(y-pi+r10)**2)/r10**2)**fac)
!          q2a(ix,iy) = 0.d0          
          
       enddo
    enddo
    call ptosp(q1)
    call ptosp(q2a)
    do kx=1,nx
       do ky=1,ny
          q2(kx,ky) = q2a(kx,ky) + keta*eta(kx,ky) ! add topography to q2a
       enddo
    enddo   
    
  end subroutine PVinit

  ! ============================================= !
  ! = Calculate the PV anomaly for bottom layer = !
  ! ============================================= !
  subroutine anom
    do kx=1,nx
       do ky=1,ny
          q2a(kx,ky) = q2(kx,ky) - keta*eta(kx,ky)
       enddo
    enddo
  end subroutine anom  
  
  ! =============== !
  ! =  Inversion  = !                      
  ! =============== !
  subroutine invert
    real(dp), parameter :: small = 1.d-13
    real(dp) :: a11,a12,a21,a22,det,deti

    call anom ! calculate PV anomaly in layer 2 (remove topography)

    do kx=1,nx
       do ky=1,ny
          a11 = -rks(kx,ky) - kd1 + kd0
          a12 = kd1
          a21 = kd2
          a22 = -rks(kx,ky) - kd2
          det = a11*a22-a12*a21
          if (det .gt. small) then
             deti = 1.d0/det
          else
             deti = 0.d0
          endif
          p1(kx,ky) = ( a22*q1(kx,ky)-a12*q2a(kx,ky))*deti
          p2(kx,ky) = (-a21*q1(kx,ky)+a11*q2a(kx,ky))*deti
       enddo
    enddo   
  end subroutine invert  

  ! =========================== !
  ! =  U,V as well as qx,qy   = !
  ! ==  output nonlinear term = !
  ! == -udq/dx-vdq/dy         = !
  ! =========================== !
  subroutine uv(iopt)
    integer :: kxc,kyc,iopt
    real(dp) :: uloc1,uloc2,uloc

    call invert ! invertion step

    do kx=1,nx
       u1(kx,1) = 0.d0
       q1y(kx,1) = 0.d0
       u1(kx,nwy+1) = 0.d0
       q1y(kx,nwy+1) = 0.d0
       u2(kx,1) = 0.d0
       q2y(kx,1) = 0.d0
       u2(kx,nwy+1) = 0.d0
       q2y(kx,nwy+1) = 0.d0
       do ky=2,nwy
          kyc = ny+2-ky
          q1y(kx,ky)  = -rkyf(ky)*q1(kx,kyc)
          u1(kx,ky)   =  rkyf(ky)*p1(kx,kyc)
          q1y(kx,kyc) =  rkyf(ky)*q1(kx,ky)
          u1(kx,kyc)  = -rkyf(ky)*p1(kx,ky)
          q2y(kx,ky)  = -rkyf(ky)*q2(kx,kyc)
          u2(kx,ky)   =  rkyf(ky)*p2(kx,kyc)
          q2y(kx,kyc) =  rkyf(ky)*q2(kx,ky)
          u2(kx,kyc)  = -rkyf(ky)*p2(kx,ky)
       enddo
    enddo
    do ky=1,ny
       v1(1,ky) = 0.d0
       q1x(1,ky) = 0.d0
       v1(nwx+1,ky) = 0.d0
       q1x(nwx+1,ky) = 0.d0
       v2(1,ky) = 0.d0
       q2x(1,ky) = 0.d0
       v2(nwx+1,ky) = 0.d0
       q2x(nwx+1,ky) = 0.d0
       do kx=2,nwx
          kxc = nx+2-kx
          q1x(kx,ky)  = -rkxf(kx)*q1(kxc,ky)
          v1(kx,ky)   = -rkxf(kx)*p1(kxc,ky)
          q1x(kxc,ky) =  rkxf(kx)*q1(kx,ky)
          v1(kxc,ky)  =  rkxf(kx)*p1(kx,ky)
          q2x(kx,ky)  = -rkxf(kx)*q2(kxc,ky)
          v2(kx,ky)   = -rkxf(kx)*p2(kxc,ky)
          q2x(kxc,ky) =  rkxf(kx)*q2(kx,ky)
          v2(kxc,ky)  =  rkxf(kx)*p2(kx,ky)
       enddo
    enddo
    
    call sptop(q1x)
    call sptop(q1y)
    call sptop(u1)
    call sptop(v1)
    call sptop(q2x)
    call sptop(q2y)
    call sptop(u2)
    call sptop(v2)

!    do ix=1,nx
!       do iy=1,ny
!          debug(ix,iy) = u2(ix,iy)**2+v2(ix,iy)**2
!       enddo   
!    enddo
    
    do ix=1,nx
       do iy=1,ny
          p1(ix,iy) = -(q1x(ix,iy)*u1(ix,iy)+q1y(ix,iy)*v1(ix,iy))
          p2(ix,iy) = -(q2x(ix,iy)*u2(ix,iy)+q2y(ix,iy)*v2(ix,iy))
       enddo
    enddo

    call ptosp(p1)
    call ptosp(p2)

    if (iopt == 0) then
       umax = 0.d0
       do ix=1,nx
          do iy=1,ny
             uloc1 = sqrt(u1(ix,iy)**2+v1(ix,iy)**2)
             uloc2 = sqrt(u2(ix,iy)**2+v2(ix,iy)**2)
             uloc = max(uloc1,uloc2)
             umax = max(umax,uloc)
          enddo
       enddo
    endif
    
  end subroutine uv  

  ! ===================================== !
  ! =        Time stepping              = !
  ! ===================================== !

  subroutine stepping

    real(dp) :: tmp1,tmp2

    do kx=1,nx
       do ky=1,ny
          q1i(kx,ky) = q1(kx,ky)
          q2i(kx,ky) = q2(kx,ky)
       enddo   
    enddo

    call uv(0)
!    open(unit=99,file='un.r4',form='unformatted',access='direct',status='replace',recl=4*(nx*ny+1))
!    write(99,rec=1) sngl(t) , sngl(debug)
!    close(99)
    
    dt = cfl*dg/umax
    dt2 = 0.5*dt
    dt3 = dt/3.d0
    dt6 = dt/6.d0
    dfac = dt2*nu
    do kx=1,nx
       do ky=1,ny
          epo(kx,ky) = exp(dfac*rks(kx,ky))
          emo(kx,ky) = 1.d0/epo(kx,ky)
       enddo
    enddo

    do kx=1,nx
       do ky=1,ny
          q1(kx,ky) = emo(kx,ky)*(q1i(kx,ky) + dt2*p1(kx,ky))
          q2(kx,ky) = emo(kx,ky)*(q2i(kx,ky) + dt2*p2(kx,ky))
          q1f(kx,ky) = dt6*p1(kx,ky)
          q2f(kx,ky) = dt6*p2(kx,ky)
       enddo
    enddo
    call uv(1)
    do kx=1,nx
       do ky=1,ny
          tmp1 = epo(kx,ky)*p1(kx,ky)
          tmp2 = epo(kx,ky)*p2(kx,ky)
          q1(kx,ky) = emo(kx,ky)*(q1i(kx,ky) + dt2*tmp1)
          q2(kx,ky) = emo(kx,ky)*(q2i(kx,ky) + dt2*tmp2)
          q1f(kx,ky) = q1f(kx,ky) + dt3*tmp1
          q2f(kx,ky) = q2f(kx,ky) + dt3*tmp2
       enddo
    enddo
    call uv(1)
    do kx=1,nx
       do ky=1,ny
          tmp1 = epo(kx,ky)**2*p1(kx,ky)
          tmp2 = epo(kx,ky)**2*p2(kx,ky)
          q1(kx,ky) = emo(kx,ky)**2*(q1i(kx,ky) + dt2*tmp1)
          q2(kx,ky) = emo(kx,ky)**2*(q2i(kx,ky) + dt2*tmp2)
          q1f(kx,ky) = q1f(kx,ky) + dt3*tmp1
          q2f(kx,ky) = q2f(kx,ky) + dt3*tmp2
       enddo
    enddo
    call uv(1)
    do kx=1,nx
       do ky=1,ny
          tmp1 = epo(kx,ky)**2*p1(kx,ky)
          tmp2 = epo(kx,ky)**2*p2(kx,ky)
          q1(kx,ky) = emo(kx,ky)**2*(q1i(kx,ky) + dt2*tmp1)
          q2(kx,ky) = emo(kx,ky)**2*(q2i(kx,ky) + dt2*tmp2)
          q1f(kx,ky) = q1f(kx,ky) + dt6*tmp1
          q2f(kx,ky) = q2f(kx,ky) + dt6*tmp2
       enddo
    enddo
    do kx=1,nx
       do ky=1,ny
          q1(kx,ky) = emo(kx,ky)**2*(q1i(kx,ky) + q1f(kx,ky))
          q2(kx,ky) = emo(kx,ky)**2*(q2i(kx,ky) + q2f(kx,ky))
       enddo
    enddo
  end subroutine stepping

  ! ===================================== !
  ! =        Time stepping (Euler!)     = !
  ! ===================================== !

  subroutine stepping_euler

    call uv(0)
    dt = 0.5d0*cfl*dg/umax
    dfac = dt*nu
    do kx=1,nx
       do ky=1,ny
          epo(kx,ky) = exp(dfac*rks(kx,ky))
          emo(kx,ky) = 1.d0/epo(kx,ky)
       enddo
    enddo

    do kx=1,nx
       do ky=1,ny
          q1(kx,ky) = emo(kx,ky)*(q1(kx,ky) + dt*p1(kx,ky))
          q2(kx,ky) = emo(kx,ky)*(q2(kx,ky) + dt*p2(kx,ky))
       enddo
    enddo
    
  end subroutine stepping_euler

  ! ======================================= !
  ! =        Time stepping (Euler back!)  = !
  ! ======================================= !

  subroutine stepping_semi

    call uv(0)
    dt = 0.5d0*cfl*dg/umax
    dfac = dt*nu

    do kx=1,nx
       do ky=1,ny
          q1(kx,ky) = (q1(kx,ky)+dt*p1(kx,ky))/(1.d0+dfac*rks(kx,ky))
          q2(kx,ky) = (q2(kx,ky)+dt*p2(kx,ky))/(1.d0+dfac*rks(kx,ky)) 
       enddo
    enddo
    
  end subroutine stepping_semi

  ! ======================================= !
  ! =        Time stepping (Leapfrog)     = !
  ! ======================================= !

  subroutine stepping_leap

    real(dp) :: ddt
    
    call uv(0)
!    dt = 0.5d0*cfl*dg/umax
    ddt = 2.d0*dt
    dfac = ddt*nu
    
    do kx=1,nx
       do ky=1,ny
          q1f(kx,ky) = q1(kx,ky)
          q2f(kx,ky) = q2(kx,ky)
       enddo
    enddo
    
    do kx=1,nx
       do ky=1,ny
          q1(kx,ky) = (q1i(kx,ky)+ddt*p1(kx,ky))/(1.d0+dfac*rks(kx,ky))
          q2(kx,ky) = (q2i(kx,ky)+ddt*p2(kx,ky))/(1.d0+dfac*rks(kx,ky)) 
       enddo
    enddo

    do kx=1,nx
       do ky=1,ny
          q1i(kx,ky) = q1f(kx,ky)
          q2i(kx,ky) = q2f(kx,ky)
       enddo
    enddo    
    
  end subroutine stepping_leap
  
  ! ======================================= !
  ! =        Time stepping (Leapfrog)     = !
  ! ======================================= !

  subroutine stepping_euler_2

    real(dp) :: ddt
    
    call uv(0)
    dt = 0.5d0*cfl*dg/umax
    dfac = dt*nu
    do kx=1,nx
       do ky=1,ny
          q1i(kx,ky) = q1(kx,ky)   
          q2i(kx,ky) = q2(kx,ky)
       enddo
    enddo
    
    do kx=1,nx
       do ky=1,ny
          q1(kx,ky) = (q1i(kx,ky)+dt*p1(kx,ky))/(1.d0+dfac*rks(kx,ky))
          q2(kx,ky) = (q2i(kx,ky)+dt*p2(kx,ky))/(1.d0+dfac*rks(kx,ky))
       enddo
    enddo
    
  end subroutine stepping_euler_2
    
  ! =============================== !
  ! =  dump results periodically  = !  
  ! =============================== !
  
  subroutine dump
    character(len=40) :: outfile1
    character(len=40) :: outfile2
    integer :: nbytes
    
    real(dp), dimension(:,:), allocatable :: localarr
    outfile1 = 'q1_000.r4'
    outfile2 = 'restart.r8'

    nbytes=4*(nx*ny+1)

    allocate(localarr(nx,ny))
    
    write(outfile1(4:6),'(i3.3)') ite
!    write(outfile2(4:6),'(i3.3)') ite

    open(unit=11,file=outfile1,form='unformatted',access='direct',status='replace',recl=nbytes)
    open(unit=12,file=outfile2,form='unformatted',access='direct',status='replace',recl=2*nbytes)
    
    do kx=1,nx
       do ky=1,ny
          localarr(kx,ky) = q1(kx,ky)
       enddo   
    enddo
    call sptop(localarr)
    write(11,rec=1) sngl(t),sngl(localarr)
    write(12,rec=1) t,localarr

    do kx=1,nx
       do ky=1,ny
          localarr(kx,ky) = q2(kx,ky) - keta*eta(kx,ky)
       enddo   
    enddo
    call sptop(localarr)
    write(11,rec=2) sngl(t),sngl(localarr)
    write(12,rec=2) t,localarr
    close(12)
    close(11)
    
    deallocate(localarr)
    
  end subroutine dump

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  READ PV: used to restart !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  
  subroutine readPV
    integer :: nbytes
    nbytes=4*(nx*ny+1)
    open(unit=12,file='restart.dat',form='unformatted',access='direct',status='replace',recl=2*nbytes)
    read(12,rec=1) t,q1
    read(12,rec=2) t,q2
    close(12)
    call ptosp(q1)
    call ptosp(q2a)
    do kx=1,nx
       do ky=1,ny
          q2(kx,ky) = q2a(kx,ky) + keta*eta(kx,ky) ! add topography to q2a
       enddo   
    enddo
  end subroutine readPV  
  
  
end program two
