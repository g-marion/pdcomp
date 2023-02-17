
  subroutine pdcomp(ni,nj,nk,wbc,ebc,sbc,nbc,dx,dy,          &
                    zh,rho0,th0,qv0,pi0,u0,v0,u,v,w,pavgin,  &
                    pb,pdn,pdl,ptdn,fpb,fptdn,fpdn,          &
                    pex,psh,psum,fpex,fpsh,fpsum,fpdl,fpdlx, &
                    fpdly,fpdlz,fpdnx,fpdny,fpdnz,fpbx,fpby, &
                    fpbz,fptdnx,fptdny,fptdnz)

  use singleton
  implicit none

  integer, intent(in) :: ni,nj,nk
  integer, intent(in) :: wbc,ebc,sbc,nbc
  real, intent(in) :: dx,dy
  real, intent(in), dimension(nk) :: zh,rho0,th0,qv0,pi0,u0,v0
  real, intent(in), dimension(-2:ni+4,-2:nj+3,nk) :: u
  real, intent(in), dimension(-2:ni+3,-2:nj+4,nk) :: v
  real, intent(in), dimension(-2:ni+3,-2:nj+3,nk+1) :: w
! JT add
  real*8, intent(in) :: pavgin
  real, intent(inout), dimension(ni,nj,nk) :: pb,pdn,pdl
! JT:  ptdn is total dynamic pressure
  real, intent(inout), dimension(ni,nj,nk) :: ptdn,fpb,fptdn,fpdn
  real, intent(inout), dimension(ni,nj,nk) :: fpdlx,fpdly,fpdlz
  real, intent(inout), dimension(ni,nj,nk) :: fpdnx,fpdny,fpdnz
  real, intent(inout), dimension(ni,nj,nk) :: fpbx,fpby,fpbz
  real, intent(inout), dimension(ni,nj,nk) :: fptdnx,fptdny,fptdnz
  real, intent(inout), dimension(ni,nj,nk) :: fpdl  
  real, intent(inout), dimension(ni,nj,nk) :: pex,psh,fpex,fpsh
  real, intent(inout), dimension(ni,nj,nk) :: psum,fpsum
! JT:  pex is pressure due to fluid extension
!      psh is pressure due to fluid shear
!      psum is the sum of these two effects
!      fpex, and fpsh are the corresponding forcing contributions
!      fpsum is the forcing contributions from the sum of shear & exten

!-----------------------------------------------------------------------
!
!  pdcomp - a fortran90 subroutine to retrieve pressure perturbations
!           from cloud-scale numerical model output.  Three terms 
!           are retrieved:
!              pb  = buoyancy pressure
!              pdn = nonlinear dynamic pressure
!              pdl = linear dynamic pressure
!              ptdn = total dynamic pressure  (JT)
!
!  Version 1.03                           Last modified:  20 February 2013
!
!  Author:  George H. Bryan
!           Mesoscale and Microscale Meteorology Division
!           National Center for Atmospheric Research
!           Boulder, Colorado, USA
!           gbryan@ucar.edu
!
!  Disclaimer:  This code is made available WITHOUT WARRANTY.
!
!  References:  Rotunno and Klemp, 1982, MWR, p. 150
!               Weisman and Rotunno, 2000, JAS, p. 1461
!
!-----------------------------------------------------------------------
!
! Input:
!   Integer variables:
!     ni  = number of grid points in x
!     nj  = number of grid points in y
!     nk  = number of grid points in z
!
!     wbc = west boundary condition (see below)
!     ebc = east boundary condition (see below)
!     sbc = south boundary condition (see below)
!     nbc = north boundary condition (see below)
!
!  for boundary conditions:   1 = periodic ;   2 = open ;   3 = rigid wall
!       
!
!   Real variables:
!     dx  = grid spacing in x (m)
!     dy  = grid spacing in y (m)  (must be same as dx, for now!)
!
!   Real one-dimensional arrays:
!     zh  (nk)      = height of model's half levels (scalar points) (m)
!     rho0(nk)      = base-state density (kg/m^3)
!     th0 (nk)      = base-state potential temperature (K)
!     qv0 (nk)      = base-state water vapor mixing ratio (kg/kg)
!     pi0 (nk)      = base-state nondimensional pressure (dimensionless)
!      u0 (nk)      = base-state wind profile in x-direction (m/s)
!      v0 (nk)      = base-state wind profile in y-direction (m/s)
!
! Output:
!
!   Real three-dimensional arrays:
!     pb (ni,nj,nk) = buoyancy pressure
!     pdn(ni,nj,nk) = nonlinear dynamic pressure
!     pdl(ni,nj,nk) = linear dynamic pressure
!
!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!            No need to modify anything below here:
!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! JT notes:
! 1.  uses 6th order approximations to 1st derivative
!-----------------------------------------------------------------------

  integer :: i,j,k,nloop,ipb,ipe,jpb,jpe,kpb,kpe,imirror,jmirror
  real :: rdx,rdy
  real*8 :: dpi,pavg,frac
  real, dimension(0:nk+1) :: thv0
  real, dimension(:,:,:), allocatable :: dum1,dum2,dum3,divx,uten,vten,wten

  real, dimension(0:nk+1) :: r1,rf0,rr0,mh,mf,zf

  complex, dimension(:,:), allocatable :: rhs,trans
  complex, dimension(:,:,:), allocatable :: deft
  complex, dimension(0:nk+1) :: cfb,lgbth,lgbph
  real, dimension(0:nk+1) :: cfa,cfc
! JT:  for buoyancy
  real, dimension(ni,nj,nk) :: bb 

!-----------------------------------------------------------------------

  real, parameter :: g     = 9.81
  real, parameter :: rd    = 287.04
  real, parameter :: rv    = 461.5
  real, parameter :: cp    = 1005.7
  real, parameter :: pi    = 3.14159265358979323
  real, parameter :: p00   = 100000.0

  real, parameter :: eps   = rd/rv
  real, parameter :: reps  = 1.0/eps
  real, parameter :: rp00  = 1.0/p00
  real, parameter :: rddcp = rd/cp

!-----------------------------------------------------------------------

    print *,'ni, nj, nz', ni, nj, nk
! JT:  get buoyancy field
    do k=1,nk
    do j=1,nj
    do i=1,ni
      bb(i,j,k) = pb(i,j,k)
    enddo
    enddo
    enddo

  if(wbc.eq.1.and.ebc.ne.1)then
    print *
    print *,'  Inconsistent boundary conditions'
    print *
    stop 1111
  endif
  if(ebc.eq.1.and.wbc.ne.1)then
    print *
    print *,'  Inconsistent boundary conditions'
    print *
    stop 1112
  endif
  if(sbc.eq.1.and.nbc.ne.1)then
    print *
    print *,'  Inconsistent boundary conditions'
    print *
    stop 1113
  endif
  if(nbc.eq.1.and.sbc.ne.1)then
    print *
    print *,'  Inconsistent boundary conditions'
    print *
    stop 1114
  endif
  if( (ebc.lt.1.or.ebc.gt.3) .or.   &
      (wbc.lt.1.or.wbc.gt.3) .or.   &
      (sbc.lt.1.or.sbc.gt.3) .or.   &
      (nbc.lt.1.or.nbc.gt.3) )then
    print *
    print *,'  Invalid setting for boundary conditions'
    print *
    stop 1115
  endif

!-----------------------------------------------------------------------

  ipb=1
  ipe=ni

  jpb=1
  jpe=nj

  imirror = 0
  jmirror = 0

  if( (wbc.eq.2.or.wbc.eq.3).or.(ebc.eq.2.or.ebc.eq.3) )then

    imirror = 1
    ipe = ni*2

  endif

  if( (sbc.eq.2.or.sbc.eq.3).or.(nbc.eq.2.or.nbc.eq.3) )then

    jmirror = 1
    jpe = nj*2

  endif

  kpb=0
  kpe=nk+1

  print *,'  ipb,ipe         = ',ipb,ipe
  print *,'  jpb,jpe         = ',jpb,jpe
  print *,'  kpb,kpe         = ',kpb,kpe
  print *,'  imirror,jmirror = ',imirror,jmirror
  allocate(  deft(ipb:ipe,jpb:jpe,kpb:kpe) )

!----- constants -----

  rdx = 1.0/dx
  rdy = 1.0/dy

  dpi = 4.0d0*datan(1.0d0)
  print *,'  dpi = ',dpi

  do k=1,nk
    rr0(k)=1.0/rho0(k)
  enddo

  do k=2,nk
    rf0(k)=0.5*(rho0(k-1)+rho0(k))
  enddo

  zf(1)=0.0
  do k=2,nk+1
    zf(k)=2.0*zh(k-1)-zf(k-1)
  enddo

  do k=1,nk
    mh(k)=1.0/(zf(k+1)-zf(k))
  enddo

  do k=2,nk
    mf(k)=1.0/(zh(k)-zh(k-1))
  enddo
  mf(1)=mf(2)
  mf(nk+1)=mf(nk)

  do k=1,nk
    thv0(k) = th0(k)*(1.0+reps*qv0(k))/(1.0+qv0(k))
!    print *, 'thv0,k=',thv0(k),k
  enddo
  thv0(0) = thv0(1)
  thv0(nk+1) = thv0(nk)

!!!  do k=1,nk+1
!!!    print *,'zf,mf,mh:',zf(k),mf(k),mh(k)
!!!  enddo

!----- coefficients -----

  cfa = 0.0
  cfc = 0.0

  do k=1,nk
    cfa(k)=mh(k)*mf(k  )*rf0(k  )*0.5*(thv0(k-1)+thv0(k))   &
          /(rho0(k)*thv0(k))
    cfc(k)=mh(k)*mf(k+1)*rf0(k+1)*0.5*(thv0(k)+thv0(k+1))   &
          /(rho0(k)*thv0(k))
  enddo
  cfa( 1) = 0.0
  cfc(nk) = 0.0

  call bcu(ni,nj,nk,wbc,ebc,sbc,nbc,u)
  call bcv(ni,nj,nk,wbc,ebc,sbc,nbc,v)
  call bcw(ni,nj,nk,wbc,ebc,sbc,nbc,w)

!----- retrieve pressure -----

!    print *,'  alloc ... '
    allocate( uten(0:ni+1,0:nj+1,0:nk+1) )
    allocate( vten(0:ni+1,0:nj+1,0:nk+1) )
    allocate( wten(0:ni+1,0:nj+1,0:nk+1) )
!    print *,'  ... ok '

  DO nloop=1,3
!  DO nloop=1,5

    print *,'  nloop = ',nloop

    uten = 0.0
    vten = 0.0
    wten = 0.0

!-----------------------------------------------------------------------
!  forcing for buoyancy pressure

  IF(nloop.eq.1)THEN

    do k=2,nk
    do j=1,nj
    do i=1,ni
      wten(i,j,k) = 0.5*(pb(i,j,k-1)+pb(i,j,k))
! trying this...doesn't seem to matter
!      wten(i,j,k) = 0.5*(bb(i,j,k-1)+bb(i,j,k))
    enddo
    enddo
    enddo

    do k=1,nk
    do j=1,nj
    do i=1,ni
      deft(i,j,k) = (rf0(k+1)*wten(i,j,k+1)-rf0(k)*wten(i,j,k))   &
                  /(zf(k+1)-zf(k))                               &
                  /(cp*rho0(k)*thv0(k))
    enddo
    enddo
    enddo

!-----------------------------------------------------------------------
!  forcing for linear dynamic pressure

  ELSEIF(nloop.eq.2)THEN

!    print *,'  alloc ... '
    allocate( dum1(0:ni+1,0:nj+1,0:nk+1) )
    dum1 = 0.0
    allocate( dum2(0:ni+1,0:nj+1,0:nk+1) )
    dum2 = 0.0
    allocate( dum3(0:ni+1,0:nj+1,0:nk+1) )
    dum3 = 0.0
!    print *,'  ... ok '

    do k=2,nk
    do j=1,nj
    do i=1,ni+1
      dum1(i,j,k) = (u0(k)-u0(k-1))/(zh(k)-zh(k-1))   &
                   *(w(i,j,k)-w(i-1,j,k))*rdx
    enddo
    enddo
    enddo

    do j=1,nj
    do i=1,ni+1
      dum1(i,j,1) = 0.0
      dum1(i,j,nk+1) = 0.0
    enddo
    enddo

    do k=2,nk
    do j=1,nj+1
    do i=1,ni
      dum2(i,j,k) = (v0(k)-v0(k-1))/(zh(k)-zh(k-1))   &
                   *(w(i,j,k)-w(i,j-1,k))*rdx
    enddo
    enddo
    enddo

    do j=1,nj+1
    do i=1,ni
      dum2(i,j,1) = 0.0
      dum2(i,j,nk+1) = 0.0
    enddo
    enddo

    do k=1,nk
    do j=1,nj
    do i=1,ni
      deft(i,j,k) = -2.0*rho0(k)*0.25*(                                        &
             ((dum1(i,j,k)+dum1(i+1,j,k))+(dum1(i,j,k+1)+dum1(i+1,j,k+1)))    &
            +((dum2(i,j,k)+dum2(i,j+1,k))+(dum2(i,j,k+1)+dum2(i,j+1,k+1))) )  &
                  /(cp*rho0(k)*thv0(k))
    enddo
    enddo
    enddo

!    print *,'  dealloc ... '
    deallocate( dum1 )
    deallocate( dum2 )
    deallocate( dum3 )
!    print *,' ... ok '

!-----------------------------------------------------------------------
!  forcing for total dynamic component
!    (nonlinear is obtained by subtracting linear dynamic pressure
!     from total dynamic pressure)

  ELSEIF(nloop.eq.3)THEN

!    print *,'  alloc ... '
    allocate( dum1(0:ni+1,0:nj+1,0:nk+1) )
    allocate( dum2(0:ni+1,0:nj+1,0:nk+1) )
    allocate( dum3(0:ni+1,0:nj+1,0:nk+1) )
    allocate( divx(0:ni+1,0:nj+1,0:nk+1) )
    divx = 0.0
!    print *,'  ... ok '

    do k=1,nk
    do j=0,nj+1
    do i=0,ni+1
      divx(i,j,k)=rho0(k)*(u(i+1,j,k)-u(i,j,k))*rdx   &
                 +rho0(k)*(v(i,j+1,k)-v(i,j,k))*rdy   &
                 +(rf0(k+1)*w(i,j,k+1)-rf0(k)*w(i,j,k))/(zf(k+1)-zf(k))
    enddo
    enddo
    enddo

    print*,'Got divx'

!--------------------
!  u-tendency

    dum1 = 0.0
    dum2 = 0.0
    dum3 = 0.0

    do k=1,nk
    do j=1,nj
    do i=0,ni+1
      dum1(i,j,k)=rho0(k)*(u(i,j,k)+u(i+1,j,k))       &
                  *( 37.0*(u(i+1,j,k)+u(i  ,j,k))     &
                     -8.0*(u(i+2,j,k)+u(i-1,j,k))     &
                         +(u(i+3,j,k)+u(i-2,j,k)) )/120.0
    enddo
    enddo
    enddo

    do k=1,nk
    do j=1,nj+1
    do i=1,ni+1
      dum2(i,j,k)=rho0(k)*(v(i,j,k)+v(i-1,j,k))       &
                  *( 37.0*(u(i,j  ,k)+u(i,j-1,k))     &
                     -8.0*(u(i,j+1,k)+u(i,j-2,k))     &
                         +(u(i,j+2,k)+u(i,j-3,k)) )/120.0
    enddo
    enddo
    enddo

    do j=1,nj
      do k=4,nk-2
      do i=1,ni+1
        dum3(i,j,k)=rf0(k)*(w(i,j,k)+w(i-1,j,k))   &
                   *( 37.0*(u(i,j,k  )+u(i,j,k-1)) &
                      -8.0*(u(i,j,k+1)+u(i,j,k-2)) &
                          +(u(i,j,k+2)+u(i,j,k-3)) )/120.0
      enddo
      enddo

      do k=3,(nk-1),(nk-4)
      do i=1,ni+1
        dum3(i,j,k)=rf0(k)*(w(i,j,k)+w(i-1,j,k))  &
                   *( 7.0*(u(i,j,k  )+u(i,j,k-1)) &
                         -(u(i,j,k+1)+u(i,j,k-2)) )/24.0
      enddo
      enddo

      do k=2,nk,(nk-2)
      do i=1,ni+1
        dum3(i,j,k)=0.25*rf0(k)*(w(i,j,k)+w(i-1,j,k))    &
                               *(u(i,j,k-1)+u(i,j,k))
      enddo
      enddo
    enddo

    do k=1,nk
    do j=1,nj
    do i=1,ni+1
      uten(i,j,k) = ( -(dum1(i,j,k)-dum1(i-1,j,k))*rdx              &
                      -(dum2(i,j+1,k)-dum2(i,j,k))*rdy              &
                      -(dum3(i,j,k+1)-dum3(i,j,k))/(zf(k+1)-zf(k))  &
                      +u(i,j,k)*0.5*(divx(i,j,k)+divx(i-1,j,k)) )*rr0(k)
    enddo
    enddo
    enddo

    print*,'uten calculated'

!--------------------
!  v-tendency

    dum1 = 0.0
    dum2 = 0.0
    dum3 = 0.0

    do k=1,nk
    do j=1,nj+1
    do i=1,ni+1
      dum1(i,j,k)=rho0(k)*(u(i,j,k)+u(i,j-1,k))       &
                  *( 37.0*(v(i  ,j,k)+v(i-1,j,k))     &
                     -8.0*(v(i+1,j,k)+v(i-2,j,k))     &
                         +(v(i+2,j,k)+v(i-3,j,k)) )/120.0
    enddo
    enddo
    enddo

    do k=1,nk
    do j=0,nj+1
    do i=1,ni
      dum2(i,j,k)=rho0(k)*(v(i,j,k)+v(i,j+1,k))       &
                  *( 37.0*(v(i,j+1,k)+v(i,j  ,k))     &
                     -8.0*(v(i,j+2,k)+v(i,j-1,k))     &
                         +(v(i,j+3,k)+v(i,j-2,k)) )/120.0
    enddo
    enddo
    enddo

    do j=1,nj+1
      do k=4,nk-2
      do i=1,ni
        dum3(i,j,k)=rf0(k)*(w(i,j,k)+w(i,j-1,k))   &
                   *( 37.0*(v(i,j,k  )+v(i,j,k-1)) &
                      -8.0*(v(i,j,k+1)+v(i,j,k-2)) &
                          +(v(i,j,k+2)+v(i,j,k-3)) )/120.0
      enddo
      enddo

      do k=3,(nk-1),(nk-4)
      do i=1,ni
        dum3(i,j,k)=rf0(k)*(w(i,j,k)+w(i,j-1,k))  &
                   *( 7.0*(v(i,j,k  )+v(i,j,k-1)) &
                         -(v(i,j,k+1)+v(i,j,k-2)) )/24.0
      enddo
      enddo

      do k=2,nk,(nk-2)
      do i=1,ni
        dum3(i,j,k)=0.25*rf0(k)*(w(i,j,k)+w(i,j-1,k))    &
                               *(v(i,j,k-1)+v(i,j,k))
      enddo
      enddo
    enddo

    do k=1,nk
    do j=1,nj+1
    do i=1,ni
      vten(i,j,k) = ( -(dum1(i+1,j,k)-dum1(i,j,k))*rdx              &
                      -(dum2(i,j,k)-dum2(i,j-1,k))*rdy              &
                      -(dum3(i,j,k+1)-dum3(i,j,k))/(zf(k+1)-zf(k))  &
                      +v(i,j,k)*0.5*(divx(i,j,k)+divx(i,j-1,k)) )*rr0(k)
    enddo
    enddo
    enddo

    print*,'vten calculated'

!--------------------
!  w-tendency

    dum1 = 0.0
    dum2 = 0.0
    dum3 = 0.0

    do k=2,nk
    do j=1,nj
    do i=1,ni+1
      dum1(i,j,k)=(rho0(k)*u(i,j,k)+rho0(k-1)*u(i,j,k-1))  &
                  *( 37.0*(w(i  ,j,k)+w(i-1,j,k))          &
                     -8.0*(w(i+1,j,k)+w(i-2,j,k))          &
                         +(w(i+2,j,k)+w(i-3,j,k)) )/120.0
    enddo
    enddo
    enddo

    do k=2,nk
    do j=1,nj+1
    do i=1,ni
      dum2(i,j,k)=(rho0(k)*v(i,j,k)+rho0(k-1)*v(i,j,k-1))  &
                  *( 37.0*(w(i,j  ,k)+w(i,j-1,k))          &
                     -8.0*(w(i,j+1,k)+w(i,j-2,k))          &
                         +(w(i,j+2,k)+w(i,j-3,k)) )/120.0
    enddo
    enddo
    enddo

    do j=1,nj
      do k=3,nk-2
      do i=1,ni
        dum3(i,j,k)=(rf0(k)*w(i,j,k)+rf0(k+1)*w(i,j,k+1))   &
                   *( 37.0*(w(i,j,k+1)+w(i,j,k  ))          &
                      -8.0*(w(i,j,k+2)+w(i,j,k-1))          &
                          +(w(i,j,k+3)+w(i,j,k-2)) )/120.0
      enddo
      enddo

      do k=2,(nk-1),(nk-3)
      do i=1,ni
        dum3(i,j,k)=(rf0(k)*w(i,j,k)+rf0(k+1)*w(i,j,k+1))   &
                   *( 7.0*(w(i,j,k+1)+w(i,j,k  ))           &
                         -(w(i,j,k+2)+w(i,j,k-1)) )/24.0
      enddo
      enddo

      do k=1,nk,(nk-1)
      do i=1,ni
        dum3(i,j,k)=0.25*(rf0(k)*w(i,j,k)+rf0(k+1)*w(i,j,k+1))  &
                        *(w(i,j,k)+w(i,j,k+1))
      enddo
      enddo
    enddo

    do k=2,nk
    do j=1,nj
    do i=1,ni
      wten(i,j,k) = ( -(dum1(i+1,j,k)-dum1(i,j,k))*rdx              &
                      -(dum2(i,j+1,k)-dum2(i,j,k))*rdy              &
                      -(dum3(i,j,k)-dum3(i,j,k-1))/(zh(k)-zh(k-1))  &
                      +w(i,j,k)*0.5*(divx(i,j,k)+divx(i,j,k-1)) )/rf0(k)
    enddo
    enddo
    enddo

    do j=1,nj
    do i=1,ni
      wten(i,j,1)=0.0
      wten(i,j,nk+1)=0.0
    enddo
    enddo

    print*,'wten calculated'

!--------------------
!  total dynamic forcing

    do k=1,nk
    do j=1,nj
    do i=1,ni
      deft(i,j,k) = (                                                         &
                rho0(k)*( (uten(i+1,j,k)-uten(i,j,k))*rdx                    &
                         +(vten(i,j+1,k)-vten(i,j,k))*rdy )                  &
               +(rf0(k+1)*wten(i,j,k+1)-rf0(k)*wten(i,j,k))/(zf(k+1)-zf(k))  &
                   )/(cp*rho0(k)*thv0(k))
    enddo
    enddo
    enddo

    print*,'total dynamic forcing calculated'

!    print *,'  dealloc ... '
    deallocate( dum1 )
    deallocate( dum2 )
    deallocate( dum3 )
    deallocate( divx )
!    print *,' ... ok '

!-----------------------------------------------------------------------
!  forcing for fluid shear contribution

  ELSEIF(nloop.eq.4)THEN

!    print *,'  alloc ... '
    allocate( dum1(0:ni+1,0:nj+1,0:nk+1) )
    allocate( dum2(0:ni+1,0:nj+1,0:nk+1) )
    allocate( dum3(0:ni+1,0:nj+1,0:nk+1) )
    allocate( divx(0:ni+1,0:nj+1,0:nk+1) )
    dum1=0.
    dum2=0.
    dum3=0.
    divx = 0.0
!    print *,'  ... ok '

!    do k=1,nk
    do k=2,nk-1
    do j=1,nj
    do i=1,ni
      dum1(i,j,k)=(0.5*((v(i+1,j,k)-v(i-1,j,k))        &
                      +(v(i+1,j+1,k)-v(i-1,j+1,k)))*rdx*0.5)   &
                 *(0.5*((u(i,j+1,k)-u(i,j-1,k))        &
                     +(u(i+1,j+1,k)-u(i+1,j-1,k)))*rdy*0.5)    
    enddo
    enddo
    enddo

    do k=2,nk-1
    do j=1,nj
    do i=1,ni
      dum2(i,j,k)=(0.5*((w(i+1,j,k)-w(i-1,j,k))               &
                 +(w(i+1,j,k+1)-w(i,j,k+1)))*rdx*0.5)    &
                 *(0.5*((u(i,j,k+1)-u(i,j,k-1))               &
                 +(u(i+1,j,k+1)-u(i+1,j,k-1)))/(zh(k+1)-zh(k-1)))
    enddo
    enddo
    enddo

    do k=2,nk-1
    do j=1,nj
    do i=1,ni
      dum3(i,j,k)=(0.5*((w(i,j+1,k)-w(i,j-1,k))               &
                 +(w(i,j+1,k+1)-w(i,j-1,k+1)))*rdy*0.5)  &
                 *(0.5*((v(i,j,k+1)-v(i,j,k-1))               &
                 +(v(i,j+1,k+1)-v(i,j+1,k-1)))/(zh(k+1)-zh(k-1)))
    enddo
    enddo
    enddo

!--------------------
!  total dynamic forcing:  fluid shear

    do k=1,nk
    do j=1,nj
    do i=1,ni
      deft(i,j,k) = (                                                    &
                -2*rho0(k)*(dum1(i,j,k)+dum2(i,j,k)+dum3(i,j,k))         &
                   )/(cp*rho0(k)*thv0(k))
    enddo
    enddo
    enddo

    do j=1,nj
    do i=1,ni
       deft(i,j,1)=0.
       deft(i,j,nk)=0.
    enddo
    enddo


!    print *,'  dealloc ... '
    deallocate( dum1 )
    deallocate( dum2 )
    deallocate( dum3 )
    deallocate( divx )
!    print *,' ... ok '

!-----------------------------------------------------------------------
!  forcing for fluid extension contribution

  ELSEIF(nloop.eq.5)THEN

!    print *,'  alloc ... '
!    allocate( dum1(ni,nj,nk) )
!    allocate( dum2(ni,nj,nk) )
!    allocate( dum3(ni,nj,nk) )
!    allocate( divx(ni,nj,nk) )
    allocate( dum1(0:ni+1,0:nj+1,0:nk+1) )
    allocate( dum2(0:ni+1,0:nj+1,0:nk+1) )
    allocate( dum3(0:ni+1,0:nj+1,0:nk+1) )
    allocate( divx(0:ni+1,0:nj+1,0:nk+1) )
!    print *,'  ... ok '

    dum1=0.
    dum2=0.
    dum3=0.
    divx=0.

!    deft=0.

!    print *,'entering dum1'
    do k=1,nk
    do j=1,nj
    do i=1,ni-1
      dum1(i,j,k)=(u(i+1,j,k)-u(i,j,k))*rdx 
    enddo
    enddo
    enddo

!    print *,'after dum1, next to dum2'

    do k=1,nk
!    do j=1,nj
    do j=1,nj+1
    do i=1,ni
      dum2(i,j,k)=(v(i,j+1,k)-v(i,j,k))*rdy 
!      vten(i,j,k)=(v(i,j+1,k)-v(i,j,k))*rdy 
    enddo
    enddo
    enddo

!    print *,'after dum2, next to dum3'

    do k=1,nk-1
    do j=1,nj
    do i=1,ni
      dum3(i,j,k)=(w(i,j,k+1)-w(i,j,k))/(zf(k+1)-zf(k))
!      wten(i,j,k)=(w(i,j,k+1)-w(i,j,k))/(zf(k+1)-zf(k))
    enddo
    enddo
    enddo

!    print *,'after dum3, next to divx'

    do k=2,nk-1
    do j=1,nj
    do i=1,ni
!    do j=0,nj+1
!    do i=0,ni+1
      divx(i,j,k)=(0.5*(w(i,j,k+1)+w(i,j,k))                  &
                  *0.5*(w(i,j,k+1)+w(i,j,k)))                 &
                  *(alog(rho0(k-1))-2.*alog(rho0(k))        &
                   +alog(rho0(k+1)))                         &
                  /(2.*(zh(k+1)-zh(k))*(zh(k+1)-zh(k)))
    enddo
    enddo
    enddo

    do j=1,nj
    do i=1,ni
       divx(i,j,1)=0.
       divx(i,j,nk)=0.
       dum1(i,j,1)=0.
       dum1(i,j,nk)=0.
       dum2(i,j,1)=0.
       dum2(i,j,nk)=0.
       dum3(i,j,1)=0.
       dum3(i,j,nk)=0.
  
    enddo
    enddo

!    print *,'final calc'

!--------------------
!  total dynamic forcing:  fluid extension

!    do k=1,nk
!    print *,'zh k-1=',zh(k-1)
    do k=2,nk-1
    do j=1,nj
    do i=1,ni
!      deft(i,j,k) = -rho0(k)*(    &
!                              (u(i+1,j,k)-u(i,j,k))*rdx              &
!                             *(u(i+1,j,k)-u(i,j,k))*rdx  )        &
!                             /(cp*rho0(k)*thv0(k))
! test only...
      deft(i,j,k) = -2.0*rho0(k)*0.25*((u0(k)-u0(k-1))/(zh(k)-zh(k-1)))*(    &
                              (v(i,j,k)-v(i-1,j,k))*rdx              &
                             +(v(i,j,k+1)-v(i-1,j,k+1))*rdx          &
                             +(v(i+1,j,k)-v(i,j,k))*rdx              &
                             +(v(i+1,j,k+1)-v(i,j,k+1))*rdx)         &
                             /(cp*rho0(k)*thv0(k))
! test only...
!                -rho0(k)*(dum1(i,j,k)**2.+dum2(i,j,k)**2.    &
!                         +dum3(i,j,k)**2.)               & 
!                   )/(cp*rho0(k)*thv0(k))
!                -rho0(k)*(dum1(i,j,k)*dum1(i,j,k)    &
!                         +dum2(i,j,k)*dum2(i,j,k)    &
!                         +dum3(i,j,k)*dum3(i,j,k)    &
!                         -divx(i,j,k))               & 
!                   )/(cp*rho0(k)*thv0(k))
    enddo
    enddo
    enddo

!    do k=1,nk
!    do j=1,nj
!       deft(1,j,k)=0.
!       deft(ni,j,k)=0.
!    enddo
!    enddo

!    do j=1,nj
!    do i=1,ni
!       deft(i,j,1)=0.
!       deft(i,j,nk)=0.
!    enddo
!    enddo

!    print *,'  dealloc ... '
    deallocate( dum1 )
    deallocate( dum2 )
    deallocate( dum3 )
    deallocate( divx )
!    print *,' ... ok '

  ENDIF  ! of the first nloop parts

!-----------------------------------------------------------------------
!  p solver

  print *,'  ipb,ipe,jpb,jpe,kpb,kpe = ',ipb,ipe,jpb,jpe,kpb,kpe
  print *,'   alloc 1 '
  allocate(   rhs(ipb:ipe,jpb:jpe) )
  print *,'   alloc 2 '
  allocate( trans(ipb:ipe,jpb:jpe) )

  DO k=1,nk
!    print *,'  k = ',k

    do j=1,nj
    do i=1,ni
!!!      rhs(i,j)=cmplx(def(i,j,k),0.0)
      rhs(i,j)=deft(i,j,k)
    enddo
    enddo

    if(imirror.eq.1)then

      do j=1,nj
      do i=1,ni
        rhs(ipe+1-i,j)=rhs(i,j)
      enddo
      enddo

    endif

    if(jmirror.eq.1)then

      do j=1,nj
      do i=1,ni
        rhs(i,jpe+1-j)=rhs(i,j)
      enddo
      enddo

    endif

    if(imirror.eq.1.and.jmirror.eq.1)then

      do j=1,nj
      do i=1,ni
        rhs(ipe+1-i,jpe+1-j)=rhs(i,j)
      enddo
      enddo

    endif

    trans=fft(rhs)

    do j=jpb,jpe
    do i=ipb,ipe
      deft(i,j,k)=trans(i,j)
    enddo
    enddo

  ENDDO

  DO j=jpb,jpe
  DO i=ipb,ipe

    do k=1,nk
      cfb(k)=2.0d0*( dcos(2.0d0*dpi*dble(i-1)/dble(ipe))          &
                    +dcos(2.0d0*dpi*dble(j-1)/dble(jpe))          &
                    -2.0d0)/(dx*dx) - cfa(k) - cfc(k)
    enddo

    if(i.eq.1.and.j.eq.1)then
      r1(nk+1)=0.0
      r1(nk)=0.0
      do k=nk,2,-1
        r1(k-1)=(deft(i,j,k)-cfc(k)*r1(k+1)-cfb(k)*r1(k))/cfa(k)
      enddo
      do k=1,nk
        deft(i,j,k)=cmplx( r1(k) , 0.0 )
      enddo
    else
      lgbth(1)=-cfc(1)/cfb(1)
      lgbph(1)= deft(i,j,1)/cfb(1)
      do k=2,nk
        lgbth(k)=-cfc(k)/(cfa(k)*lgbth(k-1)+cfb(k))
        lgbph(k)=(deft(i,j,k)-cfa(k)*lgbph(k-1))/(cfa(k)*lgbth(k-1)+cfb(k))
      enddo
      deft(i,j,nk)=lgbph(nk)
      do k=nk-1,1,-1
        deft(i,j,k)=lgbth(k)*deft(i,j,k+1)+lgbph(k)
      enddo
    endif

  ENDDO
  ENDDO

  DO k=1,nk
!    print *,'  k = ',k

    do j=jpb,jpe
    do i=ipb,ipe
      rhs(i,j)=deft(i,j,k)
    enddo
    enddo

    trans=fft(rhs,inv=.true.)

    do j=1,nj
    do i=1,ni
      deft(i,j,k)=real(trans(i,j))
    enddo
    enddo

  ENDDO

  deallocate(   rhs )
  deallocate( trans )

!----- adjust mean pressure --------------------------------------------
!  for pb, adjust mean pi along upper boundary to match total
!  for pd, adjust mean pi along upper boundary to be zero

  pavg = 0.0d0
  frac = 0.0d0

  do j=1,nj
  do i=1,ni
    frac = frac + real(deft(i,j,nk))
  enddo
  enddo

  frac = frac / (ni*nj)

  IF(nloop.eq.1)THEN

    pavg = pavgin

  ENDIF

  print *
  print *,'  frac,pavg = ',frac,pavg

  frac = pavg - frac
  pavg = 0.0d0

! this step is unclear to me...
  do k=1,nk
  do j=1,nj
  do i=1,ni
    deft(i,j,k) = deft(i,j,k) + frac
    if(k.eq.nk) pavg = pavg + deft(i,j,k)
  enddo
  enddo
  enddo

  pavg = pavg / (ni*nj)

  print *,'  pavg      = ',pavg
  print *

!---------------------------------------------------

  if(nloop.eq.1)then

    do k=1,nk
    do j=1,nj
    do i=1,ni
      pb(i,j,k)=real(deft(i,j,k))
    enddo
    enddo
    enddo

  elseif(nloop.eq.2)then

    do k=1,nk
    do j=1,nj
    do i=1,ni
      pdl(i,j,k)=real(deft(i,j,k))
    enddo
    enddo
    enddo

  elseif(nloop.eq.3)then

    do k=1,nk
    do j=1,nj
    do i=1,ni
! JT:  ptdn is the total dynamic pressure
      ptdn(i,j,k)=real(deft(i,j,k))
      pdn(i,j,k)=real(deft(i,j,k))-pdl(i,j,k)
    enddo
    enddo
    enddo

  elseif(nloop.eq.4)then

    print *, 'in nloop 4'

    do k=1,nk
    do j=1,nj
    do i=1,ni
      psh(i,j,k)=real(deft(i,j,k))
    enddo
    enddo
    enddo

  elseif(nloop.eq.5)then

    print *, 'in nloop 5'

    do k=1,nk
    do j=1,nj
    do i=1,ni
      pex(i,j,k)=real(deft(i,j,k))
      psum(i,j,k)=pex(i,j,k)+psh(i,j,k)
    enddo
    enddo
    enddo

  endif

  ENDDO   ! enddo for nloop

  deallocate( uten )
  deallocate( vten )
  deallocate( wten )

! JT:  in the following, compute the dynamic pressure forcing, and buoyancy pressure forcing, for comparison with WK90, etc.
    do k=2,nk-1
    do j=2,nj-1
    do i=2,ni-1
!!      fptdn(i,j,k)=-cp*thv0(k)*(ptdn(i,j,k+1)-ptdn(i,j,k-1))/     &
!!                   (zh(k+1)-zh(k-1))
!      fptdn(i,j,k)=fptdn(i,j,k)-cp*thv0(k)*(pi0(k+1)-pi0(k-1))/   &
!                   (zh(k+1)-zh(k-1))
!!      fpdn(i,j,k)=-cp*thv0(k)*(pdn(i,j,k+1)-pdn(i,j,k-1))/        &
!!                   (zh(k+1)-zh(k-1))
!!      fpdl(i,j,k)=-cp*thv0(k)*(pdl(i,j,k+1)-pdl(i,j,k-1))/        &
!!                   (zh(k+1)-zh(k-1))
!      fpdn(i,j,k)=fpdn(i,j,k)-cp*thv0(k)*(pi0(k+1)-pi0(k-1))/     &
!                   (zh(k+1)-zh(k-1))
!!      fpb(i,j,k)=-cp*thv0(k)*(pb(i,j,k+1)-pb(i,j,k-1))/           &
!!                   (zh(k+1)-zh(k-1)) + bb(i,j,k)
!      fpb(i,j,k)=fpb(i,j,k)-cp*thv0(k)*(pi0(k+1)-pi0(k-1))/       &
!                   (zh(k+1)-zh(k-1))
!      fpex(i,j,k)=-cp*thv0(k)*(pex(i,j,k+1)-pex(i,j,k-1))/        &
!                   (zh(k+1)-zh(k-1))
!      fpsh(i,j,k)=-cp*thv0(k)*(psh(i,j,k+1)-psh(i,j,k-1))/        &
!                   (zh(k+1)-zh(k-1))
!      fpsum(i,j,k)=fpex(i,j,k)+fpsh(i,j,k)

      fpbx(i,j,k)=-cp*thv0(k)*(pb(i+1,j,k)-pb(i-1,j,k))/(2.*dx) + bb(i,j,k)
      fpby(i,j,k)=-cp*thv0(k)*(pb(i,j+1,k)-pb(i,j-1,k))/(2.*dy) + bb(i,j,k)
      fpbz(i,j,k)=-cp*thv0(k)*(pb(i,j,k+1)-pb(i,j,k-1))/           &
                   (zh(k+1)-zh(k-1)) + bb(i,j,k)
      fpdnx(i,j,k)=-cp*thv0(k)*(pdn(i+1,j,k)-pdn(i-1,j,k))/(2.*dx)
      fpdny(i,j,k)=-cp*thv0(k)*(pdn(i,j+1,k)-pdn(i,j-1,k))/(2.*dy)
      fpdnz(i,j,k)=-cp*thv0(k)*(pdn(i,j,k+1)-pdn(i,j,k-1))/        &
                   (zh(k+1)-zh(k-1))
      fpdlx(i,j,k)=-cp*thv0(k)*(pdl(i+1,j,k)-pdl(i-1,j,k))/(2.*dx)
      fpdly(i,j,k)=-cp*thv0(k)*(pdl(i,j+1,k)-pdl(i,j-1,k))/(2.*dy)
      fpdlz(i,j,k)=-cp*thv0(k)*(pdl(i,j,k+1)-pdl(i,j,k-1))/        &
                   (zh(k+1)-zh(k-1))
      fptdnx(i,j,k)=-cp*thv0(k)*(ptdn(i+1,j,k)-ptdn(i-1,j,k))/(2.*dx)
      fptdny(i,j,k)=-cp*thv0(k)*(ptdn(i,j+1,k)-ptdn(i,j,k-1))/(2.*dy)
      fptdnz(i,j,k)=-cp*thv0(k)*(ptdn(i,j,k+1)-ptdn(i,j,k-1))/     &
                   (zh(k+1)-zh(k-1))
    enddo
    enddo
    enddo

  return
  end subroutine pdcomp

