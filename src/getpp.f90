
    program getpp
    use netcdf
    implicit none

!-----------------------------------------------------------------------
!
!  080919:  Program to read CM1 GrADS-format output and calculate
!           pressure perturbations from the subroutine "pdcomp"
!
!  Last modified:  20 February 2013
!
!-----------------------------------------------------------------------

    ! NOTE:  this code (as well as pdcomp) assumes NO TERRAIN
    ! (i.e., perfectly flat lower boundary must be used)

!    integer, parameter :: nx = 200 
!    integer, parameter :: ny = 200
!    integer, parameter :: nz = 80

    ! (NOTE:  can't use this code with horizontal grid stretching)
    ! (dx and dy must be constant)
!    real, parameter :: dx = 100.0
!    real, parameter :: dy = 100.0

    ! (NOTE:  vertical grid spacing is obtained from "cm1.print.out" file)

    ! number of moisture variables for microphysics scheme (ptype)
    ! (liquid and solid components only;  neglecting water vapor and
    !  double-moment variables)
    ! So, in this case:  qc,qr,qi,qs,qg,qhl
!    integer, parameter :: numq = 6

!  for boundary conditions:   1 = periodic ;   2 = open ;   3 = rigid wall
!       wbc = west boundary condition
!       ebc = east boundary condition
!       sbc = south boundary condition
!       nbc = north boundary condition

!    integer, parameter :: wbc = 2
!    integer, parameter :: ebc = 2
!    integer, parameter :: sbc = 2
!    integer, parameter :: nbc = 2

!-----------------------------------------------------------------------
    
    integer :: nx,ny,nz,ptype,wbc,ebc,sbc,nbc
    real    :: dx,dy
    integer :: nk

    integer :: i,j,k,imax,jmax,ntim
!    integer :: nx,ny,nz,num2d,num3d,ith,iprs,iqv,dx,dy
!!    real, dimension(nk) :: zh,prs0,u0,v0,rho0,pi0,th0,qv0
!!    real, dimension(nk+1) :: zf

!!    real*8 :: pavg,pavgin
!!    real, dimension(-2:nx+4,-2:ny+3,nz) :: u
!!    real, dimension(-2:nx+3,-2:ny+4,nz) :: v
!!    real, dimension(-2:nx+3,-2:ny+3,nz+1) :: w
!!    real, dimension(nx,ny,nz) :: bb,pb,pdn,pdl,ppi
    
!!    real, dimension(nx,ny,nz) :: ptdn,fpb,fptdn,fpdn
!!    real, dimension(nx,ny,nz) :: fpdlx,fpdly,fpdlz
!!    real, dimension(nx,ny,nz) :: fpdnx,fpdny,fpdnz
!!    real, dimension(nx,ny,nz) :: fpbx,fpby,fpbz
!!    real, dimension(nx,ny,nz) :: fptdnx,fptdny,fptdnz
!!    real, dimension(nx,ny,nz) :: fpdl
!!    real, dimension(nx,ny,nz) :: pex,psh,fpex,fpsh
!!    real, dimension(nx,ny,nz) :: psum,fpsum
!!    real, dimension(nx,ny,nz) :: qvpert,thpert,thv,prspert
!!!    real, dimension(nx,ny)     :: th,qv,q,wnc,prs
!!    real, dimension(nx+1,ny,nz)   :: unc
!!    real, dimension(nx,ny+1,nz)   :: vnc
!!    real, dimension(nx,ny,nz+1)   :: wnc

!!    real, dimension(nx,ny,nz) :: th,qv,prs
!!    real, dimension(nx,ny,nz) :: qc,qr,qi,qs,qg,qhl

    real, dimension(:),allocatable :: zh,prs0,u0,v0,rho0,pi0,th0,qv0
    real, dimension(:),allocatable :: zf

    real*8 :: pavg,pavgin
    real, dimension(:,:,:),allocatable :: u
    real, dimension(:,:,:),allocatable :: v
    real, dimension(:,:,:),allocatable :: w
    real, dimension(:,:,:),allocatable :: bb,pb,pdn,pdl,ppi

    real, dimension(:,:,:),allocatable :: ptdn,fpb,fptdn,fpdn
    real, dimension(:,:,:),allocatable :: fpdlx,fpdly,fpdlz
    real, dimension(:,:,:),allocatable :: fpdnx,fpdny,fpdnz
    real, dimension(:,:,:),allocatable :: fpbx,fpby,fpbz
    real, dimension(:,:,:),allocatable :: fptdnx,fptdny,fptdnz
    real, dimension(:,:,:),allocatable :: fpdl
    real, dimension(:,:,:),allocatable :: pex,psh,fpex,fpsh
    real, dimension(:,:,:),allocatable :: psum,fpsum
    real, dimension(:,:,:),allocatable :: qvpert,thpert,thv,prspert
!!    real, dimension(nx,ny)     :: th,qv,q,wnc,prs
    real, dimension(:,:,:),allocatable   :: unc
    real, dimension(:,:,:),allocatable   :: vnc
    real, dimension(:,:,:),allocatable   :: wnc

    real, dimension(:,:,:),allocatable :: th,qv,prs
    real, dimension(:,:,:),allocatable :: qc,qr,qi,qs,qg,qhl

!    character(len=100)  :: slash,casename,casepath,outputname
!    character(len=100)  :: filepath,ntimstr,outputpath
!    character(len=100)  :: fl,ncext
!    character(len=100)  :: outputn,output
!    character(len=100)  :: pdout,pdoutnum,flpd
    character(len=100)  :: infile,outfile

    real, dimension(:,:),allocatable    :: qtot

    real :: sumng,sumprs0,sumth0,sumqv0,sumu0,sumv0
    integer, parameter :: avgstart = 1250
    integer, parameter :: avglen = 200

!-----------------------------------------------------------------------

    real, parameter :: p00   = 100000.0
    real, parameter :: rp00  = 1.0/p00
    real, parameter :: rd    = 287.04
    real, parameter :: cp    = 1005.7
    real, parameter :: rv     = 461.5
    real, parameter :: cv     = cp-rd
    real, parameter :: g     = 9.81
    real, parameter :: reps   = rv/rd
    real, parameter :: repsm1 = rv/rd-1.0
    real, parameter :: rovcp  = rd/cp
    real, parameter :: rddcp  = rd/cp
    real, parameter :: rddcv  = rd/cv
    real, parameter :: rddrv  = rd/rv
    real, parameter :: cvdrd  = cv/rd
    real, parameter :: cpdrd  = cp/rd

    namelist /inputs/ infile,outfile,nx,ny,nz,dx,dy,ptype,wbc,ebc,sbc,nbc

!-----------------Read netCDF--------------------

    write(*,*) 'Reading in namelist'

    open(8,file='pdcomp.input')
    read(8,nml=inputs)

    write(*,*)
    write(*,*) 'Input variables:'
    write(*,*) '   infile      = ',infile
    write(*,*) '   outfile     = ',outfile
    write(*,*) '   nx          = ',nx
    write(*,*) '   ny          = ',ny
    write(*,*) '   nz          = ',nz
!    write(*,*) '   maxnx       = ',maxnx
!    write(*,*) '   maxny       = ',maxny
!    write(*,*) '   maxnz       = ',maxnz
    write(*,*) '   dx          = ',dx
    write(*,*) '   dy          = ',dy
!    write(*,*) '   istart      = ',istart
!    write(*,*) '   jstart      = ',jstart
    write(*,*) '   ptype       = ',ptype
    write(*,*) '   wbc         = ',wbc
    write(*,*) '   ebc         = ',ebc
    write(*,*) '   sbc         = ',sbc
    write(*,*) '   nbc         = ',nbc

    if(ptype.ne.27.and.ptype.ne.5)then
      print*,'Invalid ptype. Can only use Morrison/NSSLDM without modification'
      STOP
    endif

    nk = nz
!    print*,nk

    allocate( zh(nk) )
    allocate( prs0(nk) )
    allocate( u0(nk) )
    allocate( v0(nk) )
    allocate( rho0(nk) )
    allocate( pi0(nk) )
    allocate( th0(nk) )
    allocate( qv0(nk) )
    allocate( zf(nk+1) )
    allocate( u(-2:nx+4,-2:ny+3,nz) )
    allocate( v(-2:nx+3,-2:ny+4,nz) )
    allocate( w(-2:nx+3,-2:ny+3,nz+1) )
    allocate( bb(nx,ny,nz) )
    allocate( pb(nx,ny,nz) )
    allocate( pdn(nx,ny,nz) )
    allocate( pdl(nx,ny,nz) )
    allocate( ptdn(nx,ny,nz) )
    allocate( ppi(nx,ny,nz) )
    allocate( fpdlx(nx,ny,nz) )
    allocate( fpdly(nx,ny,nz) )
    allocate( fpdlz(nx,ny,nz) )
    allocate( fpdnx(nx,ny,nz) )
    allocate( fpdny(nx,ny,nz) )
    allocate( fpdnz(nx,ny,nz) )
    allocate( fpbx(nx,ny,nz) )
    allocate( fpby(nx,ny,nz) )
    allocate( fpbz(nx,ny,nz) )
    allocate( fptdnx(nx,ny,nz) )
    allocate( fptdny(nx,ny,nz) )
    allocate( fptdnz(nx,ny,nz) )
    allocate( pex(nx,ny,nz) )
    allocate( psh(nx,ny,nz) )
    allocate( fpex(nx,ny,nz) )
    allocate( fpsh(nx,ny,nz) )
    allocate( psum(nx,ny,nz) )
    allocate( fpsum(nx,ny,nz) )
    allocate( qvpert(nx,ny,nz) )
    allocate( thpert(nx,ny,nz) )
    allocate( thv(nx,ny,nz) )
    allocate( prspert(nx,ny,nz) )
    allocate( th(nx,ny,nz) )
    allocate( qv(nx,ny,nz) )
    allocate( prs(nx,ny,nz) )
    allocate( qc(nx,ny,nz) )
    allocate( qr(nx,ny,nz) )
    allocate( qi(nx,ny,nz) )
    allocate( qs(nx,ny,nz) )
    allocate( qg(nx,ny,nz) )
    allocate( qhl(nx,ny,nz) )
    allocate( unc(nx+1,ny,nz) )
    allocate( vnc(nx,ny+1,nz) )
    allocate( wnc(nx,ny,nz+1) )
    allocate( qtot(nx,ny) )

    allocate( fpb(nx,ny,nz) )
    allocate( fpdl(nx,ny,nz) )
    allocate( fpdn(nx,ny,nz) )
    allocate( fptdn(nx,ny,nz) )

    call readnc0(infile,nz,zh,zf,th0,qv0,pi0,u0,v0,prs0,rho0)

    print*,'Checkpoint 0'
!-------Calculate values needed by pdcomp--------

!    call readnc(filepath,nx,ny,nz,ntim,istart,jstart,nk,k,th,  &
!                qv,prs,wnc,unc,vnc,qc,qr,qi,qs,qg)

!    call readnc(outputpath,nx,ny,nz,ntim,nk,k,istart,jstart,wnc,unc,vnc,&
!                      prspert,thpert,qvpert)
   
    call readnc(infile,nx,ny,nz,th,qv,qc,qr,qi,qs,qg,qhl,  &
                prs,unc,vnc,wnc,thpert,qvpert,prspert)

    print*,'Checkpoint 1'

!    do k=1,nz
!      sumng = 0
!      sumprs0 = 0
!      sumqv0 = 0
!      sumth0 = 0
!      sumu0 = 0
!      sumv0 = 0
!    do j=1,ny
!    do i=1,avglen
!      sumng = sumng+1
!      sumprs0 = sumprs0 + prs(avgstart+i,j,k)
!      sumqv0 = sumqv0 + qv(avgstart+i,j,k)
!      sumth0 = sumth0 + th(avgstart+i,j,k)
!      sumu0 = sumu0 + unc(avgstart+i,j,k)
!      sumv0 = sumv0 + vnc(avgstart+i,j,k)
!    enddo
!    enddo
!      prs0(k) = sumprs0/sumng
!      qv0(k) = sumqv0/sumng
!      th0(k) = sumth0/sumng
!      u0(k) = sumu0/sumng
!      v0(k) = sumv0/sumng
!    enddo

    do k=1,nz
      qtot = 0.0
    do j=1,ny
    do i=1,nx
      if(ptype.eq.27)then
        qtot(i,j) = qtot(i,j)+qc(i,j,k)+qr(i,j,k)+qi(i,j,k) &
                    +qs(i,j,k)+qg(i,j,k)+qhl(i,j,k)
      elseif(ptype.eq.5)then
        qtot(i,j) = qtot(i,j)+qc(i,j,k)+qr(i,j,k)+qi(i,j,k) &
                    +qs(i,j,k)+qhl(i,j,k)
      endif
      ppi(i,j,k) = (prs(i,j,k)*rp00)**rovcp - pi0(k)
      pb(i,j,k)  = g*( (th(i,j,k)-th0(k))/th0(k)  &
                     +repsm1*(qv(i,j,k)-qv0(k))   &
                     -qtot(i,j))

!      pb(i,j,k)  = g*( (th(i,j,k)-th0(k))/th0(k)  &
!                     +repsm1*(qv(i,j,k)-qv0(k)))

!      prs(i,j,k) = prs0(k)+prspert(i,j,k)
!      th(i,j,k) = thpert(i,j,k)+th0(k)
!      qv(i,j,k) = qvpert(i,j,k)+qv0(k)

      thv(i,j,k) = th(i,j,k)*(1+0.61*qv(i,j,k))
!      pb(i,j,k)  = g*((thv(i,j,k)-thv0(k))/thv0(k))
      bb(i,j,k)  = pb(i,j,k)
!      w(i,j,k)   = wnc(i,j,k)
!      if(i.eq.180 .and. j.eq.180) then
!        print*, 'Thv is', thv(i,j,k)
!        print*, 'u-wind is', unc(i,j,k)
!      endif
    enddo
    enddo
    enddo

!    print*, pb(100,100,:)

    print*,'Checkpoint 2'

    do k=1,nz
    do j=1,ny+1
    do i=1,nx
      v(i,j,k) = vnc(i,j,k)
    enddo
    enddo
    do j=1,ny
    do i=1,nx+1
      u(i,j,k) = unc(i,j,k)
    enddo
    enddo
    enddo

    do k=1,nz+1
    do j=1,ny
    do i=1,nx
      w(i,j,k) = wnc(i,j,k)
    enddo
    enddo
    enddo

!      print *,'  Sample th :',th(1,1),th(nx,ny)
!      print *,'  Sample ppi:',ppi(1,1,k),ppi(nx,ny,k)
!      print *,'  Sample qv :',qv(1,1),qv(nx,ny)


!    ENDDO983040000

    print*,'Checkpoint 3'
!    print *,'  ... done reading data '
    print *

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    pavg = 0.0d0

    do j=1,ny
    do i=1,nx
      pavg = pavg + ppi(i,j,nk)
    enddo
    enddo

    pavg = pavg / (nx*ny)
    pavgin = pavg

    print *,'  pavg = ',pavg

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Get new pressure perturbation: 

    print *,'  calling pdcomp'
    print *,'    dx,dy = ',dx,dy
    print *,'    nx,ny,nz = ',nx,ny,nz

    call pdcomp(nx,ny,nz,wbc,ebc,sbc,nbc,dx,dy,          &
                zh,rho0,th0,qv0,pi0,u0,v0,u,v,w,pavgin,  &
                pb,pdn,pdl,ptdn,fpb,fptdn,fpdn,          &
                pex,psh,psum,fpex,fpsh,fpsum,fpdl,fpdlx, &
                fpdly,fpdlz,fpdnx,fpdny,fpdnz,fpbx,fpby, &
                fpbz,fptdnx,fptdny,fptdnz)

    print *,'  returned from pdcomp'

! Write data to netCDF-format file:
!!    call writenc(fl,nx,ny,nz,fptdn,fpb,fpdn,fpdl)
!    call writenc(fl,ntim,nx,ny,nz,fptdn,fpdn,fpdl)

    call writenc(outfile,nx,ny,nz,pb,pdl,pdn,fptdnx,  &
                 fptdny,fptdnz,fpbx,fpby,fpbz,fpdnx,  &
                 fpdny,fpdnz,fpdlx,fpdly,fpdlz)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Write data to GrADS-format file:
!    datfile = trim(fl) // trim(datext)
!    open(unit=45,file=datfile,status='unknown',   &
!         form='unformatted',access='direct',recl=4*nx*ny)
!    orec = 1

!    do n=1,5
!    do n=1,9
!    do n=1,4
!    do k=1,nz
!      if(n.eq.1) write(45,rec=orec) (( bb(i,j,k),i=1,nx),j=1,ny)
!      if(n.eq.2) write(45,rec=orec) ((ppi(i,j,k),i=1,nx),j=1,ny)
!      if(n.eq.3) write(45,rec=orec) (( pb(i,j,k),i=1,nx),j=1,ny)
!      if(n.eq.4) write(45,rec=orec) ((pdn(i,j,k),i=1,nx),j=1,ny)
!      if(n.eq.5) write(45,rec=orec) ((pdl(i,j,k),i=1,nx),j=1,ny)
!      if(n.eq.6) write(45,rec=orec) ((ptdn(i,j,k),i=1,nx),j=1,ny)
!      if(n.eq.7) write(45,rec=orec) ((pex(i,j,k),i=1,nx),j=1,ny)
!      if(n.eq.8) write(45,rec=orec) ((psh(i,j,k),i=1,nx),j=1,ny)
!      if(n.eq.9) write(45,rec=orec) ((psum(i,j,k),i=1,nx),j=1,ny)
!      if(n.eq.1) write(45,rec=orec) ((fptdn(i,j,k),i=1,nx),j=1,ny)
!      if(n.eq.2) write(45,rec=orec) ((fpb(i,j,k),i=1,nx),j=1,ny)
!      if(n.eq.3) write(45,rec=orec) ((fpdn(i,j,k),i=1,nx),j=1,ny)
!      if(n.eq.4) write(45,rec=orec) ((fpex(i,j,k),i=1,nx),j=1,ny)
!      if(n.eq.5) write(45,rec=orec) ((fpsh(i,j,k),i=1,nx),j=1,ny)
!      if(n.eq.4) write(45,rec=orec) ((fpdl(i,j,k),i=1,nx),j=1,ny)

!      orec=orec+1
!    enddo
!    enddo
!    close(unit=45)

!--------------------------------------------
!  grads descriptor file for diagnostic output:
!    ctlfile = trim(fl) // trim(ctlext)
!    open(unit=30,file=ctlfile,status='unknown')
!    write(30,201)
!!!    write(30,212)
!    write(30,202)
!    write(30,203)
!    write(30,204) nx,0.001*0.5*dx,0.001*dx
!    write(30,205) ny,0.001*0.5*dy,0.001*dy
!    write(30,206) nk
!    do k=1,nk
!      write(30,211) 0.001*zh(k)
!211   format(2x,f12.6)
!    enddo
!    write(30,207)
!    write(30,208) 5
!    write(30,208) 9
!    write(30,208) 4
!    write(30,209) 'b       ',nk,'buoyancy from model                               '
!    write(30,209) 'ppi     ',nk,'actual pi^prime from model                        '
!    write(30,209) 'pb      ',nk,'diagnosed pi^prime:  buoyant component            '
!    write(30,209) 'pdn     ',nk,'diagnosed pi^prime:  nonlinear dynamic component  '
!    write(30,209) 'pdl     ',nk,'diagnosed pi^prime:  linear dynamic component     '
!    write(30,209) 'ptdn    ',nk,'diagnosed pi^prime:  total dynamic component      '
!    write(30,209) 'pex     ',nk,'diagnosed pi^prime:  exten dynamic component      '
!    write(30,209) 'psh     ',nk,'diagnosed pi^prime:  shear dynamic component      '
!    write(30,209) 'psum    ',nk,'diagnosed pi^prime:  sum of shr + ext compnt      '
!    write(30,209) 'fptdn   ',nk,'total dynamic pressure force                      '
!    write(30,209) 'fpb     ',nk,'buoyancy pressure force                           '
!    write(30,209) 'fpdn    ',nk,'nonlinear dynamic pressure force                  '
!    write(30,209) 'fpex    ',nk,'nonlinear dynamic pressure ext                    '
!    write(30,209) 'fpsh    ',nk,'nonlinear dynamic pressure shr                    '
!    write(30,209) 'fpsum   ',nk,'nonlinear dynamic pressure sum                    '
!    write(30,209) 'fpdl    ',nk,'linear dynamic pressure                           '
!    write(30,210)
!    close(unit=30)

!201 format('dset ^pdcomp.dat')
!212 format('options template')
!202 format('title CM1 output')
!203 format('undef -99999999.')
!204 format('xdef ',i6,' linear ',f12.6,1x,f12.6)
!205 format('ydef ',i6,' linear ',f12.6,1x,f12.6)
!206 format('zdef ',i6,' levels')
!207 format('tdef 1000 linear 00Z03JUL0001 1YR')
!208 format('vars ',i3)
!209 format(a8,1x,i6,' 99 ',a50)
!210 format('endvars')

!------------------------------------------------------------------

    end program getpp


!:=========================================================

    subroutine readnc0(outputpath,nz,zh,zf,th0,qv0,pi0,u0,v0,prs0,rho0)
    use netcdf
    implicit none

    integer, intent(in) :: nz
    character(len=100), intent(in) :: outputpath
    real, dimension(nz), intent(out) :: th0,qv0,zh,zf,pi0,u0,v0,prs0
    real, dimension(nz), intent(out) :: rho0
    real, dimension(nz) :: t0
    integer :: ncid,status,p
    integer :: u0varid,v0varid,th0varid,qv0varid,prs0varid,zvarid,zfvarid

    real, parameter :: p00    = 100000.0
    real, parameter :: rp00   = 1.0/p00
    real, parameter :: rd     = 287.04
    real, parameter :: cp     = 1005.7
    real, parameter :: rv     = 461.5
    real, parameter :: cv     = cp-rd
    real, parameter :: g      = 9.81
    real, parameter :: reps   = rv/rd
    real, parameter :: repsm1 = rv/rd-1.0
    real, parameter :: rovcp  = rd/cp
    real, parameter :: rddcp  = rd/cp
    real, parameter :: rddcv  = rd/cv
    real, parameter :: rddrv  = rd/rv
    real, parameter :: cvdrd  = cv/rd
    real, parameter :: cpdrd  = cp/rd


    ncid = 0
!------------------Open netCDF-----------------------------
    status = nf90_open(outputpath,nf90_nowrite,ncid)
    print*,ncid
    if(status /= nf90_NoErr) print*,nf90_strerror(status)
!----------Get variables needed from netcdf----------------
    status = nf90_inq_varid(ncid,"z",zvarid)
    status = nf90_get_var(ncid,zvarid,zh,start=(/1/),count=(/nz/))

    status = nf90_inq_varid(ncid,"zf",zfvarid)
    status = nf90_get_var(ncid,zfvarid,zf,start=(/1/),count=(/nz+1/))

    status = nf90_inq_varid(ncid,"th0",th0varid)
    status = nf90_get_var(ncid,th0varid,th0,start=(/1,1,1,1/),count=(/1,1,nz,1/))

    status = nf90_inq_varid(ncid,"qv0",qv0varid)
    status = nf90_get_var(ncid,qv0varid,qv0,start=(/1,1,1,1/),count=(/1,1,nz,1/))

    status = nf90_inq_varid(ncid,"prs0",prs0varid)
    status = nf90_get_var(ncid,prs0varid,prs0,start=(/1,1,1,1/),count=(/1,1,nz,1/))

    status = nf90_inq_varid(ncid,"u0",u0varid)
    status = nf90_get_var(ncid,u0varid,u0,start=(/1,1,1,1/),count=(/1,1,nz,1/))

    status = nf90_inq_varid(ncid,"v0",v0varid)
    status = nf90_get_var(ncid,v0varid,v0,start=(/1,1,1,1/),count=(/1,1,nz,1/))

!    print*,zh(:)
!    print*,zf(:)
    do p=1,nz
!      print*,p
      pi0(p) = (prs0(p)/p00)**(rd/cp)
      t0(p) = th0(p)*(prs0(p)/p00)**(0.2854)
!      thv0(p) = th0(p)*(1+repsm1*qv0(p))
      rho0(p) = prs0(p)/(rd*th0(p)*pi0(p)*(1.0+reps*qv0(p)))
      zh(p)   = zh(p)*1000
      zf(p)   = zf(p)*1000
    enddo


!    print*,zh(:)

!------------------Close netCDF----------------------------
    status = nf90_close(ncid)

    return
    end
!=========================================================

!=================== Read NetCDF ==========================
    subroutine readnc(outputpath,nx,ny,nz,   &
                      th,qv,qc,qr,qi,qs,qg,qhl,prs,unc,vnc,wnc,thpert,qvpert,&
                      prspert)
    use netcdf
    implicit none

    integer, intent(in) :: nx,ny,nz
    character(len=100), intent(in) :: outputpath
    real, dimension(nx,ny,nz), intent(out) :: wnc,th,qv,prs,qc,qr,qi,qs,qg,qhl
    real, dimension(nx,ny,nz),   intent(out) :: thpert,qvpert,prspert
    real, dimension(nx+1,ny,nz), intent(out) :: unc
    real, dimension(nx,ny+1,nz), intent(out) :: vnc

    integer :: ncid,thpertvarid,qvpertvarid,uvarid,vvarid,wvarid,status
    integer :: prspertvarid,thvarid,qvvarid,prsvarid
    integer :: qcvarid,qrvarid,qivarid,qsvarid,qgvarid,qhlvarid

    real, parameter :: p00   = 100000.0
    real, parameter :: rp00  = 1.0/p00
    real, parameter :: rd    = 287.04
    real, parameter :: cp    = 1005.7
    real, parameter :: rv     = 461.5
    real, parameter :: cv     = cp-rd
    real, parameter :: g     = 9.81
    real, parameter :: reps   = rv/rd
    real, parameter :: repsm1 = rv/rd-1.0
    real, parameter :: rovcp  = rd/cp
    real, parameter :: rddcp  = rd/cp
    real, parameter :: rddcv  = rd/cv
    real, parameter :: rddrv  = rd/rv
    real, parameter :: cvdrd  = cv/rd
    real, parameter :: cpdrd  = cp/rd

    ncid = 0
!------------------Open netCDF-----------------------------
    status = nf90_open(outputpath,nf90_nowrite,ncid)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)
!    print*,'Problem opening netcdf...'
!    print*,ncid
!----------Get variables needed from netcdf----------------
    status = nf90_inq_varid(ncid,"th",thvarid)
    status = nf90_get_var(ncid,thvarid,th,start=(/1,1,1,1/),count=(/nx,ny,nz,1/))

    status = nf90_inq_varid(ncid,"qv",qvvarid)
    status = nf90_get_var(ncid,qvvarid,qv,start=(/1,1,1,1/),count=(/nx,ny,nz,1/))

    status = nf90_inq_varid(ncid,"qc",qcvarid)
    status = nf90_get_var(ncid,qcvarid,qc,start=(/1,1,1,1/),count=(/nx,ny,nz,1/))

    status = nf90_inq_varid(ncid,"qr",qrvarid)
    status = nf90_get_var(ncid,qrvarid,qr,start=(/1,1,1,1/),count=(/nx,ny,nz,1/))

    status = nf90_inq_varid(ncid,"qi",qivarid)
    status = nf90_get_var(ncid,qivarid,qi,start=(/1,1,1,1/),count=(/nx,ny,nz,1/))

    status = nf90_inq_varid(ncid,"qs",qsvarid)
    status = nf90_get_var(ncid,qsvarid,qs,start=(/1,1,1,1/),count=(/nx,ny,nz,1/))

    status = nf90_inq_varid(ncid,"qg",qgvarid)
    status = nf90_get_var(ncid,qgvarid,qg,start=(/1,1,1,1/),count=(/nx,ny,nz,1/))

    status = nf90_inq_varid(ncid,"qhl",qhlvarid)
    status = nf90_get_var(ncid,qhlvarid,qhl,start=(/1,1,1,1/),count=(/nx,ny,nz,1/))

    status = nf90_inq_varid(ncid,"prs",prsvarid)
    status = nf90_get_var(ncid,prsvarid,prs,start=(/1,1,1,1/),count=(/nx,ny,nz,1/))

    status = nf90_inq_varid(ncid,"u",uvarid)
    status = nf90_get_var(ncid,uvarid,unc,start=(/1,1,1,1/),count=(/nx+1,ny,nz,1/))

    status = nf90_inq_varid(ncid,"v",vvarid)
    status = nf90_get_var(ncid,vvarid,vnc,start=(/1,1,1,1/),count=(/nx,ny+1,nz,1/))

    status = nf90_inq_varid(ncid,"w",wvarid)
    status = nf90_get_var(ncid,wvarid,wnc,start=(/1,1,1,1/),count=(/nx,ny,nz+1,1/))

    status = nf90_inq_varid(ncid,"thpert",thpertvarid)
    status = nf90_get_var(ncid,thpertvarid,thpert,start=(/1,1,1,1/),count=(/nx,ny,nz,1/))

    status = nf90_inq_varid(ncid,"qvpert",qvpertvarid)
    status = nf90_get_var(ncid,qvpertvarid,qvpert,start=(/1,1,1,1/),count=(/nx,ny,nz,1/))

    status = nf90_inq_varid(ncid,"prspert",prspertvarid)
    status = nf90_get_var(ncid,prspertvarid,prspert,start=(/1,1,1,1/),count=(/nx,ny,nz,1/))

!------------------Close netCDF----------------------------
    status = nf90_close(ncid)
    print*,'Finished reading netCDF'

    return
    end
!:=========================================================

!=================== Write netCDF ==========================
!    subroutine writenc(fl,nx,ny,nz,fptdn,fpb,fpdn,fpdl)
!    subroutine writenc(fl,ntim,nx,ny,nz,fptdn,fpdn,fpdl)
    subroutine writenc(fl,nx,ny,nz,pb,pdl,pdn,fptdnx,  &
                       fptdny,fptdnz,fpbx,fpby,fpbz,   &
                       fpdnx,fpdny,fpdnz,fpdlx,fpdly,  &
                       fpdlz)

    use netcdf
    implicit none
 
    character(len=100), intent(inout) :: fl
    integer, intent(in)            :: nx,ny,nz
!    real, dimension(nx,ny,nz), intent(in) :: fptdn,fpb,fpdn,fpdl
    real, dimension(nx,ny,nz), intent(in) :: fpdlx,fpdly,fpdlz
    real, dimension(nx,ny,nz), intent(in) :: fpdnx,fpdny,fpdnz
    real, dimension(nx,ny,nz), intent(in) :: fpbx,fpby,fpbz
    real, dimension(nx,ny,nz), intent(in) :: fptdnx,fptdny,fptdnz
    real, dimension(nx,ny,nz), intent(in) :: pb,pdl,pdn

    integer :: ncid,status,niDimID,njDimID,nkDimID,timeDimID
    integer :: fptdnVarID,fpbVarID,fpdnVarID,fpdlVarID
    integer :: fptdnxVarID,fptdnyVarID,fptdnzVarID
    integer :: fpbxVarID,fpbyVarID,fpbzVarID
    integer :: fpdlxVarID,fpdlyVarID,fpdlzVarID
    integer :: fpdnxVarID,fpdnyVarID,fpdnzVarID
    integer :: pbVarID,pdlVarID,pdnVarID

!    character(len=100) :: outputpath

    fl = trim(fl)

!----------------- Create and open netCDF -----------------
    status = nf90_create(fl,NF90_64BIT_OFFSET,ncid)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)

!------------ Define dimensions and variables -------------

    status = nf90_def_dim(ncid,"ni",nx,niDimID)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)
    
    status = nf90_def_dim(ncid,"nj",ny,njDimID)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)

    status = nf90_def_dim(ncid,"nk",nz,nkDimID)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)

    status = nf90_def_dim(ncid,"time",1,timeDimID)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)

!    status = nf90_def_var(ncid,"fptdn",nf90_float, &
!                         (/niDimID,njDimID,nkDimID,timeDimID/),fptdnVarID)
!    if(status /= nf90_NoErr) print*,nf90_strerror(status)

    status = nf90_def_var(ncid,"fptdnx",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),fptdnxVarID)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)

    status = nf90_def_var(ncid,"fptdny",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),fptdnyVarID)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)

    status = nf90_def_var(ncid,"fptdnz",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),fptdnzVarID)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)

!    status = nf90_def_var(ncid,"fpb",nf90_float, &
!                         (/niDimID,njDimID,nkDimID,timeDimID/),fpbVarID)
!    if(status /= nf90_NoErr) print*,nf90_strerror(status)

    status = nf90_def_var(ncid,"fpbx",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),fpbxVarID)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)

    status = nf90_def_var(ncid,"fpby",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),fpbyVarID)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)

    status = nf90_def_var(ncid,"fpbz",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),fpbzVarID)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)

!    status = nf90_def_var(ncid,"fpdn",nf90_float, &
!                         (/niDimID,njDimID,nkDimID,timeDimID/),fpdnVarID)
!    if(status /= nf90_NoErr) print*,nf90_strerror(status)

    status = nf90_def_var(ncid,"fpdnx",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),fpdnxVarID)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)

    status = nf90_def_var(ncid,"fpdny",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),fpdnyVarID)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)

     status = nf90_def_var(ncid,"fpdnz",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),fpdnzVarID)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)

!    status = nf90_def_var(ncid,"fpdl",nf90_float, &
!                         (/niDimID,njDimID,nkDimID,timeDimID/),fpdlVarID)
!    if(status /= nf90_NoErr) print*,nf90_strerror(status)

    status = nf90_def_var(ncid,"fpdlx",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),fpdlxVarID)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)

    status = nf90_def_var(ncid,"fpdly",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),fpdlyVarID)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)

    status = nf90_def_var(ncid,"fpdlz",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),fpdlzVarID)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)

    status = nf90_def_var(ncid,"pb",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),pbVarID)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)

    status = nf90_def_var(ncid,"pdl",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),pdlVarID)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)

    status = nf90_def_var(ncid,"pdn",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),pdnVarID)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)

!    status = nf90_def_var(ncid,"fpex",nf90_float, &
!                         (/niDimID,njDimID,nkDimID,timeDimID/),fpexVarID)
!    if(status /= nf90_NoErr) print*,nf90_strerror(status)

!    status = nf90_def_var(ncid,"fpsh",nf90_float, &
!                         (/niDimID,njDimID,nkDimID,timeDimID/),fpshrVarID)
!    if(status /= nf90_NoErr) print*,nf90_strerror(status)

    print*,'Finished writing vars'

    status = nf90_enddef(ncid)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)
!------------ Write dimensions and variables -------------

!    status = nf90_put_var(ncid,niDimID,nx)
!    if(status /= nf90_NoErr) print*,nf90_strerror(status)

!    status = nf90_put_var(ncid,njDimID,ny)
!    if(status /= nf90_NoErr) print*,nf90_strerror(status)

!    status = nf90_put_var(ncid,nkDimID,nz)
!    if(status /= nf90_NoErr) print*,nf90_strerror(status)

!!    status = nf90_put_var(ncid,fptdnVarID,fptdn)
!!    if(status /= nf90_NoErr) print*,nf90_strerror(status)

    status = nf90_put_var(ncid,fptdnxVarID,fptdnx)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)

    status = nf90_put_var(ncid,fptdnyVarID,fptdny)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)

    status = nf90_put_var(ncid,fptdnzVarID,fptdnz)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)

!!    status = nf90_put_var(ncid,fpbVarID,fpb)
!!    if(status /= nf90_NoErr) print*,nf90_strerror(status)

    status = nf90_put_var(ncid,fpbxVarID,fpbx)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)

    status = nf90_put_var(ncid,fpbyVarID,fpby)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)

    status = nf90_put_var(ncid,fpbzVarID,fpbz)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)

!!    status = nf90_put_var(ncid,fpdnVarID,fpdn)
!!    if(status /= nf90_NoErr) print*,nf90_strerror(status)

    status = nf90_put_var(ncid,fpdnxVarID,fpdnx)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)

    status = nf90_put_var(ncid,fpdnyVarID,fpdny)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)

    status = nf90_put_var(ncid,fpdnzVarID,fpdnz)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)

!!    status = nf90_put_var(ncid,fpdlVarID,fpdl)
!!    if(status /= nf90_NoErr) print*,nf90_strerror(status)

    status = nf90_put_var(ncid,fpdlxVarID,fpdlx)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)

    status = nf90_put_var(ncid,fpdlyVarID,fpdly)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)

    status = nf90_put_var(ncid,fpdlzVarID,fpdlz)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)

    status = nf90_put_var(ncid,pbVarID,pb)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)

    status = nf90_put_var(ncid,pdlVarID,pdl)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)

    status = nf90_put_var(ncid,pdnVarID,pdn)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)

    status = nf90_close(ncid)
    if(status /= nf90_NoErr) print*,nf90_strerror(status)

    return
    end
