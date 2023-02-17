!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

    subroutine bcu(ni,nj,nk,wbc,ebc,sbc,nbc,u)
    implicit none

    integer :: ni,nj,nk,wbc,ebc,sbc,nbc
    real, dimension(-2:ni+4,-2:nj+3,nk) :: u

    integer :: i,j,k,ib,ie,jb,je

    ib = -2
    ie = ni+3
    jb = -2
    je = nj+3

    DO k=1,nk

!-----------------------------------

      if(wbc.eq.1)then
        do j=jb,je
          u(   0,j,k)=u(ni  ,j,k)
          u(  -1,j,k)=u(ni-1,j,k)
          u(  -2,j,k)=u(ni-2,j,k)
        enddo
      elseif(wbc.eq.2)then
        do j=0,nj+1
          u(   0,j,k)=u(1,j,k)
          u(  -1,j,k)=u(1,j,k)
          u(  -2,j,k)=u(1,j,k)
        enddo
      elseif(wbc.eq.3)then
        do j=0,nj+1
          u(   1,j,k)=0.
          u(   0,j,k)=-u(   2,j,k)
          u(  -1,j,k)=-u(   3,j,k)
          u(  -2,j,k)=-u(   4,j,k)
        enddo
      endif

!-----------------------------------

      if(ebc.eq.1)then
        do j=jb,je
          u(ni+2,j,k)=u(2,j,k)
          u(ni+3,j,k)=u(3,j,k)
          u(ni+4,j,k)=u(4,j,k)
        enddo
      elseif(ebc.eq.2)then
        do j=0,nj+1
          u(ni+2,j,k)=u(ni+1,j,k)
          u(ni+3,j,k)=u(ni+1,j,k)
          u(ni+4,j,k)=u(ni+1,j,k)
        enddo
      elseif(ebc.eq.3)then
        do j=0,nj+1
          u(ni+1,j,k)=0.
          u(ni+2,j,k)=-u(ni  ,j,k)
          u(ni+3,j,k)=-u(ni-1,j,k)
          u(ni+4,j,k)=-u(ni-2,j,k)
        enddo
      endif

!-----------------------------------

      if(sbc.eq.1)then
        do i=ib,ie+1
          u(i,   0,k)=u(i,nj  ,k)
          u(i,  -1,k)=u(i,nj-1,k)
          u(i,  -2,k)=u(i,nj-2,k)
        enddo
      elseif(sbc.eq.2)then
        do i=0,ni+2
          u(i,   0,k)=u(i, 1,k)
          u(i,  -1,k)=u(i, 1,k)
          u(i,  -2,k)=u(i, 1,k)
        enddo
      elseif(sbc.eq.3)then
        do i=0,ni+2
          u(i,   0,k)=u(i, 1,k)
          u(i,  -1,k)=u(i, 2,k)
          u(i,  -2,k)=u(i, 3,k)
        enddo
      endif

!-----------------------------------

      if(nbc.eq.1)then
        do i=ib,ie+1
          u(i,nj+1,k)=u(i,   1,k)
          u(i,nj+2,k)=u(i,   2,k)
          u(i,nj+3,k)=u(i,   3,k)
        enddo
      elseif(nbc.eq.2)then
        do i=0,ni+2
          u(i,nj+1,k)=u(i,nj,k)
          u(i,nj+2,k)=u(i,nj,k)
          u(i,nj+3,k)=u(i,nj,k)
        enddo
      elseif(nbc.eq.3)then
        do i=0,ni+2
          u(i,nj+1,k)=u(i,nj  ,k)
          u(i,nj+2,k)=u(i,nj-1,k)
          u(i,nj+3,k)=u(i,nj-2,k)
        enddo
      endif

!-----------------------------------

    ENDDO

    return
    end subroutine bcu

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

    subroutine bcv(ni,nj,nk,wbc,ebc,sbc,nbc,v)
    implicit none

    integer :: ni,nj,nk,wbc,ebc,sbc,nbc
    real, dimension(-2:ni+3,-2:nj+4,nk) :: v

    integer :: i,j,k,ib,ie,jb,je

    ib = -2
    ie = ni+3
    jb = -2
    je = nj+3

    DO k=1,nk

!-----------------------------------

      if(sbc.eq.1)then
        do i=ib,ie
          v(i,   0,k)=v(i,nj  ,k)
          v(i,  -1,k)=v(i,nj-1,k)
          v(i,  -2,k)=v(i,nj-2,k)
        enddo
      elseif(sbc.eq.2)then
        do i=0,ni+1
          v(i,   0,k)=v(i,1,k)
          v(i,  -1,k)=v(i,1,k)
          v(i,  -2,k)=v(i,1,k)
        enddo
      elseif(sbc.eq.3)then
        do i=0,ni+1
          v(i,   1,k)=0.
          v(i,   0,k)=-v(i,2,k)
          v(i,  -1,k)=-v(i,3,k)
          v(i,  -2,k)=-v(i,4,k)
        enddo
      endif

!-----------------------------------

      if(nbc.eq.1)then
        do i=ib,ie
          v(i,nj+2,k)=v(i,   2,k)
          v(i,nj+3,k)=v(i,   3,k)
          v(i,nj+4,k)=v(i,   4,k)
        enddo
      elseif(nbc.eq.2)then
        do i=0,ni+1
          v(i,nj+2,k)=v(i,nj+1,k)
          v(i,nj+3,k)=v(i,nj+1,k)
          v(i,nj+4,k)=v(i,nj+1,k)
        enddo
      elseif(nbc.eq.3)then
        do i=0,ni+1
          v(i,nj+1,k)=0.
          v(i,nj+2,k)=-v(i,nj  ,k)
          v(i,nj+3,k)=-v(i,nj-1,k)
          v(i,nj+4,k)=-v(i,nj-2,k)
        enddo
      endif

!-----------------------------------

      if(wbc.eq.1)then
        do j=jb,je+1
          v(   0,j,k)=v(ni  ,j,k)
          v(  -1,j,k)=v(ni-1,j,k)
          v(  -2,j,k)=v(ni-2,j,k)
        enddo
      elseif(wbc.eq.2)then
        do j=0,nj+2
          v(   0,j,k)=v( 1,j,k)
          v(  -1,j,k)=v( 1,j,k)
          v(  -2,j,k)=v( 1,j,k)
        enddo
      elseif(wbc.eq.3)then
        do j=0,nj+2
          v(   0,j,k)=v( 1,j,k)
          v(  -1,j,k)=v( 2,j,k)
          v(  -2,j,k)=v( 3,j,k)
        enddo
      endif

!-----------------------------------

      if(ebc.eq.1)then
        do j=jb,je+1
          v(ni+1,j,k)=v(   1,j,k)
          v(ni+2,j,k)=v(   2,j,k)
          v(ni+3,j,k)=v(   3,j,k)
        enddo
      elseif(ebc.eq.2)then
        do j=0,nj+2
          v(ni+1,j,k)=v(ni,j,k)
          v(ni+2,j,k)=v(ni,j,k)
          v(ni+3,j,k)=v(ni,j,k)
        enddo
      elseif(ebc.eq.3)then
        do j=0,nj+2
          v(ni+1,j,k)=v(ni  ,j,k)
          v(ni+2,j,k)=v(ni-1,j,k)
          v(ni+3,j,k)=v(ni-2,j,k)
        enddo
      endif

!-----------------------------------

    ENDDO

    return
    end subroutine bcv

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

    subroutine bcw(ni,nj,nk,wbc,ebc,sbc,nbc,w)
    implicit none

    integer :: ni,nj,nk,wbc,ebc,sbc,nbc
    real, dimension(-2:ni+3,-2:nj+3,nk+1) :: w

    integer :: i,j,k,ib,ie,jb,je

    ib = -2
    ie = ni+3
    jb = -2
    je = nj+3

    DO k=2,nk

!-----------------------------------

      if(wbc.eq.1)then
        do j=jb,je
          w(   0,j,k)=w(ni  ,j,k)
          w(  -1,j,k)=w(ni-1,j,k)
          w(  -2,j,k)=w(ni-2,j,k)
        enddo
      elseif(wbc.eq.2)then
        do j=0,nj+1
          w(   0,j,k)=w( 1,j,k)
          w(  -1,j,k)=w( 1,j,k)
          w(  -2,j,k)=w( 1,j,k)
        enddo
      elseif(wbc.eq.3)then
        do j=0,nj+1
          w(   0,j,k)=w( 1,j,k)
          w(  -1,j,k)=w( 2,j,k)
          w(  -2,j,k)=w( 3,j,k)
        enddo
      endif

!-----------------------------------

      if(ebc.eq.1)then
        do j=jb,je
          w(ni+1,j,k)=w(   1,j,k)
          w(ni+2,j,k)=w(   2,j,k)
          w(ni+3,j,k)=w(   3,j,k)
        enddo
      elseif(ebc.eq.2)then
        do j=0,nj+1
          w(ni+1,j,k)=w(ni,j,k)
          w(ni+2,j,k)=w(ni,j,k)
          w(ni+3,j,k)=w(ni,j,k)
        enddo
      elseif(ebc.eq.3)then
        do j=0,nj+1
          w(ni+1,j,k)=w(ni  ,j,k)
          w(ni+2,j,k)=w(ni-1,j,k)
          w(ni+3,j,k)=w(ni-2,j,k)
        enddo
      endif

!-----------------------------------

      if(sbc.eq.1)then
        do i=ib,ie
          w(i,   0,k)=w(i,nj  ,k)
          w(i,  -1,k)=w(i,nj-1,k)
          w(i,  -2,k)=w(i,nj-2,k)
        enddo
      elseif(sbc.eq.2)then
        do i=0,ni+1
          w(i,   0,k)=w(i, 1,k)
          w(i,  -1,k)=w(i, 1,k)
          w(i,  -2,k)=w(i, 1,k)
        enddo
      elseif(sbc.eq.3)then
        do i=0,ni+1
          w(i,   0,k)=w(i, 1,k)
          w(i,  -1,k)=w(i, 2,k)
          w(i,  -2,k)=w(i, 3,k)
        enddo
      endif

!-----------------------------------

      if(nbc.eq.1)then
        do i=ib,ie
          w(i,nj+1,k)=w(i,   1,k)
          w(i,nj+2,k)=w(i,   2,k)
          w(i,nj+3,k)=w(i,   3,k)
        enddo
      elseif(nbc.eq.2)then
        do i=0,ni+1
          w(i,nj+1,k)=w(i,nj,k)
          w(i,nj+2,k)=w(i,nj,k)
          w(i,nj+3,k)=w(i,nj,k)
        enddo
      elseif(nbc.eq.3)then
        do i=0,ni+1
          w(i,nj+1,k)=w(i,nj  ,k)
          w(i,nj+2,k)=w(i,nj-1,k)
          w(i,nj+3,k)=w(i,nj-2,k)
        enddo
      endif

!-----------------------------------

    ENDDO

    do j=-2,nj+3
    do i=-2,ni+3
      w(i,j,   1) = 0.0
      w(i,j,nk+1) = 0.0
    enddo
    enddo

!-----------------------------------

    return
    end subroutine bcw

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
