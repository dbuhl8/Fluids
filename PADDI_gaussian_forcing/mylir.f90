!program test
!implicit none

!integer, parameter :: nn = 105
!double precision, dimension(nn) :: x,y
!double precision :: z,z2
!integer :: i,inter

!do i=1,nn
!   x(i) = (i-1.d0)/(nn-1.d0)
!   y(i) = x(i)**2
!enddo

!z = 2.0

!call mylir(z,x,z2,y,1,nn,inter)
!write(*,*) inter,z,z2-z**2

!end program test

subroutine mylir(z,zi,y,yi,ii,nt,inter)
  implicit none

!
!
!  Linear interpolation/extrapolation routine
!
!  Input arrays: zi(nt) is input independent variable
!                yi(ii,nt) is array of input dependent variables
!                z is independent variable of new point
!  Output:       y(ii) is array of depn. vars at new point.
!     
!     inter is set to 1 for interpolation and 0 for extrapolation
! 
!
 
  double precision, dimension(nt) :: zi
  double precision, dimension(ii) :: y
  double precision, dimension(ii,nt) :: yi
  double precision :: z,z1,z2,y1,y2,diff1,diff,diff2
  integer :: ii,nt,i,inter,n,indx

!     check nt and reset il if necessary
  if(nt.lt.2) then
     write(*,*) 'Not enough datapoints'
     inter = 0
     goto 200
  endif

!     addressing constants
  inter = 1

!     determine position of z within zi
  do i=1,nt-1
     diff = zi(i+1) - zi(i) 
     if(diff.eq.0.d0) goto 199
     diff1 = z-zi(i)
     if (diff1.eq.0.d0) then
        indx = i
        goto 10
     endif
     diff2 = z-zi(i+1)
     if (diff2.eq.0.d0) then
        indx = i + 1 
        goto 10
     endif
     if(diff1*diff2.lt.0) goto 20
  enddo

! If point scanned through mesh without finding anything
  inter = 2
  diff = zi(nt)-zi(1)
  if(diff.gt.0.d0) then  ! case of increasing mesh
     if(z.gt.zi(nt)) then
        y1 = (z-zi(nt))/(zi(nt)-zi(nt-1))
        do i=1,ii
           y(i) = yi(i,nt) + y1*(yi(i,nt)-yi(i,nt-1))
        enddo
     else 
        if(z.lt.zi(1)) then
           y1 = (z-zi(1))/(zi(1)-zi(2))
           do i=1,ii
              y(i) = yi(i,1) + y1*(yi(i,1)-yi(i,2))
           enddo
        else
           write(*,*) 'Problem'
        endif
     endif
  else
     if(z.gt.zi(1)) then
        y1 = (z-zi(1))/(zi(1)-zi(2))
        do i=1,ii
           y(i) = yi(i,1) + y1*(yi(i,1)-yi(i,2))
        enddo
     else 
        if(z.lt.zi(nt)) then
           y1 = (z-zi(nt))/(zi(nt)-zi(nt-1))
           do i=1,ii
              y(i) = yi(i,nt) + y1*(yi(i,nt)-yi(i,nt-1))
           enddo
        else
           write(*,*) 'Problem'
        endif
     endif
  endif
  goto 200
          
! set y when z lies on a mesh point
10 inter = 3  
  do i=1,ii
     y(i)=yi(i,indx)
  enddo
  goto 200

! linear interpolation if z lies between meshpoints
20  n = i+1
  z1=zi(n)
  y1=(z1-z)/(z1-zi(n-1))
  y2=1.0-y1
  do i=1,ii
     y(i)=y1*yi(i,n-1)+y2*yi(i,n)
  enddo
  goto 200

! linear extrapolation if z is outside mesh

199 write(*,*) 'Mesh non-increasing at meshpoint', i

200 continue
 
end subroutine mylir
