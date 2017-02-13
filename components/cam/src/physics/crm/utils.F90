integer function lenstr (string)
      
! returns string's length ignoring the rightmost blank and null characters

implicit none
character *(*) string
integer k
lenstr = 0
do k = 1,len(string)
 if (string(k:k).ne.' '.and.string(k:k).ne.char(0)) then
     lenstr = lenstr+1
 end if
end do
111  return
end



subroutine averageXY(f,dimx1,dimx2,dimy1,dimy2,dimz,fm)
	
use grid
implicit none
integer dimx1, dimx2, dimy1, dimy2, dimz
real f(dimx1:dimx2, dimy1:dimy2, dimz),fm(nzm) 
real(8) ff,factor
integer i,j,k	
factor = 1./dble(nx*ny)
do k =1,nzm
 ff = 0.
 do j =1,ny
  do i =1,nx
    ff = ff + f(i,j,k)
  end do
 end do
 ff = ff*factor
 fm(k) = real(ff)
end do 
end


subroutine averageXY_MPI(f,dimx1,dimx2,dimy1,dimy2,dimz,fm)
	
use grid
implicit none
integer dimx1, dimx2, dimy1, dimy2, dimz
real f(dimx1:dimx2, dimy1:dimy2, dimz),fm(nzm)
real(8) fm1(nzm),fm2(nzm),factor
integer i,j,k
factor = 1./dble(nx*ny)
do k =1,nzm
 fm1(k) = 0.
 do j =1,ny
  do i =1,nx
    fm1(k) = fm1(k) + f(i,j,k)
  end do
 end do
 fm1(k) = fm1(k) * factor
end do
if(dompi) then
 do k =1,nzm
   fm2(k) = fm1(k)
 end do
 call task_sum_real8(fm2,fm1,nzm)
 do k=1,nzm
  fm(k)=real(fm1(k)/dble(nsubdomains))
 end do
else
 do k=1,nzm
  fm(k)=real(fm1(k))
 end do
endif
end
	
		
	
	
subroutine fminmax_print(name,f,dimx1,dimx2,dimy1,dimy2,dimz)

use grid
implicit none
integer dimx1, dimx2, dimy1, dimy2, dimz
real f(dimx1:dimx2, dimy1:dimy2, dimz),fmn(nz),fmx(nz)
character *(*) name
real fmin(1),fmax(1),fff(1)
integer i,j,k
	
do k=1,dimz
 if(dimx2.eq.1.and.dimy2.eq.1) then
   fmn(k) = f(1,1,k)	
   fmx(k) = f(1,1,k)	
 else
   fmn(k) = 1.e30
   fmx(k) =-1.e30
!   do j=1,ny
!    do i=1,nx
   do j=dimy1,dimy2
    do i=dimx1,dimx2
     fmn(k) = min(fmn(k),f(i,j,k))
     fmx(k) = max(fmx(k),f(i,j,k))
    end do
   enddo
 end if
enddo
fmin(1) = 1.e30
fmax(1) =-1.e30
do k=1,dimz
 fmin(1) = min(fmin(1),fmn(k))
 fmax(1) = max(fmax(1),fmx(k))
end do
	
if(dompi) then
  fff(1)=fmax(1)
  call task_max_real(fff(1),fmax(1),1)
  fff(1)=fmin(1)
  call task_min_real(fff(1),fmin(1),1)
end if
if(masterproc) print *,name,fmin,fmax
end



	
subroutine setvalue(f,n,f0)
implicit none
integer n
real f(n), f0
integer k
do k=1,n
 f(k)=f0
end do
end
