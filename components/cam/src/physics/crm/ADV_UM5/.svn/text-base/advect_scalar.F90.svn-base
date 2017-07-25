subroutine advect_scalar( f, fadv, flux, f2leadv, f2legrad, fwleadv, doit )
 	
!	5th order ultimate-macho advection scheme
! Yamaguchi, T., D. A. Randall, and M. F. Khairoutdinov, 2011:
! Cloud Modeling Tests of the ULTIMATE-MACHO Scalar Advection Scheme.
! Monthly Weather Review. 139, pp.3248-3264


!	At this point, 
!	u = u * rho * dtn / dx
!	v = v * rho * dtn / dy
!	w = w * rho * dtn / dz
!	Division by adz has not been performed.

	use grid
	use vars, only: u, v, w, rho, rhow
        use params, only: docolumn
	
	implicit none
	
	!	input
	real, dimension(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm), intent(inout) :: f
	real, dimension(nz), intent(out) :: flux, fadv
	real, dimension(nz), intent(out) :: f2leadv, f2legrad, fwleadv
	logical, intent(in) :: doit
	
	!	local
	real, dimension(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) :: df
	integer :: i, j, k
	
	
	if(docolumn) then
		flux = 0.
		return
	endif
	
	!call t_startf ('advect_scalars')
	
    df(:,:,:) = f(:,:,:)
	
	if (RUN3D) then
		call advect_scalar3D(f, u, v, w, rho, rhow, flux)
	else
		call advect_scalar2D(f, u, w, rho, rhow, flux)	  
	endif
	
	do k = 1, nzm
			fadv(k) = 0.
			do j = 1, ny
				do i = 1, nx
					fadv(k) = fadv(k) + f(i,j,k) - df(i,j,k)
				enddo
			enddo
	enddo
		
	!call t_stopf ('advect_scalars')

end subroutine advect_scalar
