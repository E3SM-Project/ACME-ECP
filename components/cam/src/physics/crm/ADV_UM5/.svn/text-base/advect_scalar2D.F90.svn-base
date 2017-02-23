subroutine advect_scalar2D( f, u, w, rho, rhow, flux )
 	
!	Two dimentional 5th order ULTIMATE-MACHO scheme

	use grid
	use advect_um_lib
        use params, only: dowallx
	implicit none

	!	input & output
	real, dimension(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm), intent(inout) :: f
	real, dimension(dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm), intent(inout) :: u
	real, dimension(dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz ), intent(in) :: w
	real, dimension(nzm), intent(in) :: rho
	real, dimension(nz), intent(in) :: rhow
	real, dimension(nz), intent(out) :: flux

	!	local
	integer, parameter :: j = 1
	integer :: macho_order, i, k
	
	!--------------------------------------------------------------------------
	if (dowallx) then
		if ( mod(rank,nsubdomains_x) == 0 ) then
			do k = 1, nzm
				do i = dimx1_u, 1
					u(i,j,k) = 0.
				enddo
			enddo
		endif
		if ( mod(rank,nsubdomains_x) == nsubdomains_x-1 ) then
			do k = 1, nzm
				do i = nx+1, dimx2_u
					u(i,j,k) = 0.
				enddo
			enddo
		endif
	endif
	!--------------------------------------------------------------------------

	!	Convert mass-weighted courant number to non-mass weighted
	!	Inverse of rho and adz
	if ( ( nstep > nstep_adv ).and.( .not.updated_cn(icycle) ) ) then
		!!if (masterproc) print*,'cn updated'
		updated_cn(icycle) = .true.	! skip for same icycle if updated
		if (icycle == ncycle) then	! skip at ncycle if updated
			nstep_adv = nstep
			updated_cn(:) = .false.
		endif

		!	Inverse of rho and adz, adzw
		do k = 1, nzm
			irho(k)  = 1. / rho(k)
			iadz(k)  = 1. / adz(k)
			iadzw(k) = 1. / adzw(k)
		enddo
	
		!	x direction
		do k = 1, nzm
			do i = -1, nxp3	
				cu(i,j,k) = u(i,j,k) * irho(k)
			enddo
		enddo
	
		!	z direction
		cw(:,:,nz) = 0.
		cw(:,:,1) = 0.	
		do k = 2, nzm
			irhow(k) = 1. / ( rhow(k) * adz(k) )
			do i = -3, nxp4
				cw(i,j,k) = w(i,j,k) * irhow(k)		
			enddo
		enddo
	endif

	!	Top and bottom boundaryies
	fz(:,:,nz) = 0.
	fz(:,:,1) = 0.

	!	Face values
	fadv(:,:,:) = f(:,:,:)
	macho_order = mod(nstep,2)
	select case (macho_order)
	case(0)
		
		!	x-direction
		call face_x_5th( 0, nxp2, 1, 1 )
		call adv_form_update_x( 0, nxp1, 1, 1 )
				
		!	z-direction
		call face_z_5th( 0, nxp1, 1, 1 )
	
	case(1)

		!	z-direction
		call face_z_5th( -3, nxp4, 1, 1 )
		call adv_form_update_z( -3, nxp4, 1, 1 )
		
		!	x-direction
		call face_x_5th( 0, nxp2, 1, 1 )
	
	end select

	!	FCT to ensure positive definite or monotone
	if (fct) then
		call fct2D( f, u, w, flux )
	else
		!	In case...
		!fz(:,:,nz) = 0.
		!fz(:,:,1) = 0.
	
		!	Flux-form update
		flux = 0.
		do k = 1, nzm
			do i = 1, nx
				f(i,j,k) = f(i,j,k) &
					+ ( u(i,j,k) * fx(i,j,k) - u(i+1,j,k) * fx(i+1,j,k) &
					+ ( w(i,j,k) * fz(i,j,k) - w(i,j,k+1) * fz(i,j,k+1) ) * iadz(k) ) * irho(k)
				flux(k) = flux(k) + w(i,j,k) * fz(i,j,k)
			enddo
		enddo
	endif

end subroutine advect_scalar2D
