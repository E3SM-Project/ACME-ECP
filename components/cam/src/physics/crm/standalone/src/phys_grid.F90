

module phys_grid
  implicit none

contains

  real(8) function get_rlon_p(lcid, col)
    integer, intent(in) :: lcid
    integer, intent(in) :: col
    !integer :: cid
    !cid = lchunks(lcid)%cid
    !get_rlon_p = clon_p(chunks(cid)%lon(col))

    !TODO: Need to store the correct value in the NetCDF file and use that here instead of "1."
    get_rlon_p = 1.
  end function get_rlon_p

  real(8) function get_rlat_p(lcid, col)
    integer, intent(in)  :: lcid          ! local chunk id
    integer, intent(in)  :: col           ! column index
    !integer :: cid                        ! global chunk id
    !cid = lchunks(lcid)%cid
    !get_rlat_p = clat_p(chunks(cid)%lat(col))

    !TODO: Need to store the correct value in the NetCDF file and use that here instead of "1."
    get_rlat_p = 1.
  end function get_rlat_p

  subroutine get_gcol_all_p(lcid, latdim, gcols)
    integer, intent(in)  :: lcid        ! local chunk id
    integer, intent(in)  :: latdim      ! declared size of output array
    integer, intent(out) :: gcols(:)    ! array of global latitude indices
    integer :: i                        ! loop index
    !gcols=-1
    !do i=1,lchunks(lcid)%ncols
    !  gcols(i) = lchunks(lcid)%gcol(i)
    !enddo

    !TODO: Need to store the correct value in the NetCDF file and use that here instead of "-1"
    gcols = -1
  end subroutine get_gcol_all_p


end module phys_grid
