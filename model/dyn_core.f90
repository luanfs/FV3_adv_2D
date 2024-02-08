module dyn_core
!========================================================================
! This module contains the routine that computes the PPM flux
!
! reference https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere/blob/main/model/dyn_core.F90
!========================================================================
use fv_arrays,  only: fv_grid_bounds_type, fv_grid_type, R_GRID
use test_cases, only: calc_winds
use tp_core,    only: fv_tp_2d
use sw_core,    only: time_averaged_cfl
use fv_duogrid, only: ext_scalar_agrid, ext_scalar_cgrid

implicit none

contains
subroutine dy_core(qa, uc, uc_old, vc_old, vc, bd, gridstruct, time, time_centered, dt, dto2, test_case, hord, lim_fac, dp)
   type(fv_grid_bounds_type), intent(INOUT) :: bd
   type(fv_grid_type), target, intent(INOUT) :: gridstruct
   real(R_GRID), intent(in) :: time, dt, dto2
   real(R_GRID), intent(in) :: time_centered
   real(R_GRID), intent(inout) :: lim_fac
   real(R_GRID), intent(inout) :: qa(bd%isd:bd%ied, bd%jsd:bd%ied)
   real(R_GRID), intent(inout) :: uc    (bd%isd:bd%ied  , bd%jsd:bd%jed+1)
   real(R_GRID), intent(inout) :: uc_old(bd%isd:bd%ied  , bd%jsd:bd%jed+1)
   real(R_GRID), intent(inout) :: vc    (bd%isd:bd%ied+1, bd%jsd:bd%jed  )
   real(R_GRID), intent(inout) :: vc_old(bd%isd:bd%ied+1, bd%jsd:bd%jed  )
 
   integer, intent(IN) :: test_case
   integer, intent(IN) :: hord
   integer, intent(IN) :: dp

   real(R_GRID) :: dx
   real(R_GRID) :: dy

   real(R_GRID) :: inner_xfx(bd%is:bd%ie+1, bd%jsd:bd%jed)
   real(R_GRID) :: outer_xfx(bd%is:bd%ie+1, bd%jsd:bd%jed)
                        
   real(R_GRID) :: inner_crx(bd%is:bd%ie+1, bd%jsd:bd%jed)
   real(R_GRID) :: outer_crx(bd%is:bd%ie+1, bd%jsd:bd%jed)
                        
   real(R_GRID) :: inner_yfx(bd%isd:bd%ied, bd%js:bd%je+1)
   real(R_GRID) :: outer_yfx(bd%isd:bd%ied, bd%js:bd%je+1)
                        
   real(R_GRID) :: inner_cry(bd%isd:bd%ied, bd%js:bd%je+1)
   real(R_GRID) :: outer_cry(bd%isd:bd%ied, bd%js:bd%je+1)
 
   real(R_GRID) :: flux_x(bd%is:bd%ie+1, bd%js:bd%je)
   real(R_GRID) :: flux_y(bd%is:bd%ie  , bd%js:bd%je+1)

   integer :: is, ie, isd, ied, ng
   integer :: js, je, jsd, jed
   integer :: i, j
   integer :: inner_dp, outer_dp

   inner_dp = dp
   outer_dp = dp

   is  = bd%is
   ie  = bd%ie
   js  = bd%js
   je  = bd%je
 
   isd = bd%isd
   ied = bd%ied

   jsd = bd%jsd
   jed = bd%jed
 
   ng  = bd%ng
   dx  = gridstruct%dx
   dy  = gridstruct%dy

   ! winds
   call calc_winds(uc_old, vc_old, bd, gridstruct, time         , test_case)
   call calc_winds(uc    , vc    , bd, gridstruct, time_centered, test_case)

   ! periodic BC
   call ext_scalar_agrid(qa, bd)
   call ext_scalar_cgrid(uc    , vc    , bd)

   ! compute time averaged cfl
   call time_averaged_cfl(gridstruct, bd, outer_crx, outer_cry, inner_crx, inner_cry, &
                             uc_old, vc_old, uc, vc, inner_dp, outer_dp, dt)
 
   ! compute adv coeffs
   inner_xfx(is:ie+1,jsd:jed) = inner_crx(is:ie+1,jsd:jed)
   outer_xfx(is:ie+1,jsd:jed) = outer_crx(is:ie+1,jsd:jed)

   inner_yfx(isd:ied,js:je+1) = inner_cry(isd:ied,js:je+1)
   outer_yfx(isd:ied,js:je+1) = outer_cry(isd:ied,js:je+1)

   call fv_tp_2d(qa, outer_crx, outer_cry, inner_crx, inner_cry, hord, flux_x, flux_y, &
                     outer_xfx, outer_yfx, inner_xfx, inner_yfx, gridstruct, bd, lim_fac)
   ! update scalar field
   !qa(is:ie,js:jed) = qa(is:ie,js:je) - (flux_x(is+1:ie+1,js:je)-flux_x(is:ie,js:je))
   qa(is:ie,js:jed) = qa(is:ie,js:je) - (flux_y(is:ie,js+1:je+1)-flux_y(is:ie,js:je))

end subroutine dy_core

end module dyn_core 
