module sw_core
!========================================================================
! This module contains the routine that computes the time averaged winds
!
! reference https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere/blob/main/model/sw_core.F90
!========================================================================
use fv_arrays, only: R_GRID, fv_grid_bounds_type, fv_grid_type
implicit none
contains


subroutine time_averaged_cfl(gridstruct, bd, outer_crx, outer_cry, inner_crx, inner_cry, &
                             uc_old, vc_old, uc, vc, inner_dp, outer_dp, dt)
    !--------------------------------------------------
    ! Compute the time average CFL needed
    ! for the departure point scheme
    !
    !--------------------------------------------------
    type(fv_grid_bounds_type), intent(IN) :: bd
    type(fv_grid_type), intent(IN), target :: gridstruct

    real(R_GRID), intent(INOUT), dimension(bd%is:bd%ie+1  , bd%jsd:bd%jed  ) :: outer_crx
    real(R_GRID), intent(INOUT), dimension(bd%isd:bd%ied  , bd%js:bd%je+1  ) :: outer_cry

    real(R_GRID), intent(INOUT), dimension(bd%is:bd%ie+1  , bd%jsd:bd%jed  ) :: inner_crx
    real(R_GRID), intent(INOUT), dimension(bd%isd:bd%ied  , bd%js:bd%je+1  ) :: inner_cry

    real(R_GRID), intent(IN)   , dimension(bd%isd:bd%ied+1, bd%jsd:bd%jed  ) :: uc_old
    real(R_GRID), intent(IN)   , dimension(bd%isd:bd%ied  , bd%jsd:bd%jed+1) :: vc_old

    real(R_GRID), intent(IN)   , dimension(bd%isd:bd%ied+1, bd%jsd:bd%jed  ) :: uc
    real(R_GRID), intent(IN)   , dimension(bd%isd:bd%ied  , bd%jsd:bd%jed+1) :: vc

    real(R_GRID), intent(IN) :: dt

    integer, intent(IN):: inner_dp ! inner departute point method !1 - Euler; 2-RK2
    integer, intent(IN):: outer_dp ! outer departute point method !1 - Euler; 2-RK2

    ! aux
    integer :: i, j
    integer :: is,  ie,  js,  je

    is  = bd%is
    ie  = bd%ie
    js  = bd%js
    je  = bd%je

    call departure_cfl(gridstruct, bd, inner_crx, inner_cry, uc_old, vc_old, uc, vc, inner_dp, dt)

    if (inner_dp==outer_dp) then
       outer_crx(is:ie+1,:) = inner_crx(is:ie+1,:)
       outer_cry(:,js:je+1) = inner_cry(:,js:je+1)
    else
       call departure_cfl(gridstruct, bd, outer_crx, outer_cry, uc_old, vc_old, uc, vc, outer_dp, dt)
    endif

end subroutine time_averaged_cfl

subroutine departure_cfl(gridstruct, bd, crx, cry, &
                             uc_old, vc_old, uc, vc, dp, dt)
    !--------------------------------------------------
    ! Compute the departure CFL for a given 
    ! departure point scheme (dp)
    !--------------------------------------------------
    type(fv_grid_bounds_type), intent(IN) :: bd
    type(fv_grid_type), intent(IN), target :: gridstruct 

    real(R_GRID), intent(INOUT), dimension(bd%is:bd%ie+1, bd%jsd:bd%jed  ) :: crx
    real(R_GRID), intent(INOUT), dimension(bd%isd:bd%ied, bd%js:bd%je+1  ) :: cry

    real(R_GRID), intent(IN)   , dimension(bd%isd:bd%ied+1, bd%jsd:bd%jed  ) :: uc_old
    real(R_GRID), intent(IN)   , dimension(bd%isd:bd%ied  , bd%jsd:bd%jed+1) :: vc_old

    real(R_GRID), intent(IN)   , dimension(bd%isd:bd%ied+1, bd%jsd:bd%jed  ) :: uc
    real(R_GRID), intent(IN)   , dimension(bd%isd:bd%ied  , bd%jsd:bd%jed+1) :: vc

    real(R_GRID), dimension(bd%is-1:bd%ie+2  , bd%jsd:bd%jed  ) :: crx_time_centered
    real(R_GRID), dimension(bd%isd:bd%ied  , bd%js-1:bd%je+2  ) :: cry_time_centered
    real(R_GRID), dimension(bd%is:bd%ie+1  , bd%jsd:bd%jed  ) :: crx_old
    real(R_GRID), dimension(bd%isd:bd%ied  , bd%js:bd%je+1  ) :: cry_old


    real(R_GRID), intent(IN) :: dt

    integer, intent(IN):: dp ! departute point method !1 - Euler; 2-RK2

    ! aux
    real(R_GRID) :: a, a1, a2, c1, c2
    integer :: i, j
    integer :: is,  ie,  js,  je
    integer :: isd, ied, jsd, jed

    is  = bd%is
    ie  = bd%ie
    js  = bd%js
    je  = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed

    select case (dp)
      case (1)
         ! CFL for RK1
         call  compute_cfl(gridstruct, bd, crx, cry, uc, vc, dt, 0)

      case (2)
         ! CFL for RK2
         call  compute_cfl(gridstruct, bd, crx_old          , cry_old          , uc_old, vc_old, dt, 0)
         call  compute_cfl(gridstruct, bd, crx_time_centered, cry_time_centered, uc    , vc    , dt, 1)

         ! RK2
         ! cfl for dp in x direction
         do j = jsd, jed
             do i = is, ie+1
                ! Linear interpolation weight
                a = crx_old(i,j)*0.5d0
                ! Upwind linear interpolation
                if (a>0.d0) then
                   c1 = crx_time_centered(i-1,j)
                   c2 = crx_time_centered(i,j)
                   a1 = a
                   a2 = 1.d0-a
                else
                   c1 = crx_time_centered(i,j)
                   c2 = crx_time_centered(i+1,j)
                   a1 = 1.d0+a
                   a2 = -a
                end if
                crx(i,j) = a1*c1 + a2*c2
             end do
          end do

          ! cfl for dp in y direction
          do i = isd, ied
             do j = js, je+1
                ! Linear interpolation weight
                a = cry_old(i,j)*0.5d0
                ! Upwind linear interpolation
                if (a>0.d0) then
                   c1 = cry_time_centered(i,j-1)
                   c2 = cry_time_centered(i,j)
                   a1 = a
                   a2 = 1.d0-a
                else
                   c1 = cry_time_centered(i,j)
                   c2 = cry_time_centered(i,j+1)
                   a1 = 1.d0+a
                   a2 = -a
                end if
                cry(i,j) = a1*c1 + a2*c2
             end do
          end do

    end select

end subroutine departure_cfl

subroutine compute_cfl(gridstruct, bd, crx, cry, uc, vc, dt, h)
   !--------------------------------------------------
   ! Compute CFL in x and y directions at C grid
   !--------------------------------------------------
   type(fv_grid_bounds_type), intent(IN) :: bd
   type(fv_grid_type), intent(IN), target :: gridstruct 
   real(R_GRID), intent(INOUT), dimension(bd%is-h:bd%ie+1+h  , bd%jsd:bd%jed  ) :: crx
   real(R_GRID), intent(INOUT), dimension(bd%isd:bd%ied  , bd%js-h:bd%je+1+h  ) :: cry
   real(R_GRID), intent(IN)   , dimension(bd%isd:bd%ied+1, bd%jsd:bd%jed  ) :: uc
   real(R_GRID), intent(IN)   , dimension(bd%isd:bd%ied  , bd%jsd:bd%jed+1) :: vc
   real(R_GRID), intent(IN) :: dt
   integer, intent(IN):: h

   ! aux
   integer :: i, j
   integer :: is,  ie,  js,  je
   integer :: isd, ied, jsd, jed
   real(R_GRID) :: dx, dy

   is  = bd%is
   ie  = bd%ie
   js  = bd%js
   je  = bd%je
   isd = bd%isd
   ied = bd%ied
   jsd = bd%jsd
   jed = bd%jed

   dx = gridstruct%dx
   dy = gridstruct%dy
 
   ! Compute CFL at timestep n
   do j=jsd,jed
      do i=is-h,ie+h+1
        crx(i,j) = uc(i,j)*dt/dx
     enddo
   enddo

   do j=js-h,je+h+1
      do i=isd,ied
        cry(i,j) = vc(i,j)*dt/dy
     enddo
   enddo
end subroutine compute_cfl

end module sw_core 
