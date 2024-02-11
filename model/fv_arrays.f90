module fv_arrays
!========================================================================
! This module contains the arrays data structure
!
! Reference:
! https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere/blob/main/model/fv_arrays.F90
!========================================================================
public
integer, public, parameter :: R_GRID = selected_real_kind(12,100)
real(R_GRID), public, parameter :: pio4  = datan(1.d0)
real(R_GRID), public, parameter :: pi    = 4.d0*pio4
real(R_GRID), public, parameter :: pio2  = pi*0.5d0
real(R_GRID), public, parameter :: twopi = pi*2.d0
real(R_GRID), public, parameter :: erad  = 6.37122e6
real(R_GRID), public, parameter :: eradi = 1.d0/erad
real(R_GRID), public, parameter :: day2sec = 86400.d0
real(R_GRID), public, parameter :: sec2day = 1.d0/day2sec
character(len=32) :: datadir = "data/"
character(len=32) :: pardir = "par/"

!-----------------------------------------------------------------------------
! domain bounds
!-----------------------------------------------------------------------------
type fv_grid_bounds_type
   integer :: is,  ie   ! x direction interior cell indexes (for agrid; interior cgrid ranges from is to ie+1)
   integer :: js,  je   ! y direction interior cell indexes (for agrid; interior dgrid ranges from js to je+1)
   integer :: isd, ied  ! x direction indexes including ghost cells (for agrid; cgrid ranges from isd to ied+1)
   integer :: jsd, jed  ! y direction indexes including ghost cells (for agrid; dgrid ranges from jsd to jed+1)
 
   integer :: ng  = 3   ! ghost cell layers
   integer :: npx       ! number of interior cells in x direction
   integer :: npy       ! number of interior cells in y direction
end type fv_grid_bounds_type

!-------------------------------------------------
! point structure
!-------------------------------------------------
type point_structure
   ! Cartesian coordinates (R^2)
   !real(kind=R_GRID), dimension(1:2) :: p

   ! Local coordinate
   real(kind=R_GRID) :: x, y

   ! Polar coordinates in radians ([-pi, pi[)
   !real(kind=R_GRID) :: theta
end type point_structure

!-------------------------------------------------
! grid structure
!-------------------------------------------------
type fv_grid_type
   type(point_structure), allocatable, dimension(:,:) :: agrid
   type(point_structure), allocatable, dimension(:,:) :: bgrid
   type(point_structure), allocatable, dimension(:,:) :: cgrid
   type(point_structure), allocatable, dimension(:,:) :: dgrid
   real(R_GRID), allocatable ::  area(:,:)
   real(R_GRID), allocatable :: rarea(:,:)
   real(R_GRID) :: dx, dy
   integer :: npx       ! number of interior cells (x direction)
   integer :: npy       ! number of interior cells (y direction)
end type fv_grid_type

!-------------------------------------------------
! atmosphere structure
!-------------------------------------------------
type fv_atmos_type
   character(len=128) :: simulation_name
   type(fv_grid_type) :: gridstruct
   type(fv_grid_bounds_type) :: bd
   integer  :: test_case
   integer  :: hord
   integer  :: adv_scheme
   integer  :: inner_dp
   integer  :: outer_dp
   integer  :: inner_adv
   integer  :: nplot = 0
   integer  :: nplots
   integer  :: plotstep
   real(R_GRID) :: lim_fac = 1.d0

   ! time vars
   real(R_GRID) :: dt            ! time step
   real(R_GRID) :: dto2          ! = (dt/2)
   real(R_GRID) :: time          ! current time
   real(R_GRID) :: time_centered ! current time + dto2
   real(R_GRID) :: Tf            ! final time
   integer  :: total_tsteps           ! total of time steps

   ! Fields
   real(R_GRID), allocatable :: qa (:,:) ! A grid scalar field
   real(R_GRID), allocatable :: qa0(:,:) ! A grid scalar field (IC)

   real(R_GRID), allocatable :: uc (:,:) ! (uc) are mostly used as the C grid winds
   real(R_GRID), allocatable :: uc0(:,:) ! (uc) are mostly used as the C grid winds (IC)
   real(R_GRID), allocatable :: uc_old(:,:) ! 

   real(R_GRID), allocatable :: vc (:,:) ! (vc) are mostly used as the C grid winds
   real(R_GRID), allocatable :: vc0(:,:) ! (vc) are mostly used as the C grid winds (IC)
   real(R_GRID), allocatable :: vc_old(:,:) ! 


   ! erros vars
   real(R_GRID), allocatable :: error_qa(:,:)
   real(R_GRID) :: linf_error_qa, l1_error_qa, l2_error_qa

   ! cfl
   real(R_GRID) :: cfl, cfl_x, cfl_y

   ! mass vars
   real(R_GRID) :: mass_qa0, mass_qa, mass_qa_var

   integer :: npx       ! number of interior cells
   integer :: npy       ! number of interior cells
end type fv_atmos_type

end module fv_arrays
