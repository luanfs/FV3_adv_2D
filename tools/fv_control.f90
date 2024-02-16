module fv_control
!========================================================================
! Module for data allocation
!
! Reference
! https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere/blob/main/tools/fv_control.F90
!========================================================================
use fv_arrays,  only: fv_grid_bounds_type, fv_grid_type, fv_atmos_type, &
                      point_structure, R_GRID, erad, pi, pio4, pio2
implicit none

contains

!-------------------------------------------------
! define bounds
!-------------------------------------------------
subroutine init_bounds(bd, npx, npy)
   type(fv_grid_bounds_type), intent(INOUT) :: bd
   integer                  , intent(IN)    :: npx, npy

   bd%is  = 1
   bd%ie  = npx
   bd%isd = bd%is - bd%ng
   bd%ied = bd%ie + bd%ng

   bd%js  = 1
   bd%je  = npy
   bd%jsd = bd%js - bd%ng
   bd%jed = bd%je + bd%ng
 
   bd%npx = npx
   bd%npy = npy
end subroutine init_bounds

!-------------------------------------------------
! allocate grid - bd must be filled
!-------------------------------------------------
subroutine init_grid(gridstruct, bd)
   type(fv_grid_type), target, intent(INOUT) :: gridstruct
   type(fv_grid_bounds_type), intent(IN)    :: bd
   type(point_structure), pointer, dimension(:,:) :: agrid
   type(point_structure), pointer, dimension(:,:) :: bgrid
   type(point_structure), pointer, dimension(:,:) :: cgrid
   type(point_structure), pointer, dimension(:,:) :: dgrid

   real(R_GRID), pointer, dimension(:,:) :: area
   real(R_GRID), pointer, dimension(:,:) :: rarea
   real(R_GRID), pointer :: dx, dy
   real(R_GRID):: L
   integer :: is, ie, isd, ied
   integer :: js, je, jsd, jed
   integer :: i, j

   is  = bd%is
   js  = bd%js
   isd = bd%isd
   jsd = bd%jsd
   ie  = bd%ie
   je  = bd%je
   ied = bd%ied
   jed = bd%jed

   gridstruct%npx = bd%npx
   gridstruct%npy = bd%npy

   allocate(gridstruct%agrid(isd:ied  , jsd:jed  ))
   allocate(gridstruct%bgrid(isd:ied+1, jsd:jed+1))
   allocate(gridstruct%cgrid(isd:ied+1, jsd:jed  ))
   allocate(gridstruct%dgrid(isd:ied  , jsd:jed+1))
   allocate(gridstruct% area(isd:ied  , jsd:jed  ))
   allocate(gridstruct%rarea(isd:ied  , jsd:jed  ))

   agrid => gridstruct%agrid
   bgrid => gridstruct%bgrid
   cgrid => gridstruct%cgrid
   dgrid => gridstruct%dgrid

   area  => gridstruct%area
   rarea => gridstruct%rarea
   dx    => gridstruct%dx
   dy    => gridstruct%dy

   L = pio2*erad
   dx = L/bd%npx
   dy = L/bd%npy
   !dx = 1.d0/bd%npx
   !dy = 1.d0/bd%npy
 
   area = dx*dy
   rarea = 1.d0/area
   ! compute c grid local coordinates
   do i = isd, ied+1
      do j = jsd, jed+1
         bgrid(i,j)%x = -L*0.5d0 + (i-1d0)*dx
         bgrid(i,j)%y = -L*0.5d0 + (j-1d0)*dy
      enddo
   enddo

   ! compute a grid local coordinates
   ! a grid
   agrid(isd:ied, jsd:jed)%x = (bgrid(isd+1:ied+1, jsd+1:jed+1)%x + bgrid(isd:ied    , jsd:jed)%x + &
                                bgrid(isd  :ied  , jsd+1:jed+1)%x + bgrid(isd+1:ied+1, jsd:jed)%x)* 0.25d0

   agrid(isd:ied, jsd:jed)%y = (bgrid(isd+1:ied+1, jsd+1:jed+1)%y + bgrid(isd:ied    , jsd:jed)%y + &
                                bgrid(isd  :ied  , jsd+1:jed+1)%y + bgrid(isd+1:ied+1, jsd:jed)%y)* 0.25d0

   cgrid(isd:ied+1, jsd:jed)%x = bgrid(isd:ied+1, jsd:jed)%x
   cgrid(isd:ied+1, jsd:jed)%y = (bgrid(isd:ied, jsd+1:jed+1)%y+bgrid(isd:ied, jsd:jed)%y)*0.5d0

   dgrid(isd:ied, jsd:jed)%x   = (bgrid(isd+1:ied+1, jsd:jed)%x + bgrid(isd:ied,jsd:jed)%x)*0.5d0
   dgrid(isd:ied, jsd:jed+1)%y =  bgrid(isd:ied, jsd:jed+1)%y

end subroutine init_grid

!-------------------------------------------------
! allocate atmos
!-------------------------------------------------
subroutine init_atmos(atm)
   type(fv_atmos_type), intent(INOUT) :: atm
   integer :: is, ie, isd, ied
   integer :: js, je, jsd, jed
   integer :: i
   character (len=60):: n, tc, hord, dp, iadv
   is  = atm%bd%is
   js  = atm%bd%js
   isd = atm%bd%isd
   jsd = atm%bd%jsd
   ie  = atm%bd%ie
   je  = atm%bd%je
   ied = atm%bd%ied
   jed = atm%bd%jed
   atm%npx = atm%bd%npx
   atm%npy = atm%bd%npy

   write(n   ,'(i8)') atm%npx
   write(tc  ,'(i8)') atm%test_case
   write(hord,'(i8)') atm%hord
   write(dp ,'(i8)') atm%dp
   write(iadv ,'(i8)') atm%inner_adv
   atm%simulation_name = "tc"//trim(adjustl(tc))//"_N"//trim(adjustl(n))//"_hord"//&
   trim(adjustl(hord))//"_iadv"//trim(adjustl(iadv))//"_dp"//trim(adjustl(dp))//"_"

   allocate(atm%qa (isd:ied,jsd:jed))
   allocate(atm%qa0(isd:ied,jsd:jed))

   allocate(atm%uc (isd:ied+1,jsd:jed))
   allocate(atm%uc0(isd:ied+1,jsd:jed))

   allocate(atm%vc (isd:ied  ,jsd:jed+1))
   allocate(atm%vc0(isd:ied,  jsd:jed+1))

   allocate(atm%error_qa(is:ie,js:je))

   atm%inner_dp  = atm%dp
   atm%outer_dp  = atm%dp

end subroutine init_atmos


!-------------------------------------------------
! init everything
!-------------------------------------------------
subroutine init_model(atm)
   type(fv_atmos_type), intent(inout) :: atm

   call init_bounds(atm%bd, atm%npx, atm%npy)
   call init_grid  (atm%gridstruct, atm%bd)
   call init_atmos (atm)

end subroutine init_model


end module fv_control 
