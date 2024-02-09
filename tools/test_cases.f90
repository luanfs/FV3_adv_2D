module test_cases
!========================================================================
! Module for ICS
!
! Reference
! https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere/blob/main/tools/test_cases.F90
!========================================================================
use fv_arrays,  only: fv_grid_bounds_type, fv_grid_type, fv_atmos_type, &
                      point_structure, R_GRID, pi, erad, eradi, day2sec, twopi, pio2
implicit none

contains

!-------------------------------------------------
! compute the ICS
!-------------------------------------------------
subroutine init_case(atm)
   type(fv_atmos_type), intent(INOUT) :: atm

   call init_scalar(atm%qa0, atm%bd, atm%gridstruct, atm%test_case)
   call init_winds (atm%uc0, atm%vc0, atm%bd, atm%gridstruct, atm%test_case)
   call calc_winds (atm%uc , atm%vc , atm%bd, atm%gridstruct, atm%time_centered, atm%test_case)
end subroutine init_case

!-------------------------------------------------
! compute initial scalar field
!-------------------------------------------------
subroutine init_scalar(qa, bd, gridstruct, test_case)
   type(fv_grid_type), target, intent(INOUT) :: gridstruct
   type(fv_grid_bounds_type), intent(IN) :: bd
   type(point_structure), pointer, dimension(:,:) :: agrid
   real(R_GRID), intent(OUT) :: qa(bd%isd:bd%ied,bd%jsd:bd%jed)
   real(R_GRID) :: x, y, c
   integer, intent(IN) :: test_case
   integer :: is, ie
   integer :: js, je
   integer :: i, j

   is = bd%is
   ie = bd%ie
   js = bd%js
   je = bd%je
 
   agrid => gridstruct%agrid

   if (test_case==1) then
      do i = is, ie
         do j = js, je
            x = agrid(i,j)%x
            y = agrid(i,j)%y
            c = 2d0*pi*erad
            c = 1.d0!2d0*pi*erad
            if (x<=0.4*c .or. x>=0.6d0*c .or. y<=0.4*c .or. y>=0.6d0*c) then
               qa(i,j) = 0.d0
            else
               qa(i,j) = 1.d0
            endif
         enddo
      enddo

   else if (test_case==2) then
      do i = is, ie
         do j = js, je
            x = agrid(i,j)%x
            y = agrid(i,j)%y
            qa(i,j) = 0.d0 + 1.d0* dexp(-10*(dcos(x*pi))**2)*dexp(-10*(dcos(y*pi))**2)
         enddo
      enddo
   else if (test_case==3 .or. test_case==4) then
      do i = is, ie
         do j = js, je
            x = agrid(i,j)%x*eradi/twopi
            y = agrid(i,j)%y*eradi/twopi
            qa(i,j) = 0.d0    + 1.d0* dexp(-10*(dsin((x-0.6)*pi))**2)*dexp(-10*(dcos(y*pi))**2)
            qa(i,j) = qa(i,j) + 1.d0* dexp(-10*(dsin((x-0.4)*pi))**2)*dexp(-10*(dcos(y*pi))**2)
         enddo
      enddo
   else


   endif


end subroutine init_scalar

!-------------------------------------------------
! compute initial winds
!-------------------------------------------------
subroutine init_winds(uc, vc, bd, gridstruct, test_case)
   type(fv_grid_type), target, intent(INOUT) :: gridstruct
   type(fv_grid_bounds_type), intent(IN) :: bd
   type(point_structure), pointer, dimension(:,:) :: cgrid
   real(R_GRID), intent(OUT) :: uc(bd%isd:bd%ied+1, bd%jsd:bd%jed  )
   real(R_GRID), intent(OUT) :: vc(bd%isd:bd%ied  , bd%jsd:bd%jed+1)
   integer, intent(IN) :: test_case

   call calc_winds(uc, vc, bd, gridstruct, 0.d0, test_case)

end subroutine init_winds


!-------------------------------------------------
! compute winds at a given time
!-------------------------------------------------
subroutine calc_winds(uc, vc, bd, gridstruct, time, test_case)
   type(fv_grid_type), target, intent(INOUT) :: gridstruct
   type(fv_grid_bounds_type), intent(IN) :: bd
   type(point_structure), pointer, dimension(:,:) :: cgrid, dgrid
   real(R_GRID), intent(OUT) :: uc(bd%isd:bd%ied+1,bd%jsd:bd%jed  )
   real(R_GRID), intent(OUT) :: vc(bd%isd:bd%ied  ,bd%jsd:bd%jed+1)
   real(R_GRID), intent(IN) :: time
   integer, intent(IN) :: test_case
   integer :: is, ie
   integer :: js, je
   integer :: i, j

   is = bd%is
   ie = bd%ie
   js = bd%js
   je = bd%je
 
   cgrid => gridstruct%cgrid
   dgrid => gridstruct%dgrid

   ! u at cgrid
   do i = is, ie+1
      do j = js, je
         call compute_wind_u(uc(i,j), cgrid(i,j)%x, cgrid(i,j)%y, time, test_case)
      enddo
   enddo

   ! v at dgrid
    do i = is, ie
      do j = js, je+1
         call compute_wind_v(vc(i,j), dgrid(i,j)%x, dgrid(i,j)%y, time, test_case)
      enddo
   enddo

  
end subroutine calc_winds

!-------------------------------------------------
! compute wind at a given time and position
!-------------------------------------------------
subroutine compute_wind_u(u, x, y, t, test_case)
   real(R_GRID), intent(OUT) :: u
   real(R_GRID), intent(IN)  :: x, y, t
   integer, intent(IN) :: test_case
   real(R_GRID) :: c, Lx, Ly, Tf, arg1, arg2, arg3, x1, y1, Ubar
   real(R_GRID) :: u0, u1

   Tf = 12.d0*day2sec
   Tf = 5.d0
   select case (test_case)
      case(1,2)
         !u = 2d0*pi*erad/Tf
         u = 0.2d0

      case(3)
         Lx = twopi
         Ly = pi

         c = (10.d0/Tf)*(Lx/twopi)**2

         x1 = x*eradi/twopi
         y1 = y*eradi/twopi

         x1 = -pi + twopi*x1
         y1 = -pio2 + pi*y1

         arg1 = twopi*(x1/Lx - t/Tf)
         arg2 = pi*y1/Ly
         arg3 = pi*t/Tf

         u = dsin(arg1)**2
         u = u*2.d0
         u = u*dcos(arg2)
         u = u*dsin(arg2)
         u = u*dcos(arg3)
         u = u*c
         u = u*pi/Ly
         u = u/twopi
         u = u-Lx/Tf
         u = -u*erad

      case(4)
         !Ubar = (2.d0*pi*erad)/Tf
         !u = -Ubar*(dsin((x+pi)/2.d0)**2)*(dsin(y))*(dcos(y*0.5d0)**2)*(dcos(pi*t/Tf))
         c = 2d0*pi*erad/Tf
         u0 = c
         u1 = c
         x1 = x/(2d0*pi*erad)
         y1 = y/(2d0*pi*erad)
         u = u0*dsin(pi*(x1-t/Tf))**2*dcos(pi*t/Tf) + u1
      case default
         print*, 'error in compute_wind: invalid testcase, ', test_case
         stop
   end select
end subroutine compute_wind_u

subroutine compute_wind_v(v, x, y, t, test_case)
   real(R_GRID), intent(OUT) :: v
   real(R_GRID), intent(IN)  :: x, y, t
   integer, intent(IN) :: test_case
   real(R_GRID) :: c, Lx, Ly, Tf, arg1, arg2, arg3, x1, y1, Ubar
   real(R_GRID) :: u0, u1

   Tf = 12.d0*day2sec
   Tf = 5.d0
   select case (test_case)
      case(1,2)
         v  = 2d0*pi*erad/Tf
         v = 0.2d0

      case(3)
         Lx = twopi
         Ly = pi

         c = (10.d0/Tf)*(Lx/twopi)**2

         x1 = x*eradi/twopi
         y1 = y*eradi/twopi

         x1 = -pi + twopi*x1
         y1 = -pio2 + pi*y1

         arg1 = twopi*(x1/Lx - t/Tf)
         arg2 = pi*y1/Ly
         arg3 = pi*t/Tf

         v = 2.d0
         v = v*dsin(arg1)
         v = v*dcos(arg1)
         v = v*dcos(arg2)**2
         v = v*dcos(arg3)
         v = v*c
         v = v*twopi/Ly
         v = -v/pi
         v = v*erad

      case(4)
         !Ubar = (2.d0*pi*erad)/Tf
         !v = (Ubar/2.d0)*(dsin((x+pi)))*(dcos(y/2.d0)**3)*(dcos(pi*t/Tf))
         c = 2d0*pi*erad/Tf
         u0 = c
         u1 = c
         x1 = x/(2d0*pi*erad)
         y1 = y/(2d0*pi*erad)
         v = u0*dsin(pi*(y1-t/Tf))**2*dcos(pi*t/Tf) + u1

      case default
         print*, 'error in compute_wind: invalid testcase, ', test_case
         stop
   end select
end subroutine compute_wind_v



end module test_cases
