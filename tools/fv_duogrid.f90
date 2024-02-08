module fv_duogrid
!========================================================================
! Module to handle with ghost cell interpolations
!
!========================================================================
use fv_arrays,  only: fv_grid_bounds_type,  &
                      R_GRID
implicit none

contains

!--------------------------------------------------------------------------
! extend agrid field
subroutine ext_scalar_agrid(qa, bd)
   type(fv_grid_bounds_type), intent(IN) :: bd
   real(R_GRID), intent(INOUT) :: qa(bd%isd:bd%ied, bd%jsd:bd%jed)
   integer :: is, ie, isd, ied, ng
   integer :: js, je, jsd, jed

   is = bd%is
   js = bd%js
   ie = bd%ie
   je = bd%je
   isd = bd%isd
   jsd = bd%jsd
   ied = bd%ied
   jed = bd%jed
   ng = bd%ng

   ! periodic BC
   qa(isd:is-1, js:je) = qa(ie-ng+1:ie, js:je)
   qa(ie+1:ied, js:je) = qa(is:is+ng-1, js:je)
   qa(isd:ied, jsd:js-1) = qa(isd:ied,je-ng+1:je)
   qa(isd:ied, je+1:jed) = qa(isd:ied,js:js+ng-1)

end subroutine ext_scalar_agrid
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! extend cgrid field
subroutine ext_scalar_cgrid(uc, vc, bd)
   type(fv_grid_bounds_type), intent(IN) :: bd
   real(R_GRID), intent(INOUT) :: uc(bd%isd:bd%ied+1, bd%jsd:bd%jed)
   real(R_GRID), intent(INOUT) :: vc(bd%isd:bd%ied  , bd%jsd:bd%jed+1)
   integer :: is, ie, isd, ied, ng
   integer :: js, je, jsd, jed

   is = bd%is
   js = bd%js

   ie = bd%ie
   je = bd%je

   isd = bd%isd
   jsd = bd%jsd

   ied = bd%ied
   jed = bd%jed
   ng = bd%ng

   ! periodic BC
   uc(isd:is-1, js:je   ) = uc(ie-ng+1:ie    , js:je)
   uc(ie+2:ied, js:je   ) = uc(is+1:is+1+ng-1, js:je)
   uc(isd:ied , jsd:js-1) = uc(isd:ied, je-ng+1:je)
   uc(isd:ied , je+1:jed) = uc(isd:ied, js:js+ng-1)

   vc(is:ie   , jsd:js-1) = vc(is:ie, je-ng+1:je)
   vc(is:ie   , je+2:jed) = vc(is:ie, js+1:js+1+ng-1)
   vc(isd:is-1, jsd:jed)  = vc(ie-ng+1:ie, jsd:jed)
   vc(ie+1:ied, jsd:jed)  = vc(is:is+ng-1, jsd:jed)


end subroutine ext_scalar_cgrid
!--------------------------------------------------------------------------


end module fv_duogrid
