! Module: background.f90

!modified by oue 20170106: added Jb_grid
subroutine wind_cost(u, v, w, ub, vb, wb, wgt_ub, wgt_vb, wgt_wb, &
                     fill_value, proc, nx, ny, nz, Jb, Jb_grid)
!                     fill_value, proc, nx, ny, nz, Jb)

! -----------------------------------------------------------------------------
! Compute the background wind field constraint cost Jb.
!
! Parameters
! ----------
! u : array, dim(z,y,x), float64
!   The eastward wind component of the analysis wind field in meters per
!   second.
! v : array, dim(z,y,x), float64
!   The northward wind component of the analysis wind field in meters per
!   second.
! w : array, dim(z,y,x), float64
!   The vertical wind component of the analysis wind field in meters per
!   second.
! ub : array, dim(z,y,x), float64
!   The eastward wind component of the background wind field in meters per
!   second.
! vb : array, dim(z,y,x), float64
!   The northward wind component of the background wind field in meters per
!   second.
! wb : array, dim(z,y,x), float64
!   The vertical wind component of the background wind field in meters per
!   second.
! wgt_ub : float64
!   Background constraint weight for the eastward wind component.
! wgt_vb : float64
!   Background constraint weight for the northward wind component.
! wgt_wb : float64
!   Background constraint weight for the vertical wind component. This should
!   usually be set to zero since vertical air motion information is typically
!   not provided by the background wind field.
!  fill_value : float64
!     Value indicating missing or bad grid points.
!  proc : int32
!     The number of parallel threads (CPUs) to use.
!
! Returns
! -------
! Jb : float64
!   Background wind field constraint cost.
!
! -----------------------------------------------------------------------------

   implicit none

   integer(kind=4), intent(in)                    :: nx, ny, nz, proc
   real(kind=8), intent(in)                       :: wgt_ub, wgt_vb, &
                                                     wgt_wb, fill_value
   real(kind=8), intent(in), dimension(nz,ny,nx)  :: u, v, w, &
                                                     ub, vb, wb
   real(kind=8), intent(out)                      :: Jb
   real(kind=8), intent(out), dimension(nz,ny,nx)  :: Jb_grid


!  Define local variables =====================================================

   integer(kind=4) :: i, j, k

!  ============================================================================


!  F2PY directives ============================================================

   !f2py integer(kind=4), optional, intent(in) :: nx, ny, nz
   !f2py integer(kind=4), intent(in)           :: proc
   !f2py real(kind=8), intent(in)              :: wgt_ub, wgt_vb, wgt_wb
   !f2py real(kind=8), intent(in)              :: fill_value
   !f2py real(kind=8), intent(in)              :: u, v, w, ub, vb, wb
   !f2py real(kind=8), intent(out)             :: Jb

!  ============================================================================


!  We are attempting to minimize a function of the form,
!
!  J = J(u1, u2, ... , uN, v1, v2, ... , vN, w1, w2, ... , wN)
!
!  which is a function of 3N variables. Note that J is typically the sum of
!  multiple different constraints, including the background cost Jb, which is
!  given by,
!
!  Jb = 0.5 * [ sum( wgt_ub * (u - ub)**2 ) + sum( wgt_vb * (v - vb)**2 ) +
!             + sum( wgt_wb * (w - wb)**2 ) ]
!
!  where the summations are over the N Cartesian grid points.

   Jb = 0.d0

   !$omp parallel num_threads(proc)

   !$omp do
   do i = 1, nx
      do j = 1, ny
         do k = 1, nz

!         Computing --> background wind field cost Jb
          Jb = Jb + 0.5d0 * (wgt_ub * (u(k,j,i) - ub(k,j,i))**2 + &
                             wgt_vb * (v(k,j,i) - vb(k,j,i))**2 + &
                             wgt_wb * (w(k,j,i) - wb(k,j,i))**2)
          Jb_grid(k,j,i)=0.5d0 * (wgt_ub * (u(k,j,i) - ub(k,j,i))**2 + &
                             wgt_vb * (v(k,j,i) - vb(k,j,i))**2 + &
                             wgt_wb * (w(k,j,i) - wb(k,j,i))**2)

         enddo
      enddo
   enddo
   !$omp end do

   !$omp end parallel

   return

end subroutine wind_cost


subroutine wind_grad(u, v, w, ub, vb, wb, wgt_ub, wgt_vb, wgt_wb, &
                     fill_value, proc, nx, ny, nz, dJbdu, dJbdv, dJbdw)

! -----------------------------------------------------------------------------
! Compute the background wind field constraint gradient.
!
! Parameters
! ----------
! u : array, dim(z,y,x), float64
!   The eastward wind component of the analysis wind field in meters per
!   second.
! v : array, dim(z,y,x), float64
!   The northward wind component of the analysis wind field in meters per
!   second.
! w : array, dim(z,y,x), float64
!   The vertical wind component of the analysis wind field in meters per
!   second.
! ub : array, dim(z,y,x), float64
!   The eastward wind component of the background wind field in meters per
!   second.
! vb : array, dim(z,y,x), float64
!   The northward wind component of the background wind field in meters per
!   second.
! wb : array, dim(z,y,x), float64
!   The vertical wind component of the background wind field in meters per
!   second.
! wgt_ub : float64
!   Background constraint weight for the eastward wind component.
! wgt_vb : float64
!   Background constraint weight for the northward wind component.
! wgt_wb : float64
!   Background constraint weight for the vertical wind component. This should
!   usually be set to zero since vertical air motion information is typically
!   not provided by the background wind field.
!  fill_value : float64
!     The value indicating missing or bad grid points.
!  proc : int32
!     The number of parallel threads (CPUs) to use.
!
! Returns
! -------
! dJbdu : array, dim(z,y,x), float64
!   Gradient of background constraint with respect to eastward wind component.
! dJbdv : array, dim(z,y,x), float64
!   Gradient of background constraint with respect to northward wind component.
! dJbdw : array, dim(z,y,x), float64
!   Gradient of background constraint with respect to vertical wind component.
!
! -----------------------------------------------------------------------------

   implicit none

   integer(kind=4), intent(in)                    :: nx, ny, nz, proc
   real(kind=8), intent(in)                       :: wgt_ub, wgt_vb, &
                                                     wgt_wb, fill_value
   real(kind=8), intent(in), dimension(nz,ny,nx)  :: u, v, w, ub, vb, wb
   real(kind=8), intent(out), dimension(nz,ny,nx) :: dJbdu, dJbdv, dJbdw


!  Define local variables ====================================================

   integer(kind=4) :: i, j, k

!  ===========================================================================


!  F2PY directives ===========================================================

   !f2py integer(kind=4), optional, intent(in) :: nx, ny, nz
   !f2py integer(kind=4), intent(in)           :: proc
   !f2py real(kind=8), intent(in)              :: wgt_ub, wgt_vb, wgt_wb
   !f2py real(kind=8), intent(in)              :: fill_value
   !f2py real(kind=8), intent(in)              :: u, v, w, ub, vb, wb
   !f2py real(kind=8), intent(out)             :: dJbdu, dJbdv, dJbdw

!  ===========================================================================


!  We are attempting to minimize a function of the form,
!
!  J = J(u1, u2, ... , uN, v1, v2, ... , vN, w1, w2, ... , wN)
!
!  which is a function of 3N variables. Note that J is typically the sum of
!  multiple different constraints, including the background cost Jb, which is
!  given by,
!
!  Jb = 0.5 * [ sum( wgt_ub * (u - ub)**2 ) + sum( wgt_vb * (v - vb)**2 ) +
!             + sum( wgt_wb * (w - wb)**2 ) ]
!
!  where the summations are over the N Cartesian grid points. We need to
!  compute dJ/du, dJ/dv, and dJ/dw, since a minimum in J corresponds with these
!  3 derivatives vanishing. Therefore, we need to compute dJb/du, dJb/dv, and
!  dJb/dw. Each of these terms will eventually need to be vectors of length N
!  since,
!
!  dJ/d(u,v,w) = (dJ/du1, ... , dJ/duN, dJ/dv1, ... ,
!                 dJ/dvN, dJ/dw1, ... , dJ/dwN)
!
!  We minimize J in the vector space (1-D) as shown above, but we will
!  initially compute the gradient of Jb in the grid space since this is the
!  most intuitive way.

   !$omp parallel num_threads(proc)
   !$omp do
   do i = 1, nx
      do j = 1, ny
         do k = 1, nz

!        Compute the gradient of the background constaint with respect to the
!        3 control variables (u, v, w), which means we need to compute dJb/du,
!        dJb/dv, and dJb/dw. These terms are easily derived from Jb,
!
!        dJb/du = wgt_ub * (u - ub) for all N
!        dJb/dv = wgt_vb * (v - vb) for all N
!        dJb/dw = wgt_wb * (w - wb) for all N

         dJbdu(k,j,i) = wgt_ub * (u(k,j,i) - ub(k,j,i))
         dJbdv(k,j,i) = wgt_vb * (v(k,j,i) - vb(k,j,i))
         dJbdw(k,j,i) = wgt_wb * (w(k,j,i) - wb(k,j,i))

         enddo
      enddo
   enddo
   !$omp end do
   !$omp end parallel

   return

end subroutine wind_grad
