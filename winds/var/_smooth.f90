! Module: smooth.f90

! modified by oue 2017/01/06: add Js_grid
subroutine potvin_cost(d2udx2, d2udy2, d2udz2, d2vdx2, d2vdy2, d2vdz2, &
                       d2wdx2, d2wdy2, d2wdz2, wgt_s1, wgt_s2, wgt_s3, &
                       wgt_s4, fill_value, nx, ny, nz, Js, Js_grid)
!                       wgt_s4, fill_value, nx, ny, nz, Js)

! -----------------------------------------------------------------------------
! Compute the spatial smoothness constraint cost Js define in Potvin et al.
! 2012, "Assessing errors in variational dual-Doppler wind syntheses of
! supercell thunderstorms observed by storm-scale mobile radars".
!
! The second order partial derivatives are none other than the vector
! Laplacian operator applied to the three wind components (u, v, w).
!
! Parameters
! ----------
! d2udx2 : array, dim(z,y,x), float64
!   Second order partial derivative of eastward wind component (u) with respect
!   to the x dimension.
!
! Returns
! -------
! Js : float64
!   Spatial smoothness constraint cost.
!
! -----------------------------------------------------------------------------

  implicit none

  integer(kind=4), intent(in)                    :: nx, ny, nz
  real(kind=8), intent(in)                       :: wgt_s1, wgt_s2, &
                                                    wgt_s3, wgt_s4, &
                                                    fill_value
  real(kind=8), intent(in), dimension(nz,ny,nx)  :: d2udx2, d2udy2, d2udz2, &
                                                    d2vdx2, d2vdy2, d2vdz2, &
                                                    d2wdx2, d2wdy2, d2wdz2
  real(kind=8), intent(out)                      :: Js
  real(kind=8), intent(out), dimension(nz,ny,nx)  :: Js_grid


! Define local variables ======================================================

  integer(kind=4) :: i, j, k

! =============================================================================


! F2PY directives =============================================================

  !f2py integer(kind=4), optional, intent(in) :: nx, ny, nz
  !f2py real(kind=8), intent(in)              :: wgt_s1, wgt_s2, wgt_s3,
  !f2py real(kind=8), intent(in)              :: wgt_s4, fill_value
  !f2py real(kind=8), intent(in)              :: d2udx2, d2udy2, d2udz2
  !f2py real(kind=8), intent(in)              :: d2vdx2, d2vdy2, d2vdz2
  !f2py real(kind=8), intent(in)              :: d2wdx2, d2wdy2, d2wdz2
  !f2py real(kind=8), intent(out)             :: Js

! =============================================================================


! We are attempting to minimize a function of the form,
!
! J = J(u1, u2, ... , uN, v1, v2, ... , vN, w1, w2, ... , wN)
!
! which is a function of 3N variables. Note that J is typically the sum of
! multiple different constraints, including the spatial smoothness constraint
! Js, which in Potvin et al. (2012) is given by,
!
! Js = 0.5 * sum( wgt_s1 * [ (d2u/dx2)**2 + (d2v/dx2)**2 +
!                           + (d2u/dy2)**2 + (d2v/dy2)**2 ] +
!                 wgt_s2 * [ (d2u/dz2)**2 + (d2v/dz2)**2 ] +
!                 wgt_s3 * [ (d2w/dx2)**2 + (d2w/dy2)**2 ] +
!                 wgt_s4 * [ (d2w/dz2)**2 ] )
!
! where the summation is over the N Cartesian grid points
  Js = 0.d0

  !$omp parallel

  !$omp do
  do i = 1, nx
  do j = 1, ny
  do k = 1, nz

! Compute the spatial smoothness constraint cost at each grid point
  Js = Js + 0.5d0 * (wgt_s1 * (d2udx2(k,j,i)**2 + d2udy2(k,j,i)**2 + &
                               d2vdx2(k,j,i)**2 + d2vdy2(k,j,i)**2) + &
                     wgt_s2 * (d2udz2(k,j,i)**2 + d2vdz2(k,j,i)**2) + &
                     wgt_s3 * (d2wdx2(k,j,i)**2 + d2wdy2(k,j,i)**2) + &
                     wgt_s4 * (d2wdz2(k,j,i)**2))
  Js_grid(k,j,i) = + 0.5d0 * (wgt_s1 * (d2udx2(k,j,i)**2 + d2udy2(k,j,i)**2 + &
                               d2vdx2(k,j,i)**2 + d2vdy2(k,j,i)**2) + &
                     wgt_s2 * (d2udz2(k,j,i)**2 + d2vdz2(k,j,i)**2) + &
                     wgt_s3 * (d2wdx2(k,j,i)**2 + d2wdy2(k,j,i)**2) + &
                     wgt_s4 * (d2wdz2(k,j,i)**2))
  enddo
  enddo
  enddo
  !$omp end do

  !$omp end parallel

  return

end subroutine potvin_cost


subroutine potvin_grad(d2udx2, d2udy2, d2udz2, d2vdx2, d2vdy2, d2vdz2, &
                       d2wdx2, d2wdy2, d2wdz2, wgt_s1, wgt_s2, wgt_s3, &
                       wgt_s4, dx, dy, dz, finite_order, fill_value, &
                       nx, ny, nz, dJsdu, dJsdv, dJsdw)

! -----------------------------------------------------------------------------
! Compute the gradient of spatial smoothness constraint Js with respect to the
! three control variables (u, v, w). The constraint Js is that defined in
! Potvin et al. 2012, "Assessing errors in variational dual-Doppler wind
! syntheses of supercell thunderstorms observed by storm-scale mobile radars".
! Currently this routine does not support grids with variable resolutions.
!
! The second order partial derivatives are none other than the vector
! Laplacian operator applied to the three wind components (u, v, w).
!
! Parameters
! ----------
! d2udx2 : array, dim(z,y,x), float64
!   Second order partial derivative of eastward wind component (u) with respect
!   to the x dimension.
!
! Returns
! -------
! dJsdu : array, dim(z,y,x), float64
!   Gradient of the spatial smoothness constraint with respect to the eastward
!   wind component (u).
! dJsdv : array, dim(z,y,x), float64
!   Gradient of the spatial smoothness constraint with respect to the northward
!   wind component (v).
! dJsdw : array, dim(z,y,x), float64
!   Gradient of the spatial smoothness constraint with respect to the vertical
!   wind component (w).
!
! -----------------------------------------------------------------------------

  implicit none

  integer(kind=4), intent(in)                    :: nx, ny, nz
  character(len=16), intent(in)                  :: finite_order
  real(kind=8), intent(in)                       :: wgt_s1, wgt_s2, &
                                                    wgt_s3, wgt_s4, &
                                                    dx, dy, dz, fill_value
  real(kind=8), intent(in), dimension(nz,ny,nx)  :: d2udx2, d2udy2, d2udz2, &
                                                    d2vdx2, d2vdy2, d2vdz2, &
                                                    d2wdx2, d2wdy2, d2wdz2
  real(kind=8), intent(out), dimension(nz,ny,nx) :: dJsdu, dJsdv, dJsdw


! Define local variables ====================================================

  real(kind=8)    :: dxx, dyy, dzz
  integer(kind=4) :: i, j, k

! ===========================================================================


! F2PY directives ===========================================================

  !f2py integer(kind=4), optional, intent(in) :: nx, ny, nz
  !f2py character(len=16), intent(in)         :: finite_order
  !f2py real(kind=8), intent(in)              :: wgt_s1, wgt_s2
  !f2py real(kind=8), intent(in)              :: wgt_s3, wgt_s4
  !f2py real(kind=8), intent(in)              :: dx, dy, dz, fill_value
  !f2py real(kind=8), intent(in)              :: d2udx2, d2udy2, d2udz2
  !f2py real(kind=8), intent(in)              :: d2vdx2, d2vdy2, d2vdz2
  !f2py real(kind=8), intent(in)              :: d2wdx2, d2wdy2, d2wdz2
  !f2py real(kind=8), intent(out)             :: dJsdu, dJsdv, dJsdw

!=============================================================================


! We are attempting to minimize a function of the form,
!
! J = J(u1,u2,...,uN,v1,v2,...,vN,w1,w2,...,wN)
!
! which is a function of 3N variables. Note that J is typically the sum of
! multiple different constraints, including the spatial smoothness constraint
! Js, which in Potvin et al. (2012) is given by,
!
! Js = 0.5 * sum( wgt_s1 * [ (d2u/dx2)**2 + (d2v/dx2)**2 +
!                           + (d2u/dy2)**2 + (d2v/dy2)**2 ] +
!                 wgt_s2 * [ (d2u/dz2)**2 + (d2v/dz2)**2 ] +
!                 wgt_s3 * [ (d2w/dx2)**2 + (d2w/dy2)**2 ] +
!                 wgt_s4 * [ (d2w/dz2)**2 ] )
!
! where the summation is over the N Cartesian grid points. We need to compute
! dJ/du, dJ/dv, and dJ/dw, since a minimum in J corresponds with these three
! derivatives vanishing. Therefore, we need to compute dJs/du, dJs/dv, and
! dJs/dw. Each of these terms will eventually need to be vectors of length N
! since,
!
! dJ/d(u,v,w) = (dJ/du1, ... , dJ/duN, dJ/dv1, ... , dJ/dvN, ...
!                dJ/dw1, ..., dJ/dwN)
!
! We minimize J in the vector space (1-D) as shown above, but we will
! initially compute the gradient of Js in the more natural grid space to
! simplify things.
!
! The partial derivatives in the definition of Js above must be approximated by
! finite differences. This means that analytical solutions to dJs/du, dJs/dv
! and dJs/dw do not exist. Each one of these terms requires an in-depth
! analysis as to how it was computed. In particular, it requires knowledge of
! the underlying finite difference schemes used to approximate the partial
! derivatives since the (u1, ... , uN, v1, ... , vN, w1, ... , wN) terms may
! have influence on surrounding grid points

  dJsdu = 0.d0
  dJsdv = 0.d0
  dJsdw = 0.d0

  dxx = dx**2
  dyy = dy**2
  dzz = dz**2

  !$omp parallel

! The first block is for when low-order finite difference schemes have been
! used to approximate the partial derivatives found in Js
  if (finite_order == 'low') then

    !$omp do
    do i = 1, nx
    do j = 1, ny
    do k = 1, nz

!   Compute the gradient of the spatial smoothness cost with respect to the
!   three control variables, which means we need to compute dJs/du, dJs/dv,
!   and dJs/dw. These 3 terms are highly dependent on the finite differences
!   used to approximate the partial derivatives, and a careful analysis must be
!   done in order to derive the correct gradient terms, especially near the
!   grid boundaries, i.e.,
!
!   i = 1, 2, 3, 4, nx-3, nx-2, nx-1, nx or
!   j = 1, 2, 3, 4, ny-3, ny-2, ny-1, ny or
!   k = 1, 2, 3, 4, nz-3, nz-2, nz-1, nz
!
!   When not near the grid boundaries there exists a general solution, but
!   closer to the grid boundaries the solution becomes highly dependent on the
!   grid point in question

!   Compute the contribution of all x derivative terms d2u/dx2, d2v/dx2, and
!   d2w/dx2 to dJs/du, dJs/dv, and dJs/dw. In other words, focus on,
!
!   Js = 0.5 * sum( wgt_s1 * [ (d2u/dx2)**2 + (d2v/dx2)**2 ] +
!                   wgt_s3 * [ (d2w/dx2)**2 ] )
    if (i > 3 .and. i < nx - 2) then
      dJsdu(k,j,i) = dJsdu(k,j,i) + wgt_s1 * (d2udx2(k,j,i-1) - 2.d0 * &
                                              d2udx2(k,j,i) + &
                                              d2udx2(k,j,i+1)) / dxx
      dJsdv(k,j,i) = dJsdv(k,j,i) + wgt_s1 * (d2vdx2(k,j,i-1) - 2.d0 * &
                                              d2vdx2(k,j,i) + &
                                              d2vdx2(k,j,i+1)) / dxx
      dJsdw(k,j,i) = dJsdw(k,j,i) + wgt_s3 * (d2wdx2(k,j,i-1) - 2.d0 * &
                                              d2wdx2(k,j,i) + &
                                              d2wdx2(k,j,i+1)) / dxx

    elseif (i == 3) then
      dJsdu(k,j,i) = dJsdu(k,j,i) + wgt_s1 * (d2udx2(k,j,i-2) + &
                                              d2udx2(k,j,i-1) - 2.d0 * &
                                              d2udx2(k,j,i) + &
                                              d2udx2(k,j,i+1)) / dxx
      dJsdv(k,j,i) = dJsdv(k,j,i) + wgt_s1 * (d2vdx2(k,j,i-2) + &
                                              d2vdx2(k,j,i-1) - 2.d0 * &
                                              d2vdx2(k,j,i) + &
                                              d2vdx2(k,j,i+1)) / dxx
      dJsdw(k,j,i) = dJsdw(k,j,i) + wgt_s3 * (d2wdx2(k,j,i-2) + &
                                              d2wdx2(k,j,i-1) - 2.d0 * &
                                              d2wdx2(k,j,i) + &
                                              d2wdx2(k,j,i+1)) / dxx

    elseif (i == 2) then
      dJsdu(k,j,i) = dJsdu(k,j,i) + wgt_s1 * (d2udx2(k,j,i+1) - 2.d0 * &
                                              d2udx2(k,j,i) - 2.d0 * &
                                              d2udx2(k,j,i-1)) / dxx
      dJsdv(k,j,i) = dJsdv(k,j,i) + wgt_s1 * (d2vdx2(k,j,i+1) - 2.d0 * &
                                              d2vdx2(k,j,i) - 2.d0 * &
                                              d2vdx2(k,j,i-1)) / dxx
      dJsdw(k,j,i) = dJsdw(k,j,i) + wgt_s3 * (d2wdx2(k,j,i+1) - 2.d0 * &
                                              d2wdx2(k,j,i) - 2.d0 * &
                                              d2wdx2(k,j,i-1)) / dxx

    elseif (i == 1) then
      dJsdu(k,j,i) = dJsdu(k,j,i) + wgt_s1 * (d2udx2(k,j,i) + &
                                              d2udx2(k,j,i+1)) / dxx
      dJsdv(k,j,i) = dJsdv(k,j,i) + wgt_s1 * (d2vdx2(k,j,i) + &
                                              d2vdx2(k,j,i+1)) / dxx
      dJsdw(k,j,i) = dJsdw(k,j,i) + wgt_s3 * (d2wdx2(k,j,i) + &
                                              d2wdx2(k,j,i+1)) / dxx

    elseif (i == nx - 2) then
      dJsdu(k,j,i) = dJsdu(k,j,i) + wgt_s1 * (d2udx2(k,j,i+2) + &
                                              d2udx2(k,j,i+1) - 2.d0 * &
                                              d2udx2(k,j,i) + &
                                              d2udx2(k,j,i-1)) / dxx
      dJsdv(k,j,i) = dJsdv(k,j,i) + wgt_s1 * (d2vdx2(k,j,i+2) + &
                                              d2vdx2(k,j,i+1) - 2.d0 * &
                                              d2vdx2(k,j,i) + &
                                              d2vdx2(k,j,i-1)) / dxx
      dJsdw(k,j,i) = dJsdw(k,j,i) + wgt_s3 * (d2wdx2(k,j,i+2) + &
                                              d2wdx2(k,j,i+1) - 2.d0 * &
                                              d2wdx2(k,j,i) + &
                                              d2wdx2(k,j,i-1)) / dxx

    elseif (i == nx - 1) then
      dJsdu(k,j,i) = dJsdu(k,j,i) + wgt_s1 * (d2udx2(k,j,i-1) - 2.d0 * &
                                              d2udx2(k,j,i) - 2.d0 * &
                                              d2udx2(k,j,i+1)) / dxx
      dJsdv(k,j,i) = dJsdv(k,j,i) + wgt_s1 * (d2vdx2(k,j,i-1) - 2.d0 * &
                                              d2vdx2(k,j,i) - 2.d0 * &
                                              d2vdx2(k,j,i+1)) / dxx
      dJsdw(k,j,i) = dJsdw(k,j,i) + wgt_s3 * (d2wdx2(k,j,i-1) - 2.d0 * &
                                              d2wdx2(k,j,i) - 2.d0 * &
                                              d2wdx2(k,j,i+1)) / dxx

    else
      dJsdu(k,j,i) = dJsdu(k,j,i) + wgt_s1 * (d2udx2(k,j,i) + &
                                              d2udx2(k,j,i-1)) / dxx
      dJsdv(k,j,i) = dJsdv(k,j,i) + wgt_s1 * (d2vdx2(k,j,i) + &
                                              d2vdx2(k,j,i-1)) / dxx
      dJsdw(k,j,i) = dJsdw(k,j,i) + wgt_s3 * (d2wdx2(k,j,i) + &
                                              d2wdx2(k,j,i-1)) / dxx
    endif

!   Compute the contribution of all the y derivative terms d2u/dy2, d2v/dy2,
!   and d2w/dy2) to dJs/du, dJs/dv, and dJs/dw. In other words, focus on,
!
!   Js = 0.5 * sum( wgt_s1 * [ (d2u/dy2)**2 + (d2v/dy2)**2 ] +
!                   wgt_s3 * [ (d2w/dy2)**2 ] )
    if (j > 3 .and. j < ny - 2) then
      dJsdu(k,j,i) = dJsdu(k,j,i) + wgt_s1 * (d2udy2(k,j-1,i) - 2.d0 * &
                                              d2udy2(k,j,i) + &
                                              d2udy2(k,j+1,i)) / dyy
      dJsdv(k,j,i) = dJsdv(k,j,i) + wgt_s1 * (d2vdy2(k,j-1,i) - 2.d0 * &
                                              d2vdy2(k,j,i) + &
                                              d2vdy2(k,j+1,i)) / dyy
      dJsdw(k,j,i) = dJsdw(k,j,i) + wgt_s3 * (d2wdy2(k,j-1,i) - 2.d0 * &
                                              d2wdy2(k,j,i) + &
                                              d2wdy2(k,j+1,i)) / dyy

    elseif (j == 3) then
      dJsdu(k,j,i) = dJsdu(k,j,i) + wgt_s1 * (d2udy2(k,j-2,i) + &
                                              d2udy2(k,j-1,i) - 2.d0 * &
                                              d2udy2(k,j,i) + &
                                              d2udy2(k,j+1,i)) / dyy
      dJsdv(k,j,i) = dJsdv(k,j,i) + wgt_s1 * (d2vdy2(k,j-2,i) + &
                                              d2vdy2(k,j-1,i) - 2.d0 * &
                                              d2vdy2(k,j,i) + &
                                              d2vdy2(k,j+1,i)) / dyy
      dJsdw(k,j,i) = dJsdw(k,j,i) + wgt_s3 * (d2wdy2(k,j-2,i) + &
                                              d2wdy2(k,j-1,i) - 2.d0 * &
                                              d2wdy2(k,j,i) + &
                                              d2wdy2(k,j+1,i)) / dyy

    elseif (j == 2) then
      dJsdu(k,j,i) = dJsdu(k,j,i) + wgt_s1 * (d2udy2(k,j+1,i) - 2.d0 * &
                                              d2udy2(k,j,i) - 2.d0 * &
                                              d2udy2(k,j-1,i)) / dyy
      dJsdv(k,j,i) = dJsdv(k,j,i) + wgt_s1 * (d2vdy2(k,j+1,i) - 2.d0 * &
                                              d2vdy2(k,j,i) - 2.d0 * &
                                              d2vdy2(k,j-1,i)) / dyy
      dJsdw(k,j,i) = dJsdw(k,j,i) + wgt_s3 * (d2wdy2(k,j+1,i) - 2.d0 * &
                                              d2wdy2(k,j,i) - 2.d0 * &
                                              d2wdy2(k,j-1,i)) / dyy

    elseif (j == 1) then
      dJsdu(k,j,i) = dJsdu(k,j,i) + wgt_s1 * (d2udy2(k,j,i) + &
                                              d2udy2(k,j+1,i)) / dyy
      dJsdv(k,j,i) = dJsdv(k,j,i) + wgt_s1 * (d2vdy2(k,j,i) + &
                                              d2vdy2(k,j+1,i)) / dyy
      dJsdw(k,j,i) = dJsdw(k,j,i) + wgt_s3 * (d2wdy2(k,j,i) + &
                                              d2wdy2(k,j+1,i)) / dyy

    elseif (j == ny - 2) then
      dJsdu(k,j,i) = dJsdu(k,j,i) + wgt_s1 * (d2udy2(k,j+2,i) + &
                                              d2udy2(k,j+1,i) - 2.d0 * &
                                              d2udy2(k,j,i) + &
                                              d2udy2(k,j-1,i)) / dyy
      dJsdv(k,j,i) = dJsdv(k,j,i) + wgt_s1 * (d2vdy2(k,j+2,i) + &
                                              d2vdy2(k,j+1,i) - 2.d0 * &
                                              d2vdy2(k,j,i) + &
                                              d2vdy2(k,j-1,i)) / dyy
      dJsdw(k,j,i) = dJsdw(k,j,i) + wgt_s3 * (d2wdy2(k,j+2,i) + &
                                              d2wdy2(k,j+1,i) - 2.d0 * &
                                              d2wdy2(k,j,i) + &
                                              d2wdy2(k,j-1,i)) / dyy

    elseif (j == ny - 1) then
      dJsdu(k,j,i) = dJsdu(k,j,i) + wgt_s1 * (d2udy2(k,j-1,i) - 2.d0 * &
                                              d2udy2(k,j,i) - 2.d0 * &
                                              d2udy2(k,j+1,i)) / dyy
      dJsdv(k,j,i) = dJsdv(k,j,i) + wgt_s1 * (d2vdy2(k,j-1,i) - 2.d0 * &
                                              d2vdy2(k,j,i) - 2.d0 * &
                                              d2vdy2(k,j+1,i)) / dyy
      dJsdw(k,j,i) = dJsdw(k,j,i) + wgt_s3 * (d2wdy2(k,j-1,i) - 2.d0 * &
                                              d2wdy2(k,j,i) - 2.d0 * &
                                              d2wdy2(k,j+1,i)) / dyy

    else
      dJsdu(k,j,i) = dJsdu(k,j,i) + wgt_s1 * (d2udy2(k,j,i) + &
                                              d2udy2(k,j-1,i)) / dyy
      dJsdv(k,j,i) = dJsdv(k,j,i) + wgt_s1 * (d2vdy2(k,j,i) + &
                                              d2vdy2(k,j-1,i)) / dyy
      dJsdw(k,j,i) = dJsdw(k,j,i) + wgt_s3 * (d2wdy2(k,j,i) + &
                                              d2wdy2(k,j-1,i)) / dyy

    endif

!   Compute the contribution of all the z derivative terms d2u/dz2, d2v/dz2,
!   and d2w/dz2 to dJs/du, dJs/dv, and dJs/dw. In other words, focus on,
!
!   Js = 0.5 * sum( wgt_s2 * [ (d2u/dz2)**2 + (d2v/dz2)**2 ] +
!                   wgt_s4 * [ (d2w/dz2)**2 ] )
    if (k > 3 .and. k < nz - 2) then
      dJsdu(k,j,i) = dJsdu(k,j,i) + wgt_s2 * (d2udz2(k-1,j,i) - 2.d0 * &
                                              d2udz2(k,j,i) + &
                                              d2udz2(k+1,j,i)) / dzz
      dJsdv(k,j,i) = dJsdv(k,j,i) + wgt_s2 * (d2vdz2(k-1,j,i) - 2.d0 * &
                                              d2vdz2(k,j,i) + &
                                              d2vdz2(k+1,j,i)) / dzz
      dJsdw(k,j,i) = dJsdw(k,j,i) + wgt_s4 * (d2wdz2(k-1,j,i) - 2.d0 * &
                                              d2wdz2(k,j,i) + &
                                              d2wdz2(k+1,j,i)) / dzz

    elseif (k == 3) then
      dJsdu(k,j,i) = dJsdu(k,j,i) + wgt_s2 * (d2udz2(k-2,j,i) + &
                                              d2udz2(k-1,j,i) - 2.d0 * &
                                              d2udz2(k,j,i) + &
                                              d2udz2(k+1,j,i)) / dzz
      dJsdv(k,j,i) = dJsdv(k,j,i) + wgt_s2 * (d2vdz2(k-2,j,i) + &
                                              d2vdz2(k-1,j,i) - 2.d0 * &
                                              d2vdz2(k,j,i) + &
                                              d2vdz2(k+1,j,i)) / dzz
      dJsdw(k,j,i) = dJsdw(k,j,i) + wgt_s4 * (d2wdz2(k-2,j,i) + &
                                              d2wdz2(k-1,j,i) - 2.d0 * &
                                              d2wdz2(k,j,i) + &
                                              d2wdz2(k+1,j,i)) / dzz

    elseif (k == 2) then
      dJsdu(k,j,i) = dJsdu(k,j,i) + wgt_s2 * (d2udz2(k+1,j,i) - 2.d0 * &
                                              d2udz2(k,j,i) - 2.d0 * &
                                              d2udz2(k-1,j,i)) / dzz
      dJsdv(k,j,i) = dJsdv(k,j,i) + wgt_s2 * (d2vdz2(k+1,j,i) - 2.d0 * &
                                              d2vdz2(k,j,i) - 2.d0 * &
                                              d2vdz2(k-1,j,i)) / dzz
      dJsdw(k,j,i) = dJsdw(k,j,i) + wgt_s4 * (d2wdz2(k+1,j,i) - 2.d0 * &
                                              d2wdz2(k,j,i) - 2.d0 * &
                                              d2wdz2(k-1,j,i)) / dzz

    elseif (k == 1) then
      dJsdu(k,j,i) = dJsdu(k,j,i) + wgt_s2 * (d2udz2(k,j,i) + &
                                              d2udz2(k+1,j,i)) / dzz
      dJsdv(k,j,i) = dJsdv(k,j,i) + wgt_s2 * (d2vdz2(k,j,i) + &
                                              d2vdz2(k+1,j,i)) / dzz
      dJsdw(k,j,i) = dJsdw(k,j,i) + wgt_s4 * (d2wdz2(k,j,i) + &
                                              d2wdz2(k+1,j,i)) / dzz

    elseif (k == nz - 2) then
      dJsdu(k,j,i) = dJsdu(k,j,i) + wgt_s2 * (d2udz2(k+2,j,i) + &
                                              d2udz2(k+1,j,i) - 2.d0 * &
                                              d2udz2(k,j,i) + &
                                              d2udz2(k-1,j,i)) / dzz
      dJsdv(k,j,i) = dJsdv(k,j,i) + wgt_s2 * (d2vdz2(k+2,j,i) + &
                                              d2vdz2(k+1,j,i) - 2.d0 * &
                                              d2vdz2(k,j,i) + &
                                              d2vdz2(k-1,j,i)) / dzz
      dJsdw(k,j,i) = dJsdw(k,j,i) + wgt_s4 * (d2wdz2(k+2,j,i) + &
                                              d2wdz2(k+1,j,i) - 2.d0 * &
                                              d2wdz2(k,j,i) + &
                                              d2wdz2(k-1,j,i)) / dzz

    elseif (k == nz - 1) then
      dJsdu(k,j,i) = dJsdu(k,j,i) + wgt_s2 * (d2udz2(k-1,j,i) - 2.d0 * &
                                              d2udz2(k,j,i) - 2.d0 * &
                                              d2udz2(k+1,j,i)) / dzz
      dJsdv(k,j,i) = dJsdv(k,j,i) + wgt_s2 * (d2vdz2(k-1,j,i) - 2.d0 * &
                                              d2vdz2(k,j,i) - 2.d0 * &
                                              d2vdz2(k+1,j,i)) / dzz
      dJsdw(k,j,i) = dJsdw(k,j,i) + wgt_s4 * (d2wdz2(k-1,j,i) - 2.d0 * &
                                              d2wdz2(k,j,i) - 2.d0 * &
                                              d2wdz2(k+1,j,i)) / dzz

    else
      dJsdu(k,j,i) = dJsdu(k,j,i) + wgt_s2 * (d2udz2(k,j,i) + &
                                              d2udz2(k-1,j,i)) / dzz
      dJsdv(k,j,i) = dJsdv(k,j,i) + wgt_s2 * (d2vdz2(k,j,i) + &
                                              d2vdz2(k-1,j,i)) / dzz
      dJsdw(k,j,i) = dJsdw(k,j,i) + wgt_s4 * (d2wdz2(k,j,i) + &
                                              d2wdz2(k-1,j,i)) / dzz

    endif

    enddo
    enddo
    enddo
    !$omp end do


! The second block is for when high-order finite difference schemes have been
! used to approximate the partial derivatives found in Js
  elseif (finite_order == 'high') then

    !$omp do
    do i = 1, nx
    do j = 1, ny
    do k = 1, nz

!   TODO: Add capabilities to compute the gradient of the spatial smoothness
!   constraint when high-order finite difference schemes have been used to
!   approximate the partial derivatives

    enddo
    enddo
    enddo
    !$omp end do

  else

    stop

  endif

  !$omp end parallel

  return

end subroutine potvin_grad
