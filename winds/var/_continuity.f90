!  ----------------------------------------------------------------------------
!  winds.var.continuity
!  ====================
!
!
!  ----------------------------------------------------------------------------

!- modified by oue 2017/01/06: add Jc_grid
subroutine potvin_cost(u, v, w, dudx, dvdy, dwdz, rho, drhodx, drhody, &
                       drhodz, wgt_c, fill_value, nx, ny, nz, Jc, Jc_grid)
!                       drhodz, wgt_c, fill_value, nx, ny, nz, Jc)

!  -----------------------------------------------------------------------------
!  Compute the anelastic mass continuity constraint cost Jc similar to that
!  defined in Potvin et al. 2012, "Assessing errors in variational dual-Doppler
!  wind syntheses of supercell thunderstorms observed by storm-scale mobile
!  radars".
!
!  Parameters
!  ----------
!  u : array, dim(z,y,x), float64
!     Eastward wind component in meters per second.
!  v : array, dim(z,y,x), float64
!     Northward wind component in meters per second.
!  w : array, dim(z,y,x), float64
!     Vertical wind component in meters per second.
!  dudx : array, dim(z,y,x), float64
!     Rate of change of eastward wind component with respect to the x dimension
!     in per seconds.
!  dvdy : array, dim(z,y,x), float64
!     Rate of change of northward wind component with respect to the y
!     dimension in per seconds.
!  dwdz : array, dim(z,y,x), float64
!     Rate of change of vertical wind component with respect to the z dimension
!     in per seconds.
!  rho : array, dim(z,y,x), float64
!     Air density in kilograms per cubic meter.
!  drhodx : array, dim(z,y,x), float64
!     Rate of change of air density with respect to the x dimension. This term
!     is usually very small and can often be neglected, e.g., set to zero
!     everywhere.
!  drhody : array, dim(z,y,x), float64
!     Rate of change of air density with respect to the y dimension. This term
!     is usually very small and can often be neglected, e.g., set to zero
!     everywhere.
!  drhodz : array, dim(z,y,x), float64
!     Rate of change of air density with respect to the z dimension.
!  wgt_c : float64
!     Anelastic mass continuity constraint weight.
!  fill_value : float64
!     Value indicating missing or bad grid points.
!
!  Returns
!  -------
!  Jc : float64
!     Anelastic mass continuity constraint cost.
!
! -----------------------------------------------------------------------------

   implicit none

   integer(kind=4), intent(in)                   :: nx, ny, nz
   real(kind=8), intent(in)                      :: wgt_c, fill_value
   real(kind=8), intent(in), dimension(nz,ny,nx) :: u, v, w, dudx, dvdy, &
                                                    dwdz, rho, drhodx, &
                                                    drhody, drhodz
   real(kind=8), intent(out)                     :: Jc
   real(kind=8), intent(out),dimension(nz,ny,nx) :: Jc_grid


!  Define local variables
   real(kind=8)    :: D
   integer(kind=4) :: i, j, k
!  F2PY directives
!  f2py integer(kind=4), optional, intent(in) :: nx, ny, nz
!  f2py real(kind=8), intent(in)              :: wgt_c, fill_value
!  f2py real(kind=8), intent(in)              :: u, v, w, dudx, dvdy, dwdz
!  f2py real(kind=8), intent(in)              :: rho, drhodx, drhody, drhodz
!  f2py real(kind=8), intent(out)             :: Jc


!  We are attempting to minimize a function of the form,
!
!  J = J(u1, u2, ... , uN, v1, v2, ... , vN, w1, w2, ... , wN)
!
!  which is a function of 3N variables. Note that J is typically the sum of
!  multiple different constraints, including the anelastic mass continuity cost
!  Jc, which is given by,
!
!  Jc = 0.5 * sum(wgt_c * [du/dx + dv/dy + dw/dz + (u * drho/dx + ...
!                          + v * drho/dy + w * drho/dz) / rho]**2)
!
!  where the summation is over the N Cartesian grid points. The first 3 terms
!  in the square brackets correspond to the mass flux due to wind divergence.
!  The remaining term in the parentheses represents the mass flux due to
!  advection.

   Jc = 0.d0

!  $omp parallel

!  $omp do
   do i = 1, nx
   do j = 1, ny
   do k = 1, nz

!     Calculate the main term in the square brackets found in Jc
      D = dudx(k,j,i) + dvdy(k,j,i) + dwdz(k,j,i) + &
          u(k,j,i) * drhodx(k,j,i) / rho(k,j,i) + &
          v(k,j,i) * drhody(k,j,i) / rho(k,j,i) + &
          w(k,j,i) * drhodz(k,j,i) / rho(k,j,i)

      Jc = Jc + 0.5d0 * wgt_c * D**2
      Jc_grid(k,j,i) = 0.5d0 * wgt_c * D**2

   enddo
   enddo
   enddo
!  $omp end do

!  $omp end parallel

  return

end subroutine potvin_cost


subroutine potvin_grad(u, v, w, dudx, dvdy, dwdz, rho, drhodx, drhody, &
                       drhodz, wgt_c, dx, dy, dz, finite_order, &
                       fill_value, nx, ny, nz, dJcdu, dJcdv, dJcdw)

!  -----------------------------------------------------------------------------
!  Compute the gradient of the anelastic mass continuity constraint Jc with
!  respect to the control variables (u, v, w). The constraint Jc is similar to
!  that defined in Potvin et al. 2012, "Assessing errors in variational
!  dual-Doppler wind syntheses of supercell thunderstorms observed by
!  storm-scale mobile radars".
!
!  Currently this routine does not support grids with variable resolutions.
!
!  Parameters
!  ----------
!  u : array, dim(z,y,x), float64
!     Eastward wind component in meters per second.
!  v : array, dim(z,y,x), float64
!     Northward wind component in meters per second.
!  w : array, dim(z,y,x), float64
!     Vertical wind component in meters per second.
!  dudx : array, dim(z,y,x), float64
!     Rate of change of eastward wind component with respect to the x dimension
!     in per seconds.
!  dvdy : array, dim(z,y,x), float64
!     Rate of change of northward wind component with respect to the y dimension
!     in per seconds.
!  dwdz : array, dim(z,y,x), float64
!     Rate of change of vertical wind component with respect to the z dimension
!     in per seconds.
!  rho : array, dim(z,y,x), float64
!     Air density in kilograms per cubic meter.
!  drhodx : array, dim(z,y,x), float64
!     Rate of change of air density with respect to the x dimension. This term
!     is usually very small and can often be neglected, e.g., set to zero
!     everywhere.
!  drhody : array, dim(z,y,x), float64
!     Rate of change of air density with respect to the y dimension. This term
!     is usually very small and can often be neglected, e.g., set to zero
!     everywhere.
!  drhodz : array, dim(z,y,x), float64
!     Rate of change of air density with respect to the z dimension.
!  wgt_c : float64
!     Anelastic mass continuity constraint weight.
!  dx : float64
!     Grid resolution in the x dimension.
!  dy : float64
!     Grid resolution in the y dimension.
!  dz : float64
!     Grid resolution in the z dimension.
!  finite_order : 'low' or 'high', character
!     The finite difference order used to compute differences. The order should
!     match the order used to compute the wind divergence and the air density
!     gradient.
!  fill_value : float64
!     Value indicating missing or bad grid points.
!
!  Returns
!  -------
!  dJcdu : array, dim(z,y,x), float64
!     Gradient of anelastic mass continuity constraint with respect to the
!     eastward wind component.
!  dJcdv : array, dim(z,y,x), float64
!     Gradient of anelastic mass continuity constraint with respect to the
!     northward wind component.
!  dJcdw : array, dim(z,y,x), float64
!     Gradient of anelastic mass continuity constraint with respect to the
!     vertical wind component.
!
! -----------------------------------------------------------------------------

   implicit none

   integer(kind=4), intent(in)                    :: nx, ny, nz
   character(len=16), intent(in)                  :: finite_order
   real(kind=8), intent(in)                       :: wgt_c, dx, dy, dz, &
                                                     fill_value
   real(kind=8), intent(in), dimension(nz,ny,nx)  :: u, v, w, dudx, dvdy, &
                                                     dwdz, rho, drhodx, &
                                                     drhody, drhodz
   real(kind=8), intent(out), dimension(nz,ny,nx) :: dJcdu, dJcdv, dJcdw


!  Define local variables
   real(kind=8), dimension(nz,ny,nx) :: D
   integer(kind=4)                   :: i, j, k

!  F2PY directives
!  f2py integer(kind=4), optional, intent(in) :: nx, ny, nz
!  f2py character(len=16), intent(in)         :: finite_order
!  f2py real(kind=8), intent(in)              :: wgt_c, dx, dy, dz, fill_value
!  f2py real(kind=8), intent(in)              :: u, v, w, dudx, dvdy, dwdz
!  f2py real(kind=8), intent(in)              :: rho, drhodx, drhody, drhodz
!  f2py real(kind=8), intent(out)             :: dJcdu, dJcdv, dJcdw

!  We are attempting to minimize a function of the form,
!
!  J = J(u1, u2, ... , uN, v1, v2, ... , vN, w1, w2, ... , wN)
!
!  which is a function of 3N variables. Note that J is typically the sum of
!  multiple different constraints, including the anelastic mass continuity
!  constraints Jc, which is given by,
!
!  Jc = 0.5 * sum(wgt_c * [du/dx + dv/dy + dw/dz + (u * drho/dx + ...
!                          + v * drho/dy + w * drho/dz) / rho]**2)
!
!  where the summation is over the N Cartesian grid points.
!
!  The first 3 terms in the square brackets correspond to the mass flux due to
!  wind divergence. The remaining term in the parentheses represents the mass
!  flux due to advection.
!
!  We need to compute dJ/du, dJ/dv, and dJ/dw, since a minimum in J corresponds
!  with these three derivatives vanishing. Therefore, we need to compute
!  dJc/du, dJc/dv, and dJc/dw. Each of these terms will eventually need to be
!  vectors of length N since,
!
!  dJ/d(u,v,w) = (dJ/du1, ... , dJ/duN, dJ/dv1, ... , dJ/dvN, ...
!                 dJ/dw1, ..., dJ/dwN)
!
!  We minimize J in the vector space as shown above, but we will initially
!  compute the gradient of Jc in its more natural grid space for simplicity.
!
!  The partial derivatives in the definition of Jc above must be approximated
!  by finite differences. This means that analytical solutions to dJc/du,
!  dJc/dv and dJc/dw do not exist. Each one of these terms requires an in-depth
!  analysis as to how it is computed. In particular, it requires knowledge of
!  the underlying finite difference schemes used to approximate the partial
!  derivatives since the (u1, ... , uN, v1, ... , vN, w1, ... , wN) terms may
!  have influence on surrounding grid points.


!  Calculate the main term found in the square brackets of Jc. This term is
!  required because Jc is proportional to the square of this term, which means
!  the gradient of Jc with respect to the 3 control variables will also depend
!  on this term via the chain rule

!  $omp parallel

!  $omp do
   do i = 1, nx
   do j = 1, ny
   do k = 1, nz

      D(k,j,i) = dudx(k,j,i) + dvdy(k,j,i) + dwdz(k,j,i) + &
                 u(k,j,i) * drhodx(k,j,i) / rho(k,j,i) + &
                 v(k,j,i) * drhody(k,j,i) / rho(k,j,i) + &
                 w(k,j,i) * drhodz(k,j,i) / rho(k,j,i)

   enddo
   enddo
   enddo
!  $omp end do

!  The first block is for when low-order finite difference schemes were used
!  to compute any of the partial derivatives found in the anelastic mass
!  continuity equation
!##   if (finite_order == 'low') then !# comment out to also work when finite_order == 'high' 

!     $omp do
      do i = 1, nx
      do j = 1, ny
      do k = 1, nz

!     Compute the gradient of the anelastic continuity constraint with respect
!     to the three control variables (u, v, w), which means we need to compute,
!     dJc/du, dJc/dv, and dJc/dw
!
!     Computing --> dJcdu
      if (i > 2 .and. i < nx - 1) then
         dJcdu(k,j,i) = wgt_c * ((D(k,j,i-1) - D(k,j,i+1)) / (2.d0 * dx) + &
                                 D(k,j,i) * drhodx(k,j,i) / rho(k,j,i))

      elseif (i == 2) then
         dJcdu(k,j,i) = wgt_c * ((2.d0 * D(k,j,i-1) - D(k,j,i+1)) / (2.d0 * dx) + &
                                 D(k,j,i) * drhodx(k,j,i) / rho(k,j,i))

      elseif (i == 1) then
         dJcdu(k,j,i) = -wgt_c * ((2.d0 * D(k,j,i) + D(k,j,i+1)) / (2.d0 * dx) - &
                                  D(k,j,i) * drhodx(k,j,i) / rho(k,j,i))

      elseif (i == nx - 1) then
         dJcdu(k,j,i) = wgt_c * ((D(k,j,i-1) - 2.d0 * D(k,j,i+1)) / (2.d0 * dx) + &
                                 D(k,j,i) * drhodx(k,j,i) / rho(k,j,i))

      else
         dJcdu(k,j,i) = wgt_c * ((D(k,j,i-1) + 2.d0 * D(k,j,i)) / (2.d0 * dx) + &
                                 D(k,j,i) * drhodx(k,j,i) / rho(k,j,i))

      endif

!     Computing --> dJcdv
      if (j > 2 .and. j < ny - 1) then
         dJcdv(k,j,i) = wgt_c * ((D(k,j-1,i) - D(k,j+1,i)) / (2.d0 * dy) + &
                                 D(k,j,i) * drhody(k,j,i) / rho(k,j,i))

      elseif (j == 2) then
         dJcdv(k,j,i) = wgt_c * ((2.d0 * D(k,j-1,i) - D(k,j+1,i)) / (2.d0 * dy) + &
                                 D(k,j,i) * drhody(k,j,i) / rho(k,j,i))

      elseif (j == 1) then
         dJcdv(k,j,i) = -wgt_c * ((2.d0 * D(k,j,i) + D(k,j+1,i)) / (2.d0 * dy) - &
                                  D(k,j,i) * drhody(k,j,i) / rho(k,j,i))

      elseif (j == ny - 1) then
         dJcdv(k,j,i) = wgt_c * ((D(k,j-1,i) - 2.d0 * D(k,j+1,i)) / (2.d0 * dy) + &
                                 D(k,j,i) * drhody(k,j,i) / rho(k,j,i))

      else
         dJcdv(k,j,i) = wgt_c * ((D(k,j-1,i) + 2.d0 * D(k,j,i)) / (2.d0 * dy) + &
                                 D(k,j,i) * drhody(k,j,i) / rho(k,j,i))

      endif

!     Computing --> dJcdw
      if (k > 2 .and. k < nz - 1) then
         dJcdw(k,j,i) = wgt_c * ((D(k-1,j,i) - D(k+1,j,i)) / (2.d0 * dz) + &
                                 D(k,j,i) * drhodz(k,j,i) / rho(k,j,i))

      elseif (k == 2) then
         dJcdw(k,j,i) = wgt_c * ((2.d0 * D(k-1,j,i) - D(k+1,j,i)) / (2.d0 * dz) + &
                                 D(k,j,i) * drhodz(k,j,i) / rho(k,j,i))

      elseif (k == 1) then
         dJcdw(k,j,i) = -wgt_c * ((2.d0 * D(k,j,i) + D(k+1,j,i)) / (2.d0 * dz) - &
                                  D(k,j,i) * drhodz(k,j,i) / rho(k,j,i))

      elseif (k == nz - 1) then
         dJcdw(k,j,i) = wgt_c * ((D(k-1,j,i) - 2.d0 * D(k+1,j,i)) / (2.d0 * dz) + &
                                 D(k,j,i) * drhodz(k,j,i) / rho(k,j,i))

      else
         dJcdw(k,j,i) = wgt_c * ((D(k-1,j,i) + 2.d0 * D(k,j,i)) / (2.d0 * dz) + &
                                 D(k,j,i) * drhodz(k,j,i) / rho(k,j,i))

      endif

      enddo
      enddo
      enddo
!     $omp end do

!  The second block is for when high-order finite difference schemes were used
!  to compute any of the partial differential equations found in the anelastic
!  mass continuity equation
!##   elseif (finite_order == 'high') then

!     $omp do
!##      do i = 1, nx
!##      do j = 1, ny
!##      do k = 1, nz

!     TODO: implement capabilites to compute gradients when high-order finite
!     difference schemes have been used

!##      enddo
!##      enddo
!##      enddo
!     $omp end do

!##   else

!##      stop

!##   endif

!  $omp end parallel

  return

end subroutine potvin_grad

!modified by oue 20170113
!subroutine anelastic_upwards(dudx, dvdy, rho, drhodx, drhody, drhodz, &
subroutine anelastic_upwards(dudx, dvdy, rho, hdiv, drhodx, drhody, drhodz, & 
                             dx, dy, dz, scheme, fill_value, nx, ny, nz, w)
!  ----------------------------------------------------------------------------
!  Integrate the anelastic mass continuity equation upwards from a lower (e.g.,
!  surface) boundary condition. The lower boundary condition is assumed to be
!  the an impermeability condition, i.e. vertical air motion vanishes at the
!  bottom boundary.
!
!  Parameters
!  ----------
!  dudx : array, dim(z,y,x), float64
!     Rate of change of eastward wind component with respect to the x dimension
!     in per seconds.
!  dvdy : array, dim(z,y,x), float64
!     Rate of change of northward wind component with respect to the y dimension
!     in per seconds.
!  rho : array, dim(z,y,x), float64
!     Air density in kilograms per cubic meter.
!  drhodx : array, dim(z,y,x), float64
!     Rate of change of air density with respect to the x dimension. This term
!     is usually very small and can often be neglected, e.g., set to zero
!     everywhere.
!  drhody : array, dim(z,y,x), float64
!     Rate of change of air density with respect to the y dimension. This term
!     is usually very small and can often be neglected, e.g., set to zero
!     everywhere.
!  drhodz : array, dim(z,y,x), float64
!     Rate of change of air density with respect to the z dimension.
!  dx : float64
!     Grid resolution in the x dimension.
!  dy : float64
!     Grid resolution in the y dimension.
!  dz : float64
!     Grid resolution in the z dimension.
! scheme : 'finite' or 'integral', character
!  Scheme to use when solving the anelastic mass continuity equation. If
!  'finite', the solution for w is derived from expressing the anelastic mass
!  continuity in finite  difference form and solving for w. If 'integral', the
!  solution for w is derived from integrating between consecutive height levels
!  (e.g., between k and k-1) and solving for w at k.
! fill_value : float64
!   Value indication missing or bad grid points.
!
! Returns
! -------
! w : array, dim(z,y,x), float64
!   Vertical wind component in meters per second derived from integrating the
!   anelastic mass continuity equation upwards from a lower (e.g., surface)
!   boundary condition.
!
!  ----------------------------------------------------------------------------

   implicit none

   integer(kind=4), intent(in)                    :: nx, ny, nz
   character(len=16), intent(in)                  :: scheme
   real(kind=8), intent(in)                       :: dx, dy, dz, fill_value
   real(kind=8), intent(in), dimension(nz,ny,nx)  :: dudx, dvdy, rho
   real(kind=8), intent(in), dimension(nz,ny,nx)  :: hdiv
   real(kind=8), intent(in), dimension(nz,ny,nx)  :: drhodx, drhody, drhodz
   real(kind=8), intent(out), dimension(nz,ny,nx) :: w

!  Define local variables
   integer(kind=4)             :: i, j, k

!  F2PY derivatives
!  f2py integer(kind=4), optional, intent(in) :: nx, ny, nz
!  f2py character(len=16), intent(in)         :: scheme
!  f2py real(kind=8), intent(in)              :: dx, dy, dz, fill_value
!  f2py real(kind=8), intent(in)              :: dudx, dvdy, rho
!  f2py real(kind=8), intent(in)              :: drhodx, drhody, drhodz
!  f2py real(kind=8), intent(out)             :: w

!  Anelastic mass continuity is expressed as,
!
!  du/dx + dv/dy + dw/dz + (u * drho/dx + v * drho/dy + w * drho/dz) / rho = 0
!
!  Note that the rate of change of density with respect to the x and y
!  dimensions is usually small compared to the z dimension and is therefore
!  often ignored
!
!  When the partial derivatives are expressed as finite differences, this
!  equation can be rearranged to give a solution for vertical air motion.
!
!  d(rho * w)/dz = - rho * (du/dx + dv/dy),
!
!  whereby taking the integral between consecutive height levels of both sides
!  of this equation provides an estimate for w

!  $omp parallel

!  $omp do
   do i = 1, nx
   do j = 1, ny
   do k = 1, nz

!  This block is for grid points which are not at the bottom boundary
   if (k > 1) then

   if (scheme == 'finite') then
!      w(k,j,i) = rho(k) * (w(k-1,j,i) - hdiv(k,j,i) * dz) / (rho(k) + drho(k))
      w(k,j,i) = rho(k,j,i) * (w(k-1,j,i) - hdiv(k,j,i) * dz) / (rho(k,j,i) + drhodz(k,j,i)) !modifiedd by oue 20170112

   elseif (scheme == 'integral') then
       w(k,j,i) = (rho(k-1,j,i) * w(k-1,j,i) - &
                  drhodx(k,j,i) + drhodx(k-1,j,i) - & !modifiedd by oue 20170112
                  drhodx(k,j,i) + drhody(k-1,j,i)) / rho(k,j,i) !modifiedd by oue 20170112
!                  drhoudx(k,j,i) + drhoudx(k-1,j,i) - & 
!                  drhovdx(k,j,i) + drhovdy(k-1,j,i)) / rho(k,j,i)

   else

      stop

   endif

!  This block is for grid points which are at the bottom boundary
   else
      w(k,j,i) = 0.d0

   endif

   enddo
   enddo
   enddo
!  $omp end do

!  $omp end parallel

  return

end subroutine anelastic_upwards


!- modified by oue 2017/01/06: add Jc_grid
subroutine anelastic_cost_trad(w, wc, wgt_c, fill_value, nx, ny, nz, Jc, Jc_grid)
!subroutine anelastic_cost_trad(w, wc, wgt_c, fill_value, nx, ny, nz, Jc)

! -----------------------------------------------------------------------------
!
!
! Parameters
! ----------
! w : array, dim(z,y,x), float64
!   Vertical wind component in meters per second.
! wc : array, dim(z,y,x), float64
!   Vertical wind component in meters per second estimated through integration
!   of the anelastic mass continuity equation.
! wgt_c : float64
!   Anelastic mass continuity constraint weight.
! fill_value : float64
!   The value indicating missing or bad grid points.
! nx : int32, optional
!   Size of the x dimension. If not provided, value is derived from the input
!   array size.
! ny : int32, optional
!   Size of the y dimension. If not provided, value is derived from the input
!   array size.
! nz : int32, optional
!   Size of the z dimension. If not provided, value is derived from the input
!   array size.
!
! Returns
! -------
! Jc : float64
!   Traditional anelastic mass continuity constraint cost.
!
! -----------------------------------------------------------------------------

  implicit none

  integer(kind=4), intent(in)                    :: nx, ny, nz
  real(kind=8), intent(in)                       :: wgt_c, fill_value
  real(kind=8), intent(in), dimension(nz,ny,nx)  :: w, wc
  real(kind=8), intent(out)                      :: Jc
  real(kind=8), intent(out) , dimension(nz,ny,nx) :: Jc_grid


! Define local variables ======================================================

  integer(kind=4) :: i, j, k

! =============================================================================


! F2PY directives ===========================================================

  !f2py integer(kind=4), optional, intent(in) :: nx, ny, nz
  !f2py real(kind=8), intent(in)              :: w, wc, wgt_c, fill_value
  !f2py real(kind=8), intent(out)             :: Jc

! =============================================================================


! We are attempting to minimize a function of the form,
!
! J = J(u1, u2, ... , uN, v1, v2, ... , vN, w1, w2, ... , wN)
!
! which is a function of 3N variables. Note that J is typically the sum of
! multiple different constraints, including mass continuity constraint Jc,
! which in this case is given by,
!
! Jc = 0.5 * sum( wgt_c * (w - wc)**2 )
!
! where the summation is over the N Cartesian grid points

  Jc = 0.d0

  !$omp parallel

  !$omp do
  do i = 1, nx
    do j = 1, ny
      do k = 1, nz

        Jc = Jc + 0.5d0 * wgt_c * (w(k,j,i) - wc(k,j,i))**2
        Jc_grid(k,j,i) = 0.5d0 * wgt_c * (w(k,j,i) - wc(k,j,i))**2

        enddo
    enddo
  enddo
  !$omp end do

  !$omp end parallel

  return

end subroutine anelastic_cost_trad


subroutine anelastic_grad_trad(w, wc, wgt_c, fill_value, nx, ny, nz, &
                               dJcdu, dJcdv, dJcdw)

! -----------------------------------------------------------------------------
!
!
!
! -----------------------------------------------------------------------------

  implicit none

  integer(kind=4), intent(in)                    :: nx, ny, nz
  real(kind=8), intent(in)                       :: wgt_c, fill_value
  real(kind=8), intent(in), dimension(nz,ny,nx)  :: w, wc
  real(kind=8), intent(out), dimension(nz,ny,nx) :: dJcdu, dJcdv, dJcdw


! Define local variables ======================================================

  integer(kind=4) :: i, j, k

! =============================================================================


! F2PY directives ===========================================================

  !f2py integer(kind=4), optional, intent(in) :: nx, ny, nz
  !f2py real(kind=8), intent(in)              :: w, wc, wgt_c, fill_value
  !f2py real(kind=8), intent(out)             :: dJcdu, dJcdv, dJcdw

! =============================================================================


! We are attempting to minimize a function of the form,
!
! J = J(u1, u2, ... , uN, v1, v2, ... , vN, w1, w2, ... , wN)
!
! which is a function of 3N variables. Note that J is typically the sum of
! multiple different constraints, including the anelastic mass continuity
! constraint Jc, which in this case is given by,
!
! Jc = 0.5 * sum( wgt_c * (w - wc)**2 )
!
! where the summation is over the N Cartesian grid points.
!
! We need to compute dJ/du, dJ/dv, and dJ/dw, since a minimum in J
! corresponds with these three derivatives vanishing. Therefore, we need to
! compute dJc/du, dJc/dv, and dJc/dw. Each of these terms will eventually
! need to be vectors of length N since,
!
! dJ/d(u,v,w) = (dJ/du1, ... , dJ/duN, dJ/dv1, ... , ...
!                dJ/dvN, dJ/dw1, ... , dJ/dwN)
!
! We minimize J in the vector space (1-D) as shown above, but we will
! initially compute the gradient of Jc in the natural grid space (3-D) to
! simplify matters
!
! Given the definition of Jc above, we quickly see that it has no
! dependence on u or v, i.e.,
!
! dJc/du = (0, ... , 0) for all (u1, u2, ... , uN)
! dJc/dv = (0, ... , 0) for all (v1, v2, ... , vN)

  !$omp parallel

  !$omp do
  do i = 1, nx
    do j = 1, ny
      do k = 1, nz

!     Compute the gradient of the continuity cost with respect to the three
!     control variables (u, v, w), which means we need to compute dJc/du,
!     dJc/dv, and dJc/dw. However, since the continuity cost has no dependence
!     on u or v, we have
!
!     dJc/du = 0 for all N
!     dJc/dv = 0 for all N
!
!     Furthemore, note how dJc/dw can be easily derived from Jc, and is
!     equivalent to,
!
!     dJc/dw = wgt_c * (w - wc) for all N
      dJcdu(k,j,i) = 0.d0
      dJcdv(k,j,i) = 0.d0
      dJcdw(k,j,i) = wgt_c * (w(k,j,i) - wc(k,j,i))

      enddo
    enddo
  enddo
  !$omp end do

  !$omp end parallel

  return

end subroutine anelastic_grad_trad
