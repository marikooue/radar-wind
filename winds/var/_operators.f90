!  ----------------------------------------------------------------------------
!  winds.var.operators
!  ===================
!
!
!  ----------------------------------------------------------------------------


subroutine gradient(f, dx, dy, dz, finite_order, fill_value, &
                    proc, nx, ny, nz, dfdx, dfdy, dfdz)
!  ----------------------------------------------------------------------------
!  Compute the gradient of a scalar function defined in three-dimensional
!  space.
!
!  Currently this routine does not support grids with variable resolutions.
!
!  Parameters
!  ----------
!  f : array, dim(z,y,x), float64
!     The samples of a scalar function in three-dimensional space.
!  dx : float64
!     Grid resolution in the x dimension in meters.
!  dy : float64
!     Grid resolution in the y dimension in meters.
!  dz : float64
!     Grid resolution in the z dimension in meters.
!  finite_order : 'low' or 'high', character
!     The finite difference order to use when computing differences. A 'low'
!     order will use 1st and 2nd order finite differences, while a 'high' order
!     will use 6th order differences where possible.
!  fill_value : float64
!     The value indicating missing or bad grid points.
!  proc : int32
!     The number of parallel threads (CPUs) to use. If this value is larger
!     than the maximum number of available threads, the maximum is used.
!
!  Returns
!  -------
!  dfdx : array, dim(z,y,x), float64
!     The derivative of the input samples with respect to the x dimension.
!  dfdy : array, dim(z,y,x), float64
!     The derivative of the input samples with respect to the y dimension.
!  dfdz : array, dim(z,y,x), float64
!     The derivative of the input samples with respect to the z dimension.
!
!  ----------------------------------------------------------------------------

   implicit none

!  Define inputs and outputs
   integer(kind=4), intent(in)                    :: nx, ny, nz, proc
   character(len=16), intent(in)                  :: finite_order
   real(kind=8), intent(in)                       :: dx, dy, dz, fill_value
   real(kind=8), intent(in), dimension(nz,ny,nx)  :: f
   real(kind=8), intent(out), dimension(nz,ny,nx) :: dfdx, dfdy, dfdz

!  Define local variables
   integer(kind=4)         :: i, j, k
   real(kind=8), parameter :: a1=2.d0/3.d0, a2=1.d0/12.d0, &
                              b1=3.d0/4.d0, b2=3.d0/20.d0, &
                              b3=-1.d0/60.d0, c0=-49.d0/20.d0, &
                              c1=6.d0, c2=-15.d0/2.d0, &
                              c3=20.d0/3.d0, c4=-15.d0/4.d0, &
                              c5=6.d0/5.d0, c6=-1.d0/6.d0

!  F2PY directives
!  f2py integer(kind=4), optional, intent(in) :: nx, ny, nz
!  f2py integer(kind=4), intent(in)           :: proc
!  f2py character(len=16), intent(in)         :: finite_order
!  f2py real(kind=8), intent(in)              :: dx, dy, dz, fill_value
!  f2py real(kind=8), intent(in)              :: f
!  f2py real(kind=8), intent(out)             :: dfdx, dfdy, dfdz

!  $omp parallel num_threads(proc)

   if (finite_order == 'low') then

!     $omp do
      do i = 1, nx
      do j = 1, ny
      do k = 1, nz

!     For interior grid points, i.e.,
!
!     i = [2, nx-1] or
!     j = [2, ny-1] or
!     k = [2, nz-1]
!
!     use a centered difference scheme with p = 2. When at the grid boundaries,
!     i.e.,
!
!     i = 1, nx or
!     j = 1, ny or
!     k = 1, nz
!
!     use either a forward or backward difference scheme, both with p = 1.

!     Computing --> df/dx
      if (i > 1 .and. i < nx) then
         dfdx(k,j,i) = (f(k,j,i+1) - f(k,j,i-1)) / (2.d0 * dx)

      elseif (i == 1) then
         dfdx(k,j,i) = (f(k,j,i+1) - f(k,j,i)) / dx

      else
         dfdx(k,j,i) = (f(k,j,i) - f(k,j,i-1)) / dx

      endif

!     Computing --> df/dy
      if (j > 1 .and. j < ny) then
         dfdy(k,j,i) = (f(k,j+1,i) - f(k,j-1,i)) / (2.d0 * dy)

      elseif (j == 1) then
         dfdy(k,j,i) = (f(k,j+1,i) - f(k,j,i)) / dy

      else
         dfdy(k,j,i) = (f(k,j,i) - f(k,j-1,i)) / dy

      endif

!     Computing --> df/dz
      if (k > 1 .and. k < nz) then
         dfdz(k,j,i) = (f(k+1,j,i) - f(k-1,j,i)) / (2.d0 * dz)

      elseif (k == 1) then
         dfdz(k,j,i) = (f(k+1,j,i) - f(k,j,i)) / dz

      else
         dfdz(k,j,i) = (f(k,j,i) - f(k-1,j,i)) / dz

      endif

      enddo
      enddo
      enddo
!     $omp end do

!  The second block is for high-order finite difference schemes
   elseif (finite_order == 'high') then

      do i = 1, nx
      do j = 1, ny
      do k = 1, nz

!     TODO: implement the higher-order schemes

      enddo
      enddo
      enddo

   else

      stop

   endif

!  $omp end parallel

   return

end subroutine gradient


subroutine div_wind(u, v, w, dx, dy, dz, finite_order, fill_value, &
                    proc, nx, ny, nz, div, dudx, dvdy, dwdz)
! -----------------------------------------------------------------------------
!  Compute the divergence of the wind field (u, v, w). Currently this routine
!  does not support grids with variable resolutions.
!
!  Parameters
!  ----------
!  u : array, dim(z,y,x), float64
!     Estward wind component in meters per second.
!  v : array, dim(z,y,x), float64
!     Northward wind component in meters per second.
!  w : array, dim(z,y,x), float64
!     Vertical wind component in meters per second.
!  dx : float64
!     Grid resolution in the x dimension in meters.
!  dy : float64
!     Grid resolution in the y dimension in meters.
!  dz : float64
!     Grid resolution in the z dimension in meters.
!  finite_order : 'low' or 'high', character
!     The finite difference order to use when computing divergence field. A
!     'low' order will use 1st and 2nd order finite differences, while a 'high'
!     order will use 6th order differences where possible.
!  fill_value : float64
!     The value indicating missing or bad grid points.
!  proc : int32
!     The number of parallel threads (CPUs) to use. If this value is larger
!     than the maximum number of available threads, the maximum is used.
!
!  Returns
!  -------
!  div : array, dim(z,y,x), float64
!     The wind divergence in per second.
!  dudx : array, dim(z,y,x), float64
!     The rate of change of the eastward wind component with respect to the x
!     dimension.
!  dvdy : array, dim(z,y,x), float64
!     The rate of change of the northward wind component with respect to the y
!     dimension.
!  dwdz : array, dim(z,y,x), float64
!     The rate of change of the vertical wind component with respect to the z
!     dimension.
!
! -----------------------------------------------------------------------------

   implicit none

   integer(kind=4), intent(in)                    :: nx, ny, nz, proc
   character(len=16), intent(in)                  :: finite_order
   real(kind=8), intent(in)                       :: dx, dy, dz, fill_value
   real(kind=8), intent(in), dimension(nz,ny,nx)  :: u, v, w
   real(kind=8), intent(out), dimension(nz,ny,nx) :: div, dudx, dvdy, dwdz

!  Define local variables
   integer(kind=4)         :: i, j, k
   real(kind=8), parameter :: a1=2.d0/3.d0, a2=1.d0/12.d0, &
                              b1=3.d0/4.d0, b2=3.d0/20.d0, &
                              b3=-1.d0/60.d0, c0=-49.d0/20.d0, &
                              c1=6.d0, c2=-15.d0/2.d0, &
                              c3=20.d0/3.d0, c4=-15.d0/4.d0, &
                              c5=6.d0/5.d0, c6=-1.d0/6.d0

!  F2PY directives
!  f2py integer(kind=4), optional, intent(in) :: nx, ny, nz
!  f2py integer(kind=4), intent(in)           :: proc
!  f2py character(len=16), intent(in)         :: finite_order
!  f2py real(kind=8), intent(in)              :: dx, dy, dz, fill_value
!  f2py real(kind=8), intent(in)              :: u, v, w
!  f2py real(kind=8), intent(out)             :: div, dudx, dvdy, dwdz


!  The wind divergence is given by the sum of the 3 terms du/dx, dv/dy
!  and dw/dz. These partial derivatives should be computed in their natural
!  three-dimensional space rather than in a vector space for simplicity

!  $omp parallel num_threads(proc)

!  The first block is for low-order finite difference schemes
   if (finite_order == 'low') then

!     $omp do
      do i = 1, nx
      do j = 1, ny
      do k = 1, nz

!     For interior grid points, i.e.,
!
!     i = [2, nx-1] or
!     j = [2, ny-1] or
!     k = [2, nz-1]
!
!     use a centered difference scheme with p = 2. When at the grid boundaries,
!     i.e.,
!
!     i = 1, nx or
!     j = 1, ny or
!     k = 1, nz
!
!     use either a forward or backward difference scheme, both with p = 1.

!     Computing --> du/dx
      if (i > 1 .and. i < nx) then
         dudx(k,j,i) = (u(k,j,i+1) - u(k,j,i-1)) / (2.d0 * dx)

      elseif (i == 1) then
         dudx(k,j,i) = (u(k,j,i+1) - u(k,j,i)) / dx

      else
         dudx(k,j,i) = (u(k,j,i) - u(k,j,i-1)) / dx

      endif

!     Computing --> dv/dy
      if (j > 1 .and. j < ny) then
         dvdy(k,j,i) = (v(k,j+1,i) - v(k,j-1,i)) / (2.d0 * dy)

      elseif (j == 1) then
         dvdy(k,j,i) = (v(k,j+1,i) - v(k,j,i)) / dy

      else
         dvdy(k,j,i) = (v(k,j,i) - v(k,j-1,i)) / dy

      endif

!     Computing --> dw/dz
      if (k > 1 .and. k < nz) then
         dwdz(k,j,i) = (w(k+1,j,i) - w(k-1,j,i)) / (2.d0 * dz)

      elseif (k == 1) then
         dwdz(k,j,i) = (w(k+1,j,i) - w(k,j,i)) / dz

      else
         dwdz(k,j,i) = (w(k,j,i) - w(k-1,j,i)) / dz

      endif

      enddo
      enddo
      enddo
!     $omp end do


!  The second block is for high-order finite difference schemes
   elseif (finite_order == 'high') then

      !$omp do
      do i = 1, nx
         do j = 1, ny
            do k = 1, nx

!           For the interior grid points, i.e.,
!
!           i = [4, nx-3] or
!           j = [4, ny-3] or
!           k = [4, nz-3]
!
!           use a centered difference scheme with p = 6. When closer to the
!           grid boundaries, i.e.,
!
!           i = 2, 3, nx-2, nx-1 or
!           j = 2, 3, ny-2, ny-1 or
!           k = 2, 3, nz-2, nz-1
!
!           still use centered difference schemes, but of lower accuracy,
!           i.e. p = 2 or p = 4. When at the grid boundaries, i.e.,
!
!           i = 1, nx or
!           j = 1, ny or
!           k = 1, nz
!
!           use either a forward or backward difference scheme, both with
!           p = 6.

!           Computing --> du/dx

            if (i > 3 .and. i < nx - 2) then
               dudx(k,j,i) = (b3 * u(k,j,i-3) + b2 * u(k,j,i-2) - &
                              b1 * u(k,j,i-1) + b1 * u(k,j,i+1) - &
                              b2 * u(k,j,i+2) - b3 * u(k,j,i+3)) / dx

            elseif (i > 2 .and. i < nx - 1) then
               dudx(k,j,i) = (a2 * u(k,j,i-2) - a1 * u(k,j,i-1) + &
                              a1 * u(k,j,i+1) - a2 * u(k,j,i+2)) / dx

            elseif (i > 1 .and. i < nx) then
               dudx(k,j,i) = (u(k,j,i+1) - u(k,j,i-1)) / (2.d0 * dx)

            elseif (i == 1) then
               dudx(k,j,i) = (c0 * u(k,j,i) + c1 * u(k,j,i+1) + &
                              c2 * u(k,j,i+2) + c3 * u(k,j,i+3) + &
                              c4 * u(k,j,i+4) + c5 * u(k,j,i+5) + &
                              c6 * u(k,j,i+6)) / dx

            else
               dudx(k,j,i) = (-c0 * u(k,j,i) - c1 * u(k,j,i-1) - &
                               c2 * u(k,j,i-2) - c3 * u(k,j,i-3) - &
                               c4 * u(k,j,i-4) - c5 * u(k,j,i-5) - &
                               c6 * u(k,j,i-6)) / dx

            endif

!           Computing --> dv/dy
            if (j > 3 .and. j < ny - 2) then
               dvdy(k,j,i) = (b3 * v(k,j-3,i) + b2 * v(k,j-2,i) - &
                              b1 * v(k,j-1,i) + b1 * v(k,j+1,i) - &
                              b2 * v(k,j+2,i) - b3 * v(k,j+3,i)) / dy

            elseif (j > 2 .and. j < ny - 1) then
               dvdy(k,j,i) = (a2 * v(k,j-2,i) - a1 * v(k,j-1,i) + &
                              a1 * v(k,j+1,i) - a2 * v(k,j+2,i)) / dy

            elseif (j > 1 .and. j < ny) then
               dvdy(k,j,i) = (v(k,j+1,i) - v(k,j-1,i)) / (2.d0 * dy)

            elseif (j == 1) then
               dvdy(k,j,i) = (c0 * v(k,j,i) + c1 * v(k,j+1,i) + &
                              c2 * v(k,j+2,i) + c3 * v(k,j+3,i) + &
                              c4 * v(k,j+4,i) + c5 * v(k,j+5,i) + &
                              c6 * v(k,j+6,i)) / dy

            else
               dvdy(k,j,i) = (-c0 * v(k,j,i) - c1 * v(k,j-1,i) - &
                               c2 * v(k,j-2,i) - c3 * v(k,j-3,i) - &
                               c4 * v(k,j-4,i) - c5 * v(k,j-5,i) - &
                               c6 * v(k,j-6,i)) / dy

            endif

!           Computing --> dw/dz
            if (k > 3 .and. k < nz - 2) then
               dwdz(k,j,i) = (b3 * w(k-3,j,i) + b2 * w(k-2,j,i) - &
                              b1 * w(k-1,j,i) + b1 * w(k+1,j,i) - &
                              b2 * w(k+2,j,i) - b3 * w(k+3,j,i)) / dz

            elseif (k > 2 .and. k < nz - 1) then
               dwdz(k,j,i) = (a2 * w(k-2,j,i) - a1 * w(k-1,j,i) + &
                              a1 * w(k+1,j,i) - a2 * w(k+2,j,i)) / dz

            elseif (k > 1 .and. k < nz) then
               dwdz(k,j,i) = (w(k+1,j,i) - w(k-1,j,i)) / (2.d0 * dz)

            elseif (k == 1) then
               dwdz(k,j,i) = (c0 * w(k,j,i) + c1 * w(k+1,j,i) + &
                              c2 * w(k+2,j,i) + c3 * w(k+3,j,i) + &
                              c4 * w(k+4,j,i) + c5 * w(k+5,j,i) + &
                              c6 * w(k+6,j,i)) / dz

            else
               dwdz(k,j,i) = (-c0 * w(k,j,i) - c1 * w(k-1,j,i) - &
                               c2 * w(k-2,j,i) - c3 * w(k-3,j,i) - &
                               c4 * w(k-4,j,i) - c5 * w(k-5,j,i) - &
                               c6 * w(k-6,j,i)) / dz

            endif

            enddo
         enddo
      enddo
      !$omp end do

   else

      stop

   endif

!  Computing --> wind divergence
   div = dudx + dvdy + dwdz

   !$omp end parallel

   return

end subroutine div_wind


subroutine div_horiz_wind(u, v, dx, dy, finite_order, fill_value, &
                          proc, nx, ny, nz, div, dudx, dvdy)
! -----------------------------------------------------------------------------
!  Compute the divergence of the horizontal wind field (u, v). Currently this
!  routine does not support grids with variable resolutions.
!
!  Parameters
!  ----------
!  u : array, dim(z,y,x), float64
!     The eastward wind component in meters per second.
!  v : array, dim(z,y,x), float64
!     The northward wind component in meters per second.
!  dx : float64
!     Grid resolution in the x dimension.
!  dy : float64
!     Grid resolution in the y dimension.
!  finite_order : 'low' or 'high', character
!     The finite difference order to use when computing the divergence. A low
!     order will use 1st and 2nd order finite differences, while a high order
!     will use 6th order differences when possible.
!  fill_value : float64
!     The value indicating missing or bad grid points.
!  proc : int32
!     The number of parallel threads (CPUs) to use.
!
!  Returns
!  -------
!  div : array, dim(z,y,x), float64
!     The horizontal wind divergence in per second.
!  dudx : array, dim(z,y,x), float64
!     The rate of change of the eastward wind component (u) with respect to the
!     x dimension.
!  dvdy : array, dim(z,y,x), float64
!     The rate of change of the northward wind component (v) with respect to
!     the y dimension.
!
!------------------------------------------------------------------------------

   implicit none

   integer(kind=4), intent(in)                    :: nx, ny, nz, proc
   character(len=16), intent(in)                  :: finite_order
   real(kind=8), intent(in)                       :: dx, dy, fill_value
   real(kind=8), intent(in), dimension(nz,ny,nx)  :: u, v
   real(kind=8), intent(out), dimension(nz,ny,nx) :: div, dudx, dvdy


!  Define local variables =====================================================

   real(kind=8), parameter :: a1=2.d0/3.d0, a2=1.d0/12.d0, &
                              b1=3.d0/4.d0, b2=3.d0/20.d0, &
                              b3=-1.d0/60.d0, c0=-49.d0/20.d0, &
                              c1=6.d0, c2=-15.d0/2.d0, &
                              c3=20.d0/3.d0, c4=-15.d0/4.d0, &
                              c5=6.d0/5.d0, c6=-1.d0/6.d0

   integer(kind=4)         :: i, j, k

!  ============================================================================


!  F2PY directives ============================================================

   !f2py integer(kind=4), optional, intent(in) :: nx, ny, nz
   !f2py integer(kind=4), intent(in)           :: proc
   !f2py character(len=16), intent(in)         :: finite_order
   !f2py real(kind=8), intent(in)              :: u, v, dx, dy, fill_value
   !f2py real(kind=8), intent(out)             :: div, dudx, dvdy

!  ============================================================================


!  The horizontal wind divergence is given by the sum of du/dx and dv/dy. These
!  partial derivatives should be computed in the natural grid space rather than
!  in a vector space simply to make it easier to interpret.
!
!  The first block is for low-order finite difference schemes.

   !$omp parallel num_threads(proc)
   if (finite_order == 'low') then

      !$omp do
      do i = 1, nx
         do j = 1, ny
            do k = 1, nz

!           For interior grid points, i.e.,
!
!           i = [2, nx-1] or
!           j = [2, ny-1]
!
!           use a centered difference scheme with p = 2, and at the grid
!           boundaries, i.e.,
!
!           i = 1, nx or
!           j = 1, ny
!
!           use either a forward or backward difference scheme, both with
!           p = 1.

!           Computing --> du/dx
            if (i > 1 .and. i < nx) then
               dudx(k,j,i) = (u(k,j,i+1) - u(k,j,i-1)) / (2.d0 * dx)

            elseif (i == 1) then
               dudx(k,j,i) = (u(k,j,i+1) - u(k,j,i)) / dx

            else
               dudx(k,j,i) = (u(k,j,i) - u(k,j,i-1)) / dx
            endif

!           Computing --> dv/dy
            if (j == 1 .and. ny < 2) then  !added by oue for RHI grid
               dvdy(k,j,i) = 0.0
            else
            if (j > 1 .and. j < ny) then
               dvdy(k,j,i) = (v(k,j+1,i) - v(k,j-1,i)) / (2.d0 * dy)

            elseif (j == 1) then
               dvdy(k,j,i) = (v(k,j+1,i) - v(k,j,i)) / dy

            else
               dvdy(k,j,i) = (v(k,j,i) - v(k,j-1,i)) / dy
            endif
            endif!added by oue for RHI grid

            enddo
         enddo
      enddo
      !$omp end do


!  The second block is for high-order finite difference schemes

   elseif (finite_order == 'high') then

      !$omp do
      do i =  1, nx
         do j = 1, ny
            do k = 1, nz

!           For interior points, i.e.,
!
!           i = [4, nx-3] or
!           j = [4, ny-3]
!
!           use a centered difference scheme with p = 6. Closer to the grid
!           boundaries, i.e.,
!
!           i = 2, 3, nx-2, nx-1 or
!           j = 2, 3, ny-2, ny-1
!
!           still use centered difference schemes, but of lower accuracy, i.e.,
!           p = 2 or p = 4. When at the grid boundaries, i.e.,
!
!           i = 1, nx or
!           j = 1, ny
!
!           use either a forward or backward difference scheme, both with
!           p = 6.

!           Computing --> du/dx
            if (i > 3 .and. i < nx - 2) then
               dudx(k,j,i) = (b3 * u(k,j,i-3) + b2 * u(k,j,i-2) - &
                              b1 * u(k,j,i-1) + b1 * u(k,j,i+1) - &
                              b2 * u(k,j,i+2) - b3 * u(k,j,i+3)) / dx

            elseif (i > 2 .and. i < nx - 1) then
               dudx(k,j,i) = (a2 * u(k,j,i-2) - a1 * u(k,j,i-1) + &
                              a1 * u(k,j,i+1) - a2 * u(k,j,i+2)) / dx

            elseif (i > 1 .and. i < nx) then
               dudx(k,j,i) = (u(k,j,i+1) - u(k,j,i-1)) / (2.d0 * dx)

            elseif (i == 1) then
               dudx(k,j,i) = (c0 * u(k,j,i) + c1 * u(k,j,i+1) + &
                              c2 * u(k,j,i+2) + c3 * u(k,j,i+3) + &
                              c4 * u(k,j,i+4) + c5 * u(k,j,i+5) + &
                              c6 * u(k,j,i+6)) / dx

            else
               dudx(k,j,i) = (-c0 * u(k,j,i) - c1 * u(k,j,i-1) - &
                               c2 * u(k,j,i-2) - c3 * u(k,j,i-3) - &
                               c4 * u(k,j,i-4) - c5 * u(k,j,i-5) - &
                               c6 * u(k,j,i-6)) / dx
            endif

!           Computing --> dv/dy
            if (j > 3 .and. j < ny - 2) then
               dvdy(k,j,i) = (b3 * v(k,j-3,i) + b2 * v(k,j-2,i) - &
                              b1 * v(k,j-1,i) + b1 * v(k,j+1,i) - &
                              b2 * v(k,j+2,i) - b3 * v(k,j+3,i)) / dy

            elseif (j > 2 .and. j < ny - 1) then
               dvdy(k,j,i) = (a2 * v(k,j-2,i) - a1 * v(k,j-1,i) + &
                              a1 * v(k,j+1,i) - a2 * v(k,j+2,i)) / dy

            elseif (j > 1 .and. j < ny) then
               dvdy(k,j,i) = (v(k,j+1,i) - v(k,j-1,i)) / (2.d0 * dy)

            elseif (j == 1) then
               dvdy(k,j,i) = (c0 * v(k,j,i) + c1 * v(k,j+1,i) + &
                              c2 * v(k,j+2,i) + c3 * v(k,j+3,i) + &
                              c4 * v(k,j+4,i) + c5 * v(k,j+5,i) + &
                              c6 * v(k,j+6,i)) / dy

            else
               dvdy(k,j,i) = (-c0 * v(k,j,i) - c1 * v(k,j-1,i) - &
                               c2 * v(k,j-2,i) - c3 * v(k,j-3,i) - &
                               c4 * v(k,j-4,i) - c5 * v(k,j-5,i) - &
                               c6 * v(k,j-6,i)) / dy
            endif

            enddo
         enddo
      enddo
      !$omp end do

   else

      stop

   endif

!  Computing --> horizontal wind divergence
   div = dudx + dvdy

   !$omp end parallel

   return

end subroutine div_horiz_wind


subroutine laplace_wind(u, v, w, dx, dy, dz, finite_order, fill_value, &
                        proc, nx, ny, nz, d2udx2, d2udy2, d2udz2, d2vdx2, &
                        d2vdy2, d2vdz2, d2wdx2, d2wdy2, d2wdz2)

!  ----------------------------------------------------------------------------
!  Compute the vector Laplacian of the wind field (u, v, w). Whereas the scalar
!  Laplacian applies to scalar fields and returns a scalar quantity, the vector
!  Laplacian applies to vector fields and returns a vector quantity.
!
!  Parameters
!  ----------
!  u : array, dim(z,y,x), float64
!     The eastward wind component in meters per second.
!  v : array, dim(z,y,x), float64
!     The northward wind component in meters per second.
!  w : array, dim(z,y,x), float64
!     The vertical wind component in meters per second.
!  dx : float64
!     Grid resolution in the x dimension in meters.
!  dy : float64
!     Grid resolution in the y dimension in meters.
!  dz : float64
!     Grid resolution in the z dimension in meters.
!  finite_order : character
!     The finite difference order to use when computing the divergence. Options
!     are 'low' or 'high'. A low order will use 1st and 2nd order finite
!     differences, while a high order will use 6th order differences when
!     possible.
!  fill_value : float64
!     The value indicating missing or bad grid points.
!  proc : int32
!     The number of parallel threads (CPUs) to use.
!
!  Results
!  -------
!  d2udx2 : array, dim(z,y,x), float64
!   Second order partial derivative of eastward wind component (u) with respect
!   to the x dimension.
!
!  ----------------------------------------------------------------------------

   implicit none

   integer(kind=4), intent(in)                    :: nx, ny, nz, proc
   character(len=16), intent(in)                  :: finite_order
   real(kind=8), intent(in)                       :: dx, dy, dz, fill_value
   real(kind=8), intent(in), dimension(nz,ny,nx)  :: u, v, w
   real(kind=8), intent(out), dimension(nz,ny,nx) :: d2udx2, d2udy2, d2udz2, &
                                                     d2vdx2, d2vdy2, d2vdz2, &
                                                     d2wdx2, d2wdy2, d2wdz2


!  Define local variables =====================================================

   real(kind=8)            :: dxx, dyy, dzz
   integer(kind=4)         :: i, j, k

   real(kind=8), parameter :: a0=-5.d0/4.d0, a1=2.d0/3.d0, &
                              a2=-1.d0/24.d0, b0=-49.d0/36.d0, &
                              b1=3.d0/4.d0, b2=-3.d0/40.d0, &
                              b3=1.d0/180.d0, c0=469.d0/180.d0, &
                              c1=-223.d0/20.d0, c2=879.d0/40.d0, &
                              c3=-949.d0/36.d0, c4=41.d0/2.d0, &
                              c5=-201.d0/20.d0, c6=1019.d0/360.d0, &
                              c7=-7.d0/20.d0

!  ============================================================================


!  F2PY directives ============================================================

   !f2py integer(kind=4), optional, intent(in) :: nx, ny, nz
   !f2py integer(kind=4), intent(in)           :: proc
   !f2py character(len=16), intent(in)         :: finite_order
   !f2py real(kind=8), intent(in)              :: dx, dy, dz, fill_value
   !f2py real(kind=8), intent(in)              :: u, v, w
   !f2py real(kind=8), intent(out)             :: d2udx2, d2udy2, d2udz2
   !f2py real(kind=8), intent(out)             :: d2vdx2, d2vdy2, d2vdz2
   !f2py real(kind=8), intent(out)             :: d2wdx2, d2wdy2, d2wdz2

!  ============================================================================


!  Compute the vector Laplacian of the 3 wind components (u, v, w),
!  which means we have a total of 9 terms to compute at each grid point,
!
!  d2u/dx2, d2u/dy2, d2u/dz2
!  d2v/dx2, d2v/dy2, d2v/dz2
!  d2w/dx2, d2w/dy2, d2w/dz2
!
!  We compute these partial derivatives in their natural grid space rather than
!  the vector space to simplify interpretation.

   dxx = dx**2
   dyy = dy**2
   dzz = dz**2

   !$omp parallel num_threads(proc)

!  The first block is for low-order finite difference schemes
   if (finite_order == 'low') then

      !$omp do
      do i = 1, nx
      do j = 1, ny
      do k = 1, nz

!     For interior grid points, i.e.,
!
!     i = [2, nx-1] or
!     j = [2, ny-1] or
!     k = [2, nz-1]
!
!     use a centered difference scheme with p = 2. When at the grid boundaries,
!     i.e.,
!
!     i = 1, nx or
!     j = 1, ny or
!     k = 1, nz
!
!     use either a forward or backward difference scheme, both with p = 1.

!     Computing --> d2u/dx2, d2v/dx2, d2w/dx2
      if (i > 1 .and. i < nx) then
         d2udx2(k,j,i) = (u(k,j,i+1) - 2.d0 * u(k,j,i) + u(k,j,i-1)) / dxx
         d2vdx2(k,j,i) = (v(k,j,i+1) - 2.d0 * v(k,j,i) + v(k,j,i-1)) / dxx
         d2wdx2(k,j,i) = (w(k,j,i+1) - 2.d0 * w(k,j,i) + w(k,j,i-1)) / dxx

      elseif (i == 1) then
         d2udx2(k,j,i) = (u(k,j,i) - 2.d0 * u(k,j,i+1) + u(k,j,i+2)) / dxx
         d2vdx2(k,j,i) = (v(k,j,i) - 2.d0 * v(k,j,i+1) + v(k,j,i+2)) / dxx
         d2wdx2(k,j,i) = (w(k,j,i) - 2.d0 * w(k,j,i+1) + w(k,j,i+2)) / dxx

      else
         d2udx2(k,j,i) = (u(k,j,i) - 2.d0 * u(k,j,i-1) + u(k,j,i-2)) / dxx
         d2vdx2(k,j,i) = (v(k,j,i) - 2.d0 * v(k,j,i-1) + v(k,j,i-2)) / dxx
         d2wdx2(k,j,i) = (w(k,j,i) - 2.d0 * w(k,j,i-1) + w(k,j,i-2)) / dxx

      endif

!     Computing --> d2u/dy2, d2v/dy2, d2w/dy2
      if (j > 1 .and. j < ny) then
         d2udy2(k,j,i) = (u(k,j+1,i) - 2.d0 * u(k,j,i) + u(k,j-1,i)) / dyy
         d2vdy2(k,j,i) = (v(k,j+1,i) - 2.d0 * v(k,j,i) + v(k,j-1,i)) / dyy
         d2wdy2(k,j,i) = (w(k,j+1,i) - 2.d0 * w(k,j,i) + w(k,j-1,i)) / dyy

      elseif (j == 1) then
         d2udy2(k,j,i) = (u(k,j,i) - 2.d0 * u(k,j+1,i) + u(k,j+2,i)) / dyy
         d2vdy2(k,j,i) = (v(k,j,i) - 2.d0 * v(k,j+1,i) + v(k,j+2,i)) / dyy
         d2wdy2(k,j,i) = (w(k,j,i) - 2.d0 * w(k,j+1,i) + w(k,j+2,i)) / dyy

      else
         d2udy2(k,j,i) = (u(k,j,i) - 2.d0 * u(k,j-1,i) + u(k,j-2,i)) / dyy
         d2vdy2(k,j,i) = (v(k,j,i) - 2.d0 * v(k,j-1,i) + v(k,j-2,i)) / dyy
         d2wdy2(k,j,i) = (w(k,j,i) - 2.d0 * w(k,j-1,i) + w(k,j-2,i)) / dyy

      endif

!     Computing --> d2u/dz2, d2v/dz2, d2w/dz2
      if (k > 1 .and. k < nz) then
         d2udz2(k,j,i) = (u(k+1,j,i) - 2.d0 * u(k,j,i) + u(k-1,j,i)) / dzz
         d2vdz2(k,j,i) = (v(k+1,j,i) - 2.d0 * v(k,j,i) + v(k-1,j,i)) / dzz
         d2wdz2(k,j,i) = (w(k+1,j,i) - 2.d0 * w(k,j,i) + w(k-1,j,i)) / dzz

      elseif (k == 1) then
         d2udz2(k,j,i) = (u(k,j,i) - 2.d0 * u(k+1,j,i) + u(k+2,j,i)) / dzz
         d2vdz2(k,j,i) = (v(k,j,i) - 2.d0 * v(k+1,j,i) + v(k+2,j,i)) / dzz
         d2wdz2(k,j,i) = (w(k,j,i) - 2.d0 * w(k+1,j,i) + w(k+2,j,i)) / dzz

      else
         d2udz2(k,j,i) = (u(k,j,i) - 2.d0 * u(k-1,j,i) + u(k-2,j,i)) / dzz
         d2vdz2(k,j,i) = (v(k,j,i) - 2.d0 * v(k-1,j,i) + v(k-2,j,i)) / dzz
         d2wdz2(k,j,i) = (w(k,j,i) - 2.d0 * w(k-1,j,i) + w(k-2,j,i)) / dzz

      endif

      enddo
      enddo
      enddo
      !$omp end do


!  The second block is for high-order finite difference schemes
   elseif (finite_order == 'high') then

      !$omp do
      do i = 1, nx
      do j = 1, ny
      do k = 1, nz

!     For the interior grid points, i.e.,
!
!     i = [4, nx-3] or
!     j = [4, ny-3] or
!     k = [4, nz-3]
!
!     use a centered difference scheme with p = 6. When closer to the grid
!     boundaries, i.e.,
!
!     i = 2, 3, nx-2, nx-1 or
!     j = 2, 3, ny-2, ny-1 or
!     k = 2, 3, nz-2, nz-1
!
!     still use a centered difference schemes, but of lower accuracy, i.e.,
!     p = 2 or p = 4. When at the grid boundaries, i.e.,
!
!     i = 1, nx or
!     j = 1, ny or
!     k = 1, nz
!
!     use either a forward of backward difference scheme, both with p = 6.

!     Computing --> d2u/dx2, d2v/dx2, d2w/dx2
      if (i > 3 .and. i < nx - 2) then
         d2udx2(k,j,i) = (b3 * u(k,j,i-3) + b2 * u(k,j,i-2) + &
                          b1 * u(k,j,i-1) + b0 * u(k,j,i) + &
                          b1 * u(k,j,i+1) + b2 * u(k,j,i+2) + &
                          b3 * u(k,j,i+3)) / dxx
         d2vdx2(k,j,i) = (b3 * v(k,j,i-3) + b2 * v(k,j,i-2) + &
                          b1 * v(k,j,i-1) + b0 * v(k,j,i) + &
                          b1 * v(k,j,i+1) + b2 * v(k,j,i+2) + &
                          b3 * v(k,j,i+3)) / dxx
         d2wdx2(k,j,i) = (b3 * w(k,j,i-3) + b2 * w(k,j,i-2) + &
                          b1 * w(k,j,i-1) + b0 * w(k,j,i) + &
                          b1 * w(k,j,i+1) + b2 * w(k,j,i+2) + &
                          b3 * w(k,j,i+3)) / dxx

      elseif (i > 2 .and. i < nx - 1) then
         d2udx2(k,j,i) = (a2 * u(k,j,i-2) + a1 * u(k,j,i-1) + &
                          a0 * u(k,j,i) + a1 * u(k,j,i+1) + &
                          a2 * u(k,j,i+2)) / dxx
         d2vdx2(k,j,i) = (a2 * v(k,j,i-2) + a1 * v(k,j,i-1) + &
                          a0 * v(k,j,i) + a1 * v(k,j,i+1) + &
                          a2 * v(k,j,i+2)) / dxx
         d2wdx2(k,j,i) = (a2 * w(k,j,i-2) + a1 * w(k,j,i-1) + &
                          a0 * w(k,j,i) + a1 * w(k,j,i+1) + &
                          a2 * w(k,j,i+2)) / dxx

      elseif (i > 1 .and. i < nx) then
         d2udx2(k,j,i) = (u(k,j,i+1) - 2.d0 * u(k,j,i) + u(k,j,i-1)) / dxx
         d2vdx2(k,j,i) = (v(k,j,i+1) - 2.d0 * v(k,j,i) + v(k,j,i-1)) / dxx
         d2wdx2(k,j,i) = (w(k,j,i+1) - 2.d0 * w(k,j,i) + w(k,j,i-1)) / dxx

      elseif (i == 1) then
         d2udx2(k,j,i) = (c0 * u(k,j,i) + c1 * u(k,j,i+1) + &
                          c2 * u(k,j,i+2) + c3 * u(k,j,i+3) + &
                          c4 * u(k,j,i+4) + c5 * u(k,j,i+5) + &
                          c6 * u(k,j,i+6) + c7 * u(k,j,i+7)) / dxx
         d2vdx2(k,j,i) = (c0 * v(k,j,i) + c1 * v(k,j,i+1) + &
                          c2 * v(k,j,i+2) + c3 * v(k,j,i+3) + &
                          c4 * v(k,j,i+4) + c5 * v(k,j,i+5) + &
                          c6 * v(k,j,i+6) + c7 * v(k,j,i+7)) / dxx
         d2wdx2(k,j,i) = (c0 * w(k,j,i) + c1 * w(k,j,i+1) + &
                          c2 * w(k,j,i+2) + c3 * w(k,j,i+3) + &
                          c4 * w(k,j,i+4) + c5 * w(k,j,i+5) + &
                          c6 * w(k,j,i+6) + c7 * w(k,j,i+7)) / dxx

      else
         d2udx2(k,j,i) = (c0 * u(k,j,i) + c1 * u(k,j,i-1) + &
                          c2 * u(k,j,i-2) + c3 * u(k,j,i-3) + &
                          c4 * u(k,j,i-4) + c5 * u(k,j,i-5) + &
                          c6 * u(k,j,i-6) + c7 * u(k,j,i-7)) / dxx
         d2vdx2(k,j,i) = (c0 * v(k,j,i) + c1 * v(k,j,i-1) + &
                          c2 * v(k,j,i-2) + c3 * v(k,j,i-3) + &
                          c4 * v(k,j,i-4) + c5 * v(k,j,i-5) + &
                          c6 * v(k,j,i-6) + c7 * v(k,j,i-7)) / dxx
         d2wdx2(k,j,i) = (c0 * w(k,j,i) + c1 * w(k,j,i-1) + &
                          c2 * w(k,j,i-2) + c3 * w(k,j,i-3) + &
                          c4 * w(k,j,i-4) + c5 * w(k,j,i-5) + &
                          c6 * w(k,j,i-6) + c7 * w(k,j,i-7)) / dxx

      endif

!     Computing --> d2u/dy2, d2v/dy2, d2w/dy2
      if (j > 3 .and. j < nx - 2) then
         d2udy2(k,j,i) = (b3 * u(k,j-3,i) + b2 * u(k,j-2,i) + &
                          b1 * u(k,j-1,i) + b0 * u(k,j,i) + &
                          b1 * u(k,j+1,i) + b2 * u(k,j+2,i) + &
                          b3 * u(k,j+3,i)) / dyy
         d2vdy2(k,j,i) = (b3 * v(k,j-3,i) + b2 * v(k,j-2,i) + &
                          b1 * v(k,j-1,i) + b0 * v(k,j,i) + &
                          b1 * v(k,j+1,i) + b2 * v(k,j+2,i) + &
                          b3 * v(k,j+3,i)) / dyy
         d2wdy2(k,j,i) = (b3 * w(k,j-3,i) + b2 * w(k,j-2,i) + &
                          b1 * w(k,j-1,i) + b0 * w(k,j,i) + &
                          b1 * w(k,j+1,i) + b2 * w(k,j+2,i) + &
                          b3 * w(k,j+3,i)) / dyy

      elseif (j > 2 .and. j < nx - 1) then
         d2udy2(k,j,i) = (a2 * u(k,j-2,i) + a1 * u(k,j-1,i) + &
                          a0 * u(k,j,i) + a1 * u(k,j+1,i) + &
                          a2 * u(k,j+2,i)) / dyy
         d2vdy2(k,j,i) = (a2 * v(k,j-2,i) + a1 * v(k,j-1,i) + &
                          a0 * v(k,j,i) + a1 * v(k,j+1,i) + &
                          a2 * v(k,j+2,i)) / dyy
         d2wdy2(k,j,i) = (a2 * w(k,j-2,i) + a1 * w(k,j-1,i) + &
                          a0 * w(k,j,i) + a1 * w(k,j+1,i) + &
                          a2 * w(k,j+2,i)) / dyy

      elseif (j > 1 .and. j < ny) then
         d2udy2(k,j,i) = (u(k,j+1,i) - 2.d0 * u(k,j,i) + u(k,j-1,i)) / dyy
         d2vdy2(k,j,i) = (v(k,j+1,i) - 2.d0 * v(k,j,i) + v(k,j-1,i)) / dyy
         d2wdy2(k,j,i) = (w(k,j+1,i) - 2.d0 * w(k,j,i) + w(k,j-1,i)) / dyy

      elseif (j == 1) then
         d2udy2(k,j,i) = (c0 * u(k,j,i) + c1 * u(k,j+1,i) + &
                          c2 * u(k,j+2,i) + c3 * u(k,j+3,i) + &
                          c4 * u(k,j+4,i) + c5 * u(k,j+5,i) + &
                          c6 * u(k,j+6,i) + c7 * u(k,j+7,i)) / dyy
         d2vdy2(k,j,i) = (c0 * v(k,j,i) + c1 * v(k,j+1,i) + &
                          c2 * v(k,j+2,i) + c3 * v(k,j+3,i) + &
                          c4 * v(k,j+4,i) + c5 * v(k,j+5,i) + &
                          c6 * v(k,j+6,i) + c7 * v(k,j+7,i)) / dyy
         d2wdy2(k,j,i) = (c0 * w(k,j,i) + c1 * w(k,j+1,i) + &
                          c2 * w(k,j+2,i) + c3 * w(k,j+3,i) + &
                          c4 * w(k,j+4,i) + c5 * w(k,j+5,i) + &
                          c6 * w(k,j+6,i) + c7 * w(k,j+7,i)) / dyy

      else
         d2udy2(k,j,i) = (c0 * u(k,j,i) + c1 * u(k,j-1,i) + &
                          c2 * u(k,j-2,i) + c3 * u(k,j-3,i) + &
                          c4 * u(k,j-4,i) + c5 * u(k,j-5,i) + &
                          c6 * u(k,j-6,i) + c7 * u(k,j-7,i)) / dyy
         d2vdy2(k,j,i) = (c0 * v(k,j,i) + c1 * v(k,j-1,i) + &
                          c2 * v(k,j-2,i) + c3 * v(k,j-3,i) + &
                          c4 * v(k,j-4,i) + c5 * v(k,j-5,i) + &
                          c6 * v(k,j-6,i) + c7 * v(k,j-7,i)) / dyy
         d2wdy2(k,j,i) = (c0 * w(k,j,i) + c1 * w(k,j-1,i) + &
                          c2 * w(k,j-2,i) + c3 * w(k,j-3,i) + &
                          c4 * w(k,j-4,i) + c5 * w(k,j-5,i) + &
                          c6 * w(k,j-6,i) + c7 * w(k,j-7,i)) / dyy

      endif

!     Computing --> d2u/dz2, d2v/dz2, d2w/dz2
      if (k > 3 .and. k < nz - 2) then
         d2udz2(k,j,i) = (b3 * u(k-3,j,i) + b2 * u(k-2,j,i) + &
                          b1 * u(k-1,j,i) + b0 * u(k,j,i) + &
                          b1 * u(k+1,j,i) + b2 * u(k+2,j,i) + &
                          b3 * u(k+3,j,i)) / dzz
         d2vdz2(k,j,i) = (b3 * v(k-3,j,i) + b2 * v(k-2,j,i) + &
                          b1 * v(k-1,j,i) + b0 * v(k,j,i) + &
                          b1 * v(k+1,j,i) + b2 * v(k+2,j,i) + &
                          b3 * v(k+3,j,i)) / dzz
         d2wdz2(k,j,i) = (b3 * w(k-3,j,i) + b2 * w(k-2,j,i) + &
                          b1 * w(k-1,j,i) + b0 * w(k,j,i) + &
                          b1 * w(k+1,j,i) + b2 * w(k+2,j,i) + &
                          b3 * w(k+3,j,i)) / dzz

      elseif (k > 2 .and. k < nz - 1) then
         d2udz2(k,j,i) = (a2 * u(k-2,j,i) + a1 * u(k-1,j,i) + &
                          a0 * u(k,j,i) + a1 * u(k+1,j,i) + &
                          a2 * u(k+2,j,i)) / dzz
         d2vdz2(k,j,i) = (a2 * v(k-2,j,i) + a1 * v(k-1,j,i) + &
                          a0 * v(k,j,i) + a1 * v(k+1,j,i) + &
                          a2 * v(k+2,j,i)) / dzz
         d2wdz2(k,j,i) = (a2 * w(k-2,j,i) + a1 * w(k-1,j,i) + &
                          a0 * w(k,j,i) + a1 * w(k+1,j,i) + &
                          a2 * w(k+2,j,i)) / dzz

      elseif (k > 1 .and. k < nz) then
         d2udz2(k,j,i) = (u(k+1,j,i) - 2.d0 * u(k,j,i) + u(k-1,j,i)) / dzz
         d2vdz2(k,j,i) = (v(k+1,j,i) - 2.d0 * v(k,j,i) + v(k-1,j,i)) / dzz
         d2wdz2(k,j,i) = (w(k+1,j,i) - 2.d0 * w(k,j,i) + w(k-1,j,i)) / dzz

      elseif (k == 1) then
         d2udz2(k,j,i) = (c0 * u(k,j,i) + c1 * u(k+1,j,i) + &
                          c2 * u(k+2,j,i) + c3 * u(k+3,j,i) + &
                          c4 * u(k+4,j,i) + c5 * u(k+5,j,i) + &
                          c6 * u(k+6,j,i) + c7 * u(k+7,j,i)) / dzz
         d2vdz2(k,j,i) = (c0 * v(k,j,i) + c1 * v(k+1,j,i) + &
                          c2 * v(k+2,j,i) + c3 * v(k+3,j,i) + &
                          c4 * v(k+4,j,i) + c5 * v(k+5,j,i) + &
                          c6 * v(k+6,j,i) + c7 * v(k+7,j,i)) / dzz
         d2wdz2(k,j,i) = (c0 * w(k,j,i) + c1 * w(k+1,j,i) + &
                          c2 * w(k+2,j,i) + c3 * w(k+3,j,i) + &
                          c4 * w(k+4,j,i) + c5 * w(k+5,j,i) + &
                          c6 * w(k+6,j,i) + c7 * w(k+7,j,i)) / dzz

      else
         d2udz2(k,j,i) = (c0 * u(k,j,i) + c1 * u(k-1,j,i) + &
                          c2 * u(k-2,j,i) + c3 * u(k-3,j,i) + &
                          c4 * u(k-4,j,i) + c5 * u(k-5,j,i) + &
                          c6 * u(k-6,j,i) + c7 * u(k-7,j,i)) / dzz
         d2vdz2(k,j,i) = (c0 * v(k,j,i) + c1 * v(k-1,j,i) + &
                          c2 * v(k-2,j,i) + c3 * v(k-3,j,i) + &
                          c4 * v(k-4,j,i) + c5 * v(k-5,j,i) + &
                          c6 * v(k-6,j,i) + c7 * v(k-7,j,i)) / dzz
         d2wdz2(k,j,i) = (c0 * w(k,j,i) + c1 * w(k-1,j,i) + &
                          c2 * w(k-2,j,i) + c3 * w(k-3,j,i) + &
                          c4 * w(k-4,j,i) + c5 * w(k-5,j,i) + &
                          c6 * w(k-6,j,i) + c7 * w(k-7,j,i)) / dzz

      endif

      enddo
      enddo
      enddo
      !$omp end do

   else

      stop

   endif

   !$omp end parallel

   return

end subroutine laplace_wind
