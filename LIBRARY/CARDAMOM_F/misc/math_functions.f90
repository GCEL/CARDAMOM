!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CARbon DAta MOdel fraMework (CARDAMOM) and DALEC terrestrial ecosystem model suite
! CARDAMOM is a Bayesian model-data fusion software framework. CARDAMOM is used to
! assimilate observations and ecological theory to retrieve parameters for the
! DALEC suite of intermediate complexity terrestrial ecosystem models. DALEC can be
! used as a fully integrated component of CARDAMOM or independently.
! Copyright (C) 2024  University of Edinburgh,
!                     Mathew Williams (mat.williams@ed.ac.uk),
!                     T. Luke Smallman (t.l.smallman@ed.ac.uk)
! UoE = University of Edinburgh

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see < https://www.gnu.org/licenses/>.

!!!!!!!!!!!! File specific description !!!!!!!!!!
! Module contains functions needed to mathematical calculations in CARDAMOM
!
! This code is based on the original C verion of the University of Edinburgh
! CARDAMOM framework created by A. A. Bloom (now at the Jet Propulsion Laboratory).
! All code translation into Fortran, integration into the University of
! Edinburgh CARDAMOM code and subsequent modifications by:
! T. L. Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)
! J. F. Exbrayat (University of Edinburgh)
! See function/subroutine specific comments for exceptions and contributors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module math_functions

   ! assume default private
   private

   ! make explicit bits we want others to see
   public:: linear_model_gradient

contains
   !
   !--------------------------------------------------------------------
   !
   subroutine matrix_vector_func(uplo, n, alpha, A, lda, X, incx, beta, Y, incy)

      ! Performs the matrix-vector operation
      ! y := alpha*A*x+beta*y,
      ! where alpha and beta are scalars, x and y are n element vectors and
      ! A is an n by n symmetric matrix.

      !  Arguments:
      !  ==========
      !
      ! intent(in):: UPLO
      !          UPLO is CHARACTER*1
      !          On entry, UPLO specifies whether the upper or lower
      !          triangular part of the array A is to be referenced as
      !          follows:
      !
      !              UPLO = 'U' or 'u'   Only the upper triangular part of A
      !                                  is to be referenced.
      !
      !              UPLO = 'L' or 'l'   Only the lower triangular part of A
      !                                  is to be referenced.!
      !
      ! intent(in):: N
      !           N is INTEGER
      !           On entry, N specifies the order of the matrix A.
      !           N must be at least zero.
      !
      ! intent(in):: ALPHA
      !           ALPHA is DOUBLE PRECISION.
      !           On entry, ALPHA specifies the scalar alpha.
      !
      ! intent(in):: A
      !           A is DOUBLE PRECISION array, dimension ( LDA, N )
      !           Before entry with  UPLO = 'U' or 'u', the leading n by n
      !           upper triangular part of the array A must contain the upper
      !           triangular part of the symmetric matrix and the strictly
      !           lower triangular part of A is not referenced.
      !           Before entry with UPLO = 'L' or 'l', the leading n by n
      !           lower triangular part of the array A must contain the lower
      !           triangular part of the symmetric matrix and the strictly
      !           upper triangular part of A is not referenced.
      !
      ! intent(in):: LDA
      !           LDA is INTEGER
      !           On entry, LDA specifies the first dimension of A as declared
      !           in the calling (sub) program. LDA must be at least
      !           max( 1, n ).
      !
      ! intent(in):: X
      !           X is DOUBLE PRECISION array, dimension at least
      !           ( 1 + ( n-1 )*abs( INCX ) ).
      !           Before entry, the incremented array X must contain the n
      !           element vector x.
      !
      ! intent(in):: INCX
      !           INCX is INTEGER
      !           On entry, INCX specifies the increment for the elements of
      !           X. INCX must not be zero.
      !
      ! intent(in):: BETA
      !           BETA is DOUBLE PRECISION.
      !           On entry, BETA specifies the scalar beta. When BETA is
      !           supplied as zero then Y need not be set on input.
      !
      ! intent(inout):: Y
      !           Y is DOUBLE PRECISION array, dimension at least
      !           ( 1 + ( n-1 )*abs( INCY ) ).
      !           Before entry, the incremented array Y must contain the n
      !           element vector y. On exit, Y is overwritten by the updated
      !           vector y.
      !
      ! intent(in):: INCY
      !           INCY is INTEGER
      !           On entry, INCY specifies the increment for the elements of
      !           Y. INCY must not be zero.
      !
      !  Authors:
      !  ========
      !
      ! author Univ. of Tennessee
      ! author Univ. of California Berkeley
      ! author Univ. of Colorado Denver
      ! author NAG Ltd.
      !
      ! Date: December 2016
      !
      ! Further Details:
      ! =====================
      !
      !  Level 2 Basic Linear Algebra Subprograms (BLAS) routines.
      !  The vector and matrix arguments are not referenced when N = 0, or M = 0
      !
      !  -- Written on 22-October-1986.
      !     Jack Dongarra, Argonne National Lab.
      !     Jeremy Du Croz, Nag Central Office.
      !     Sven Hammarling, Nag Central Office.
      !     Richard Hanson, Sandia National Labs.
      !
      !  -- Code modifed for inclusion into CARDAMOM 27th June 2019
      !     T. Luke Smallman (UoE; t.l.smallman@ed.ac.uk)
      !
      !  =====================================================================
      !
      !  -- Reference BLAS level2 routine (version 3.7.0) --
      !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
      !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      !     December 2016

      ! Arguments
      double precision, intent(in):: alpha, beta
      integer, intent(in):: incx, incy, lda, n
      character, intent(in):: uplo
      double precision, intent(in):: X(n)
      double precision, intent(in):: A(lda, n)
      double precision, intent(inout):: Y(n)

      !
      !  =====================================================================
      !

      ! local variables
      double precision, parameter:: one = 1d0, zero = 0d0
      double precision:: tmp1, tmp2
      integer:: i, info, ix, iy, j, jx, jy, kx, ky

      ! Test the input parameters.

      info = 0
      if (trim(uplo) /= 'U' .and. trim(uplo) /= 'L') THEN
         info = 1
      else if (n <= 0) then
         info = 2
      else if (lda < max(1, n)) then
         info = 5
      else if (incx == 0) then
         info = 7
      else if (incy == 0) then
         info = 10
      end if
      if (info /= 0) then
         print *, "Inputs to matrix_vector_func not correct-error code = ", info
         stop
      end if

      !     Quick return if possible.
      if (alpha == zero .and. beta == one) return

      !
      !     Set up the start points in  X  and  Y.
      !
      if (incx > 0) then
         kx = 1
      else
         kx = 1 - (n - 1)*incx
      end if
      if (incy > 0) then
         ky = 1
      else
         ky = 1 - (n - 1)*incy
      end if

      !     Start the operations. In this version the elements of A are
      !     accessed sequentially with one pass through the triangular part
      !     of A.
      !
      !     First form  y := beta*y.

      if (beta /= one) then
         if (incy == 1) then
            if (beta == zero) then
               do i = 1, n
                  y(i) = zero
               end do
            else
               DO i = 1, n
                  y(i) = beta*y(i)
               end do
            end if
         else
            iy = ky
            if (beta == zero) then
               do i = 1, n
                  y(iy) = zero
                  iy = iy + incy
               end do
            else
               do i = 1, n
                  y(iy) = beta*y(iy)
                  iy = iy + incy
               end do
            end if
         end if
      end if
      ! have done beta effect, no more change will occur if alpha is zero therefore return
      if (alpha == zero) return

      ! Now apply alpha component on calculation

      if (trim(uplo) == 'U') then

         ! Form  y  when A is stored in upper triangle.

         if ((incx .eq. 1) .and. (incy .eq. 1)) then

            do j = 1, n
               tmp1 = alpha*x(j)
               tmp2 = zero
               do i = 1, j - 1
                  y(i) = y(i) + tmp1*A(i, j)
                  tmp2 = tmp2 + a(i, j)*x(i)
               end do
               y(j) = y(j) + tmp1*A(j, j) + alpha*tmp2
            end do

         else
            jx = kx
            jy = ky
            do j = 1, n
               tmp1 = alpha*x(jx)
               tmp2 = zero
               ix = kx
               iy = ky
               do i = 1, j - 1
                  y(iy) = y(iy) + tmp1*a(i, j)
                  tmp2 = tmp2 + A(i, j)*x(ix)
                  ix = ix + incx
                  iy = iy + incy
               end do
               y(jy) = y(jy) + tmp1*A(j, j) + alpha*tmp2
               jx = jx + incx
               jy = jy + incy
            end do
         end if

      else  ! assume we must be in lower triangle

         ! Form  y  when A is stored in lower triangle.

         if (incx .eq. 1 .and. incy .eq. 1) then
            do j = 1, n
               tmp1 = alpha*x(j)
               tmp2 = zero
               y(j) = y(j) + tmp1*A(j, j)
               do i = j + 1, n
                  y(i) = y(i) + tmp1*A(i, j)
                  tmp2 = tmp2 + A(i, j)*x(i)
               end do
               y(j) = y(j) + alpha*tmp2
            end do
         else  ! incx .eq. 1 .and. incy .eq. 1
            jx = kx
            jy = ky
            do j = 1, n
               tmp1 = alpha*x(jx)
               tmp2 = zero
               y(jy) = y(jy) + tmp1*A(j, j)
               ix = jx
               iy = jy
               do i = j + 1, n
                  ix = ix + incx
                  iy = iy + incy
                  y(iy) = y(iy) + tmp1*A(i, j)
                  tmp2 = tmp2 + A(i, j)*x(ix)
               end do
               y(jy) = y(jy) + alpha*tmp2
               jx = jx + incx
               jy = jy + incy
            end do

         end if  ! incx .eq. 1 .and. incy .eq. 1

      end if ! (trim(uplo) == 'U')

      ! now return finally to the user
      return

   end subroutine matrix_vector_func
   !
   !--------------------------------------------------------------------
   !
   double precision function linear_model_gradient(x, y, interval)

      ! Function to calculate the gradient of a linear model for a given depentent
      ! variable (y) based on predictive variable (x). The typical use of this
      ! function will in fact be to assume that x is time.

      implicit none

      ! declare input variables
      integer:: interval  ! the total number of variables being regressed over
      double precision, dimension(interval):: x, y

      ! declare local variables
      double precision:: sum_x, sum_y

      ! calculate the sum of x
      sum_x = sum(x)
      ! calculate the sum of y
      sum_y = sum(y)
      ! calculate the sum of squares of x
      !sumsq_x = sum(x*x)
      ! calculate the sum of the product of xy
      !sum_product_xy = sum(x*y)
      ! calculate the gradient
      !linear_model_gradient = ( (dble(interval)*sum_product_xy) - (sum_x*sum_y) ) &
      !                      / ( (dble(interval)*sumsq_x) - (sum_x*sum_x) )
      ! Linear regression done as single line to reduce assignment requirements
      linear_model_gradient = ((dble(interval)*sum(x*y)) - (sum_x*sum_y)) &
                              /((dble(interval)*sum(x*x)) - (sum_x*sum_x))

      ! for future reference here is how to calculate the intercept
!    intercept = ( (sum_y*sumsq_x) - (sum_x*sum_product_xy) ) &
!              / ( (dble(interval)*sumsq_x) - (sum_x*sum_x) )

      ! don't forget to return to the user
      return

   end function linear_model_gradient
end module math_functions
