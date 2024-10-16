
module math_functions

  !!!!!!!!!!!
  ! Module contains functions needed to mathematical calculations in CARDAMOM
  !!!!!!!!!!!

  implicit none

  ! assume default private
  private

  ! make explicit bits we want others to see
  public:: randn, std, idum, covariance_matrix, &
            random_normal, random_uniform, rnstrt, &
            random_multivariate, increment_covariance_matrix, &
            par2nor, nor2par, log_par2nor, log_nor2par, &
            cholesky_factor, inverse_matrix, matrix_vector_func, &
            calculate_variance, increment_variance

  !!!!!!!!!!!
  ! Authorship contributions
  !
  ! This code is based on the original C verion of the University of Edinburgh
  ! CARDAMOM framework created by A. A. Bloom (now at the Jet Propulsion Laboratory).
  ! All code translation into Fortran, integration into the University of
  ! Edinburgh CARDAMOM code and subsequent modifications by:
  ! T. L. Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)
  ! J. F. Exbrayat (University of Edinburgh)
  ! See function/subroutine specific comments for exceptions and contributors
  !!!!!!!!!!!

  !!!!!!!!!!!
  ! Subroutines rand(), narray() and rnstrt() are from:
  !!!!!!!!!!!

  ! Code converted using TO_F90 by Alan Miller
  ! Date: 2000-09-10  Time: 16:37:48
  ! Latest revision-16 January 2003

  ! FORTRAN 77 version of "ran_array"
  ! from Seminumerical Algorithms by D E Knuth, 3rd edition (1997)
  !       including the MODIFICATIONS made in the 9th printing (2002)
  ! ********* see the book for explanations and caveats! *********
  ! Author: Steve Kifowit
  ! http://ourworld.compuserve.com/homepages/steve_kifowit
  ! with modifications by Alan Miller to rnarry and rnstrt based upon
  ! Knuth's code.

  ! For Donald Knuth's Fortran 77 versions, go to:
  ! http://www-cs-faculty.stanford.edu/~knuth/programs
  ! Look for frng.f and frngdb.f

  ! NOTE: that minimum number of values to be returned is 100
  !!!!!!!!!!!

  !!!!!!!!!!!
  ! Function randn()
  !!!!!!!!!!!

  ! Code is a modified version of that found in the uniform distribution generator from Numerical receipes
  ! Modified to give 0-1 unform on input of value 0 and normal distribution (mean = 0 sd = 1) on input of 1
  ! Modified by TLS

  !!!!!!!!!!!

  ! randn() related seed value
  double precision:: idum

  ! rand(), narray(), rnstrt() related
  integer, parameter  :: kk = 100, ll = 37, mm = 2**30, tt = 70, kkk = kk+kk-1
  integer, save       :: ranx(kk)

  contains
  !
  !--------------------------------------------------------------------
  !
  subroutine calculate_variance(sample, mean_par, naccepted, variance)

    ! Subroutine to estimate the sample variance
    ! Var = Σ ( Xi-X )*2 / (N-1)
    ! X = mean for parameter
    ! Xi = ith member of the vector
    ! N = number of parameters accepted so far
    ! This code was based on CARDAMOM routines provided by A. A. Bloom, 
    ! available at github.com/CARDAMOM-framework/CARDAMOM_2.1.6c
    ! (contact abloom@jpl.nasa.gov for access)

    implicit none

    ! Arguments
    integer, intent(in):: naccepted
    double precision, intent(in):: sample(naccepted)
    double precision, intent(out):: mean_par, variance

    ! local variables
    integer:: i
    double precision, dimension(:), allocatable:: deviances

    ! allocate memory to local variable
    allocate(deviances(naccepted))

    ! calculate components needed for variance
    mean_par = sum(sample) / dble(naccepted)
    ! estimate deviance
    deviances = sample-mean_par
    ! estimate the variance
    ! NOTE: that naccepted-1 makes this the sample variance
    variance = sum(deviances*deviances) * dble(naccepted-1)**(-1)

    ! tidy up
    deallocate(deviances)

    ! return to user
    return

  end subroutine calculate_variance
  !
  !--------------------------------------------------------------------
  !
  subroutine increment_variance(sample, mean_par, cur, new, variance)

    ! Subroutine for incremental update of the variance
    ! CMOUT = CM*(N-1)/(N-1+ar) + (N*M'*M-(N+ar)*Mi'*Mi+x'*x*ar)/(N-1+ar)
    ! M  = mean vector for parameters
    ! Mi = new mean vector for updated variance_matrix
    ! ar = number of new parameters to be added
    ! N = number of parameters accepted so far
    ! This code was based on CARDAMOM routines provided by A. A. Bloom, 
    ! available at github.com/CARDAMOM-framework/CARDAMOM_2.1.6c
    ! (contact abloom@jpl.nasa.gov for access)
    ! Translation to fortran and subsequent modifications by T. L. Smallman
    ! University of Edinburgh, t.l.smallman@ed.ac.uk

    implicit none

    ! Arguments
    integer, intent(in):: new
    double precision, intent(in):: sample(new)
    double precision, intent(inout):: cur, mean_par, variance

    ! local variables
    integer:: n, i, j
    double precision:: new_mean_par, nnew

    nnew = 1d0
    ! loop through each accepted parameter set...
    do n = 1, new
       ! ...estimate the new mean value for each parameter...
       new_mean_par = ((mean_par*cur) + (sample(n) * nnew)) &
                    / (cur+nnew)
       ! ...update the variance with each new parameter vector in turn
       variance = variance*(cur-1d0)/(cur-1d0+nnew) &
                + (cur*mean_par*mean_par- &
                  (cur+nnew)*new_mean_par*new_mean_par + &
                   nnew*sample(n)*sample(n))/(cur-1d0+nnew)
       ! update running totals and mean for the next iteration
       cur = cur+1; mean_par = new_mean_par
    end do  ! new_accepted

    ! return to user
    return

  end subroutine increment_variance
  !
  !--------------------------------------------------------------------
  !
  subroutine covariance_matrix(PARSALL, mean_par, npars, naccepted, covariance)

    ! Subroutine to estimate the covariance matrix
    ! Cov(X, Y) = Σ ( Xi-X ) ( Yi-Y ) / (N-1)
    ! X = mean for parameter 1
    ! Y = mean parameter 2
    ! Xi = ith member of the vector
    ! N = number of parameters accepted so far
    ! This code was based on CARDAMOM routines provided by A. A. Bloom, 
    ! available at github.com/CARDAMOM-framework/CARDAMOM_2.1.6c
    ! (contact abloom@jpl.nasa.gov for access)
    ! Translation to fortran and subsequent modifications by T. L. Smallman
    ! University of Edinburgh, t.l.smallman@ed.ac.uk

    implicit none

    ! Arguments
    integer, intent(in):: npars, naccepted
    double precision, intent(in):: PARSALL(npars, naccepted)
    double precision, intent(out):: mean_par(npars), covariance(npars, npars)

    ! local variables
    integer:: i
    double precision, dimension(:,:), allocatable:: deviances

    ! allocate memory to local variable
    allocate(deviances(npars, naccepted))

    ! calculate components needed for covariance
    do i = 1, npars
       mean_par(i) = sum(PARSALL(i, :)) / dble(naccepted)
       ! estimate deviance
       deviances(i, :) = PARSALL(i, :) - mean_par(i)
    end do

    ! use matrix multiplication to estimate covariance
    ! NOTE: that naccepted-1 makes this the sample covariance
    covariance = matmul(deviances, transpose(deviances)) * dble(naccepted-1)**(-1)

    ! tidy up
    deallocate(deviances)

    ! return to user
    return

  end subroutine covariance_matrix
  !
  !--------------------------------------------------------------------
  !
  subroutine increment_covariance_matrix(PARSALL, mean_par, npars, cur, new, covariance)

    ! Subroutine for incremental update of a covariance matrix
    ! CMOUT = CM*(N-1)/(N-1+ar) + (N*M'*M-(N+ar)*Mi'*Mi+x'*x*ar)/(N-1+ar)
    ! M  = mean vector for parameters
    ! Mi = new mean vector for updated covariance_matrix
    ! ar = number of new parameters to be added
    ! N = number of parameters accepted so far
    ! This code was based on CARDAMOM routines provided by A. A. Bloom, 
    ! available at github.com/CARDAMOM-framework/CARDAMOM_2.1.6c
    ! (contact abloom@jpl.nasa.gov for access)

    implicit none

    ! Arguments
    integer, intent(in):: npars, new
    double precision, intent(in):: PARSALL(npars, new)
    double precision, intent(inout):: cur, mean_par(npars), covariance(npars, npars)

    ! local variables
    integer:: n, i, j
    double precision:: new_mean_par(npars), nnew

    nnew = 1d0
    ! loop through each accepted parameter set...
    do n = 1, new
       ! ...estimate the new mean value for each parameter...
       new_mean_par = ((mean_par*cur) + (PARSALL(:,n) * nnew)) &
                    / (cur+nnew)
       ! ...update the covariance matrix with each new parameter vector in turn
       do i = 1, npars
          do j = 1, npars
             covariance(i, j) = covariance(i, j)*(cur-1d0)/(cur-1d0+nnew) &
                             + (cur*mean_par(i)*mean_par(j)- &
                               (cur+nnew)*new_mean_par(i)*new_mean_par(j) + &
                                nnew*PARSALL(i, n)*PARSALL(j, n))/(cur-1d0+nnew)
          end do  ! j = 1, npars
       end do  ! i = 1, npars
       ! update running totals and mean for the next iteration
       cur = cur+1; mean_par = new_mean_par
    end do  ! new_accepted

    ! return to user
    return

  end subroutine increment_covariance_matrix
  !
  !--------------------------------------------------------------------
  !
  subroutine inverse_matrix(n, a, c)

     !============================================================
     ! Inverse for positive definite symmetric matrix
     ! Method: Based on Doolittle LU factorization for Ax = b
     ! Alex L. Godunov December 2009
     ! Modifed for CARDAMOM: T. Luke Smallman (June 2019)
     !                       t.l.smallman@ed.ac.uk
     ! Warning if matrix not positive definite this function will fail
     !-----------------------------------------------------------
     ! input ...
     ! a(n, n) - array of coefficients for matrix A
     ! n      - dimension
     ! output ...
     ! c(n, n) - inverse matrix of A
     ! comments ...
     ! the original matrix a(n, n) will be destroyed
     ! during the calculation
     !===========================================================

     implicit none

     ! arguments
     integer, intent(in):: n
     double precision, dimension(n, n), intent(in):: a
     double precision, dimension(n, n), intent(out):: c

     ! local arguments
     double precision, dimension(n, n):: L, U, a_local
     double precision, dimension(n):: b(n), d(n), x(n)
     double precision:: coeff
     integer:: i, j, k

     ! Step 0: initialization for matrices L and U and b
     ! Fortran 90/95 allows such operations on matrices
     L = 0d0; U = 0d0; b = 0d0; a_local = a

     ! Step 1: forward elimination
     do k = 1, n-1
        do i = k+1, n
           coeff = a_local(i, k)/a_local(k, k)
           L(i, k) = coeff
           do j = k+1, n
              a_local(i, j) = a_local(i, j)-coeff*a_local(k, j)
           end do
        end do
     end do

     ! Step 2: prepare L and U matrices

     ! L matrix is a matrix of the elimination coefficient
     ! + the diagonal elements are 1.0

     do i = 1, n
        L(i, i) = 1d0
     end do

     ! U matrix is the upper triangular part of A
     do j = 1, n
        do i = 1, j
           U(i, j) = a_local(i, j)
        end do
     end do

     ! Step 3: compute columns of the inverse matrix C
     do k = 1, n
        b(k) = 1d0
        d(1) = b(1)
        ! Step 3a: Solve Ld = b using the forward substitution
        do i = 2, n
           d(i) = b(i)
           do j = 1, i-1
              d(i) = d(i) - L(i, j)*d(j)
           end do
        end do
        ! Step 3b: Solve Ux = d using the back substitution
        x(n) = d(n)/U(n, n)
        do i = n-1, 1, -1
           x(i) = d(i)
           do j = n, i+1, -1
              x(i) = x(i)-U(i, j)*x(j)
           end do
           x(i) = x(i)/U(i, i)
        end do
        ! Step 3c: fill the solutions x(n) into column k of C
        do i = 1, n
           c(i, k) = x(i)
        end do
        b(k) = 0d0
     end do

     ! return to user
     return

  end subroutine inverse_matrix
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
    double precision, intent (in):: A(lda, n)
    double precision, intent (inout):: Y(n)

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
        print*,"Inputs to matrix_vector_func not correct-error code = ",info
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
        kx = 1 - (n-1)*incx
    end if
    if (incy > 0) then
        ky = 1
    else
        ky = 1 - (n-1)*incy
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
                    iy = iy+incy
                end do
            else
                do i = 1, n
                    y(iy) = beta*y(iy)
                    iy = iy+incy
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
              do i = 1, j-1
                 y(i) = y(i) + tmp1*A(i, j)
                 tmp2 = tmp2+a(i, j)*x(i)
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
              do i = 1, j-1
                 y(iy) = y(iy) + tmp1*a(i, j)
                 tmp2 = tmp2+A(i, j)*x(ix)
                 ix = ix+incx
                 iy = iy+incy
              end do
              y(jy) = y(jy) + tmp1*A(j, j) + alpha*tmp2
              jx = jx+incx
              jy = jy+incy
           end do
       end if

    else  ! assume we must be in lower triangle

    ! Form  y  when A is stored in lower triangle.

       if (incx .eq. 1 .and. incy .eq. 1) then
           do j = 1, n
              tmp1 = alpha*x(j)
              tmp2 = zero
              y(j) = y(j) + tmp1*A(j, j)
              do i = j+1, n
                 y(i) = y(i) + tmp1*A(i, j)
                 tmp2 = tmp2+A(i, j)*x(i)
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
              do i = j+1, n
                 ix = ix+incx
                 iy = iy+incy
                 y(iy) = y(iy) + tmp1*A(i, j)
                 tmp2 = tmp2+A(i, j)*x(ix)
              end do
              y(jy) = y(jy) + alpha*tmp2
              jx = jx+incx
              jy = jy+incy
           end do

       end if  ! incx .eq. 1 .and. incy .eq. 1

    end if ! (trim(uplo) == 'U')

    ! now return finally to the user
    return

  end subroutine matrix_vector_func
  !
  !--------------------------------------------------------------------
  !
  double precision function std(a, n)

    ! Function to determine the standard deviation
    ! inputs are the vector of values and number of values included

    implicit none

    ! declare inputs
    integer, intent(in):: n  ! number of values in vector
    double precision, dimension(n), intent(in):: a

    ! declare local variables
    integer:: i
    double precision:: mean, sq_diff_sum, diff, variance, sample

    ! multiple use variable
    sample = dble(n)
    ! first calculate the mean
    mean = sum(a) / sample
    ! ensure zero values
    sq_diff_sum = 0d0
    ! calculate cumulative square difference
    do i = 1, n
       diff = a(i)-mean
       sq_diff_sum = sq_diff_sum + (diff*diff)
    end do
    ! calculate the variance
    variance = sq_diff_sum / (sample-1d0)
    ! return the standard deviation
    std = sqrt(variance)

    return

  end function std
  !
  !------------------------------------------------------------------
  !
  subroutine par2nor(niter, initial_par, min_par, max_par, out_par)

    ! functions to normalised parameter values and return them back to
    ! un-normalised value.

    ! converting parameters on log scale between 0-1 for min/max values
    implicit none
    integer, intent(in):: niter     ! number of iterations in current vector
    double precision, intent(in):: min_par, max_par
    double precision, dimension(niter), intent(in):: initial_par
    double precision, dimension(niter), intent(out):: out_par

    ! then normalise
    out_par = (initial_par-min_par)/(max_par-min_par)

    ! explicit return
    return

  end subroutine par2nor
  !
  !---------------------and vise versa------------------------------
  !
  subroutine nor2par(niter, initial_par, min_par, max_par, out_par)

    ! Converting values back from normalised (0-1) to 'real' numbers

    implicit none
    integer, intent(in):: niter     ! number of iterations in current vector
    double precision, intent(in):: min_par, max_par
    double precision, dimension(niter), intent(in):: initial_par
    double precision, dimension(niter), intent(out):: out_par

    ! ...then un-normalise without logs as we cross zero and logs wont work
    out_par = min_par+(max_par-min_par)*initial_par

    ! explicit return
    return

  end subroutine nor2par
  !
  !------------------------------------------------------------------
  !
  subroutine log_par2nor(niter, initial_par, min_par, max_par, par_adj, out_par)

    ! Functions to normalised-log parameter values.

    ! Converting parameters on log scale between 0-1 for min/max values
    implicit none
    integer, intent(in):: niter     ! number of iterations in current vector
    double precision, intent(in):: min_par, max_par, par_adj
    double precision, dimension(niter), intent(in):: initial_par
    double precision, dimension(niter), intent(out):: out_par
    ! local values
    double precision:: invar(niter), minvar, maxvar

    ! Assign inputs to the local variables
    invar = initial_par+par_adj
    minvar = min_par+par_adj
    maxvar = max_par+par_adj

    ! Then normalise
    !out_par = log(initial_par/min_par)/log(max_par/min_par)
    out_par = log(invar/minvar)/log(maxvar/minvar)

    ! explicit return
    return

  end subroutine log_par2nor
  !
  !---------------------and vise versa------------------------------
  !
  subroutine log_nor2par(niter, initial_par, min_par, max_par, par_adj, out_par)

    ! Converting values back from log-normalised (0-1) to 'real' numbers

    implicit none
    integer, intent(in):: niter     ! number of iterations in current vector
    double precision, intent(in):: min_par, max_par, par_adj  ! adjustment prevents negative values being fed into the analysis
    double precision, dimension(niter), intent(in):: initial_par
    double precision, dimension(niter), intent(out):: out_par
    ! local values
    double precision:: minvar, maxvar

    ! Assign inputs to the local variables
    minvar = min_par+par_adj
    maxvar = max_par+par_adj

    ! ...then un-normalise without logs as we cross zero and logs wont work
    !out_par = min_par*(max_par/min_par)**initial_par
    out_par = (minvar*(maxvar/minvar)**initial_par) - par_adj

    ! explicit return
    return

  end subroutine log_nor2par
  !
  !------------------------------------------------------------------
  !
  double precision function randn(option)

    ! From Numerical Recipes p271 Press et al., 1986 2nd Edition Chapter 7, 
    ! Random Numbers function returns real random number between 0-1 based on an initial start
    ! point (ran1). The start point (default = -1) is reinitialised every time the model runs
    ! providing the same distribution each run. To ensure random numbers each time use the
    ! sum of the current system time.

    ! This code was modified based on CARDAMOM routines provided by A. A. Bloom, 
    ! available at github.com/CARDAMOM-framework/CARDAMOM_2.1.6c
    ! (contact abloom@jpl.nasa.gov for access)

    implicit none
    integer:: IA, IM, IQ, IR, NTAB, NDIV, option
    double precision:: AM, EPS, RNMX, const, r1, r2, pi
    parameter(IA = 16807, IM = 2147483647, AM = 1d0/dble(IM), IQ = 127773, IR = 2836, NTAB = 32, NDIV = 1+(IM-1)/NTAB, EPS = 1.2d-30, RNMX = 1d0-EPS)
    integer:: j, k, iv(NTAB), iy
    SAVE iv, iy
    DATA iv/NTAB*0/
    DATA iy/0/

    const = 1d0
    pi = 3.141592653589793d0

    if (option == 0) then
        if (idum < 0d0 .or. iy == 0) then
            idum = max(-idum, const)
            do j = (NTAB+8), 1, -1
               k = nint(idum/dble(IQ))
               idum = dble(IA)*(idum-dble(k*IQ))-dble(IR*k)
               if (idum < 0d0) idum = idum+dble(IM)
               if (j < NTAB) iv(j) = nint(idum)
            enddo
            iy = iv(1)
        endif
        k = nint(idum/dble(IQ))
        idum = dble(IA)*(idum-dble(k*IQ))-dble(IR*k)
        if (idum < 0d0) idum = idum+dble(IM)
        j = 1+iy/NDIV
        iy = iv(j)
        iv(j) = nint(idum)

        ! output now
        randn = min(AM*dble(iy), RNMX)
        return

    else

        if (idum < 0d0 .or. iy == 0) then
            idum = max(-idum, const)
            do j = (NTAB+8), 1, -1
               k = nint(idum/dble(IQ))
               idum = dble(IA)*(idum-dble(k*IQ))-dble(IR*k)
               if (idum < 0d0) idum = idum+dble(IM)
               if (j < NTAB) iv(j) = nint(idum)
            enddo
            iy = iv(1)
        endif
        k = nint(idum)/IQ
        idum = dble(IA)*(idum-dble(k*IQ))-dble(IR*k)
        if (idum < 0d0) idum = idum+dble(IM)
        j = 1+iy/NDIV
        iy = iv(j)
        iv(j) = nint(idum)
        r1 = max(min(AM*dble(iy), RNMX), 1d-30)

        if (idum < 0d0 .or. iy == 0) then
            idum = max(-idum, const)
            do j = (NTAB+8), 1, -1
               k = nint(idum)/IQ
               idum = dble(IA)*(idum-dble(k*IQ))-dble(IR*k)
               if (idum < 0d0) idum = idum+dble(IM)
               if (j < NTAB) iv(j) = nint(idum)
            enddo
            iy = iv(1)
        endif
        k = nint(idum)/IQ
        idum = dble(IA)*(idum-dble(k*IQ))-dble(IR*k)
        if (idum < 0d0) idum = idum+dble(IM)
        j = 1+iy/NDIV
        iy = iv(j)
        iv(j) = nint(idum)
        r2 = max(min(AM*dble(iy), RNMX), 1d-30)

        ! output now
        randn = sqrt(-2d0*log(r1)) * cos(2d0*pi*r2)
        return

    endif

  end function randn
  !
  !--------------------------------------------------------------------
  !
  subroutine random_normal(uniform, random_length, uniform_random_vector, fn_val)

    ! Generate a random normal deviate using the polar method.
    ! Reference: Marsaglia, G. & Bray, T.A. 'A convenient method for
    ! generating normal variables',
    ! Siam Rev., vol.6, 260-264, 1964.
    ! Code modified from that created by Alan Miller
    ! (https://jblevins.org/mirror/amiller/rnorm.f90, last updated February 2004)

    implicit none

    ! arguments
    integer, intent(in):: random_length
    integer, intent(inout):: uniform
    double precision, intent(out):: fn_val
    double precision, dimension(random_length), intent(inout):: uniform_random_vector

    ! Local variables
    double precision:: u, sumsq
    double precision, save:: v, sln
    logical, save:: second = .false.
    double precision, parameter:: one = 1d0, vsmall = tiny( one )
    double precision, parameter:: sample_mean = 0d0, sample_std = 1d0

    if (second) then

        ! If second, use the second random number generated on last call
        second = .false.
        fn_val = sample_mean+sample_std*v * sln

    else

        ! First call; generate a pair of random normals
        second = .true. ; sumsq = 0d0
        do while (sumsq > one .or. sumsq == 0d0)
           u = uniform_random_vector(uniform); uniform = uniform+1
           v = uniform_random_vector(uniform); uniform = uniform+1
           u = u*2d0-one
           v = v*2d0-one
           sumsq = u*u + v*v
           if (uniform >= random_length) then
               call random_uniform(uniform_random_vector, random_length)
               uniform = 1
           endif
        end do
        sln = sqrt(-2d0*log(sumsq) / sumsq)
        fn_val = sample_mean+sample_std*u * sln
    end if

    return

  end subroutine random_normal
  !
  !--------------------------------------------------------------------
  !
  subroutine random_uniform(u, n)

    ! Generate an array of n double precision values between 0 and 1.
    ! from Seminumerical Algorithms by D E Knuth, 3rd edition (1997)
    !       including the MODIFICATIONS made in the 9th printing (2002)
    ! ********* see the book for explanations and caveats! *********
    ! Author: Steve Kifowit
    ! http://ourworld.compuserve.com/homepages/steve_kifowit
    ! with modifications by Alan Miller to rnarry and rnstrt based upon
    ! Knuth's code.
    ! Code converted using TO_F90 by Alan Miller
    ! Date: 2000-09-10, last update 16 January 2003
    ! Modified for integration into CARDAMOM by T. Luke Smallman (t.l.smallman@ed.ac.uk)
    ! 03/05/2019

    integer, intent(in)  :: n  ! number of random values wanted
    double precision, intent(out):: u(n)  ! output vector

    ! Local array
    integer, allocatable, dimension(:)  :: aa

    ! allocate memory
    allocate(aa(n))

    call rnarry(aa, n)
    u(1:n) = scale( dble(aa), -30)

    ! tidy
    deallocate(aa)

    return

  end subroutine random_uniform
  !
  !--------------------------------------------------------------------
  !
  subroutine random_multivariate ( m, n, u, u_size, uniform_random_vector, a, mu, x )

    !
    !  Discussion:
    !
    !    The multivariate normal distribution for the M dimensional vector X
    !    has the form:
    !
    !      pdf(X) = (2*pi*det(A))**(-M/2) * exp(-0.5*(X-MU)'*inverse(A)*(X-MU))
    !
    !    where MU is the mean vector, and A is a positive definite symmetric
    !    matrix called the variance-covariance matrix.
    !
    !  Licensing: This code is distributed under the GNU LGPL license.
    !
    !  Last Modified: Tue 08 Oct 2024 10:20:49 BST
    !
    !  Original Author: John Burkardt (07 December 2009)
    !
    !  Modified for coupling into CARDAMOM (03 May 2019):
    !    T. L. Smallman (t.l.smallman@ed.ac.uk)
    !
    !  Parameters:
    !
    !    Input, M, the dimension of the space.
    !    Input, N, the number of points.
    !    Input, A(M, M), the variance-covariance
    !    Matrix. A must be positive definite symmetric.
    !    Input,  MU(M), the mean vector.
    !
    !    Output, X(M), the points.
    !

    implicit none

    ! arguments
    integer, intent(in):: m, & ! number of parameters
                           n, & ! number of multivariate samples wanted per parameter
                      u_size    ! vector length for the uniform_random_vector
    integer, intent(inout):: u  ! current position within the random uniform vector
    double precision, intent(in):: a(m, m), mu(m)
    double precision, intent(out):: x(m, n)
    double precision, dimension(u_size), intent(inout):: uniform_random_vector

    ! local variables
    integer:: info, i, j
    double precision:: r(m, m)

    !
    !  Compute the upper triangular Cholesky factor R of the variance-covariance
    !  matrix.
    !

    ! Requires variance-covariance matrix, however as the matrix is over written, 
    ! make a duplicate
    r(1:m, 1:m) = a(1:m, 1:m)
    call cholesky_factor ( m, r, info )

    ! error checking from Cholesky Factor calculation
    if ( info /= 0 ) then

        ! Normal_multivariate-Fatal error!
        ! The variance-covariance matrix is not positive definite symmetric.
        ! Non-multivariate normal sample returned instead as default.'

        do i = 1, m  ! number of parameters searching for
           do j = 1, n  ! number of samples needed per parameter
              ! Draw random number from normal distribution
              call random_normal(u, u_size, uniform_random_vector, x(i, j))
              ! Updated variance based on variance from covariance matrix
              x(i, j) = x(i, j) * a(i, i)
           end do  ! j
        end do  ! i

        ! return early
        return

    end if

    !
    !  Get an MxN matrix of samples of the 1D normal distribution with mean 0
    !  and variance 1.
    !

    do i = 1, m
       do j = 1, n
          ! draw random number from normal distribution
          call random_normal(u, u_size, uniform_random_vector, x(i, j))
       end do  ! j
    end do  ! i

    !
    !  Compute R' * X.
    !  We actually carry out this computation in the equivalent form X' * R.
    !

    ! Whole multivariate estimated via matrix multplication
    ! Each desires set of combinations wanted it iterated
    do j = 1, n
       x(1:m, j) = mu(1:m) + matmul ( x(1:m, j), r(1:m, 1:m) )
    end do

    return

  end subroutine random_multivariate
  !
  !--------------------------------------------------------------------
  !
  subroutine cholesky_factor ( n, a, info )

    !
    !  Discussion:
    !
    !    Subroutine calculates the symmetric positive definite
    !    matrix and its inverse. The Cholesky factor is an
    !    upper triangular matrix.
    !
    !    Only the diagonal and upper triangle of the square array are used.
    !    For clarity, the lower triangle is set to zero.
    !
    !    The positive definite symmetric matrix A has a Cholesky factorization
    !    of the form:
    !
    !      A = R' * R
    !
    !    where R is an upper triangular matrix with positive elements on
    !    its diagonal.  This routine overwrites the matrix A with its
    !    factor R.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Last Modified: Tue 08 Oct 2024 10:20:49 BST
    !
    !    03/05/2019
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
    !    FORTRAN90 version by John Burkardt.
    !    Modified for coupling to CARDAMOM by T. L. Smallman (t.l.smallman@ed.ac.uk)
    !
    !  Reference:
    !
    !    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart, 
    !    LINPACK User's Guide, 
    !    SIAM, 1979, 
    !    ISBN13: 978-0-898711-72-1, 
    !    LC: QA214.L56.
    !
    !  Parameters:
    !
    !    Input, integer N, the order of the matrix.
    !
    !    Input/output, double precision A(N, N).
    !    On output, the Cholesky factor R.
    !
    !    Output, integer INFO, error flag.
    !    0, normal return.
    !    K, error condition. The principal minor of order K is not
    !    positive definite, and the factorization was not completed.
    !

    implicit none

    ! arguments
    integer, intent(in):: n
    integer, intent(out):: info
    double precision, intent(inout):: a(n, n)

    ! local arguments
    integer:: i, j, k
    double precision:: s

    ! Loop through the matrix along one dimension
    do j = 1, n

       ! doing the upper triangle only
       do k = 1, j-1
          a(k, j) = ( a(k, j) - sum ( a(1:k-1, k) * a(1:k-1, j) ) ) / a(k, k)
       end do

       s = a(j, j) - sum ( a(1:j-1, j)**2 )

       ! error checking
       if ( s <= 0.0D+00 ) then
           info = j
           return
       end if

       ! final calculation of Cholesky
       a(j, j) = sqrt ( s )

    end do  ! j = 1, n

    info = 0

    !
    !  Since the Cholesky factor is upper right corner only, be sure to
    !  zero out the lower triangle.
    !

    do i = 1, n
       do j = 1, i-1
          a(i, j) = 0.0D+00
       end do
    end do

    return

  end subroutine cholesky_factor
  !
  !--------------------------------------------------------------------
  !
  subroutine rnarry(aa, n)

    ! Generate an array of n integers between 0 and 2^30-1.
    ! Part of process to an array of n double precision values between 0 and 1.
    ! from Seminumerical Algorithms by D E Knuth, 3rd edition (1997)
    !       including the MODIFICATIONS made in the 9th printing (2002)
    ! ********* see the book for explanations and caveats! *********
    ! Author: Steve Kifowit
    ! http://ourworld.compuserve.com/homepages/steve_kifowit
    ! with modifications by Alan Miller to rnarry and rnstrt based upon
    ! Knuth's code.
    ! Code converted using TO_F90 by Alan Miller
    ! Date: 2000-09-10, last update 16 January 2003
    ! Modified for integration into CARDAMOM by T. Luke Smallman (t.l.smallman@ed.ac.uk)
    ! 03/05/2019

    integer, intent(in)   :: n
    integer, intent(out)  :: aa(n)

    ! Local variables
    integer  :: j

    aa(1:kk) = ranx(1:kk)
    do j = kk+1, n
       aa(j) = aa(j-kk) - aa(j-ll)
       if (aa(j) < 0) aa(j) = aa(j) + mm
    end do
    do j = 1, ll
       ranx(j) = aa(n+j-kk) - aa(n+j-ll)
       if (ranx(j) < 0) ranx(j) = ranx(j) + mm
    end do
    do j = ll+1, kk
       ranx(j) = aa(n+j-kk) - ranx(j-ll)
       if (ranx(j) < 0) ranx(j) = ranx(j) + mm
    end do

    return

  end subroutine rnarry
  !
  !--------------------------------------------------------------------
  !
  subroutine rnstrt(seed)

    ! Initialize integer array ranx using the input seed.
    ! Part of process to an array of n double precision values between 0 and 1.
    ! from Seminumerical Algorithms by D E Knuth, 3rd edition (1997)
    !       including the MODIFICATIONS made in the 9th printing (2002)
    ! ********* see the book for explanations and caveats! *********
    ! Author: Steve Kifowit
    ! http://ourworld.compuserve.com/homepages/steve_kifowit
    ! with modifications by Alan Miller to rnarry and rnstrt based upon
    ! Knuth's code.
    ! Code converted using TO_F90 by Alan Miller
    ! Date: 2000-09-10, last update 16 January 2003
    ! Modified for integration into CARDAMOM by T. Luke Smallman (t.l.smallman@ed.ac.uk)
    ! 03/05/2019

    integer, intent(in)  :: seed

    ! Local variables
    integer  :: x(kkk), j, ss, sseed, t

    if (seed < 0) then
        sseed = mm-1 - mod(-1-seed, mm)
    else
        sseed = mod(seed, mm)
    end if
    ss = sseed-mod(sseed, 2) + 2
    do j = 1, kk
       x(j) = ss
       ss = ishft(ss, 1)
       if (ss >= mm) ss = ss-mm+2
    end do
    x(kk+1:kkk) = 0
    x(2) = x(2)+1
    ss = sseed
    t = tt-1
10  do j = kk, 2, -1
       x(j+j-1) = x(j)
    end do
    do j = kkk, kk+1, -1
       x(j-(kk-ll)) = x(j-(kk-ll)) - x(j)
       if (x(j-(kk-ll)) < 0) x(j-(kk-ll)) = x(j-(kk-ll)) + mm
       x(j-kk) = x(j-kk) - x(j)
       if (x(j-kk) < 0) x(j-kk) = x(j-kk) + mm
    end do
    if (mod(ss, 2) == 1) then
        do j = kk, 1, -1
           x(j+1) = x(j)
        end do
        x(1) = x(kk+1)
        x(ll+1) = x(ll+1) - x(kk+1)
        if (x(ll+1) < 0) x(ll+1) = x(ll+1) + mm
    end if
    if (ss /= 0) THEN
        ss = ishft(ss, -1)
    else
        t = t-1
    end if
    if (t > 0) GO TO 10

    do j = 1, ll
       ranx(j+kk-ll) = x(j)
    end do
    do j = ll+1, kk
       ranx(j-ll) = x(j)
    end do

    do j = 1, 10
       call rnarry(x, kkk)
    end do

    return
  end subroutine rnstrt
  !
  !--------------------------------------------------------------------
  !
end module math_functions
