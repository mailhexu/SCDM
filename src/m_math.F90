module m_math
  implicit none
  integer, parameter, public :: dp = 8
  real(dp), parameter, public :: PI=3.1415926535
  real(dp), parameter, public :: tpi=2*PI
  complex(dp), parameter, public :: tpi_im = cmplx(0.0_dp, tpi, kind=dp)

  public :: complex_QRCP_piv_only
  public :: real_svd
  public :: complex_svd
  public :: gaussian
  public :: fermi
  private

contains

  subroutine complex_QRCP_Piv_only(A, Piv)
    complex(dp), intent(in) :: A(:, :)
    integer, intent(inout) :: Piv(:)
    real(dp) :: rwork(size(A,2)*2)
    complex(dp) :: tau(min(size(A,1), size(A,2)))
    integer ::m, n
    complex(dp), allocatable :: work(:)
    complex(dp) :: tmp_work(2)
    integer :: lwork
    integer :: info
    EXTERNAL ZGEQP3

    m=size(A, 1)
    n=size(A, 2)
    rwork(:)=0.0_dp
    call ZGEQP3(m,n,A,m,piv,tau,tmp_work,-1,rwork,info)
    if(info/=0) then
       print *, "Error in doing QRCP"
    endif
    Piv(:) = 0
    rwork(:) = 0.0_DP
    lwork = INT(AINT(REAL(tmp_work(1))))

    ALLOCATE(work(lwork))
    work(:) = (0.0_DP,0.0_DP)
    print *, "work allocated"

    CALL ZGEQP3(m,n,A,m,piv,tau,work,lwork,rwork,info)
    if(info/=0) then
       print *, "Error in doing QRCP"
    endif
    DEALLOCATE(work)
  end subroutine complex_QRCP_Piv_only



  subroutine real_svd(A,  U, S, VT)
    real(dp), intent(in) :: A(:, :)
    real(dp), intent(inout) :: U(:, :), S(:), VT(:,:)
    integer:: LWMAX
    real(dp) :: tmp(5)
    real(dp), allocatable::  WORK(:)
    integer  :: M, N
    integer::         LDA, LDU, LDVT
    integer::         INFO, LWORK
    EXTERNAL         DGESVD
    M=size(A,1)
    N=size(A, 2)
    LDA = M
    LDU = M
    LDVT = N
    LWORK = -1
    LWMAX= max(size(A,1), size(A,2))*10
    CALL DGESVD( 'All', 'All', M, N, A, LDA, S, U, LDU, VT, LDVT,&
         & tmp, LWORK, INFO )
    LWORK = MIN( LWMAX, INT( tmp( 1 ) ) )

    allocate(work(lwork))
    CALL DGESVD( 'All', 'All', M, N, A, LDA, S, U, LDU, VT, LDVT, &
         WORK, LWORK, INFO )
    IF( INFO.GT.0 ) THEN
       WRITE(*,*)'The algorithm computing SVD failed to converge.'
       STOP
    END IF
    deallocate(work)
  end subroutine real_svd

  subroutine complex_svd(A,  U, S, VT)
    complex(dp), intent(in) :: A(:, :)
    complex(dp), intent(inout) :: U(:, :),  VT(:,:)
    real(dp), intent(inout) :: S(:)
    integer:: LWMAX
    complex(dp) :: tmp(2)
    complex(dp), allocatable::  WORK(:)
    real(dp), allocatable :: rwork(:)
    integer  :: M, N
    integer::         LDA, LDU, LDVT
    integer::         INFO, LWORK
    EXTERNAL         ZGESVD

    M=size(A,1)
    N=size(A,2)
    LWMAX= max(size(A,1), size(A,2))*10
    LDA = M
    LDU = M
    LDVT = N
    LWORK = -1
    CALL ZGESVD( 'All', 'All', M, N, A, LDA, S, U, LDU, VT, LDVT,&
         & tmp, LWORK, rwork, INFO )
    LWORK = MIN( LWMAX, INT( tmp( 1 ) ) )
    allocate(work(lwork))
    allocate(rwork(min(M, N)*6))
    CALL ZGESVD( 'All', 'All', M, N, A, LDA, S, U, LDU, VT, LDVT, &
         WORK, LWORK, rwork,INFO )
    IF( INFO.GT.0 ) THEN
       WRITE(*,*)'The algorithm computing SVD failed to converge.'
       STOP
    END IF
    deallocate(work)
    deallocate(rwork)
  end subroutine complex_svd


  function gaussian(x, mu, sigma) result(y)
    real(dp), intent(in) :: x, mu, sigma
    real(dp) :: y
    y=exp(-1.0 * (x - mu)**2 / sigma**2)
  end function gaussian

  !===============================================================
  ! Complementary error function.
  !===============================================================
  pure function erfcd(x) result(y)
    real(dp), intent(in):: x
    real(dp) :: y
    real(dp):: t,z
    z = abs(x)
    t = 1.0 / ( 1.0 + 0.5 * z )

    y = t * exp( -z * z - 1.26551223 + t *   &
         ( 1.00002368 + t * ( 0.37409196 + t * &
         ( 0.09678418 + t * (-0.18628806 + t * &
         ( 0.27886807 + t * (-1.13520398 + t * &
         ( 1.48851587 + t * (-0.82215223 + t * 0.17087277 )))))))))
    if ( x.lt.0.0 ) y= 2.0 - y

  end function erfcd

  function fermi(x, mu, sigma) result(y)
    real(dp), intent(in) :: x, mu, sigma
    real(dp) :: y
    y=0.5 * erfc((x - mu) / sigma)
  end function fermi

end module m_math
