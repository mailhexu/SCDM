module m_math
  implicit none
  integer, parameter :: dp = 8
  real(dp), parameter :: PI=3.1415926535
  real(dp), parameter :: tpi=2*PI
  complex(dp), parameter :: tpi_im = cmplx(0.0_dp, tpi, kind=dp)

  public :: real_svd
  public :: complex_svd
  public ::
  private

contains

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


end module m_math
