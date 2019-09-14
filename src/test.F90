program test
  use m_scdm , only: real_svd, complex_svd
  integer, parameter :: M=6, N=5, LWMAX=1000
  complex(8):: A( M, N ), U( M, M ), VT( N, N )
  real(8) :: S( N )
  DATA     A/                     &
    8.79, 6.11,-9.15, 9.57,-3.49, 9.84,  &
    9.93, 6.91,-7.93, 1.64, 4.02, 0.15, &
    9.83, 5.04, 4.86, 8.83, 9.80,-8.99, &
    5.45,-0.27, 4.85, 0.74,10.00,-6.02, &
    3.16, 7.98, 3.01, 5.80, 4.27,-5.31 &
                    /

  call complex_svd(A, U, S, VT)
  print *, "U", U
  print *, "S", S

end program test
