#include "abi_common.h"
module utests
    implicit none

contains
    subroutine test_complex_svd()
        use m_scdm_math, only: real_svd, complex_svd
        integer, parameter :: M = 6, N = 5, LWMAX = 1000
        complex(8):: A(M, N), U(M, M), VT(N, N)
        real(8) :: S(N)
        DATA A/ &
            8.79, 6.11, -9.15, 9.57, -3.49, 9.84, &
            9.93, 6.91, -7.93, 1.64, 4.02, 0.15, &
            9.83, 5.04, 4.86, 8.83, 9.80, -8.99, &
            5.45, -0.27, 4.85, 0.74, 10.00, -6.02, &
            3.16, 7.98, 3.01, 5.80, 4.27, -5.31 &
            /

        call complex_svd(A, U, S, VT, 'A')
        print *, "U", U
        print *, "S", S
    end subroutine test_complex_svd

    subroutine test_QRPiv()
        use m_scdm_math, only: complex_QRCP_piv_only
        integer, parameter :: M = 6, N = 5
        complex(8):: A(5, 6)
        complex(8):: AT(6, 5)
        integer::  Piv(5)
        integer :: d
        DATA A/ &
            8.79, 6.11, -9.15, 9.57, -3.49, 9.84, &
            9.93, 6.91, -7.93, 1.64, 4.02, 0.15, &
            9.83, 5.04, 4.86, 8.83, 9.80, -8.99, &
            5.45, -0.27, 4.85, 0.74, 10.00, -6.02, &
            3.16, 7.98, 3.01, 5.80, 4.27, -5.31 &
            /

        AT = transpose(A)
        call complex_QRCP_piv_only(AT, Piv)

        d = sum((Piv - [3, 1, 4, 5, 2])**2)

        if (d /= 0) then
            ABI_ERROR("complex QRCP piv only result is wrong")
        end if

    end subroutine test_QRPiv
end module utests

program utest
    use utests
    use m_tb_scdm, only: test_tb_scdm
    use m_kpoints, only: test_kpoints

    print *, "testing complex QRpiv"
    call test_QRPiv()


    print *, "testing kpoints"
    call test_kpoints()

    print *, "testing tight binding"
    call test_tb_scdm()

end program utest

