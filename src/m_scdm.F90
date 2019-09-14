

module m_scdm
  implicit none
  integer, parameter :: dp = 8
  real(dp), parameter :: PI=3.1415926535
  real(dp), parameter :: tpi=2*PI
  complex(dp), parameter :: tpi_im = cmplx(0.0_dp, tpi, kind=dp)
  private

  type scdmk
     real(dp), allocatable :: evals(:, :)   !(iband, ikpt)
     complex(dp), allocatable :: wfn(:,:,:) ! (ibasis, iband, ikpt)
     real(dp), allocatable :: kpts(:,:) !(idim, ikpt)
     real(dp), allocatable :: positions_red(:,:) !(idim, ibasis)
     real(dp), allocatable :: weight(:,:) !(iband, ikpt)
     integer :: nwann, nkpt, nband, nbasis
     integer :: dim ! dimension of position
     integer :: ik_anchor
   contains
     procedure :: initialize
     procedure :: finalize
     procedure :: compute_Amn ! Amn for one k-point
     procedure :: compute_Amnk ! Amnk for all kpoints
     procedure :: dephase_all  ! 
     procedure :: get_anchor_projections
     procedure :: get_weight
     procedure :: select_column
  end type scdmk


contains

  subroutine initialize(self, evals, wfn, kpts, positions_red)
    class(scdmk), intent(inout) :: self
    real(dp), allocatable :: evals(:, :)   !(iband, ikpt)
    complex(dp), allocatable :: wfn(:,:,:) ! (ibasis, iband, ikpt)
    real(dp), allocatable :: kpts(:,:) !(idim, ikpt)
    real(dp), allocatable :: positions_red(:,:) !(idim, ibasis)

    ! allocate memories
  end subroutine initialize

  subroutine finalize(self)
    class (scdmk), intent(inout) :: self
    ! deallocate
  end subroutine finalize

  subroutine compute_Amn(self, Amn)
    class(scdmk), intent(inout) :: self
    complex(dp), intent(inout) :: Amn(self%nbasis, self%nwann)
  end subroutine compute_Amn


  subroutine compute_Amnk(self, Amnk)
    class(scdmk), intent(inout) :: self
    complex(dp), intent(inout) :: Amnk(self%nbasis, self%nwann, self%nkpt)
    
  end subroutine compute_Amnk

  subroutine scdm(psi, cols, Amn)
    complex(dp), intent(in) :: psi(:,:)
    integer :: cols(:)
    complex(dp), intent(inout) :: Amn
    complex(dp), allocatable :: U(:,:), S(:,:), V(:,:)
    !!call real_svd(conjg(transpose(psi)), U, S, V)
    !Amn(:,:) = matmul(U, conjg(transpose(V)))
  end subroutine scdm

  !-------------------------------------------------------------------!
  !psi->psi_k * exp(-2pi k.r)
  !-------------------------------------------------------------------!

  subroutine dephase_all(self, dephased_psi)
    class(scdmk), intent(in) :: self
    integer :: ik, iband, ipos
    real(dp) :: pos(self%dim)
    complex(dp) :: phase(self%nbasis)
    complex(dp), intent(inout) :: dephased_psi(self%nbasis, self%nband, self%nkpt)
    do ik =1, self%nkpt
       do ipos =1, self%nbasis
          phase(ipos) = exp(tpi_im * dot_product(self%kpts(:,ik), self%positions_red(:,ipos)))
       end do
       do iband=1, self%nband
         do ipos=1, self%nbasis 
            dephased_psi(ipos,iband, ik)= self%wfn(ipos,iband, ik)*phase(ipos)
         end do
      end do
   end do
  end subroutine dephase_all


  subroutine get_projection_to_anchors(self, anchor_ikpts, anchor_ibands)
    class(scdmk), intent(inout) :: self
    integer, intent(in) :: anchor_ikpts(:), anchor_ibands(:)
    integer :: i, n, ikpt, iband
    complex(dp) :: dephased_psi(self%nbasis, self%nband, self%nkpt)
    n=size(anchor_ikpts)
    call self%dephase_all(dephased_psi)
    do i=1, n
       do ikpt=1, self%nkpt
          do iband=1, self%nband
             self%weight(iband, ikpt) =abs(self%weight(iband, ikpt) + &
                  & dot_product(conjg(dephased_psi(:,iband, ikpt)), dephased_psi(:, anchor_ibands(i), anchor_ikpts(i))))
          end do
       end do
    end do
  end subroutine get_projection_to_anchors

  subroutine get_weight(self)
    class(scdmk), intent(inout) :: self
  end subroutine get_weight

  subroutine Amn_to_H(Amn, psi, H0, H, nbasis, nwann)
    complex(dp), intent(in) :: Amn(:,:), psi(:,:), H0(:,:)
    complex(dp), intent(inout) :: H(:,:)
    integer, intent(in) :: nbasis, nwann
    complex(dp) :: tmp(nbasis, nwann)
    tmp(:,:)=matmul(psi, Amn)
    H=matmul(matmul(transpose(conjg(tmp)), H0), tmp)
  end subroutine Amn_to_H

end module m_scdm

