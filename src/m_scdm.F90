

!===============================================================
! SCDM-k
!> @description: Select Column Density matrix method for
!   generating Wannier functions.
!===============================================================
module m_scdm
  use m_math, only: complex_QRCP_piv_only, complex_svd, tpi_im, gaussian, fermi
  implicit none
  integer, parameter :: dp = 8
  real(dp), parameter :: PI=3.1415926535
  real(dp), parameter :: tpi=2*PI
  public :: scdmk
  public :: Amn_to_H
  private

  !===============================================================
  ! scdmk type:
  !> @ description: the class for scdmk method.
  !===============================================================
  type::  scdmk
     real(dp), pointer :: evals(:, :)=> null()   !(iband, ikpt)
     complex(dp), pointer :: psi(:,:,:)=>null() ! (ibasis, iband, ikpt)
     real(dp), allocatable :: kpts(:,:) !(idim, ikpt)
     real(dp), allocatable :: positions_red(:,:) !(idim, ibasis)
     real(dp), allocatable :: weight(:,:) !(iband, ikpt)
     integer, allocatable :: cols(:)
     real(dp), allocatable :: anchor_kpt(:)
     integer :: anchor_ikpt
     integer, allocatable :: anchor_ibands(:)
     integer :: nwann, nkpt, nband, nbasis, nkdim
     integer :: dim ! dimension of position
     integer :: ik_anchor
     integer :: disentangle_func_type
     real(dp) :: mu, sigma
     complex(dp), allocatable :: Amnk(:, :, :) !(nband, nwann, nkpt)
     complex(dp), allocatable :: Hwannk(:,:,:) !(nwann, nwann, nkpt)
     complex(dp), allocatable :: psi_wann_k(:,:,:)  !(nbasis, nwann, nkpt)
     integer, allocatable :: exclude_bands(:)
   contains
     procedure :: initialize
     procedure :: finalize
     procedure :: run_all
     procedure :: find_kpoint
     procedure :: remove_phase
     procedure :: get_columns
     procedure :: get_Amnk ! Amnk for all kpoints
     !procedure :: get_anchor_projections
     procedure :: get_weight
     !procedure :: select_column
     procedure :: get_HwannR
     procedure :: get_WannR
     procedure :: write_Amnk
     procedure :: write_Hwann
  end type scdmk



contains

  !===============================================================
  !
  !> @
  !> evals: pointer to eigen values. 2D real matrix. The indices are (iband, ikpt).
  !> psi: pointer to wavefunctions. complex matrix. The indices are (ibasis, iband, ikpt)
  !> kpts: kpoints. indices: (idim, ikpt)
  !> positions_red: reduced coordinate of each basis.
  !> nwann: number of Wannier functions to be calcualted.
  !> anchor_kpt: anchor kpoint (optional).
  !> anchor_ibands: the indices of mode used as anchor points (optional).
  !===============================================================
  subroutine initialize(self, evals, psi, kpts, positions_red, nwann, anchor_kpt, &
       &anchor_ibands, psi_phase, disentangle_func_type, mu, sigma, exclude_bands)
    class(scdmk), intent(inout) :: self
    integer, intent(in) :: nwann,  anchor_ibands(:)
    real(dp), intent(in), target :: evals(:, :)   !(iband, ikpt)
    complex(dp), intent(in), target :: psi(:,:,:) ! (ibasis, iband, ikpt)
    real(dp), intent(in) :: kpts(:,:)  !(idim, ikpt)
    real(dp), intent(in) :: positions_red(:,:) !(idim, ibasis)
    real(dp), intent(in) :: anchor_kpt(:)
    logical, optional, intent(in) :: psi_phase
    integer, optional, intent(in) :: disentangle_func_type
    real, optional, intent(in) :: mu, sigma
    integer, optional, intent(in) :: exclude_bands(:)
    ! allocate memories

    self%nkdim=size(kpts,1)
    self%nkpt=size(kpts,2)
    self%nbasis=size(psi, 1)
    self%nband=size(psi, 2)
    self%nwann=nwann
    self%anchor_kpt=anchor_kpt
    if(present(psi_phase)) then
       if (psi_phase) then
          allocate(self%psi(self%nbasis, self%nband, self%nkpt))
          call self%remove_phase(psi)
           else
          self%psi => psi
       end if
    else
       self%psi => psi
    endif
    self%evals => evals
    allocate(self%cols(self%nwann))
    self%cols(:)=0
    allocate(self%kpts(self%nkdim, self%nkpt))
    self%kpts=kpts
    allocate(self%weight(self%nband, self%nkpt))
    self%weight(:, :)=0.0_dp
    allocate(self%positions_red(3, self%nbasis))
    self%positions_red=positions_red
    allocate(self%anchor_kpt(size(anchor_kpt)))
    self%anchor_kpt=anchor_kpt
    allocate(self%anchor_ibands(size(anchor_ibands)))
    self%anchor_ibands = anchor_ibands
    self%anchor_ikpt=self%find_kpoint(anchor_kpt)

    if(present(disentangle_func_type)) then
       self%disentangle_func_type=disentangle_func_type
    else
       self%disentangle_func_type=0
    end if

    if(present(mu)) then
       self%mu=mu
    else
       self%mu=0
    end if

    if(present(sigma)) then
       self%sigma=sigma
    else
       self%sigma=sigma
    end if

    allocate(self%Amnk(self%nband, self%nwann, self%nkpt))
    allocate(self%psi_wann_k(self%nbasis, self%nwann, self%nkpt))
    allocate(self%Hwannk(self%nwann, self%nwann, self%nkpt))

    if(present(exclude_bands) then
        allocate(self%exclude_bands(size(exclude_bands, 1))
    end if 
  end subroutine initialize

  subroutine finalize(self)
    class (scdmk), intent(inout) :: self
    ! deallocate
    deallocate(self%cols)
    deallocate(self%kpts)
    deallocate(self%positions_red)
    deallocate(self%anchor_kpt)
    deallocate(self%anchor_ibands)
    deallocate(self%Amnk)
    deallocate(self%psi_wann_k)
    deallocate(self%Hwannk)
  end subroutine finalize

  subroutine run_all(self)
    class(scdmk), intent(inout) :: self
    integer :: ikpt, iband, iwann
    complex(dp) :: psi_dagger(self%nband, self%nbasis)
    complex(dp) :: tmp(self%nband, self%nwann)
    ! find anchor points, by default gamma
    if(size(self%anchor_ibands)/=0) then
       self%anchor_ikpt=self%find_kpoint(self%anchor_kpt)
    end if
    ! calculate weight matrix for each kpoint
    do ikpt=1, self%nkpt
       call self%get_weight(ikpt, self%disentangle_func_type, self%mu, self%sigma, self%weight(:, ikpt))
    end do
    ! at anchor-kpoint, find cols
    !psi: (ibasis, iband)
    psi_dagger= transpose(conjg(self%psi(:, :, self%anchor_ikpt)))
    do iband =1, self%nband
       psi_dagger(iband, :) = psi_dagger(iband, :) * self%weight(:, self%anchor_ikpt)
    end do
    call self%get_columns(psi_dagger, self%cols )

    ! For each kpoint, calculate Amn matrix, wannier function, and Hk at k
    do ikpt=1, self%nkpt
       psi_dagger(:, iband) = psi_dagger(:, iband) * self%weight(:, self%anchor_ikpt)

       !subroutine get_Amnk(self, psik, cols, weightk,  Amnk)
       !Amnk (nband, nwann, nkpt)
       call self%get_Amnk(psi_dagger, self%cols, self%Amnk(:, :, ikpt))
       self%psi_wann_k(:,:, ikpt)=matmul(self%psi(:,:, ikpt), self%Amnk(:,:, ikpt))

       do iwann=1, self%nwann
          tmp(:,iwann)=self%Amnk(:,iwann, ikpt) *self%evals(:,ikpt)
       end do
       ! Calculate Hamiltonian at k
       self%Hwannk(:,:,ikpt) = matmul(transpose(conjg(self%Amnk(:,:, ikpt))), tmp)
    end do

    ! Fourier transform of wannier function to real space
  end subroutine run_all

  subroutine remove_phase(self, psip)
    class(scdmk), intent(inout) :: self
    complex(dp), intent(in) :: psip(:,:,:) ! (ibasis, iband, ikpt)
    !complex(dp), intent(out) :: psi(:,:,:) ! (ibasis, iband, ikpt)
    integer :: ikpt, ibasis
    complex(dp) :: phase
    do ikpt=1, self%nkpt
       do ibasis=1, self%nbasis
          phase = exp(-tpi_im * dot_product(self%kpts(:, ikpt), self%positions_red(:,ibasis)) )
          self%psi(ibasis, :, ikpt) = psip(ibasis, :, ikpt) * phase
       end do
    end do
  end subroutine remove_phase

  !===============================================================
  ! Find one kpoint in a list of kpoints.
  !> @
  !===============================================================
  function find_kpoint(self, kpoint) result(ik)
    class (scdmk), intent(inout) :: self
    real(dp), intent(in) :: kpoint(:)
    integer :: ik
    integer :: i
    real(dp) :: a(size(self%kpts,1))
    ! should transfer back to 1st BZ?
    do i =1, size(self%kpts, 2)
       a(:)= sum((self%kpts(:, i)-kpoint) **2)
    end do
    ik=minloc(a, dim=1)
    if( a(ik)>0.001) then
       print *, "Error in finding gamma point from kpoints, gamma not found. "
    end if
  end function find_kpoint

  !===============================================================
  ! Calculate the weight function for each mode described by iband and ikpt
  ! The
  !> @
  !===============================================================
  subroutine get_weight(self, ikpt, type, mu, sigma, weight)
    class(scdmk), intent(inout) :: self
    integer, intent(in) :: ikpt, type
    real(dp), intent(in) :: mu, sigma
    real(dp), intent(inout) :: weight(self%nband)
    integer :: iband, ianchor
    select case(type)
    case(0)
       weight(:) = 1.0
    case(1)
       do iband=1, self%nband
          weight(iband) = fermi( self%evals(iband, ikpt), mu, sigma)
       end do
    case(2)
       do iband=1, self%nband
          weight(iband) = gaussian( self%evals(iband, ikpt), mu, sigma)
       end do
    end select

    if(size(self%anchor_ibands) /= 0) then
       do iband=1, self%nband
          do ianchor=1, size(self%anchor_ibands)
             weight(iband)=weight(iband)* real(dot_product(  &
                  & conjg(self%psi(:, ianchor, self%anchor_ikpt)), &
                  & self%psi(:, iband, ikpt) ))
          end do
       end do
    end if
  end subroutine get_weight

  subroutine get_columns(self, psi_dagger, cols)
    class(scdmk), intent(inout) :: self
    complex(dp), intent(in) :: psi_dagger(:,:)
    integer :: piv(size(psi_dagger, 2))
    integer, intent(inout) :: cols(self%nwann)
    !complex(dp), allocatable :: U(:,:), S(:,:), VT(:,:)
    call complex_QRCP_piv_only(psi_dagger, piv)
    cols=piv(:self%nwann)
  end subroutine get_columns


  subroutine get_Amnk(self, psi_dagger, cols,  Amnk)
    class(scdmk), intent(inout) :: self
    complex(dp), intent(in) :: psi_dagger(self%nband,self%nbasis)
    integer :: cols(:)
    complex(dp), intent(inout) :: Amnk(self%nbasis, self%nwann)
    complex(dp), allocatable :: U(:,:), VT(:,:)
    real(dp), allocatable :: S(:)
    call complex_svd(psi_dagger(:, cols), U, S, VT)
    Amnk(:,:) = matmul(U, conjg(transpose(VT)))
    deallocate(U)
    deallocate(S)
    deallocate(VT)
  end subroutine get_Amnk

  subroutine get_WannR(self, Rgrid, Wann)
    class(scdmk), intent(inout) :: self
    integer, intent(in) :: Rgrid(:, :)
    real(dp) :: Wann(self%nbasis, self%nwann, size(Rgrid, 2))
    complex(dp) :: Wann_complex(self%nbasis, self%nwann, size(Rgrid, 2))
    integer :: ik, iR, nR
    nR=size(Rgrid, 2)
    Wann_complex(:, :,:)=cmplx(0.0,0.0, dp)
    do ik=1, self%nkpt
       do iR=1, nR
          Wann_complex(:,:, iR) = self%psi_wann_k(:,:, ik) * exp(tpi_Im*dot_product(self%kpts(:, ik), Rgrid(:, iR)))
       end do
    end do
    Wann=real(real(Wann_complex))
  end subroutine get_WannR

  subroutine get_HwannR(self, Rgrid, HR)
    class(scdmk), intent(inout) :: self
    integer, intent(in) :: Rgrid(:, :)
    complex(dp), intent(out) :: HR(self%nwann, self%nwann, size(Rgrid,2))
   ! H(R)= \sum_k H(k) * exp(i2pi k.R)
    integer :: ik, iR, nR
    nR=size(Rgrid, 2)
    HR(:, :,:)=cmplx(0.0,0.0, dp)
    do ik=1, self%nkpt
       do iR=1, nR
          HR(:,:, iR) = self%Hwannk(:,:, ik) * exp(tpi_Im*dot_product(self%kpts(:, ik), Rgrid(:, iR)))
       end do
    end do
  end subroutine get_HwannR

  subroutine Amn_to_H(Amn, psi, H0, H, nbasis, nwann)
    complex(dp), intent(in) :: Amn(:,:), psi(:,:), H0(:,:)
    complex(dp), intent(inout) :: H(:,:)
    integer, intent(in) :: nbasis, nwann
    complex(dp) :: tmp(nbasis, nwann)
    tmp(:,:)=matmul(psi, Amn)
    H=matmul(matmul(transpose(conjg(tmp)), H0), tmp)
  end subroutine Amn_to_H

  subroutine write_Amnk(self, fname)
    class(scdmk), intent(inout) :: self
    character(len=*), intent(in) :: fname
    integer :: iwann, iband, ikpt, locibnd
    integer :: iun_amn

    iun_amn=101
    OPEN (unit=iun_amn, file=trim(fname)//".amn",form='formatted')

    !IF (wan_mode=='standalone') THEN
    !   iun_amn = find_free_unit()
    !   IF (ionode) OPEN (unit=iun_amn, file=trim(seedname)//".amn",form='formatted')
    !ENDIF

    !WRITE(stdout,'(a,i8)') '  AMN: iknum = ',iknum
    !
    !IF (wan_mode=='standalone') THEN
    !   CALL date_and_tim( cdate, ctime )
    !   header='Created on '//cdate//' at '//ctime//' with SCDM '
    !   IF (ionode) THEN
    !      WRITE (iun_amn,*) header
    !      WRITE (iun_amn,'(3i8,xxx,2f10.6)') numbands,  iknum, n_wannier, scdm_mu, scdm_sigma
    !   ENDIF
    !ENDIF

    do ikpt =1, self%nkpt
       do iwann = 1, self%nwann
          locibnd = 0
          do iband = 1,self%nband
             !IF (excluded_band(iband)) CYCLE
             locibnd = locibnd + 1
             WRITE(iun_amn,'(3i5,2f18.12)') locibnd, iwann, ikpt, &
                  & REAL(self%Amnk(locibnd,iwann,ikpt)), &
                  & AIMAG(self%Amnk(locibnd,iwann, ikpt))
          enddo
       enddo
    enddo

    close(iun_amn)
  end subroutine write_Amnk

  subroutine write_Hwann(self, HR, Rgrid, fname)
    class(scdmk), intent(inout) :: self
    complex(dp), intent(in) :: HR(:, :, :)
    integer, intent(in) :: Rgrid(:, :)
    character(len=*), intent(in) :: fname
    integer :: iR ,ifile, iwann1, iwann2
    ifile=103
    OPEN (unit=ifile, file=trim(fname)//".amn",form='formatted')
    do iR=1, size(Rgrid, 2)
       WRITE(ifile, '(3i5)') Rgrid(:, iR)
       do iwann1 =1, self%nwann
          do iwann2 =1, self%nwann
             WRITE(ifile, '(f18.12)') HR(iwann1, iwann2, iR)
          end do
       end do
    end do
    close(ifile)
  end subroutine write_Hwann

end module m_scdm

