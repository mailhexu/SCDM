#include "abi_common.h"

module m_tight_binding
  use defs_basis , only : dp,dpc, j_dpc, two_pi
  use m_mathfuncs, only : eigensolver
  use m_supercell_maker, only : supercell_maker_t
  use m_test_utils, only: assertion, all_close_real
  use m_hashtable_strval, only: hash_table_t
  use m_io_netcdf, only: read_tb_netcdf
  !use read_inputs, only: HamR
  !use generate_ham, only: tran, generate_hamr_from_TB

  implicit none
  complex(dpc), parameter ::  itwopi = two_pi *j_dpc

  private

  !===============================================================
  ! Tight binding Hamiltonian
  !===============================================================
  type, public :: TBHam
     integer :: nR=0 ! number of R points 
     integer :: nwann=0 ! number of Wannier functions
     ! Note: natom, ipeci....to iorb_list are for labeling the wannier functions,
     ! which are useful for finding the DMFT correlated subspace.
     integer :: natom=0
     integer, allocatable :: ispecie_list(:), iatom_list(:), iorb_list(:) ! index of specie, atom and orbital for each atom
     complex(dp), allocatable:: HamR(:,:,:) ! Hamiltonian. index: iwann, jwann, iR
     integer, allocatable:: Rlist(:, :) ! Rpoints, 3*nR
     type(hash_table_t) :: Rdict ! for finding the index of R in Rlist.
     type(eigensolver) :: esolver ! eigen solver
   contains
     procedure :: initialize
     procedure :: add_new_R
     procedure :: set_Rlist
     procedure :: set_Rlist_and_HamR
     procedure :: set_basis_indices
     procedure :: set_HijR
     procedure :: set_Hij_indR
     procedure :: read_from_netcdf
     !procedure :: read_from_wannier90_hr
     procedure :: finalize
     procedure :: get_Hk
     procedure :: solve
     procedure :: solve_all
     procedure :: fill_supercell
     procedure :: fill_supercell_with_matrix
  end type TBHam

  public :: build_1d_chain
  public :: test_tight_binding

contains

  !===============================================================
  ! Initialize Hamiltonian.
  ! The number of orbitals are fixed while the number of R is dynamic
  !> @ nwann : number of orbitals
  !===============================================================
  subroutine initialize(self,  nwann)
    class(TBHam),intent(inout) :: self
    integer, intent(in) ::  nwann
    self%nR=0
    self%nwann=nwann
    ABI_MALLOC(self%ispecie_list, (nwann))
    ABI_MALLOC(self%iatom_list, (nwann))
    ABI_MALLOC(self%iorb_list, (nwann))
    call self%Rdict%init()
  end subroutine initialize

  !===============================================================
  ! Free memory
  !===============================================================
  subroutine finalize(self)
    class(TBHam),intent(inout) :: self
    self%nR=0
    self%nwann=0
    self%natom=0
    call self%esolver%finalize()
    call self%Rdict%free()
    ABI_SFREE(self%HamR)
    ABI_SFREE(self%Rlist)
    ABI_SFREE(self%ispecie_list)
    ABI_SFREE(self%iatom_list)
    ABI_SFREE(self%iorb_list)
  end subroutine finalize

  !===============================================================
  ! Put a list of R into the Hamiltonian.
  !> @ Rlist (3, nR)
  !===============================================================
  subroutine set_Rlist(self, Rlist)
    class(TBHam),intent(inout) :: self
    integer, intent(in) :: Rlist(:,:)
    integer :: i
    self%nR=size(Rlist, 2)
    ABI_MALLOC(self%Rlist, (3, self%nR))
    ABI_MALLOC(self%hamr,(self%nwann, self%nwann, self%nR))
    self%Rlist=Rlist
    self%hamr=cmplx(0.0_dp, 0.0_dp, kind=dp)
    do i=1, self%nR
       call self%Rdict%put_int3(Rlist(:,i), i)
    end do
  end subroutine set_Rlist


  subroutine set_Rlist_and_HamR(self, Rlist, HamR)
    class(TBHam),intent(inout) :: self
    integer, intent(in) :: Rlist(:,:)
    complex(dp), intent(in) :: HamR(:,:,:)
    call self%set_Rlist(Rlist)
    self%HamR=HamR
  end subroutine set_Rlist_and_HamR

  !===============================================================
  ! set the indices for all basis
  !> @ iatom_list: the index of atom for each wannier function
  !> @ iorb_list: the index of orbital in an atom for each wannier function
  !===============================================================
  subroutine set_basis_indices(self, iatom_list, iorb_list)
    class(TBHam),intent(inout) :: self
    integer, intent(in) :: iatom_list(:), iorb_list(:)
    self%iatom_list(:)=iatom_list(:)
    self%iorb_list(:)=iorb_list(:)
  end subroutine set_basis_indices

  subroutine add_new_R(self, R)
    ! add an H(R) to the Hamiltonian.
    class(TBHam),intent(inout) :: self
    integer, intent(in) ::  R(3)
    integer, allocatable :: tmp_Rlist(:,:)
    complex(dp), allocatable :: tmp_Hamr(:,:, :)

    if (.not. self%Rdict%has_key_int3(R)) then
       self%nR=self%nR+1
       call self%Rdict%put_int3(R, self%nR)
       ABI_MALLOC(tmp_Rlist, (3,self%nR))
       if (self%nR>1) then
          tmp_Rlist(:, :self%nR-1)= self%Rlist
       end if
       ABI_MOVE_ALLOC(tmp_Rlist, self%Rlist)
       self%Rlist(:, self%nR) = R
       ABI_MALLOC(tmp_Hamr, (self%nwann, self%nwann,self%nR))
       if (self%nR>1) then
          tmp_Hamr(:,:, :self%nR-1)= self%hamr
       end if
       ABI_MOVE_ALLOC(tmp_Hamr, self%hamr)
       self%hamr(:,:,self%nR)= 0.0_dp
    end if
  end subroutine add_new_R


  ! set one Hij(R)
  subroutine set_HijR(self, i,j, R, val)
    class(TBHam),intent(inout) :: self
    integer, intent(in) :: i, j, R(3)
    complex(dp), intent(in) :: val
    integer :: iR
    call self%add_new_R(R)
    iR=self%Rdict%get_int3(R)
    call self%set_Hij_indR(i,j, iR, val)
  end subroutine set_HijR

  ! set one Hij(R) but using the index of R in the Rlist
  subroutine set_Hij_indR(self, i,j, iR, val)
    class(TBHam),intent(inout) :: self
    integer, intent(in) :: i, j, iR
    complex(dp) :: val
    self%hamR(i, j, iR) = val
  end subroutine set_Hij_indR


  !===============================================================
  ! Read tight-binding model from netcdf file generated by BandDownfolder
  ! No need to call initialize before this subroutine.
  !> @ fname: File name of the netcdf file
  !===============================================================
  subroutine read_from_netcdf(self, fname)
    class(TBHam),intent(inout) :: self
    character(len=*), intent(in) :: fname
    integer, allocatable :: Rlist(:,:)
    real(dp), allocatable :: HamR_real(:,:,:), HamR_imag(:,:,:)
    integer :: nwann
    call read_tb_netcdf(fname, Rlist, HamR_real, HamR_imag)
    nwann=size(HamR_real,1)
    call self%initialize(nwann)
    call self%set_Rlist_and_HamR(Rlist, cmplx(HamR_real, HamR_imag, kind=dp))
    ABI_SFREE(Rlist)
    ABI_SFREE(HamR_real)
    ABI_SFREE(HamR_imag)

    ! TODO: should read the index of atoms, species, orbitals.
    ! Not yet implmented.
  end subroutine read_from_netcdf

  !===============================================================
  ! Read tight-binding model from Wannier90 hr.dat file
  ! No need to call initialize before this subroutine.
  !===============================================================
  !subroutine read_from_wannier90_hr(self, hr_fname)
  !  class(TBHam),intent(inout) :: self
  !  character(len=*), intent(in) :: hr_fname
  !  integer :: iR, iwann, nwann
  !  call generate_hamr_from_TB(hr_fname=hr_fname)
  !  nwann=size(HamR, 2)
  !  print *, "nwann:", nwann
  !  call self%initialize(nwann)
  !  call self%set_Rlist(tran)
  !  print *, "nR:", self%nR
  !  do iR=1, self%nR
  !     ! TODO: should it be transposed? In non-SOC calculation it should be
  !     ! the same, but need to be checked.
  !     self%HamR(:,:, iR) = HamR(iR, :,:)
  !     !if(dot_product(tran(:,iR), tran(:,iR))==0) then
  !     !   do iwann=1, nwann
  !     !      self%HamR(iwann, iwann, iR)=HamR(iR, iwann, iwann)
  !     !   end do
  !     !end if
  !  end do

  !  ! TODO: should read the index of species, atoms, orbitals.
  !  ! This is not a feature of W90, should it come from another file?
  !end subroutine read_from_wannier90_hr

  !===============================================================
  ! Calculate Hamiltonian at k-point (Hk)
  !> @ k: the kpoint in the reciprocal lattice unit
  !> @ Hk: the outputed Hamiltonian
  !===============================================================
  subroutine get_Hk(self, k, Hk)
    class(TBHam), intent(in) :: self
    real(dp), intent(in) :: k(3)
    complex(dp), intent(inout) :: Hk(self%nwann, self%nwann)
    integer :: iR
    Hk(:,:) = 0.0_dp
    do iR =1, self%nR
       Hk=Hk+ self%HamR(:,:, iR)* exp(-itwopi* dot_product(k, self%Rlist(:,iR)))
    end do
  end subroutine Get_Hk

  ! get the eigenvalues and eigenvectors for one kpoint.
  subroutine solve(self, k, evals, evecs)
    class(TBHam), intent(inout) :: self
    real(dp), intent(in) :: k(:)
    real(dp), intent(inout) :: evals(self%nwann)
    complex(dp), intent(inout) ::   evecs(self%nwann, self%nwann)
    call self%get_Hk(k, evecs)
    call self%esolver%run(evals, evecs)
  end subroutine solve

  !===============================================================
  ! Diagonize the Hamiltonian for a set of kpoints
  !> @ kpts: dim: (3, nk)
  !> @ evals : eigen values. dim:(norb , nk)
  !> @ evecs : eigen vectors. dim: (norb, nband, nk). The first dimension
  ! is the number of basis, and gives one eigenvector.
  !===============================================================
  subroutine solve_all(self, kpts, evals, evecs)
    class(TBHam), intent(inout) :: self
    real(dp), intent(in) :: kpts(:,:)
    real(dp), intent(inout) :: evals(:,:)
    complex(dp), intent(inout) :: evecs(self%nwann,self%nwann, size(kpts, 2))
    integer :: ik
    do ik=1, size(kpts, 2)
       call self%solve(kpts(:, ik), evals(:, ik), evecs(:, :, ik))
    end do
  end subroutine solve_all


  !===============================================================
  !
  !> @ scmaker
  !> @ sc_tb: tight binding model hamiltonian ins supercell
  !===============================================================
  subroutine fill_supercell(self, scmaker, sc_tb)
    class(TBHam), intent(inout) :: self, sc_tb
    class(supercell_maker_t), intent(inout) :: scmaker
    integer :: sc_nwann, iR, ii, ij,icell,  jj
    integer :: i_sc, j_sc, Rj_sc(3)
    !complex(dp)

    sc_nwann=self%nwann*scmaker%ncells
    call sc_tb%initialize(sc_nwann)

    do iR =1, self%nR
       do icell =1, scmaker%ncells
          call scmaker%R_to_sc(self%Rlist(:, iR) + scmaker%rvecs(:,icell), Rj_sc, jj)
          do ij=1, self%nwann
             j_sc=self%nwann*(jj-1)+ij
             do ii=1, self%nwann
                i_sc=self%nwann*(icell-1)+ii
                associate(val=> self%hamr(ii, ij, iR))
                  call sc_tb%set_HijR(i_sc, j_sc, Rj_sc,val)
                end associate
             end do
          end do
       end do
    end do

    sc_tb%natom=self%natom * scmaker%ncells
    call scmaker%trans_ilist(self%natom, self%iatom_list, sc_tb%iatom_list )
    call scmaker%repeat_int1d_noalloc(self%ispecie_list, sc_tb%ispecie_list)
    call scmaker%repeat_int1d_noalloc(self%iorb_list, sc_tb%iorb_list)
  end subroutine fill_supercell

  !===============================================================
  ! Fill supercell with input supercell matrix
  !> @ sc_mat: a 3*3 integer matrix
  !> @ sc_tb: tight binding model hamiltonian ins supercell
  !===============================================================
  subroutine fill_supercell_with_matrix(self, sc_mat, sc_tb)
    class(TBHam), intent(inout) :: self, sc_tb
    integer, intent(in) :: sc_mat(3,3)
    type(supercell_maker_t)  :: scmaker
    call scmaker%initialize(sc_mat)
    call self%fill_supercell(scmaker, sc_tb)
    call scmaker%finalize()
  end subroutine fill_supercell_with_matrix


!=============================================================== For test:

  subroutine build_1d_chain(tb, on_site, hopping)
    type(TBHam) :: tb
    complex(dp) :: on_site, hopping
    integer :: Rlist(3,3)
    call tb%initialize(nwann=1)
    Rlist(:,1)=[-1,0,0]
    Rlist(:,2)=[0,0,0]
    Rlist(:,3)=[1,0,0]
    call tb%set_Rlist(Rlist)
    call tb%set_Hij_indR(i=1, j=1, iR=1, val=hopping)
    call tb%set_Hij_indR(i=1, j=1, iR=2, val=on_site)
    call tb%set_Hij_indR(i=1, j=1, iR=3, val=hopping)
  end subroutine build_1d_chain


  subroutine test_tight_binding()
    type(TBHam) :: tb, tb_sc
    integer, parameter :: n_sc=5
    integer :: scmat(3,3)
    real(dp) :: evals(1), sc_evals(n_sc)
    complex(dp) :: evecs(1, 1), sc_evecs(n_sc,n_sc)
    type(supercell_maker_t) :: scmaker
    ! single band 
    call build_1d_chain(tb, on_site=cmplx(1.0_dp,0.0_dp, kind=dp), &
         & hopping=cmplx(1.0_dp,0.0_dp, kind=dp))
    call tb%solve(k=[0.0_dp, 0.0_dp, 0.0_dp], evals=evals, evecs=evecs)
    call assertion(abs(evals(1)-3.0_dp)<1d-5, "Wrong eigen value")
    call assertion(abs(evecs(1,1)-1.0_dp)<1d-5, "Wrong eigen vector")

    call tb%solve(k=[0.5_dp, 0.0_dp, 0.0_dp], evals=evals, evecs=evecs)
    call assertion(abs(evals(1)+1.0_dp)<1d-5, "Wrong eigen value")

    scmat=0
    scmat(1,1)=n_sc
    scmat(2,2)=1
    scmat(3,3)=1
    call scmaker%initialize(scmat)
    call tb%fill_supercell(scmaker, tb_sc)

    call tb_sc%solve(k=[0.0_dp, 0.0_dp, 0.0_dp], evals=sc_evals, evecs=sc_evecs)
    !print *, sc_evals
    !print *, [-1.61803399_dp, -1.61803399_dp,  0.61803399_dp,  0.61803399_dp,  2._dp]+1.0
    call assertion(all_close_real(sc_evals,[-1.61803399_dp, -1.61803399_dp,  0.61803399_dp,  0.61803399_dp,  2._dp]+1.0, 1d-6), &
         & "wrong eigen value for supercell"  )
    call tb%finalize()
    call tb_sc%finalize()

  end subroutine test_tight_binding


end module m_tight_binding
