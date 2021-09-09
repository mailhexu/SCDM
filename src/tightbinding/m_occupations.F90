module m_occ
  use defs_basis, only: dp
  use m_bisect, only: bisect_solver
  use m_mathfuncs, only: outer_product
  implicit none
  private
  public:: find_efermi
  public:: get_density_matrix
  public:: get_density


contains

 pure function fermi(e, mu,  beta, nspin) result(y)
    real(dp), intent(in) :: e, mu, beta
    integer, intent(in) :: nspin
    real(dp) :: x, y
    x=(e-mu)/beta
    if (x< -10) then
       y=1.0/nspin
    else if (x>10) then
       y=0.0_dp
    else
       y=2.0_dp/nspin/(exp(x)+1.0)
    end if
  end function fermi

 !===============================================================
  !Fermi function. Vectorized version
  !> @ es: energy list ()
  !> @ mu: fermi energy
  !> @ beta: kb T
  !> @ nspin: number of spin in the model. 1: non-polarized. 2: polarized
  !> @ occs: occupations of each point in the band
  !===============================================================
  subroutine fermi_vec(es, mu,  beta, nspin, occs)
    real(dp), intent(in) :: es(:, :), mu, beta
    integer, intent(in) :: nspin
    real(dp), intent(out) :: occs(:, :)
    occs=(es-mu)/beta
    where (occs>10)
       occs=0.0_dp
    elsewhere
       occs=2.0_dp/nspin/(exp(occs)+1.0)
    end where
  end subroutine fermi_vec

  !===============================================================
  ! Get number of electrons from 
  !> @
  !===============================================================
  function get_ne(f, kweights) result (ne)
    real(dp) , intent(in) :: f(:,:), kweights(:)
    real(dp) :: ne
    integer :: i
    ne=0.0_dp
    do i=1, size(kweights)
       ne=ne+ sum(f(:, i)*kweights(i))
    end do
  end function

  subroutine find_efermi(evals, kweights, beta, occs, nspin, nel, thr, mu)
    real(dp), intent(in) :: evals(:,:), kweights(:), nel, thr, beta
    integer, intent(in) :: nspin
    real(dp), intent(inout) :: mu, occs(:,:)
    real(dp) :: a, b, fa, fb, fx
    type(bisect_solver) ::  solver

    a = minval(evals)-0.1
    b = maxval(evals)+0.1

    if (nel< 1e-3) then
       a=a-1000.0_dp
    else if(nel > size(evals, 1)*nspin-1.0d-3) then
       b=b+1000.0_dp
    end if

    if (nel<0 .or. nel>size(evals, 1)*2/nspin) then
       print *, "The number of electron is not in (0, nband*nspin)"
       stop 
    endif

    mu=0.5*(a+b)

    call fermi_vec(evals,a, beta, nspin, occs)
    fa=get_ne(occs, kweights)-nel

    call fermi_vec(evals,b, beta, nspin, occs)
    fb=get_ne(occs, kweights)-nel

    call solver%initialize(a, b, fa, fb)

    fx=1.0_dp
    do while(abs(fx)>thr .and. (.not. solver%too_many_steps()))
    !do while(abs(fx)>thr)

       call fermi_vec(evals,mu, beta, nspin, occs)
       fx=get_ne(occs, kweights)-nel
       call solver%nextx(mu, fx)
    end do
    !call fermi_vec(evals,mu, beta, nspin, occs)
  end subroutine find_efermi

  subroutine get_density_matrix(evecs, kweights, occs, dm)
    complex(dp) , intent(in) :: evecs(:,:,:)
    real(dp) :: occs(:,:), kweights(:)
    complex(dp) :: dm(:,:)
    integer :: ik, iband
    dm(:,:)=0.0_dp
    do ik=1, size(kweights)
       do iband=1, size(occs,1)
          dm=dm+ outer_product(conjg(evecs(:, iband, ik)), &
               &evecs(:, iband, ik))*(occs(iband, ik)*kweights(ik))
       end do
    end do
  end subroutine get_density_matrix

  subroutine get_density(evecs, kweights, occs, rho)
    complex(dp) , intent(in) :: evecs(:,:,:)
    real(dp), intent(in) :: occs(:,:), kweights(:)
    real(dp), intent(inout) :: rho(:)
    real(dp) :: tmp
    integer :: ik, iband
    rho(:)=0.0_dp
    do ik=1, size(kweights)
       do iband=1, size(occs,1)
          tmp=occs(iband, ik)*kweights(ik)
          rho(:)=rho(:)+ real(conjg(evecs(:, iband, ik))* &
           & evecs(:, iband, ik))*tmp
       end do
    end do
  end subroutine get_density

end module m_occ
