!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! User-routine for defining velocity/attenuation structure
!!
!! Model for NIED balldrop test
!! 2020.12.16- Kurama Okubo (kokubo@bosai.go.jp)
!!
!! @copyright
!!   Copyright 2013-2020 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!!
!! 2021.04.26 Kurama Okubo (kokubo@bosai.go.jp) developped m_vmodel_balldropseg.F90
!!
!<
!! ----
module m_vmodel_balldropseg_sidecoord

  use m_std
  use m_global
  use m_readini
  implicit none
  private
  save

  public :: vmodel_balldropseg_sidecoord

contains

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Define meidum velocity, density and attenuation
  !!
  !! This is a user-specific routine to define original veloicty model.
  !!
  !! Input:
  !!    io_prm                          :: I/O number of parameter file (file has been opened already)
  !!    i0,i1, j0,j1, k0,k1             :: model area by indices in i-,j-,and k- directions
  !!    xc(i0:i1), yc(j0:j1), zc(k0:k1) :: Cartesian coordinate location
  !!    vcut                            :: cut-off velocity specified by input parameter
  !!
  !! Output:
  !!    rho(k0:k1, i0:i1, j0:j1)        :: mass density (usually in g/cm^3)
  !!    lam(k0:k1, i0:i1, j0:j1)        :: Lame's parameter (usually in (g/cm^3) * (km/s)^2)
  !!    mu (k0:k1, i0:i1, j0:j1)        :: Lame's parameter (usually in (g/cm^3) * (km/s)^2)
  !!    qp (k0:k1, i0:i1, j0:j1)        :: Attenuation QP
  !!    qs (k0:k1, i0:i1, j0:j1)        :: Attenuation QS
  !!    bd (i0:i1, j0:j1, 0:NBD)        :: Boundary depths
  !!
  !! Note:
  !! bd(:,:,0) are treated as topography shape for output.
  !!    this is only for output and visualization. topography in the simulation will be automatically detected by medium params.
  !! bd(:,:,1:NBD) may contain internal boundary depths. The boundary number can be specified as source depth or station depth.
  !!
  !<
  !! ----
  subroutine vmodel_balldropseg_sidecoord( io_prm, i0, i1, j0, j1, k0, k1, xc, yc, zc, vcut, rho, lam, mu, Qp, Qs, bd )

    !! -- Arguments
    integer,  intent(in)  :: io_prm
    integer,  intent(in)  :: i0, i1                         !< i-region
    integer,  intent(in)  :: j0, j1                         !< j-region
    integer,  intent(in)  :: k0, k1                         !< k-region
    real(SP), intent(in)  :: xc  ( i0:i1 )                  !< x-coordinate location
    real(SP), intent(in)  :: yc  ( j0:j1 )                  !< y-coordinate location
    real(SP), intent(in)  :: zc  ( k0:k1 )                  !< z-coordinate location
    real(SP), intent(in)  :: vcut                           !< cut-off minimum velocity
    real(SP), intent(out) :: rho ( k0:k1, i0:i1, j0:j1 )    !< mass density [g/cm^3]
    real(SP), intent(out) :: lam ( k0:k1, i0:i1, j0:j1 )    !< Lame's parameter lambda [ (g/cm^3) * (km/s)**2 ]
    real(SP), intent(out) :: mu  ( k0:k1, i0:i1, j0:j1 )    !< Lame's parameter mu     [ (g/cm^3) * (km/s)**2 ]
    real(SP), intent(out) :: qp  ( k0:k1, i0:i1, j0:j1 )    !< P-wave attenuation
    real(SP), intent(out) :: qs  ( k0:k1, i0:i1, j0:j1 )    !< S-wave attenuation
    real(SP), intent(out) :: bd  ( i0:i1, j0:j1, 0:NBD )    !< Boundary depths
    !! --

    integer  :: i, j, k
    ! real(SP) :: vp0, vs0, rho0, qp0, qs0, topo0
    real(SP) :: vp_rock, vs_rock, rho_rock, qp0_rock, qs0_rock
    real(SP) :: vp_metal, vs_metal, rho_metal, qp0_metal, qs0_metal
    real(SP) :: balldrop_w, balldrop_h, balldrop_hm, topo0 !,balldrop_l, balldrop_lm,

    ! real(SP) :: vp1, vs1
    real(SP) :: dum
    !! ----

    !!
    !! The following dummy code is an example how to discribe the routine.
    !!

    !!
    !! subroutine readini() can access parameters defined in the input file.
    !! Any original parameters can be added in the input file.
    !!

    !! Read parameters for balldrop model
    call readini( io_prm, 'vp_rock',   vp_rock, 6.919 )
    call readini( io_prm, 'vs_rock',   vs_rock, 3.631 )
    call readini( io_prm, 'rho_rock',  rho_rock, 2.98 )
    call readini( io_prm, 'qp0_rock',   qp0_rock,  1000000.0 )
    call readini( io_prm, 'qs0_rock',   qs0_rock,  1000000.0 )

    call readini( io_prm, 'vp_metal',  vp_metal, 5.5 )
    call readini( io_prm, 'vs_metal',  vs_metal, 3.2 )
    call readini( io_prm, 'rho_metal', rho_metal, 8.0 )
    call readini( io_prm, 'qp0_metal',  qp0_metal, 1000000.0 )
    call readini( io_prm, 'qs0_metal',  qs0_metal, 1000000.0 )

    ! call readini( io_prm, 'balldrop_l',   balldrop_l, 1000.0e-6 )
    call readini( io_prm, 'balldrop_w',   balldrop_w,  100.0e-6 )
    call readini( io_prm, 'balldrop_h',   balldrop_h,  200.0e-6 )
    ! call readini( io_prm, 'balldrop_lm',   balldrop_lm,  50.0e-6 )
    call readini( io_prm, 'balldrop_hm',   balldrop_hm,  20.0e-6 )

    !! bd is not used in this model, so topo0 is fixed at 0.
    topo0 = 0.0

    !! Way to define velocity model

    !!
    !! The medium parameter must be set from given region (i0:i1, j0:j1, k0:k1)
    !! Note that the order of indices is k->i->j, for improving performance
    !!
    do j = j0, j1
      do i = i0, i1

        !! define topography shape here
        bd(i,j,0) = topo0

        do k = k0, k1
          !! Flowchart to define the velocity model for ball drop test
          !! 1. if (i, j) is outside of medium, assign the air column
          !! 2. if (i, j) is inside the medium, select if it is rock, metal or air.
          !! ---------------------- ----------------------------------------------

          if( zc(k) < 0 .or. zc(k) > balldrop_w ) then
            !! this is outside of medium, assign the air column
            rho(k,i,j) = 0.001
            mu (k,i,j) = rho(k,i,j) * vs1 * vs1
            lam(k,i,j) = rho(k,i,j) * ( vp1*vp1 - 2*vs1*vs1 )
            qp (k,i,j) = 10.0 ! artificially strong attenuation in air-column
            qs (k,i,j) = 10.0 ! artificially strong attenuation in air-column

          else if (yc( j ) > balldrop_h+balldrop_hm) then
            !! this is air under the bottom metal plate
            rho(k,i,j) = 0.001
            mu (k,i,j) = rho(k,i,j) * vs1 * vs1
            lam(k,i,j) = rho(k,i,j) * ( vp1*vp1 - 2*vs1*vs1 )
            qp (k,i,j) = 10.0 ! artificially strong attenuation in air-column
            qs (k,i,j) = 10.0 ! artificially strong attenuation in air-column

          else if (yc( j ) > balldrop_h) then
            !! this is bottom metal plate
            rho(k,i,j) = rho_metal
            mu (k,i,j) = rho(k,i,j) * vs_metal * vs_metal ! rho(k,i,j) * vs1 * vs1
            lam(k,i,j) = rho(k,i,j) * ( vp_metal*vp_metal - 2*vs_metal*vs_metal ) ! rho(k,i,j) * ( vp1*vp1 - 2*vs1*vs1 )
            qp (k,i,j) = qp0_metal ! artificially strong attenuation in air-column
            qs (k,i,j) = qs0_metal ! artificially strong attenuation in air-column

            ! if (xc(i) < -balldrop_lm .or. xc(i) > balldrop_l + balldrop_lm) then
            !   !! this is air outside of the side metal plates
            !   rho(k,i,j) = 0.001
            !   mu (k,i,j) = rho(k,i,j) * vs1 * vs1
            !   lam(k,i,j) = rho(k,i,j) * ( vp1*vp1 - 2*vs1*vs1 )
            !   qp (k,i,j) = 10.0 ! artificially strong attenuation in air-column
            !   qs (k,i,j) = 10.0 ! artificially strong attenuation in air-column
            ! else
            !   !! this is bottom metal plate
            !   rho(k,i,j) = rho_metal
            !   mu (k,i,j) = rho(k,i,j) * vs_metal * vs_metal ! rho(k,i,j) * vs1 * vs1
            !   lam(k,i,j) = rho(k,i,j) * ( vp_metal*vp_metal - 2*vs_metal*vs_metal ) ! rho(k,i,j) * ( vp1*vp1 - 2*vs1*vs1 )
            !   qp (k,i,j) = qp0_metal ! artificially strong attenuation in air-column
            !   qs (k,i,j) = qs0_metal ! artificially strong attenuation in air-column
            ! end if

          else if (yc( j ) <= balldrop_h .and. yc( j ) >= 0.0) then

            !! this is rock sample
            rho(k,i,j) = rho_rock
            mu (k,i,j) = rho(k,i,j) * vs_rock * vs_rock ! rho(k,i,j) * vs1 * vs1
            lam(k,i,j) = rho(k,i,j) * ( vp_rock*vp_rock - 2*vs_rock*vs_rock ) ! rho(k,i,j) * ( vp1*vp1 - 2*vs1*vs1 )
            qp (k,i,j) = qp0_rock ! artificially strong attenuation in air-column
            qs (k,i,j) = qs0_rock ! artificially strong attenuation in air-column

            ! if (xc(i) < -balldrop_lm  .or. xc(i) > balldrop_l + balldrop_lm) then
            !   !! this is air outside of the side metal plates
            !   rho(k,i,j) = 0.001
            !   mu (k,i,j) = rho(k,i,j) * vs1 * vs1
            !   lam(k,i,j) = rho(k,i,j) * ( vp1*vp1 - 2*vs1*vs1 )
            !   qp (k,i,j) = 10.0 ! artificially strong attenuation in air-column
            !   qs (k,i,j) = 10.0 ! artificially strong attenuation in air-column

            !! this is middle rock and side metal plates
            ! else if (xc(i) < 0.0 .or. xc(i) > balldrop_l) then
            !   !! side metal plates
            !   rho(k,i,j) = rho_metal
            !   mu (k,i,j) = rho(k,i,j) * vs_metal * vs_metal ! rho(k,i,j) * vs1 * vs1
            !   lam(k,i,j) = rho(k,i,j) * ( vp_metal*vp_metal - 2*vs_metal*vs_metal ) ! rho(k,i,j) * ( vp1*vp1 - 2*vs1*vs1 )
            !   qp (k,i,j) = qp0_metal ! artificially strong attenuation in air-column
            !   qs (k,i,j) = qs0_metal ! artificially strong attenuation in air-column
            ! else
            !   !! this is rock sample
            !   rho(k,i,j) = rho_rock
            !   mu (k,i,j) = rho(k,i,j) * vs_rock * vs_rock ! rho(k,i,j) * vs1 * vs1
            !   lam(k,i,j) = rho(k,i,j) * ( vp_rock*vp_rock - 2*vs_rock*vs_rock ) ! rho(k,i,j) * ( vp1*vp1 - 2*vs1*vs1 )
            !   qp (k,i,j) = qp0_rock ! artificially strong attenuation in air-column
            !   qs (k,i,j) = qs0_rock ! artificially strong attenuation in air-column
            ! end if

          else
            ! this is middle top air column
            rho(k,i,j) = 0.001
            mu (k,i,j) = rho(k,i,j) * vs1 * vs1
            lam(k,i,j) = rho(k,i,j) * ( vp1*vp1 - 2*vs1*vs1 )
            qp (k,i,j) = 10.0 ! artificially strong attenuation in air-column
            qs (k,i,j) = 10.0 ! artificially strong attenuation in air-column
          endif
          !! --------------------------------------------------------------------
        end do
      end do
    end do

    !! dummy value
    bd(:,:,1:NBD) = -9999

    ! substitute to a dummy variable for avoiding compiler warnings
    dum = xc(i0)
    dum = yc(j0)
    dum = zc(k0)
    dum = vcut

  end subroutine vmodel_balldropseg_sidecoord
  !! --------------------------------------------------------------------------------------------------------------------------- !!


end module m_vmodel_balldropseg_sidecoord
!! ----------------------------------------------------------------------------------------------------------------------------- !!
