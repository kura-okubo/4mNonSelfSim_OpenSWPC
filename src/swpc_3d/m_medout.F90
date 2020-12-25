!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Model parameter output
!!
!! @copyright
!!   Copyright 2013-2020 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!!
!! 2020.12.24 Kurama Okubo (kokubo@bosai.go.jp) developped m_outputmodel.F90
!<
!! ----
#include "m_debug.h"
module m_medout

  !! -- Dependency
  use m_std
  ! use m_debug
  ! use m_global
  use m_readini
  ! use m_system
  use netcdf

  !! -- Declarations
  implicit none
  private
  save

  public :: output__model

  contains


    !! --------------------------------------------------------------------------------------------------------------------------- !!
    !>
    !! output model parameters in 3d
    !<
    !! ----
    subroutine output__model( io_prm, nx, ny, nz, xc, yc, zc, rho, lam, mu, qp, qs )

      !! -- Arguments
      character(256)        :: fn_out
      integer,  intent(in)  :: io_prm
      real(SP), intent(in)  :: xc(nx)            !< x-coordinate location
      real(SP), intent(in)  :: yc(ny)            !< y-coordinate location
      real(SP), intent(in)  :: zc(nz)            !< z-coordinate location
      real(SP), intent(in)  :: rho(nz,nx,ny)      !< mass density [g/cm^3]
      real(SP), intent(in)  :: lam(nz,nx,ny)      !< Lame's parameter lambda [ (g/cm^3) * (km/s) ]
      real(SP), intent(in)  :: mu(nz,nx,ny)       !< Lame's parameter mu     [ (g/cm^3) * (km/s) ]
      real(SP), intent(in)  :: qp(nz,nx,ny)       !< P-wave attenuation
      real(SP), intent(in)  :: qs(nz,nx,ny)       !< S-wave attenuation
      integer               :: ncid
      integer               :: nx, ny, nz
      integer               :: dimid_x, dimid_y, dimid_z
      integer               :: varid_x, varid_y, varid_z
      integer               :: varid_r, varid_l, varid_m, varid_qp, varid_qs
      integer               :: dimids(2)
      !! --

      ! read file name for output model parameters
      call readini( io_prm, 'fn_model', fn_out, './out/modelparameter.nc')

      ! open netcdf file
      call nc_chk(nf90_create(fn_out, NF90_CLOBBER, ncid))
      call nc_chk(nf90_def_dim(ncid, 'x', nx, dimid_x))
      call nc_chk(nf90_def_dim(ncid, 'y', ny, dimid_y))
      call nc_chk(nf90_def_dim(ncid, 'z', nz, dimid_z))
      call nc_chk(nf90_def_var(ncid, 'x', NF90_REAL, dimid_x, varid_x))
      call nc_chk(nf90_def_var(ncid, 'y', NF90_REAL, dimid_y, varid_y))
      call nc_chk(nf90_def_var(ncid, 'z', NF90_REAL, dimid_z, varid_z))
      call nc_chk(nf90_put_att(ncid, NF90_GLOBAL, 'title', 'model parameters'))

      call nc_chk(nf90_put_att(ncid, varid_x, 'long_name', 'x'))
      call nc_chk(nf90_put_att(ncid, varid_y, 'long_name', 'y'))
      call nc_chk(nf90_put_att(ncid, varid_z, 'long_name', 'z'))

      call nc_chk(nf90_put_att(ncid, varid_x, 'units', 'km'))
      call nc_chk(nf90_put_att(ncid, varid_y, 'units', 'km'))
      call nc_chk(nf90_put_att(ncid, varid_z, 'units', 'km'))

      call nc_chk(nf90_put_att(ncid, varid_x, 'actual_range', (/xc(1), xc(nx)/)))
      call nc_chk(nf90_put_att(ncid, varid_y, 'actual_range', (/yc(1), yc(nx)/)))
      call nc_chk(nf90_put_att(ncid, varid_z, 'actual_range', (/zc(1), zc(nz)/)))

      call nc_chk(nf90_def_var(ncid, "rho", NF90_REAL, (/dimid_z, dimid_x, dimid_y/),  varid_r))
      call nc_chk(nf90_def_var(ncid, "lam", NF90_REAL, (/dimid_z, dimid_x, dimid_y/),  varid_l))
      call nc_chk(nf90_def_var(ncid,  "mu", NF90_REAL, (/dimid_z, dimid_x, dimid_y/),  varid_m))
      call nc_chk(nf90_def_var(ncid,  "Qp", NF90_REAL, (/dimid_z, dimid_x, dimid_y/),  varid_qp))
      call nc_chk(nf90_def_var(ncid,  "Qs", NF90_REAL, (/dimid_z, dimid_x, dimid_y/),  varid_qs))

      call nc_chk(nf90_put_att(ncid, varid_r, 'units', 'g/cm^3'))
      call nc_chk(nf90_put_att(ncid, varid_l, 'units', 'GPa'))
      call nc_chk(nf90_put_att(ncid, varid_m, 'units', 'GPa'))
      call nc_chk(nf90_put_att(ncid, varid_qp, 'units', ''))
      call nc_chk(nf90_put_att(ncid, varid_qs, 'units', ''))

      call nc_chk(nf90_enddef(ncid))

      call nc_chk(nf90_put_var(ncid, varid_x, xc(1:nx)))
      call nc_chk(nf90_put_var(ncid, varid_y, yc(1:ny)))
      call nc_chk(nf90_put_var(ncid, varid_z, zc(1:nz)))

      call nc_chk(nf90_put_var(ncid, varid_r, rho(1:nz,1:nx,1:ny)))
      call nc_chk(nf90_put_var(ncid, varid_l, lam(1:nz,1:nx,1:ny)))
      call nc_chk(nf90_put_var(ncid, varid_m,  mu(1:nz,1:nx,1:ny)))
      call nc_chk(nf90_put_var(ncid, varid_qp, qp(1:nz,1:nx,1:ny)))
      call nc_chk(nf90_put_var(ncid, varid_qs, qs(1:nz,1:nx,1:ny)))

      call nc_chk(nf90_close(ncid))

    contains

      !! --------------------------------------------------------------------------------------------------------------------------- !!
      subroutine nc_chk(ierr)

        integer, intent(in) :: ierr

        if(ierr /= NF90_NOERR) then
          write(STDERR,*) NF90_STRERROR(ierr)
          stop
        end if

      end subroutine nc_chk
      !! --------------------------------------------------------------------------------------------------------------------------- !!

    end subroutine output__model
    !! --------------------------------------------------------------------------------------------------------------------------- !!

end module m_medout
!! ----------------------------------------------------------------------------------------------------------------------------- !!
