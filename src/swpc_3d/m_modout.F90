!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Model parameter output
!!
!! @copyright
!!   Copyright 2013-2020 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!!
!! 2020.12.24 Kurama Okubo developped m_outputmodel.F90
!<
!! ----
#include "m_debug.h"
module m_modout

  !! -- Dependency
  use m_std
  use m_debug
  use m_global
  use m_readini
  use m_system
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
    subroutine output__model( io_prm, xc, yc, zc, rho, lam, mu, Qp, Qs )

      !! -- Arguments
      character(256)        :: fn_out
      integer,  intent(in)  :: io_prm
      real(SP), intent(in)  :: xc              !< x-coordinate location
      real(SP), intent(in)  :: yc              !< y-coordinate location
      real(SP), intent(in)  :: zc              !< z-coordinate location
      real(SP), intent(in)  :: rho             !< mass density [g/cm^3]
      real(SP), intent(in)  :: lam             !< Lame's parameter lambda [ (g/cm^3) * (km/s) ]
      real(SP), intent(in)  :: mu              !< Lame's parameter mu     [ (g/cm^3) * (km/s) ]
      real(SP), intent(in)  :: qp              !< P-wave attenuation
      real(SP), intent(in)  :: qs              !< S-wave attenuation
      integer               :: ncid
      integer               :: nx, ny, nz
      integer               :: dimid_x, dimid_y, dimid_z
      integer               :: varid_x, varid_y, varid_z, varid_r
      integer               :: dimids(2)
      !! --

      ! read file name for output model parameters
      call readini( io_prm, 'fn_model', fn_out, './out/modelparameter.nc')

      ! open netcdf file
      call nc_chk(nf90_create(fn_out, NF90_CLOBBER, ncid))

      ! write parameters

      call nc_chk(nf90_create(fn_out, NF90_CLOBBER, ncid))
      call nc_chk(nf90_def_dim(ncid, 'x', nx, dimid_x))
      call nc_chk(nf90_def_dim(ncid, 'y', ny, dimid_y))
      call nc_chk(nf90_def_dim(ncid, 'z', nz, dimid_z))
      call nc_chk(nf90_def_var(ncid, 'x', NF90_REAL, dimid_x, varid_x))
      call nc_chk(nf90_def_var(ncid, 'y', NF90_REAL, dimid_y, varid_y))
      call nc_chk(nf90_def_var(ncid, 'z', NF90_REAL, dimid_z, varid_z))
      call nc_chk(nf90_put_att(ncid, NF90_GLOBAL, 'title', 'random media'))
      call nc_chk(nf90_put_att(ncid, NF90_GLOBAL, 'ax',ax))
      call nc_chk(nf90_put_att(ncid, NF90_GLOBAL, 'ay',ay))
      call nc_chk(nf90_put_att(ncid, NF90_GLOBAL, 'az',az))
      call nc_chk(nf90_put_att(ncid, NF90_GLOBAL, 'epsilon', epsil))
      select case (ptype)
      case(1)
        call nc_chk(nf90_put_att(ncid, NF90_GLOBAL, 'PSDF type', 'Gaussian'))
      case(2)
        call nc_chk(nf90_put_att(ncid, NF90_GLOBAL, 'PSDF type', 'Exponential'))
      case(3)
        call nc_chk(nf90_put_att(ncid, NF90_GLOBAL, 'PSDF type', 'von Karman'))
        call nc_chk(nf90_put_att(ncid, NF90_GLOBAL, 'kappa',kappa))
      end select

      call nc_chk(nf90_put_att(ncid, NF90_GLOBAL, 'dx', dx))
      call nc_chk(nf90_put_att(ncid, NF90_GLOBAL, 'dy', dy))
      call nc_chk(nf90_put_att(ncid, NF90_GLOBAL, 'dz', dz))

      call nc_chk(nf90_put_att(ncid, varid_x, 'long_name', 'x'))
      call nc_chk(nf90_put_att(ncid, varid_y, 'long_name', 'y'))
      call nc_chk(nf90_put_att(ncid, varid_z, 'long_name', 'z'))
      call nc_chk(nf90_put_att(ncid, varid_x, 'units', 'km'))
      call nc_chk(nf90_put_att(ncid, varid_y, 'units', 'km'))
      call nc_chk(nf90_put_att(ncid, varid_z, 'units', 'km'))
      call nc_chk(nf90_put_att(ncid, varid_x, 'actual_range', (/xx(1), xx(nx)/)))
      call nc_chk(nf90_put_att(ncid, varid_y, 'actual_range', (/yy(1), yy(nx)/)))
      call nc_chk(nf90_put_att(ncid, varid_z, 'actual_range', (/zz(1), zz(nz)/)))
      call nc_chk(nf90_def_var(ncid, "random media", NF90_REAL, (/dimid_x, dimid_y, dimid_z/),  varid_r))
      call nc_chk(nf90_put_att(ncid, varid_r, 'long_name', 'random media'))
      call nc_chk(nf90_put_att(ncid, varid_r, 'units', ''))

      call nc_chk(nf90_enddef(ncid))




      call nc_chk(nf90_create(fn_out, NF90_CLOBBER, ncid))
      call nc_chk(nf90_def_dim(ncid, 'x', nx, dimid_x))
      call nc_chk(nf90_def_dim(ncid, 'y', ny, dimid_y))
      call nc_chk(nf90_def_dim(ncid, 'z', nz, dimid_z))
      call nc_chk(nf90_def_var(ncid, 'x', NF90_REAL, dimid_x, varid_x))
      call nc_chk(nf90_def_var(ncid, 'y', NF90_REAL, dimid_y, varid_y))
      call nc_chk(nf90_def_var(ncid, 'z', NF90_REAL, dimid_z, varid_z))
      call nc_chk(nf90_put_att(ncid, NF90_GLOBAL, 'title', 'random media'))
      call nc_chk(nf90_put_att(ncid, NF90_GLOBAL, 'ax',ax))
      call nc_chk(nf90_put_att(ncid, NF90_GLOBAL, 'ay',ay))
      call nc_chk(nf90_put_att(ncid, NF90_GLOBAL, 'az',az))
      call nc_chk(nf90_put_att(ncid, NF90_GLOBAL, 'epsilon', epsil))



      call nc_chk(nf90_put_att(ncid, NF90_GLOBAL, 'dx', dx))
      call nc_chk(nf90_put_att(ncid, NF90_GLOBAL, 'dy', dy))
      call nc_chk(nf90_put_att(ncid, NF90_GLOBAL, 'dz', dz))

      call nc_chk(nf90_put_att(ncid, varid_x, 'long_name', 'x'))
      call nc_chk(nf90_put_att(ncid, varid_y, 'long_name', 'y'))
      call nc_chk(nf90_put_att(ncid, varid_z, 'long_name', 'z'))
      call nc_chk(nf90_put_att(ncid, varid_x, 'units', 'km'))
      call nc_chk(nf90_put_att(ncid, varid_y, 'units', 'km'))
      call nc_chk(nf90_put_att(ncid, varid_z, 'units', 'km'))
      call nc_chk(nf90_put_att(ncid, varid_x, 'actual_range', (/xx(1), xx(nx)/)))
      call nc_chk(nf90_put_att(ncid, varid_y, 'actual_range', (/yy(1), yy(nx)/)))
      call nc_chk(nf90_put_att(ncid, varid_z, 'actual_range', (/zz(1), zz(nz)/)))
      call nc_chk(nf90_def_var(ncid, "random media", NF90_REAL, (/dimid_x, dimid_y, dimid_z/),  varid_r))
      call nc_chk(nf90_put_att(ncid, varid_r, 'long_name', 'random media'))
      call nc_chk(nf90_put_att(ncid, varid_r, 'units', ''))

      call nc_chk(nf90_enddef(ncid))

      call nc_chk(nf90_put_var(ncid, varid_x, xx(1:nx)))
      call nc_chk(nf90_put_var(ncid, varid_y, yy(1:ny)))
      call nc_chk(nf90_put_var(ncid, varid_z, zz(1:nz)))
      call nc_chk(nf90_put_var(ncid, varid_r, med(1:nx,1:ny,1:nz)))
      call nc_chk(nf90_close(ncid))


    end subroutine output__model
    !! --------------------------------------------------------------------------------------------------------------------------- !!

end module m_modout
!! ----------------------------------------------------------------------------------------------------------------------------- !!
