!     
!     rst_hdf5.f
!     restart from a HDF5 snapshot
!     

#ifdef HDF

      SUBROUTINE rst_hdf5( rst_file, rst_rec, rstfil,                    &
      basenm, mnx1, mnx5, mny1, mny5, mnz1, mnz5,                    &
      H_openr, H_opend, r_filer, r_filed,                            &
      gamma_ad, xzn, yzn, zzn, timej, dt,                            &
      bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,                &
      gridlx, gridly, gridlz, nx, ny, ny0, nyold, nz, nsdim,         &
      nxo, nzo, nxtot, nztot, nstep, igeomx, igeomy, igeomz,         &
      mpi_coordx, mpi_coordy, mpi_coordz,                            &
      pb, rhob, rhoeb, velb, temb, machb, csb2, dk1, dk2,            &      
      psh1, psh2, mm1, mm2, nn, tem0, tem1, alpn, betn, gamn, rbn,   &
      radius, cdenst, shift, fg, rm, brt, me, mp,                    &
      m_0, R_b, c, pc, rho_0a, q0, oang, ipr, ms, nmod, A0,          &
      ncl, ic, jc, kc, mclouds, uwi, vwi, vk, rorb,                  &
      t0, t1, vptbd, rptbd, pptbd, grav, rhoa, pa, physic)
     
      USE hdf5
      USE constantes
      USE parallel
      USE tipos
      USE workarrays
    
#if defined(PARALELO) 
      USE mpi
#endif
      IMPLICIT NONE
#if defined(PARALELO) && defined(THIN)
      include "mpif.h"
#endif

# include "types.h"
      
      integer(is) :: rst_file, rst_rec

#ifdef PARALELO
      integer ierror
#endif
      
      integer :: error, curn, dec, i, j, k, H_openr, H_opend, r_filer, r_filed, curR_R
      integer :: ierr
      character :: c1, c2, c3, c4, c5, gc1, gc2, gc3, gc4, gc5, citers(10), mc1, mc2, mc3, mc4, mc5, mc6, mc7, mc8, mc9

!-------------------------------------------------------Input Variables
      character, intent(in) :: rstfil*(*)
      character, intent(in) :: basenm*(*)
!-------------------------------------------------------Output Variables

      integer(is):: ird, mnx1, mnx5, mny1, mny5, mnz1, mnz5, mn6, mn5, mn1, freqpl,  &
      iextrax, iextraz, nx2, nz2, ixsym, izsym, ny0, igodu, mpi_coordx, mpi_coordy, mpi_coordz

  integer(is) :: nrst, nout1, nx, ny, nyold, nyol, nyh, nz, nyt,        &
                 nsdim, igeomx, igeomy, igeomz, nstep,            &
                 bndmxz, bndmny, bndmxy, bndmnz, bndmxx, bndmnx,  &
                 ipr, ms, nmod, brt, mclouds, ncl, nx0, nz0

  real(rs) :: dt, dk1, dk2, psh1, psh2, temb, machb, csb2, nn,    &
              mm1, mm2, dely, delz, tout1, trst,                  &
              gamma_ad, gridlx, gridly, gridlz, pb, rhob, rhoeb,  & 
              velb, timej, tem0, tem1, &
              uwi, vwi, vk, rorb, oang, t0, t1, vptbd, rptbd, pptbd,     &
              A0, fra(25), frs(25), pha(25)

  real(rs) :: alpn, betn, gamn,        &
              rbn, m_0, R_b, c, pc, shift, rho_0a, q0,     &
              radius, cdenst, fg, rm, me, mp, me_cgs,kb,mp_cgs   

  integer(is):: ic(2000), jc(2000), kc(2000)  
 
  real(rs) :: rhoa(-4:mny5), pa(-4:mny5), grav(-4:mny5)

      real(rs) :: xznl(mnx1), xzn(mnx1), xznr(mnx1), &
      yznl(mny1), yzn(mny1), yznr(mny1), &
      zznl(mnz1), zzn(mnz1), zznr(mnz1)  

      type(Tphysic)::  physic(-4:mnx5,-4:mny5,-4:mnz5)

      INTEGER(HID_T) :: file_id

!     dataspaces
      integer(HID_T) :: d3space, d3filespace, d3memspace
      integer(HID_T) :: dxspace, dxfilespace, dxmemspace
      integer(HID_T) :: dyspace, dyfilespace, dymemspace
      integer(HID_T) :: dzspace, dzfilespace, dzmemspace
      integer(HID_T) :: dblspace, intspace

!     groups
      integer(HID_T) :: group_id

!     datasets
      integer(HID_T) :: d3set, dxset, dyset, dzset, dblset, intset

!     dimensions
      INTEGER(HSIZE_T), DIMENSION(3) :: d3dims, d3totdims
      INTEGER(HSIZE_T), DIMENSION(1) :: dxdims, dxtotdims
      INTEGER(HSIZE_T), DIMENSION(1) :: dydims, dytotdims
      INTEGER(HSIZE_T), DIMENSION(1) :: dzdims, dztotdims
      INTEGER(HSIZE_T), DIMENSION(1) :: dbldims = (/ 1 /), intdims = (/ 1 /)

!     offsets
      INTEGER(HSIZE_T), DIMENSION(3) :: d3offset
      INTEGER(HSIZE_T), DIMENSION(1) :: dxoffset, dyoffset, dzoffset

!     attributes and properties
      integer(HID_T) :: data_transfer, file_access

!     physical grid
      integer(is) :: imin, imax, jmin, jmax, kmin, kmax

!     integer and double data
      integer, dimension(1) :: intdata
      doubleprecision, dimension(1) :: dbldata

!     MPI stuff
      integer :: hdcomm, hdinfo

!     buffers  
      real(rs) :: dbuf(-4:mnx5,-4:mny5,-4:mnz5), dbufx(mnx1), &
                  dbufy(mny1), dbufy5(mny5), dbufz(mnz1)

!-------------------------------------------------------Local Variables
      character(256) :: path_out

      integer(is) :: nxo, nyo, nzo, mpi_coordxo, mpi_coordyo, mpi_coordzo, tid, &
                     nxtot, nztot

      real(rs) :: gridlxo, gridlyo, gridlzo

 tid = 1


!     initialize dimensions
      d3dims = (/nx,ny,nz/)
      dxdims = (/nx/)
      dydims = (/ny/)
      dzdims = (/nz/)
      
!     compute physical domain indices
      imin = 1
      imax = nx
      jmin = 1
      jmax = ny
      kmin = 1
      kmax = nz

!     create filename strings
#ifdef PARALELO

!     MPI stuff
      hdcomm = CART_WORLD
      hdinfo = MPI_INFO_NULL
      call MPI_INFO_CREATE(hdinfo, ierror)

!     CPU coordinates string
      curn=mpi_coords(1)
      mc1=achar(48+(curn-mod(curn,100))/100)
      curn=mod(curn,100)
      mc2=achar(48+(curn-mod(curn,10))/10)
      curn=mod(curn,10)
      mc3=achar(48+curn)
      
      curn=mpi_coords(2)
      mc4=achar(48+(curn-mod(curn,100))/100)
      curn=mod(curn,100)
      mc5=achar(48+(curn-mod(curn,10))/10)
      curn=mod(curn,10)
      mc6=achar(48+curn)
      
      curn=mpi_coords(3)
      mc7=achar(48+(curn-mod(curn,100))/100)
      curn=mod(curn,100)
      mc8=achar(48+(curn-mod(curn,10))/10)
      curn=mod(curn,10)
      mc9=achar(48+curn)
#endif

!     compute offsets and total sizes
#ifdef PARALELO

!      if (nx == nxo .and. nz == nzo) then
         dxoffset(1)  = mpi_coords(1) * nx      
         dxtotdims(1) = mpi_dims(1) * nx
         dyoffset(1)  = mpi_coords(2) * ny       
         dytotdims(1) = mpi_dims(2) * ny !- (ny-nyold)
         dzoffset(1)  = mpi_coords(3) * nz     
         dztotdims(1) = mpi_dims(3) * nz 
!      else
!         dxoffset(1)  = mpi_coords(1) * nx - (nx-nxo)/2
!         dxtotdims(1) = mpi_dims(1) * nx
!         if (dxoffset(1) < 0) then
!            dxoffset(1) = dxoffset(1) + dxtotdims(1)
!         endif
!         dyoffset(1)  = mpi_coords(2) * ny
!         dytotdims(1) = mpi_dims(2) * ny !- (ny-nyold)
!         dzoffset(1)  = mpi_coords(3) * nz - (nz-nzo)/2
!         dztotdims(1) = mpi_dims(3) * nz
!          if (dzoffset(1) < 0) then
!             dzoffset(1) = dzoffset(1) + dztotdims(1)
!         endif
!      endif

      d3offset(1)  = dxoffset(1)
      d3offset(2)  = dyoffset(1)
      d3offset(3)  = dzoffset(1)
      d3totdims(1) = dxtotdims(1)
      d3totdims(2) = dytotdims(1)
      d3totdims(3) = dztotdims(1)

#else

      dxtotdims = dxdims
      dytotdims = dydims
      dztotdims = dzdims
      d3totdims = d3dims

#endif /* offsets and dimensions*/

!     create dataspaces for header
      call h5screate_simple_f(1, intdims, intspace, error)
      call h5screate_simple_f(1, dbldims, dblspace, error)


!     create total dataspaces
      call h5screate_simple_f(3, d3totdims, d3space, error)
      call h5screate_simple_f(1, dxtotdims, dxspace, error)
      if (nsdim.gt.1) then
         call h5screate_simple_f(1, dytotdims, dyspace, error)
         call h5screate_simple_f(1, dztotdims, dzspace, error)
      endif

#ifdef PARALELO
!     create memory spaces for MPI collective output
      call h5screate_simple_f(3, d3dims, d3memspace, error)
      call h5screate_simple_f(1, dxdims, dxmemspace, error)
      if (nsdim.gt.1) then
         call h5screate_simple_f(1, dydims, dymemspace, error)
         call h5screate_simple_f(1, dzdims, dzmemspace, error)
      endif

#endif /* buffers */

!     create transfer properties
#ifdef PARALELO
      call h5pcreate_f(H5P_DATASET_XFER_F, data_transfer, error)
      call h5pset_dxpl_mpio_f(data_transfer, H5FD_MPIO_COLLECTIVE_F, error)
      call h5pcreate_f(H5P_FILE_ACCESS_F, file_access, error)
      call h5pset_fapl_mpio_f(file_access, hdcomm, hdinfo, error)
#endif


#ifdef PARALELO
      call MPI_INFO_FREE(hdinfo, ierror)
#endif


!     create filename
      R_filer = rst_file
      curn=R_filer
      c1=achar(48+(curn-mod(curn,10000))/10000)
      curn=mod(curn,10000)
      c2=achar(48+(curn-mod(curn,1000))/1000)
      curn=mod(curn,1000)
      c3=achar(48+(curn-mod(curn,100))/100)
      curn=mod(curn,100)
      c4=achar(48+(curn-mod(curn,10))/10)
      curn=mod(curn,10)
      c5=achar(48+curn)
      
!if (nuproc == 0) then
      if (cartrank == 0) then
         print*,'[rst_hdf5] Restarting from ',&
         TRIM(outRoot)//rstfil// &
         '-'//c1//c2//c3//c4//c5//'.h5'
         print*,'[rst_hdf5] Record ',rst_rec
      endif

      call h5fopen_f(TRIM(outRoot)//rstfil// &
      '-'//c1//c2//c3//c4//c5//'.h5', &
      H5F_ACC_RDWR_F, file_id, error &
#ifdef PARALELO
      ,access_prp = file_access &
#endif
      )
      
!     open group for the desired iteration
      curR_R = rst_rec
      curn = curR_R
      gc1=achar(48+(curn-mod(curn,10000))/10000)
      curn=mod(curn,10000)
      gc2=achar(48+(curn-mod(curn,1000))/1000)
      curn=mod(curn,1000)
      gc3=achar(48+(curn-mod(curn,100))/100)
      curn=mod(curn,100)
      gc4=achar(48+(curn-mod(curn,10))/10)
      curn=mod(curn,10)
      gc5=achar(48+curn)

      call h5gopen_f(file_id, gc1//gc2//gc3//gc4//gc5, group_id, error)

!     read datasets for different vars

!     headers are read


!     #ifdef PARALELO
!     #define DBLPRM(prm) call h5dopen_f(group_id, #prm, dblset, error); call h5dread_f(dblset, H5T_NATIVE_DOUBLE, dbldata, dbldims, error, xfer_prp = data_transfer); call h5dclose_f(dblset, error); prm = dbldata(1)
!     #define INTPRM(prm) call h5dopen_f(group_id, #prm, intset, error); call h5dread_f(intset, H5T_NATIVE_INTEGER, intdata, intdims, error, xfer_prp = data_transfer); call h5dclose_f(intset, error); prm = intdata(1)
!     #else
!     #define DBLPRM(prm) call h5dopen_f(group_id, #prm, dblset, error); call h5dread_f(dblset, H5T_NATIVE_DOUBLE, dbldata, dbldims, error); call h5dclose_f(dblset, error); prm = dbldata(1)
!     #define INTPRM(prm) call h5dopen_f(group_id, #prm, intset, error); call h5dread_f(intset, H5T_NATIVE_INTEGER, intdata, intdims, error); call h5dclose_f(intset, error); prm = intdata(1)
!     #endif

!     #define prm timej
!     #include "DBLPRM_def.F"   
      call h5dopen_f(group_id, 'timej', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, timej, dbldims, error)
      call h5dclose_f(dblset, error)     

!     #define prm dt
!     #include "DBLPRM_def.F"   
     call h5dopen_f(group_id, 'dt', dblset, error) 
     call h5dread_f(dblset, H5T_NATIVE_DOUBLE, dt, dbldims, error)
     call h5dclose_f(dblset, error)          
      
!     #define prm nsdim
!     #include "INTPRM_def.F"   
      call h5dopen_f(group_id, 'nsdim', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, nsdim, intdims, error)
      call h5dclose_f(intset, error)          

!     #define prm nx
!     #include "INTPRM_def.F"
      call h5dopen_f(group_id, 'nx', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, nx0, intdims, error)
      call h5dclose_f(intset, error)  

!     #define prm nyold
!     #include "INTPRM_def.F"
      call h5dopen_f(group_id, 'nyold', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, nyo, intdims, error)
      call h5dclose_f(intset, error)  

!     #define prm ny
!     #include "INTPRM_def.F"
      call h5dopen_f(group_id, 'ny', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, nyt, intdims, error)
      call h5dclose_f(intset, error)  

      
!     #define prm nz
!     #include "INTPRM_def.F"   
      call h5dopen_f(group_id, 'nz', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, nz0, intdims, error)
      call h5dclose_f(intset, error)  
      
!     #define prm nstep
!     #include "INTPRM_def.F"
      call h5dopen_f(group_id, 'nstep', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, nstep, intdims, error)
      call h5dclose_f(intset, error)  

      if (cartrank == 0) then
         print*,'[rst_hdf5] At iteration ',nstep
         print*,'[rst_hdf5] Time = ',timej     
      endif
!     #define prm bndmnx
!     #include "INTPRM_def.F"
      call h5dopen_f(group_id, 'bndmnx', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, bndmnx, intdims, error)
      call h5dclose_f(intset, error)  
      
!     #define prm bndmxx
!     #include "INTPRM_def.F"
      call h5dopen_f(group_id, 'bndmxx', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, bndmxx, intdims, error)
      call h5dclose_f(intset, error)  
      
!     #define prm bndmny
!     #include "INTPRM_def.F"
      call h5dopen_f(group_id, 'bndmny', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, bndmny, intdims, error)
      call h5dclose_f(intset, error)  
      
!     #define prm bndmxy
!     #include "INTPRM_def.F"
      call h5dopen_f(group_id, 'bndmxy', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, bndmxy, intdims, error)
      call h5dclose_f(intset, error)  

!     #define prm bndmnz
!     #include "INTPRM_def.F"
      call h5dopen_f(group_id, 'bndmnz', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, bndmnz, intdims, error)
      call h5dclose_f(intset, error)  
      
!     #define prm bndmxz
!     #include "INTPRM_def.F"
      call h5dopen_f(group_id, 'bndmxz', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, bndmxz, intdims, error)
      call h5dclose_f(intset, error)  
      
!     #define prm gridlx
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'gridlx', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, gridlxo, dbldims, error)
      call h5dclose_f(dblset, error)     
      
!     #define prm gridly
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'gridly', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, gridlyo, dbldims, error)
      call h5dclose_f(dblset, error)     
      
!     #define prm gridlz
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'gridlz', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, gridlzo, dbldims, error)
      call h5dclose_f(dblset, error)     
      
!     #define prm igeomx
!     #include "INTPRM_def.F"
      call h5dopen_f(group_id, 'igeomx', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, igeomx, intdims, error)
      call h5dclose_f(intset, error)  
      
!     #define prm igeomy
!     #include "INTPRM_def.F"
      call h5dopen_f(group_id, 'igeomy', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, igeomy, intdims, error)
      call h5dclose_f(intset, error)  

!     #define prm igeomz
!     #include "INTPRM_def.F"
      call h5dopen_f(group_id, 'igeomz', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, igeomz, intdims, error)
      call h5dclose_f(intset, error)  
      
!     !#define prm nbproco
!     !#include "INTPRM_def.F"
!     call h5dopen_f(group_id, 'nbproco', intset, error) 
!     call h5dread_f(intset, H5T_NATIVE_INTEGER, nbproco, intdims, error)
!     call h5dclose_f(intset, error)  

!     #define prm mpi_coordx
!     #include "INTPRM_def.F"
      call h5dopen_f(group_id, 'mpi_coordx', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, mpi_coordxo, intdims, error)
      call h5dclose_f(intset, error)  

!     #define prm mpi_coordy
!     #include "INTPRM_def.F"
      call h5dopen_f(group_id, 'mpi_coordy', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, mpi_coordyo, intdims, error)
      call h5dclose_f(intset, error)  

!     #define prm mpi_coordz
!     #include "INTPRM_def.F"
      call h5dopen_f(group_id, 'mpi_coordz', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, mpi_coordzo, intdims, error)
      call h5dclose_f(intset, error)    


!     #define prm gamma_ad
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'gamma_ad', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, gamma_ad, dbldims, error)
      call h5dclose_f(dblset, error) 

!     #define prm pb
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'pb', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, pb, dbldims, error)
      call h5dclose_f(dblset, error) 
      
!     #define prm temb
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'temb', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, temb, dbldims, error)
      call h5dclose_f(dblset, error)

!     #define prm rhob
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'rhob', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, rhob, dbldims, error)
      call h5dclose_f(dblset, error) 
      
!     #define prm rhoeb
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'rhoeb', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, rhoeb, dbldims, error)
      call h5dclose_f(dblset, error) 

!     #define prm velb
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'velb', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, velb, dbldims, error)
      call h5dclose_f(dblset, error) 
      
!     #define prm machb
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'machb', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, machb, dbldims, error)
      call h5dclose_f(dblset, error) 
      
!     #define prm csb2
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'csb2', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, csb2, dbldims, error)
      call h5dclose_f(dblset, error) 
      
!     #define prm dk1
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'dk1', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, dk1, dbldims, error)
      call h5dclose_f(dblset, error) 

!     #define prm dk2
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'dk2', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, dk2, dbldims, error)
      call h5dclose_f(dblset, error)
      
!     #define prm psh1
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'psh1', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, psh1, dbldims, error)
      call h5dclose_f(dblset, error) 

!     #define prm psh2
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'psh2', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, psh2, dbldims, error)
      call h5dclose_f(dblset, error) 
            
!     #define prm mm1
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'mm1', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, mm1, dbldims, error)
      call h5dclose_f(dblset, error) 

!     #define prm mm2
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'mm2', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, mm2, dbldims, error)
      call h5dclose_f(dblset, error) 

!     #define prm nn
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'nn', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, nn, dbldims, error)
      call h5dclose_f(dblset, error) 
      
!     #define prm tem0
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'tem0', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, tem0, dbldims, error)
      call h5dclose_f(dblset, error) 

!     #define prm tem1
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'tem1', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, tem1, dbldims, error)
      call h5dclose_f(dblset, error) 

!     #define prm alpn
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'alpn', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, alpn, dbldims, error)
      call h5dclose_f(dblset, error) 
      
!     #define prm betn
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'betn', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, betn, dbldims, error)
      call h5dclose_f(dblset, error) 
      
!     #define prm gamn
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'gamn', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, gamn, dbldims, error)
      call h5dclose_f(dblset, error) 
      
!     #define prm rbn
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'rbn', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, rbn, dbldims, error)
      call h5dclose_f(dblset, error) 
      
!     #define prm m_0
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'm_0', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, m_0, dbldims, error)
      call h5dclose_f(dblset, error) 
      
!     #define prm me
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'me', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, me, dbldims, error)
      call h5dclose_f(dblset, error) 

!     #define prm mp
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'mp', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, mp, dbldims, error)
      call h5dclose_f(dblset, error) 

!     #define prm R_b
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'R_b', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, R_b, dbldims, error)
      call h5dclose_f(dblset, error) 
      
!     #define prm c
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'c', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, c, dbldims, error)
      call h5dclose_f(dblset, error) 
      
!     #define prm pc
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'pc', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, pc, dbldims, error)
      call h5dclose_f(dblset, error) 
      
!     #define prm shift
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'shift', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, shift, dbldims, error)
      call h5dclose_f(dblset, error) 
      
!     #define prm rho_0a
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'rho_0a', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, rho_0a, dbldims, error)
      call h5dclose_f(dblset, error) 
      
!     #define prm q0
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'q0', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, q0, dbldims, error)
      call h5dclose_f(dblset, error) 

!     #define prm radius
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'radius', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, radius, dbldims, error)
      call h5dclose_f(dblset, error) 

!     #define prm cdenst
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'cdenst', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, cdenst, dbldims, error)
      call h5dclose_f(dblset, error) 

!     #define prm fg
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'fg', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, fg, dbldims, error)
      call h5dclose_f(dblset, error) 

!     #define prm rm
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'rm', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, rm, dbldims, error)
      call h5dclose_f(dblset, error) 

!     #define prm oang
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'oang', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, oang, dbldims, error)
      call h5dclose_f(dblset, error) 
      
!     #define prm ipr
!     #include "INTPRM_def.F"
      call h5dopen_f(group_id, 'ipr', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, ipr, intdims, error)
      call h5dclose_f(intset, error)  
      
!     #define prm ms
!     #include "INTPRM_def.F"
      call h5dopen_f(group_id, 'ms', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, ms, intdims, error)
      call h5dclose_f(intset, error)  

!     #define prm brt
!     #include "INTPRM_def.F"
      call h5dopen_f(group_id, 'brt', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, brt, intdims, error)
      call h5dclose_f(intset, error)  

!     #define prm mclouds
!     #include "INTPRM_def.F"
      call h5dopen_f(group_id, 'mclouds', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, mclouds, intdims, error)
      call h5dclose_f(intset, error)  

!     #define prm ncl
!     #include "INTPRM_def.F"
      call h5dopen_f(group_id, 'ncl', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, ncl, intdims, error)
      call h5dclose_f(intset, error)  

!     #define prm nmod
!     #include "INTPRM_def.F"
      call h5dopen_f(group_id, 'nmod', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, nmod, intdims, error)
      call h5dclose_f(intset, error)  
      
!     #define prm A0
!     #include "DBLPRM_def.F"        
      call h5dopen_f(group_id, 'A0', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, A0, dbldims, error)
      call h5dclose_f(dblset, error) 

      if (basenm(1:1)=='M') then
!     #define prm uwi
!     #include "DBLPRM_def.F"
         call h5dopen_f(group_id, 'uwi', dblset, error) 
         call h5dread_f(dblset, H5T_NATIVE_DOUBLE, uwi, dbldims, error)
         call h5dclose_f(dblset, error) 
         
!     #define prm vwi
!     #include "DBLPRM_def.F"
         call h5dopen_f(group_id, 'vwi', dblset, error) 
         call h5dread_f(dblset, H5T_NATIVE_DOUBLE, vwi, dbldims, error)
         call h5dclose_f(dblset, error) 
         
!     #define prm rorb
!     #include "DBLPRM_def.F"
         call h5dopen_f(group_id, 'rorb', dblset, error) 
         call h5dread_f(dblset, H5T_NATIVE_DOUBLE, rorb, dbldims, error)
         call h5dclose_f(dblset, error) 

!     #define prm vk
!     #include "DBLPRM_def.F"
         call h5dopen_f(group_id, 'vk', dblset, error) 
         call h5dread_f(dblset, H5T_NATIVE_DOUBLE, vk, dbldims, error)
         call h5dclose_f(dblset, error) 

      endif

!     #define prm t0
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 't0', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, t0, dbldims, error)
      call h5dclose_f(dblset, error) 

!     #define prm t1
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 't1', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, t1, dbldims, error)
      call h5dclose_f(dblset, error) 

!     #define prm vptbd
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'vptbd', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, vptbd, dbldims, error)
      call h5dclose_f(dblset, error) 

!     #define prm rptbd
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'rptbd', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, rptbd, dbldims, error)
      call h5dclose_f(dblset, error) 

!     #define prm pptbd
!     #include "DBLPRM_def.F"
      call h5dopen_f(group_id, 'pptbd', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, pptbd, dbldims, error)
      call h5dclose_f(dblset, error) 

      call h5dopen_f(group_id, 'r_filed', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, r_filed, intdims, error)
      call h5dclose_f(intset, error)  
      


!!      print*,mpi_coords(1),mpi_coords(2),mpi_coords(3),'before split',dxoffset(1),dyoffset(1),dzoffset(1)

!      print*,dxtotdims(1),dytotdims(1),dztotdims(1),nxo,nyo,nzo

             
!          dxoffset(1)  = dxoffset(1) - (dxtotdims(1)-nxo)/2
!          dzoffset(1)  = dzoffset(1) - (dztotdims(1)-nzo)/2

!          print*,d3totdims(1),d3totdims(2),d3totdims(3),nxo,nzo,nyo

!    print*,mpi_coords(1),mpi_coords(2),mpi_coords(3),dxoffset(1),dyoffset(1),dzoffset(1)

!          imin=imin+dxoffset(1)
!          imax=imax+dxoffset(1)



!!     if (dxoffset(1) >= (nx-nxo)/2 .and. dxoffset(1) < (nx-nxo)/2 + nxo ) then


!     data is read
!     X
      call h5dopen_f(group_id, 'xzn', dxset, error)
!$$$
!$$$  call flush(6)
!$$$
!$$$
#ifdef PARALELO
      call h5dget_space_f(dxset, dxfilespace, error)
      call h5sselect_hyperslab_f(dxfilespace, H5S_SELECT_SET_F, dxoffset, dxdims, error)
#endif

!$$ print*,mpi_coords(1),mpi_coords(2),mpi_coords(3),'got space for x',imin,imax,dxfilespace,dxoffset,dxdims
!$  call flush(6)
!$$$
!$$$
#ifdef PARALELO
      call h5dread_f(dxset, H5T_NATIVE_DOUBLE, dbufx(imin:imax), dxdims, error &
      ,file_space_id = dxfilespace, mem_space_id = dxmemspace &
      ,xfer_prp = data_transfer)
#else
      call h5dread_f(dxset, H5T_NATIVE_DOUBLE, dbufx(imin:imax), dxdims, error)
#endif
!$$$      
!$$$ print*,mpi_coords(1),mpi_coords(2),mpi_coords(3),'read x'
!$$$  call flush(6)
!$$$
      call h5dclose_f(dxset, error)
#ifdef PARALELO
      call h5sclose_f(dxfilespace, error)
#endif
!$$$ print*,mpi_coords(1),mpi_coords(2),mpi_coords(3),'closed x'
!$$$
!$$$  call flush(6)
!$$$
!$$$
      do i = imin,imax
         xzn(i) = dbufx(i)
      enddo

!      print*,mpi_coords(1),mpi_coords(2),mpi_coords(3),xzn(1),xzn(imax)
!!    endif

!!     if (dzoffset(1) >= (nz-nzo)/2 .and. dzoffset(1) < (nz-nzo)/2 + nzo ) then


!$$$
!$$$
      if (nsdim.gt.2) then
!     Z
         call h5dopen_f(group_id, 'zzn', dzset, error)
!$$$
#ifdef PARALELO
         call h5dget_space_f(dzset, dzfilespace, error)
         call h5sselect_hyperslab_f(dzfilespace, H5S_SELECT_SET_F, dzoffset, dzdims, error)
#endif
#ifdef PARALELO
         call h5dread_f(dzset, H5T_NATIVE_DOUBLE, dbufz(kmin:kmax), dzdims, error &
         ,file_space_id = dzfilespace, mem_space_id = dzmemspace &
         ,xfer_prp = data_transfer)
#else
         call h5dread_f(dzset, H5T_NATIVE_DOUBLE, dbufz(kmin:kmax), dzdims, error)
#endif
         call h5dclose_f(dzset, error)
#ifdef PARALELO
         call h5sclose_f(dzfilespace, error)
#endif
      endif

      do k = kmin,kmax
         zzn(k) = dbufz(k)
      enddo
!$$$    
!$$$
!$$$
!$$$ print*,mpi_coords(1),mpi_coords(2),mpi_coords(3),'read z'
!$$  call flush(6)

!!      endif 

 

!      if (nsdim.gt.1) then

!!      if (dyoffset(1) < nyo ) then

!     Y
         call h5dopen_f(group_id, 'yzn', dyset, error)
!$$$
#ifdef PARALELO
         call h5dget_space_f(dyset, dyfilespace, error)
         call h5sselect_hyperslab_f(dyfilespace, H5S_SELECT_SET_F, dyoffset, dydims, error)
!$$$         
#endif

!      print*,mpi_coords(1),mpi_coords(2),mpi_coords(3),yzn(1),yzn(jmax)
!      print*,dydims,dymemspace,dyfilespace,data_transfer

#ifdef PARALELO
         call h5dread_f(dyset, H5T_NATIVE_DOUBLE, dbufy(jmin:jmax), dydims, error &
         ,file_space_id = dyfilespace, mem_space_id = dymemspace &
         ,xfer_prp = data_transfer)
#else
         call h5dread_f(dyset, H5T_NATIVE_DOUBLE, dbufy(jmin:jmax), dydims, error)
#endif        
         call h5dclose_f(dyset, error)
#ifdef PARALELO
         call h5sclose_f(dyfilespace, error)
#endif
!$$$
!$$$
!      print*,mpi_coords(1),mpi_coords(2),mpi_coords(3),yzn(1),yzn(jmax)  
!      print*,mpi_coords(1),mpi_coords(2),mpi_coords(3),jmin,jmax
!      print*,mpi_coords(1),mpi_coords(2),mpi_coords(3),dyoffset(1)
!      call flush(6)
!$$$

      do j = jmin,jmax
         yzn(j) = dbufy(j)
      enddo

!!      endif

!$$$
!$$$
!$$$!     rhoa
!$$$         call h5dopen_f(group_id, 'rhoa', dyset, error)
!$$$
!$$$#ifdef PARALELO
!$$$         call h5dget_space_f(dyset, dyfilespace, error)
!$$$         call h5sselect_hyperslab_f(dyfilespace, H5S_SELECT_SET_F, dyoffset, dydims, error)
!$$$         
!$$$#endif
!$$$#ifdef PARALELO
!$$$         call h5dread_f(dyset, H5T_NATIVE_DOUBLE, dbufy5(jmin:jmax+5), dydims, error &
!$$$         ,file_space_id = dyfilespace, mem_space_id = dymemspace &
!$$$         ,xfer_prp = data_transfer)
!$$$#else
!$$$         call h5dread_f(dyset, H5T_NATIVE_DOUBLE, dbufy5(jmin:jmax+5), dydims, error)
!$$$#endif        
!$$$         call h5dclose_f(dyset, error)
!$$$#ifdef PARALELO
!$$$         call h5sclose_f(dyfilespace, error)
!$$$#endif
!$$$
!$$$      do j = jmin,jmax+5
!$$$         rhoa(j) = dbufy5(j)
!$$$      enddo
!$$$
!$$$!     pa
!$$$         call h5dopen_f(group_id, 'pa', dyset, error)
!$$$
!$$$#ifdef PARALELO
!$$$         call h5dget_space_f(dyset, dyfilespace, error)
!$$$         call h5sselect_hyperslab_f(dyfilespace, H5S_SELECT_SET_F, dyoffset, dydims, error)
!$$$         
!$$$#endif
!$$$#ifdef PARALELO
!$$$         call h5dread_f(dyset, H5T_NATIVE_DOUBLE, dbufy(jmin:jmax), dydims, error &
!$$$         ,file_space_id = dyfilespace, mem_space_id = dymemspace &
!$$$         ,xfer_prp = data_transfer)
!$$$#else
!$$$         call h5dread_f(dyset, H5T_NATIVE_DOUBLE, dbufy(jmin:jmax), dydims, error)
!$$$#endif        
!$$$         call h5dclose_f(dyset, error)
!$$$#ifdef PARALELO
!$$$         call h5sclose_f(dyfilespace, error)
!$$$#endif
!$$$
!$$$      do j = jmin,jmax
!$$$         pa(j) = dbufy(j)
!$$$      enddo
!$$$
!$$$!     grav
!$$$         call h5dopen_f(group_id, 'grav', dyset, error)
!$$$
!$$$#ifdef PARALELO
!$$$         call h5dget_space_f(dyset, dyfilespace, error)
!$$$         call h5sselect_hyperslab_f(dyfilespace, H5S_SELECT_SET_F, dyoffset, dydims, error)
!$$$         
!$$$#endif
!$$$#ifdef PARALELO
!$$$         call h5dread_f(dyset, H5T_NATIVE_DOUBLE, dbufy(jmin:jmax), dydims, error &
!$$$         ,file_space_id = dyfilespace, mem_space_id = dymemspace &
!$$$         ,xfer_prp = data_transfer)
!$$$#else
!$$$         call h5dread_f(dyset, H5T_NATIVE_DOUBLE, dbufy(jmin:jmax), dydims, error)
!$$$#endif        
!$$$         call h5dclose_f(dyset, error)
!$$$#ifdef PARALELO
!$$$         call h5sclose_f(dyfilespace, error)
!$$$#endif
!      endif
!$$$
!$$$      do j = jmin,jmax
!$$$         grav(j) = dbufy(j)
!$$$      enddo
!
! print*,mpi_coords(1),mpi_coords(2),mpi_coords(3),'read vectors'
!  
!!      print*,nxo,nyo,nzo, 'before vectors'
!!     call MPI_BARRIER( MPI_COMM_WORLD, ierr )

!!      if (dyoffset(1) < nyo .and. dxoffset(1) >= (nx-nxo)/2 &
!!                            .and. dxoffset(1) < (nx-nxo)/2 + nxo  &
!!                            .and. dzoffset(1) >= (nz-nzo)/2 &
!!                            .and. dzoffset(1) < (nz-nzo)/2 + nzo ) then

! print*,mpi_coords(1),mpi_coords(2),mpi_coords(3),'read vectors'
! call flush(6)

!     densty
      call h5dopen_f(group_id, 'densty', d3set, error)
#ifdef PARALELO
      call h5dget_space_f(d3set, d3filespace, error)
      call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)      

      call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error &
      ,file_space_id = d3filespace, mem_space_id = d3memspace &
      ,xfer_prp = data_transfer)
#else
      call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error)
#endif
      call h5dclose_f(d3set, error)
#ifdef PARALELO
      call h5sclose_f(d3filespace, error)
#endif
      
      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax
               physic(i, j, k)%densty = dbuf(i,j,k)
            enddo
         enddo
      enddo

! print*,mpi_coords(1),mpi_coords(2),mpi_coords(3),physic(1,1,1)%densty,'densty'
!  call flush(6)

!     denstye
      call h5dopen_f(group_id, 'denstye', d3set, error)
#ifdef PARALELO
      call h5dget_space_f(d3set, d3filespace, error)
      call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)      

      call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error &
      ,file_space_id = d3filespace, mem_space_id = d3memspace &
      ,xfer_prp = data_transfer)
#else
      call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error)
#endif
      call h5dclose_f(d3set, error)
#ifdef PARALELO
      call h5sclose_f(d3filespace, error)
#endif
      
      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax
               physic(i, j, k)%denstye = dbuf(i,j,k)
            enddo
         enddo
      enddo

! print*,mpi_coords(1),mpi_coords(2),mpi_coords(3),physic(1,1,1)%denstye,'denstye'
!  call flush(6)

!     pres
      call h5dopen_f(group_id, 'pres', d3set, error)

#ifdef PARALELO
      call h5dget_space_f(d3set, d3filespace, error)
      call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)
      
#endif
#ifdef PARALELO
      call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error &
      ,file_space_id = d3filespace, mem_space_id = d3memspace &
      ,xfer_prp = data_transfer)
#else
      call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error)
#endif
      call h5dclose_f(d3set, error)
#ifdef PARALELO
      call h5sclose_f(d3filespace, error)
#endif

      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax
               physic(i, j, k)%pres = dbuf(i,j,k)
            enddo
         enddo
      enddo

! print*,mpi_coords(1),mpi_coords(2),mpi_coords(3),physic(1,1,1)%pres,'pres'
!  call flush(6)

!     velx
      call h5dopen_f(group_id, 'velx', d3set, error)

#ifdef PARALELO
      call h5dget_space_f(d3set, d3filespace, error)
      call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)
      
#endif
#ifdef PARALELO
      call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error &
      ,file_space_id = d3filespace, mem_space_id = d3memspace &
      ,xfer_prp = data_transfer)
#else
      call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error)
#endif
      call h5dclose_f(d3set, error)
#ifdef PARALELO
      call h5sclose_f(d3filespace, error)
#endif
      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax
               physic(i, j, k)%velx = dbuf(i,j,k)
            enddo
         enddo
      enddo

! print*,mpi_coords(1),mpi_coords(2),mpi_coords(3),physic(1,1,1)%velx,'velx'
!  call flush(6)

      if (nsdim.gt.1) then
!     vely
         call h5dopen_f(group_id, 'vely', d3set, error)

#ifdef PARALELO
         call h5dget_space_f(d3set, d3filespace, error)
         call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)
         
#endif
#ifdef PARALELO
         call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error &
         ,file_space_id = d3filespace, mem_space_id = d3memspace &
         ,xfer_prp = data_transfer)
#else
         call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error)
#endif
         call h5dclose_f(d3set, error)
#ifdef PARALELO
         call h5sclose_f(d3filespace, error)
#endif
      endif
      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax
               physic(i, j, k)%vely = dbuf(i,j,k)
            enddo
         enddo
      enddo

! print*,mpi_coords(1),mpi_coords(2),mpi_coords(3),physic(1,1,1)%vely,'vely'
!  call flush(6)

      if (nsdim.gt.2) then
!     velz
         call h5dopen_f(group_id, 'velz', d3set, error)

#ifdef PARALELO
         call h5dget_space_f(d3set, d3filespace, error)
         call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)
         
#endif
#ifdef PARALELO
         call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error &
         ,file_space_id = d3filespace, mem_space_id = d3memspace &
         ,xfer_prp = data_transfer)
#else
         call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error)
#endif
         call h5dclose_f(d3set, error)
#ifdef PARALELO
         call h5sclose_f(d3filespace, error)
#endif
      endif
      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax
               physic(i, j, k)%velz = dbuf(i,j,k)
            enddo
         enddo
      enddo

! print*,mpi_coords(1),mpi_coords(2),mpi_coords(3),physic(1,1,1)%velz,'velz'
!  call flush(6)

!     tracer
      call h5dopen_f(group_id, 'tracer', d3set, error)

#ifdef PARALELO
      call h5dget_space_f(d3set, d3filespace, error)
      call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)
      
#endif
#ifdef PARALELO
      call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error &
      ,file_space_id = d3filespace, mem_space_id = d3memspace &
      ,xfer_prp = data_transfer)
#else
      call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error)
#endif
      call h5dclose_f(d3set, error)
#ifdef PARALELO
      call h5sclose_f(d3filespace, error)
#endif
      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax
               physic(i, j, k)%tracer = dbuf(i,j,k)
            enddo
         enddo
      enddo

! print*,mpi_coords(1),mpi_coords(2),mpi_coords(3),physic(1,1,1)%tracer,'tracer'
!  call flush(6)



!     gravx
      call h5dopen_f(group_id, 'gravx', d3set, error)

#ifdef PARALELO
      call h5dget_space_f(d3set, d3filespace, error)
      call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)
      
#endif
#ifdef PARALELO
      call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error &
      ,file_space_id = d3filespace, mem_space_id = d3memspace &
      ,xfer_prp = data_transfer)
#else
      call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error)
#endif
      call h5dclose_f(d3set, error)
#ifdef PARALELO
      call h5sclose_f(d3filespace, error)
#endif
      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax
               physic(i, j, k)%gravx = dbuf(i,j,k)
            enddo
         enddo
      enddo

  

! print*,mpi_coords(1),mpi_coords(2),mpi_coords(3),physic(1,1,1)%gravx,'gravx'
!  call flush(6)

!     gravy
      call h5dopen_f(group_id, 'gravy', d3set, error)

#ifdef PARALELO
      call h5dget_space_f(d3set, d3filespace, error)
      call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)
      
#endif
#ifdef PARALELO
      call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error &
      ,file_space_id = d3filespace, mem_space_id = d3memspace &
      ,xfer_prp = data_transfer)
#else
      call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error)
#endif
      call h5dclose_f(d3set, error)
#ifdef PARALELO
      call h5sclose_f(d3filespace, error)
#endif
      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax
               physic(i, j, k)%gravy = dbuf(i,j,k)
            enddo
         enddo
      enddo

! print*,mpi_coords(1),mpi_coords(2),mpi_coords(3),physic(1,1,1)%gravy,'gravy'
!  call flush(6)

!  print*,mpi_coords(1),mpi_coords(2),mpi_coords(3),d3offset,d3dims
! print*,mpi_coords(1),mpi_coords(2),mpi_coords(3),imin,imax,jmin,jmax,kmin,kmax

!     gravz
      call h5dopen_f(group_id, 'gravz', d3set, error)

#ifdef PARALELO
      call h5dget_space_f(d3set, d3filespace, error)
      call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)
      
#endif
#ifdef PARALELO
      call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error &
      ,file_space_id = d3filespace, mem_space_id = d3memspace &
      ,xfer_prp = data_transfer)
#else
      call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error)
#endif
      call h5dclose_f(d3set, error)
#ifdef PARALELO
      call h5sclose_f(d3filespace, error)
#endif
      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax
               physic(i, j, k)%gravz = dbuf(i,j,k)
            enddo
         enddo
      enddo

! print*,mpi_coords(1),mpi_coords(2),mpi_coords(3),physic(1,1,1)%gravz,'gravz'
!  call flush(6)

!!      endif

!     #ifdef PARALELO
!     #define DDBLDAT(dat) call h5dopen_f(group_id, #dat, d3set, error); call h5dget_space_f(d3set, d3filespace, error); call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error); call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf, d3dims, error,file_space_id = d3filespace, mem_space_id = d3memspace ,xfer_prp = data_transfer); call h5dclose_f(d3set, error); call h5sclose_f(d3filespace, error); 
!     #else
!     #define DDBLDAT(dat) call h5dopen_f(group_id, #dat, d3set, error); call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf, d3dims, error); call h5dclose_f(d3set, error);
!     #endif


!     #define DDBLDAT(dat) DDBLDAT(dat); do k = kmin, kmax; do j = jmin, jmax; do i = imin, imax; dat(i, j, k) = dbuf(i, j, k); enddo; enddo; enddo

!     DBLDAT(densty)
!     DBLDAT(pres)
!     DBLDAT(velx)
!     DBLDAT(vely)
!     DBLDAT(velz)
!     DBLDAT(tracer)


!      print*,mpi_coords(1),mpi_coords(2),mpi_coords(3),physic(1,1,1)%densty,'before closing'

!     close group
      call h5gclose_f(group_id, error)

#ifdef PARALELO
!     close transfer properties
      call h5pclose_f(data_transfer, error)

!     close file properties
      call h5pclose_f(file_access, error)
#endif


#ifdef PARALELO
!     close memory spaces for MPI collective output
      if (nsdim.gt.1) then
         call h5sclose_f(dymemspace, error)
      endif
      call h5sclose_f(dxmemspace, error)
      call h5sclose_f(d3memspace, error)
#endif /* PARALELO */

!     close total memory spaces
      if (nsdim.gt.1) then
         call h5sclose_f(dyspace, error)
      endif
      call h5sclose_f(dxspace, error)
      call h5sclose_f(d3space, error)

!     create memory spaces for header
      call h5sclose_f(intspace, error)
      call h5sclose_f(dblspace, error)

!     WE WILL USE ONLY ONE GROUP FOR THE TIME BEING. PERUCHO      
!     curR_R = curR_R + 1

!     close file         
      call h5fclose_f(file_id, error)
      
!     prepare for writing a new file
!     if (curR_R.gt.RECORDER_R) then
      if (cartrank == 0) then 
         print*, &
         "Closing file "//TRIM(outRoot)//rstfil//"-" &
         //c1//c2//c3//c4//c5//".h5"
      endif
      R_filer = R_filer + 1
!     R_filed = R_filed + 1
      H_openr = 0
      H_opend = 0
!     else
!     H_openr = 1
!     endif
            
     call MPI_BARRIER( MPI_COMM_WORLD, ierr )

      RETURN
      END

#endif

