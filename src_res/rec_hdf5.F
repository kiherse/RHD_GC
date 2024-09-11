!
!     rec_hdf5.f
!     records snapshots into HDF5 files
!

#ifdef HDF

      SUBROUTINE rrhdf5_init(H_openr, H_opend)

      use hdf5
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

#ifdef PARALELO
      integer  ierror
#endif

      integer :: error, H_openr, H_opend

      call h5open_f(error)

      H_openr = 0
      H_opend = 0

      RETURN
      END SUBROUTINE rrhdf5_init 


      SUBROUTINE rrhdf5_record( ird, plane, datfil, rstfil, outfil,           &
               basenm, mnx1, mnx5, mny1, mny5, mnz1, mnz5,                    & 
               H_openr, H_opend, curR_R, r_filer, r_filed,                    &
               gamma_ad, xzn, yzn, zzn, timej,                                &
               bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,                &
               gridlx, gridly, gridlz, nx, ny, ny0, nyold, nz, nsdim,         &
               nstep, igeomx, igeomy, igeomz,                                 &
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
      
#ifdef PARALELO
      integer ierror
#endif

!------------------------------------------------------- Local Variables   
      integer :: error, curn, dec, i, j, k, r_filer, r_filed, curR_R, nx_rec, ny_rec, nz_rec
      character :: c1, c2, c3, c4, c5, gc1, gc2, gc3, gc4, gc5, citers(10), mc1, mc2, mc3, mc4, mc5, mc6, mc7, mc8, mc9
      character(256) :: path2,fich, path_out

!------------------------------------------------------- Input Variables
  character*(*), intent(in) :: datfil
  character*(*), intent(in) :: rstfil, outfil
  character(4), intent(in) :: basenm

  integer:: ird, plane, mnx1, mnx5, mny1, mny5, mnz1, mnz5, mn6, mn5, mn1, freqpl,  &
                iextrax, iextraz, nx2, nz2, ixsym, izsym, ny0, igodu, mpi_coordx,   &
                mpi_coordy, mpi_coordz, uno

  integer(is) :: nrst, nout1, nx, ny, nyold, nyh, nz, nyt,        &
                 nsdim, igeomx, igeomy, igeomz, nstep,            &
                 bndmxz, bndmny, bndmxy, bndmnz, bndmxx, bndmnx,  &
                 ipr, ms, nmod, H_openr, H_opend, brt, mclouds, ncl

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

  type(Tphysic):: physic(-4:mnx5,-4:mny5,-4:mnz5)
  
   
  !real(rs) :: dbufa(nx,ny),dbuft(nx,nz,nbproc)
   real(rs) :: dbufa(nx,ny),dbuft(nx,ny,mpi_dims(3))

      INTEGER(HID_T) :: file_id
      
!     dataspaces
      integer(HID_T) :: d3space, d3filespace, d3memspace
      integer(HID_T) :: d2space, d2filespace, d2memspace
      integer(HID_T) :: d2spacea, d2filespacea, d2memspacea
      integer(HID_T) :: dxspace, dxfilespace, dxmemspace
      integer(HID_T) :: dyspace, dyfilespace, dymemspace
      integer(HID_T) :: dzspace, dzfilespace, dzmemspace
      integer(HID_T) :: dblspace, intspace

!     groups
      integer(HID_T) :: group_id

!     datasets
      integer(HID_T) :: d3set, d2seta, d2set, dxset, dyset, dzset, dblset, intset

!     dimensions
      INTEGER(HSIZE_T), DIMENSION(3) :: d3dims, d3totdims
      INTEGER(HSIZE_T), DIMENSION(2) :: d2dimsa, d2totdimsa
      INTEGER(HSIZE_T), DIMENSION(3) :: d2dims, d2totdims
      INTEGER(HSIZE_T), DIMENSION(1) :: dxdims, dxtotdims
      INTEGER(HSIZE_T), DIMENSION(1) :: dydims, dytotdims
      INTEGER(HSIZE_T), DIMENSION(1) :: dzdims, dztotdims
      INTEGER(HSIZE_T), DIMENSION(1) :: dbldims = (/ 1 /), intdims = (/ 1 /)
    
!     offsets
      INTEGER(HSIZE_T), DIMENSION(3) :: d3offset, d2offset
      INTEGER(HSIZE_T), DIMENSION(2) :: d2offseta
      INTEGER(HSIZE_T), DIMENSION(1) :: dxoffset, dyoffset, dzoffset

!     attributes and properties
      integer(HID_T) :: data_transfer, file_access

!     physical grid
      integer :: imin, imax, jmin, jmax, kmin, kmax

!     integer and double data
      integer, dimension(1) :: intdata
      doubleprecision, dimension(1) :: dbldata

!     MPI stuff
      integer :: hdcomm, hdinfo
  
!     Buffers 3D
      real(rs), allocatable :: physicbuff(:,:,:)

!     initialize dimensions
      d3dims = (/nx,ny,nz/)
      d2dims = (/nx,ny,1/)
      dxdims = (/nx/)
      dydims = (/ny/)
      dzdims = (/nz/)
      uno = 1

      print*, 'ENTERING REC_HDF5, before allocation'
    
!     compute physical domain indices
      if (ird == 0) then 
            plane = 0
            imin = 1
            imax = nx
            jmin = 1
            jmax = ny
            kmin = 1
            kmax = nz
            allocate (physicbuff(imin:imax,jmin:jmax,kmin:kmax))
      endif
      print*, 'ENTERING REC_HDF5, after allocation'
      

      if (cartrank == 0) then
            write(6,*) '---------------------------------------'
            write(6,*) 'REC_HDF, 0: ', ird, imin, imax, kmin,kmax, jmin, jmax, mpi_coordx, mpi_coordy, mpi_coordz 
            !call flush(6)
         endif     
     
!     create filename strings
#ifdef PARALELO
!     MPI stuff
      hdcomm = MPI_COMM_WORLD
      hdinfo = MPI_INFO_NULL
      call MPI_INFO_CREATE(hdinfo, ierror)

#endif

!     compute offsets and total sizes
#ifdef PARALELO

      dxoffset(1)  = mpi_coordx * nx
      dxtotdims(1) = mpi_dims(1)   * nx
     
      dyoffset(1)  = mpi_coordy * ny
      dytotdims(1) = mpi_dims(2)   * ny
      
      dzoffset(1)  = mpi_coordz * nz
      dztotdims(1) = mpi_dims(3)   * nz
           
      d3offset(1)  = dxoffset(1)
      d3offset(2)  = dyoffset(1)
      d3offset(3)  = dzoffset(1)
     
      d2offset(1)  = dxoffset(1)
      d2offset(2)  = dyoffset(1)
      d2offset(3)  = mpi_coordz

      d3totdims(1) = dxtotdims(1)
      d3totdims(2) = dytotdims(1)
      d3totdims(3) = dztotdims(1)

      d2totdims(1) =  mpi_dims(1) * nx
      d2totdims(2) =  mpi_dims(2) * ny
      d2totdims(3) =  mpi_dims(3)

      
#else
      dxtotdims  = dxdims
      dytotdims  = dydims
      dztotdims  = dzdims
      d3totdims  = d3dims

#endif /* offsets and dimensions*/

!     create total dataspaces
      if (ird == 0) then
         call h5screate_simple_f(3, d3totdims, d3space, error)     
      endif

      call h5screate_simple_f(3, d2totdims, d2space, error)    
      call h5screate_simple_f(1, dxtotdims, dxspace, error)
      if (nsdim.gt.1) then
         call h5screate_simple_f(1, dytotdims, dyspace, error)
         call h5screate_simple_f(1, dztotdims, dzspace, error)
      endif

       if (cartrank == 0) then
            write(6,*) 'REC_HDF, 1: ', dxtotdims,dytotdims, dztotdims, d3totdims
            !call flush(6)
       endif 

#ifdef PARALELO
!     create memory spaces for MPI collective output
      if (ird == 0) then
         call h5screate_simple_f(3, d3dims, d3memspace, error)     
         call h5screate_simple_f(3, d2dims, d2memspace, error)
      endif

      call h5screate_simple_f(1, dxdims, dxmemspace, error)
      if (nsdim.gt.1) then
         call h5screate_simple_f(1, dydims, dymemspace, error)
         call h5screate_simple_f(1, dzdims, dzmemspace, error)
      endif

#endif /* buffers */

!     create dataspaces for header
      call h5screate_simple_f(1, intdims, intspace, error)
      call h5screate_simple_f(1, dbldims, dblspace, error)

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
      if (ird == 0) then
         curn=R_filer
      else
         curn=R_filed
      endif
      c1=achar(48+(curn-mod(curn,10000))/10000)
      curn=mod(curn,10000)
      c2=achar(48+(curn-mod(curn,1000))/1000)
      curn=mod(curn,1000)
      c3=achar(48+(curn-mod(curn,100))/100)
      curn=mod(curn,100)
      c4=achar(48+(curn-mod(curn,10))/10)
      curn=mod(curn,10)
      c5=achar(48+curn)

!      print*,H_openr,H_opend

    if (ird == 0) then
!     create new file if necessary
      if (H_openr == 0) then
          if (cartrank == 0) then
             call h5fcreate_f(TRIM(outRoot)//rstfil// &
             '-'//c1//c2//c3//c4//c5//'.h5', &
             H5F_ACC_TRUNC_F, file_id, error &
             )
          endif
    
         curR_R = 1
         
         H_openr = 1
         if (cartrank == 0) then
            call h5fclose_f(file_id, error)
         endif
      endif
      
      if (cartrank == 0) then 
            call h5fopen_f(TRIM(outRoot)//rstfil// &
            '-'//c1//c2//c3//c4//c5//'.h5', &
            H5F_ACC_RDWR_F, file_id, error &
            )
      endif
  
   else
!     create new file if necessary
      if (H_opend == 0) then
          if (cartrank == 0) then
               call h5fcreate_f(TRIM(outRoot)//datfil// &
               '-'//c1//c2//c3//c4//c5//'.h5', &
               H5F_ACC_TRUNC_F, file_id, error &
               )
          endif
      
         curR_R = 1
         
         H_opend = 1
         if (cartrank == 0) then
            call h5fclose_f(file_id, error)
         endif
      endif      
      if (cartrank == 0) then 
        call h5fopen_f(TRIM(outRoot)//datfil// &
        '-'//c1//c2//c3//c4//c5//'.h5', &
        H5F_ACC_RDWR_F, file_id, error &
        )
      endif
   endif

!     add new group for the current iteration
     

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

      if (cartrank == 0) &
      call h5gcreate_f(file_id, gc1//gc2//gc3//gc4//gc5, group_id, error)

     
!      print*,nuproc
      if (cartrank == 0) then
        if (ird == 0) then
!#define prm timej
!#include "DBLPRM_def.F"   
        call h5dcreate_f(group_id, 'timej', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,timej, dbldims, error )
        call h5dclose_f(dblset, error)

!#define prm dt
!#include "DBLPRM_def.F"   
        call h5dcreate_f(group_id, 'dt', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,dt,dbldims, error )
        call h5dclose_f(dblset, error)

!#define prm nsdim
!#include "INTPRM_def.F"   
         call h5dcreate_f(group_id,'nsdim', H5T_NATIVE_INTEGER, intspace, intset, error)
         call h5dwrite_f(intset, H5T_NATIVE_INTEGER,nsdim, intdims, error )
         call h5dclose_f(intset, error)    

         nx_rec = nx * mpi_dims(1)
         ny_rec = ny * mpi_dims(2)
         nz_rec = nz * mpi_dims(3)
             
!     ny = nbproc*ny

!#define prm ny
!#include "INTPRM_def.F"
         call h5dcreate_f(group_id,'nx', H5T_NATIVE_INTEGER, intspace, intset, error)
         call h5dwrite_f(intset, H5T_NATIVE_INTEGER,nx_rec,intdims, error )
         call h5dclose_f(intset, error)
        ! ny = ny/nbproc
!#define prm nyold
!#include "INTPRM_def.F"
         call h5dcreate_f(group_id,'nyold', H5T_NATIVE_INTEGER, intspace, intset, error)
         call h5dwrite_f(intset, H5T_NATIVE_INTEGER,nyold,intdims, error )
         call h5dclose_f(intset, error)        
                
!#define prm ny
!#include "INTPRM_def.F"
         call h5dcreate_f(group_id,'ny', H5T_NATIVE_INTEGER, intspace, intset, error)
         call h5dwrite_f(intset, H5T_NATIVE_INTEGER,ny_rec,intdims, error )
         call h5dclose_f(intset, error)


!#define prm nz
!#include "INTPRM_def.F"
         call h5dcreate_f(group_id,'nz', H5T_NATIVE_INTEGER, intspace, intset, error)
         call h5dwrite_f(intset, H5T_NATIVE_INTEGER,nz_rec,intdims, error )
         call h5dclose_f(intset, error)
 
     
!#define prm nstep
!#include "INTPRM_def.F"
         call h5dcreate_f(group_id,'nstep', H5T_NATIVE_INTEGER, intspace, intset, error)
         call h5dwrite_f(intset, H5T_NATIVE_INTEGER,nstep,intdims, error )
         call h5dclose_f(intset, error)
        
!#define prm bndmnx
!#include "INTPRM_def.F"
         call h5dcreate_f(group_id,'bndmnx', H5T_NATIVE_INTEGER, intspace, intset, error)
         call h5dwrite_f(intset, H5T_NATIVE_INTEGER,bndmnx,intdims, error )
         call h5dclose_f(intset, error)
       
!#define prm bndmxx
!#include "INTPRM_def.F"
         call h5dcreate_f(group_id,'bndmxx', H5T_NATIVE_INTEGER, intspace, intset, error)
         call h5dwrite_f(intset, H5T_NATIVE_INTEGER,bndmxx,intdims, error  )
         call h5dclose_f(intset, error)
        
!#define prm bndmny
!#include "INTPRM_def.F"
         call h5dcreate_f(group_id,'bndmny', H5T_NATIVE_INTEGER, intspace, intset, error)
         call h5dwrite_f(intset, H5T_NATIVE_INTEGER,bndmny,intdims, error )
         call h5dclose_f(intset, error)
       
!#define prm bndmxy
!#include "INTPRM_def.F"
         call h5dcreate_f(group_id,'bndmxy', H5T_NATIVE_INTEGER, intspace, intset, error)
         call h5dwrite_f(intset, H5T_NATIVE_INTEGER,bndmxy,intdims, error )
         call h5dclose_f(intset, error)
        
!#define prm bndmnz
!#include "INTPRM_def.F"
         call h5dcreate_f(group_id,'bndmnz', H5T_NATIVE_INTEGER, intspace, intset, error)
         call h5dwrite_f(intset, H5T_NATIVE_INTEGER,bndmnz,intdims, error )
         call h5dclose_f(intset, error)
       
!#define prm bndmxz
!#include "INTPRM_def.F"
         call h5dcreate_f(group_id,'bndmxz', H5T_NATIVE_INTEGER, intspace, intset, error)
         call h5dwrite_f(intset, H5T_NATIVE_INTEGER,bndmxz,intdims, error )
         call h5dclose_f(intset, error)
       
!#define prm gridlx
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'gridlx', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,gridlx,dbldims, error )
        call h5dclose_f(dblset, error) 

!#define prm gridly
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'gridly', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,gridly,dbldims, error )
        call h5dclose_f(dblset, error) 
       
!#define prm gridlz
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'gridlz', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,gridlz,dbldims, error )
        call h5dclose_f(dblset, error) 
   
!#define prm igeomx
!#include "INTPRM_def.F"
        call h5dcreate_f(group_id,'igeomx', H5T_NATIVE_INTEGER, intspace, intset, error)
        call h5dwrite_f(intset, H5T_NATIVE_INTEGER,igeomx,intdims, error )
        call h5dclose_f(intset, error)
        
!#define prm igeomy
!#include "INTPRM_def.F"
        call h5dcreate_f(group_id,'igeomy', H5T_NATIVE_INTEGER, intspace, intset, error)
        call h5dwrite_f(intset, H5T_NATIVE_INTEGER,igeomy,intdims, error )
        call h5dclose_f(intset, error)

!#define prm igeomz
!#include "INTPRM_def.F"
        call h5dcreate_f(group_id,'igeomz', H5T_NATIVE_INTEGER, intspace, intset, error)
        call h5dwrite_f(intset, H5T_NATIVE_INTEGER,igeomz,intdims, error )
        call h5dclose_f(intset, error)
       

!#define prm mpi_coordx
!#include "INTPRM_def.F"
        call h5dcreate_f(group_id,'mpi_coordx', H5T_NATIVE_INTEGER, intspace, intset, error)
        call h5dwrite_f(intset, H5T_NATIVE_INTEGER,mpi_coordx,intdims, error )
        call h5dclose_f(intset, error)

!#define prm mpi_coordy
!#include "INTPRM_def.F"
        call h5dcreate_f(group_id,'mpi_coordy', H5T_NATIVE_INTEGER, intspace, intset, error)
        call h5dwrite_f(intset, H5T_NATIVE_INTEGER,mpi_coordy,intdims, error )
        call h5dclose_f(intset, error)

!#define prm mpi_coordz
!#include "INTPRM_def.F"
        call h5dcreate_f(group_id,'mpi_coordz', H5T_NATIVE_INTEGER, intspace, intset, error)
        call h5dwrite_f(intset, H5T_NATIVE_INTEGER,mpi_coordz,intdims, error )
        call h5dclose_f(intset, error)


!#define prm gamma_ad
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'gamma_ad', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,gamma_ad,dbldims, error )
        call h5dclose_f(dblset, error) 
   

!#define prm pb
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'pb', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,pb,dbldims, error )
        call h5dclose_f(dblset, error) 
   
!#define prm temb
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'temb', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,temb,dbldims, error )
        call h5dclose_f(dblset, error) 

!#define prm rhob
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'rhob', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,rhob,dbldims, error )
        call h5dclose_f(dblset, error) 

!#define prm rhoeb
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'rhoeb', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,rhoeb,dbldims, error )
        call h5dclose_f(dblset, error) 
    
!#define prm velb
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'velb', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,velb,dbldims, error )
        call h5dclose_f(dblset, error) 
       
!#define prm machb
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'machb', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,machb,dbldims, error )
        call h5dclose_f(dblset, error) 

!#define prm csb2
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'csb2', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,csb2,dbldims, error )
        call h5dclose_f(dblset, error) 
    
!#define prm dk1
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'dk1', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,dk1,dbldims, error )
        call h5dclose_f(dblset, error) 
       

!#define prm dk2
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'dk2', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,dk2,dbldims, error )
        call h5dclose_f(dblset, error) 

!#define prm psh1
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'psh1', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,psh1,dbldims, error )
        call h5dclose_f(dblset, error) 
       

!#define prm psh2
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'psh2', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,psh2,dbldims, error )
        call h5dclose_f(dblset, error) 


!#define prm mm1
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'mm1', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,mm1,dbldims, error )
        call h5dclose_f(dblset, error) 
   
!#define prm mm2
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'mm2', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,mm2,dbldims, error )
        call h5dclose_f(dblset, error) 

!#define prm nn
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'nn', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,nn,dbldims, error )
        call h5dclose_f(dblset, error) 


!#define prm tem0
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'tem0', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,tem0,dbldims, error )
        call h5dclose_f(dblset, error) 

!#define prm tem1
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'tem1', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,tem1,dbldims, error )
        call h5dclose_f(dblset, error) 
      
!#define prm alpn
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'alpn', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,alpn,dbldims, error )
        call h5dclose_f(dblset, error) 

!#define prm betn
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'betn', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,betn,dbldims, error )
        call h5dclose_f(dblset, error) 
        
!#define prm gamn
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'gamn', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,gamn,dbldims, error )
        call h5dclose_f(dblset, error) 
     
!#define prm rbn
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'rbn', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,rbn,dbldims, error )
        call h5dclose_f(dblset, error) 
       
!#define prm m_0
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'm_0', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,m_0,dbldims, error )
        call h5dclose_f(dblset, error) 

!#define prm me
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'me', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,me,dbldims, error )
        call h5dclose_f(dblset, error) 

!#define prm mp
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'mp', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,mp,dbldims, error )
        call h5dclose_f(dblset, error) 

!#define prm R_b
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'R_b', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,R_b,dbldims, error )
        call h5dclose_f(dblset, error) 
        
!#define prm c
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'c', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,c,dbldims, error )
        call h5dclose_f(dblset, error) 
    
!#define prm pc
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'pc', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,pc,dbldims, error )
        call h5dclose_f(dblset, error) 
       
!#define prm shift
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'shift', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,shift,dbldims, error )
        call h5dclose_f(dblset, error) 
        
!#define prm rho_0a
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'rho_0a', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,rho_0a,dbldims, error )
        call h5dclose_f(dblset, error) 
        
!#define prm q0
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'q0', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,q0,dbldims, error )
        call h5dclose_f(dblset, error) 

!#define prm radius
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'radius', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,radius,dbldims, error )
        call h5dclose_f(dblset, error) 

!#define prm cdenst
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'cdenst', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,cdenst,dbldims, error )
        call h5dclose_f(dblset, error) 

!#define prm fg
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'fg', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,fg,dbldims, error )
        call h5dclose_f(dblset, error) 

!#define prm rm
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'rm', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,rm,dbldims, error )
        call h5dclose_f(dblset, error) 

!#define prm oang
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'oang', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,oang,dbldims, error )
        call h5dclose_f(dblset, error) 

!#define prm ipr
!#include "INTPRM_def.F"
         call h5dcreate_f(group_id,'ipr', H5T_NATIVE_INTEGER, intspace, intset, error)
         call h5dwrite_f(intset, H5T_NATIVE_INTEGER,ipr,intdims, error )
         call h5dclose_f(intset, error)
     
!#define prm ms
!#include "INTPRM_def.F"
         call h5dcreate_f(group_id,'ms', H5T_NATIVE_INTEGER, intspace, intset, error)
         call h5dwrite_f(intset, H5T_NATIVE_INTEGER,ms,intdims, error )
         call h5dclose_f(intset, error)

!#define prm brt
!#include "INTPRM_def.F"
         call h5dcreate_f(group_id,'brt', H5T_NATIVE_INTEGER, intspace, intset, error)
         call h5dwrite_f(intset, H5T_NATIVE_INTEGER,brt,intdims, error )
         call h5dclose_f(intset, error)
       
!#define prm mclouds
!#include "INTPRM_def.F"
         call h5dcreate_f(group_id,'mclouds', H5T_NATIVE_INTEGER, intspace, intset, error)
         call h5dwrite_f(intset, H5T_NATIVE_INTEGER,mclouds,intdims, error )
         call h5dclose_f(intset, error)

!#define prm ncl
!#include "INTPRM_def.F"
         call h5dcreate_f(group_id,'ncl', H5T_NATIVE_INTEGER, intspace, intset, error)
         call h5dwrite_f(intset, H5T_NATIVE_INTEGER,ncl,intdims, error )
         call h5dclose_f(intset, error)

!#define prm nmod
!#include "INTPRM_def.F"
         call h5dcreate_f(group_id,'nmod', H5T_NATIVE_INTEGER, intspace, intset, error)
         call h5dwrite_f(intset, H5T_NATIVE_INTEGER,nmod,intdims, error )
         call h5dclose_f(intset, error)
     
!#define prm A0
!#include "DBLPRM_def.F"        
        call h5dcreate_f(group_id, 'A0', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,A0,dbldims, error )
        call h5dclose_f(dblset, error) 
     

        if (basenm(1:1)=='M') then
!#define prm uwi
!#include "DBLPRM_def.F"
          call h5dcreate_f(group_id, 'uwi', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
          call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,uwi,dbldims, error )
          call h5dclose_f(dblset, error) 
        
!#define prm vwi
!#include "DBLPRM_def.F"
          call h5dcreate_f(group_id, 'vwi', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
          call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,vwi,dbldims, error  )
          call h5dclose_f(dblset, error) 
        
!#define prm rorb
!#include "DBLPRM_def.F"
          call h5dcreate_f(group_id, 'rorb', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
          call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,rorb,dbldims, error )
          call h5dclose_f(dblset, error) 
        endif

!#define prm vk
!#include "DBLPRM_def.F"
          call h5dcreate_f(group_id, 'vk', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
          call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,vk,dbldims, error  )
          call h5dclose_f(dblset, error) 
       endif


!#define prm t0
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 't0', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,t0,dbldims, error)
        call h5dclose_f(dblset, error) 

!#define prm t1
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 't1', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,t1,dbldims, error  )
        call h5dclose_f(dblset, error) 

!#define prm vptbd
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'vptbd', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,vptbd,dbldims, error &
             )
        call h5dclose_f(dblset, error) 

!#define prm rptbd
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id,'rptbd', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,rptbd,dbldims, error &
             )
        call h5dclose_f(dblset, error) 

!#define prm pptbd
!#include "DBLPRM_def.F"
        call h5dcreate_f(group_id, 'pptbd', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
        call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,pptbd,dbldims, error  )
        call h5dclose_f(dblset, error) 

        call h5dcreate_f(group_id,'r_filed', H5T_NATIVE_INTEGER, intspace, intset, error)
        call h5dwrite_f(intset, H5T_NATIVE_INTEGER,r_filed,intdims, error )
        call h5dclose_f(intset, error)



        

!     close group
        call h5gclose_f(group_id, error)
!     close file         
        call h5fclose_f(file_id, error)
      endif
 
 
!      print*,'after pars',nuproc

      if (ird == 0) then
#ifdef PARALELO
        call h5fopen_f(TRIM(outRoot)//rstfil// &
        '-'//c1//c2//c3//c4//c5//'.h5', &
        H5F_ACC_RDWR_F, file_id, error &
        ,access_prp = file_access )
#else
        call h5fopen_f(TRIM(outRoot)//rstfil// &
        '-'//c1//c2//c3//c4//c5//'.h5', &
        H5F_ACC_RDWR_F, file_id, error )
#endif
      else
#ifdef PARALELO
        call h5fopen_f(TRIM(outRoot)//datfil// &
        '-'//c1//c2//c3//c4//c5//'.h5', &
        H5F_ACC_RDWR_F, file_id, error &
        ,access_prp = file_access )
#else
        call h5fopen_f(TRIM(outRoot)//datfil// &
        '-'//c1//c2//c3//c4//c5//'.h5', &
        H5F_ACC_RDWR_F, file_id, error )
#endif
      endif
      
      call h5gOPEN_f(file_id, gc1//gc2//gc3//gc4//gc5, group_id, error)

 if (cartrank == 0) then
            write(6,*) 'REC_HDF, 2: ',imin, imax, kmin,kmax, jmin, jmax 
            !!call flush(6)
         endif 

!     data is written
      IF (ird == 0) THEN

#ifdef PARALELO
 
!   X
         !call h5dcreate_f(group_id, 'xzn', H5T_NATIVE_DOUBLE, dxspace, dxset, error)    
         !call h5dwrite_f(dxset, H5T_NATIVE_DOUBLE, xzn(imin:imax), dxdims, error)
         !call h5dclose_f(dxset, error)
         
         call h5dcreate_f(group_id, 'xzn', H5T_NATIVE_DOUBLE, dxspace, dxset, error)
        
         write(6,*) 'REC_HDF, 2.2: ',cartrank,  dxoffset, dxdims,dxset
         !!call flush(6)
           
         call h5dget_space_f(dxset, dxfilespace, error)
         call h5sselect_hyperslab_f(dxfilespace, H5S_SELECT_SET_F, dxoffset, dxdims, error)
         call h5dwrite_f(dxset, H5T_NATIVE_DOUBLE, xzn(imin:imax), dxdims, error &
         ,file_space_id = dxfilespace, mem_space_id = dxmemspace &
         ,xfer_prp = data_transfer )

!   Z
         !call h5dcreate_f(group_id, 'zzn', H5T_NATIVE_DOUBLE, dzspace, dzset, error)
         !call h5dwrite_f(dzset, H5T_NATIVE_DOUBLE, zzn(kmin:kmax), dzdims, error)
         !call h5dclose_f(dzset, error)  
         call h5dcreate_f(group_id, 'zzn', H5T_NATIVE_DOUBLE, dzspace, dzset, error)
         call h5dget_space_f(dzset, dzfilespace, error)
         call h5sselect_hyperslab_f(dzfilespace, H5S_SELECT_SET_F, dzoffset, dzdims, error)
         call h5dwrite_f(dzset, H5T_NATIVE_DOUBLE, zzn(kmin:kmax), dzdims, error &
         ,file_space_id = dzfilespace, mem_space_id = dzmemspace &
         ,xfer_prp = data_transfer )

!   Y
         call h5dcreate_f(group_id, 'yzn', H5T_NATIVE_DOUBLE, dyspace, dyset, error)
         call h5dget_space_f(dyset, dyfilespace, error)
         call h5sselect_hyperslab_f(dyfilespace, H5S_SELECT_SET_F, dyoffset, dydims, error)
         call h5dwrite_f(dyset, H5T_NATIVE_DOUBLE, yzn(jmin:jmax), dydims, error &
         ,file_space_id = dyfilespace, mem_space_id = dymemspace &
         ,xfer_prp = data_transfer )

#else
         call h5dwrite_f(dxset, H5T_NATIVE_DOUBLE, xzn(imin:imax), dxdims, error )
         call h5dwrite_f(dyset, H5T_NATIVE_DOUBLE, yzn(jmin:jmax), dydims, error )
         call h5dwrite_f(dzset, H5T_NATIVE_DOUBLE, zzn(kmin:kmax), dzdims, error )
#endif
         call h5dclose_f(dxset, error)
         call h5dclose_f(dyset, error)
         call h5dclose_f(dzset, error)
#ifdef PARALELO
         call h5sclose_f(dxfilespace, error)
         call h5sclose_f(dyfilespace, error)
         call h5sclose_f(dzfilespace, error)
#endif

!   rhoa
         call h5dcreate_f(group_id, 'rhoa', H5T_NATIVE_DOUBLE, dyspace, dyset, error)

#ifdef PARALELO
         call h5dget_space_f(dyset, dyfilespace, error)
         call h5sselect_hyperslab_f(dyfilespace, H5S_SELECT_SET_F, dyoffset, dydims, error)
      
#endif
         call h5dwrite_f(dyset, H5T_NATIVE_DOUBLE, rhoa(jmin:jmax), dydims, error &
#ifdef PARALELO
         ,file_space_id = dyfilespace, mem_space_id = dymemspace &
         ,xfer_prp = data_transfer &
#endif
              )
         call h5dclose_f(dyset, error)
#ifdef PARALELO
         call h5sclose_f(dyfilespace, error)
#endif

!   pa
         call h5dcreate_f(group_id, 'pa', H5T_NATIVE_DOUBLE, dyspace, dyset, error)

#ifdef PARALELO
         call h5dget_space_f(dyset, dyfilespace, error)
         call h5sselect_hyperslab_f(dyfilespace, H5S_SELECT_SET_F, dyoffset, dydims, error)
      
#endif
         call h5dwrite_f(dyset, H5T_NATIVE_DOUBLE, pa(jmin:jmax), dydims, error &
#ifdef PARALELO
         ,file_space_id = dyfilespace, mem_space_id = dymemspace &
         ,xfer_prp = data_transfer &
#endif
              )
         call h5dclose_f(dyset, error)
#ifdef PARALELO
         call h5sclose_f(dyfilespace, error)
#endif      


!   grav
         call h5dcreate_f(group_id, 'grav', H5T_NATIVE_DOUBLE, dyspace, dyset, error)

#ifdef PARALELO
         call h5dget_space_f(dyset, dyfilespace, error)
         call h5sselect_hyperslab_f(dyfilespace, H5S_SELECT_SET_F, dyoffset, dydims, error)
      
#endif
         call h5dwrite_f(dyset, H5T_NATIVE_DOUBLE, grav(jmin:jmax), dydims, error &
#ifdef PARALELO
         ,file_space_id = dyfilespace, mem_space_id = dymemspace &
         ,xfer_prp = data_transfer &
#endif
              )
         call h5dclose_f(dyset, error)
#ifdef PARALELO
         call h5sclose_f(dyfilespace, error)
#endif
      !endif      
  

!   denstya
      call h5dcreate_f(group_id, 'denstya', H5T_NATIVE_DOUBLE, d2space, d2set, error)       

!      if (mpi_coordz == mpi_dims(3)/2) then
!         dbuft(imin:imax, jmin:jmax,mpi_coordz) =  physic(imin:imax, jmin:jmax, kmax)%densty
!      else if (mpi_coordz == mpi_dims(3)/2 + 1) then

      dbuft(imin:imax, jmin:jmax, 1) =  physic(imin:imax, jmin:jmax, 1)%densty

      print*,cartrank,imin,imax,jmin,jmax,d2offset(1),d2offset(2),d2offset(3),d2space,d2set,d2dims
!      else
!         dbuft(imin:imax, jmin:jmax,mpi_coordz) =  physic(imin:imax, jmin:jmax, kmax)%densty
!      endif

#ifdef PARALELO
      call h5dget_space_f(d2set, d2filespace, error)
      call h5sselect_hyperslab_f(d2filespace, H5S_SELECT_SET_F, d2offset, d2dims, error)
      
#endif
      call h5dwrite_f(d2set, H5T_NATIVE_DOUBLE, dbuft, d2dims, error &
#ifdef PARALELO
           ,file_space_id = d2filespace, mem_space_id = d2memspace &
           ,xfer_prp = data_transfer &
#endif
           )
      call h5dclose_f(d2set, error)
#ifdef PARALELO
      call h5sclose_f(d2filespace, error)
#endif




 if (cartrank == 0) then
            write(6,*) 'REC_HDF, 3: ', imin, imax, kmin,kmax, jmin, jmax  
            !!call flush(6)
         endif 

!   densty
      call h5dcreate_f(group_id, 'densty', H5T_NATIVE_DOUBLE, d3space, d3set, error)

#ifdef PARALELO
      call h5dget_space_f(d3set, d3filespace, error)
      call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)
      
#endif
      physicbuff(imin:imax, jmin:jmax, kmin:kmax) =  physic(imin:imax, jmin:jmax, kmin:kmax)%densty
      call h5dwrite_f(d3set, H5T_NATIVE_DOUBLE, physicbuff(imin:imax, jmin:jmax, kmin:kmax), d3dims, error &
#ifdef PARALELO
           ,file_space_id = d3filespace, mem_space_id = d3memspace &
           ,xfer_prp = data_transfer &
#endif
           )
      call h5dclose_f(d3set, error)
#ifdef PARALELO
      call h5sclose_f(d3filespace, error)
#endif


!   denstye
      call h5dcreate_f(group_id, 'denstye', H5T_NATIVE_DOUBLE, d3space, d3set, error)

#ifdef PARALELO
      call h5dget_space_f(d3set, d3filespace, error)
      call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)
      
#endif
      physicbuff(imin:imax, jmin:jmax, kmin:kmax) =  physic(imin:imax, jmin:jmax, kmin:kmax)%denstye
      call h5dwrite_f(d3set, H5T_NATIVE_DOUBLE, physicbuff(imin:imax, jmin:jmax, kmin:kmax), d3dims, error &
#ifdef PARALELO
           ,file_space_id = d3filespace, mem_space_id = d3memspace &
           ,xfer_prp = data_transfer &
#endif
           )
      call h5dclose_f(d3set, error)
#ifdef PARALELO
      call h5sclose_f(d3filespace, error)
#endif

!     pres
      call h5dcreate_f(group_id, 'pres', H5T_NATIVE_DOUBLE, d3space, d3set, error)

#ifdef PARALELO
      call h5dget_space_f(d3set, d3filespace, error)
      call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)
      
#endif
      physicbuff(imin:imax, jmin:jmax, kmin:kmax) =  physic(imin:imax, jmin:jmax, kmin:kmax)%pres
      call h5dwrite_f(d3set, H5T_NATIVE_DOUBLE, physicbuff(imin:imax, jmin:jmax, kmin:kmax), d3dims, error &
#ifdef PARALELO
           ,file_space_id = d3filespace, mem_space_id = d3memspace &
           ,xfer_prp = data_transfer &
#endif
          )
      call h5dclose_f(d3set, error)
#ifdef PARALELO
      call h5sclose_f(d3filespace, error)
#endif


!     velx
      call h5dcreate_f(group_id, 'velx', H5T_NATIVE_DOUBLE, d3space, d3set, error)

#ifdef PARALELO
      call h5dget_space_f(d3set, d3filespace, error)
      call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)
      
#endif
       physicbuff(imin:imax, jmin:jmax, kmin:kmax) =  physic(imin:imax, jmin:jmax, kmin:kmax)%velx
      call h5dwrite_f(d3set, H5T_NATIVE_DOUBLE, physicbuff(imin:imax, jmin:jmax, kmin:kmax), d3dims, error &
#ifdef PARALELO
           ,file_space_id = d3filespace, mem_space_id = d3memspace &
           ,xfer_prp = data_transfer &
#endif
           )
      call h5dclose_f(d3set, error)
#ifdef PARALELO
      call h5sclose_f(d3filespace, error)
#endif

      if (nsdim.gt.1) then
!     vely
         call h5dcreate_f(group_id, 'vely', H5T_NATIVE_DOUBLE, d3space, d3set, error)

#ifdef PARALELO
        call h5dget_space_f(d3set, d3filespace, error)
        call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)
      
#endif
         physicbuff(imin:imax, jmin:jmax, kmin:kmax) =  physic(imin:imax, jmin:jmax, kmin:kmax)%vely
      call h5dwrite_f(d3set, H5T_NATIVE_DOUBLE, physicbuff(imin:imax, jmin:jmax, kmin:kmax), d3dims, error &
#ifdef PARALELO
           ,file_space_id = d3filespace, mem_space_id = d3memspace &
           ,xfer_prp = data_transfer &
#endif
           )
        call h5dclose_f(d3set, error)
#ifdef PARALELO
        call h5sclose_f(d3filespace, error)
#endif

      endif

      if (nsdim.gt.2) then
!     velz
        call h5dcreate_f(group_id, 'velz', H5T_NATIVE_DOUBLE, d3space, d3set, error)

#ifdef PARALELO
        call h5dget_space_f(d3set, d3filespace, error)
        call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error) 
#endif
         physicbuff(imin:imax, jmin:jmax, kmin:kmax) =  physic(imin:imax, jmin:jmax, kmin:kmax)%velz
      call h5dwrite_f(d3set, H5T_NATIVE_DOUBLE, physicbuff(imin:imax, jmin:jmax, kmin:kmax), d3dims, error &
#ifdef PARALELO
           ,file_space_id = d3filespace, mem_space_id = d3memspace &
           ,xfer_prp = data_transfer &
#endif
           )
        call h5dclose_f(d3set, error)
#ifdef PARALELO
        call h5sclose_f(d3filespace, error)
#endif
      endif

!     tracer
      call h5dcreate_f(group_id, 'tracer', H5T_NATIVE_DOUBLE, d3space, d3set, error)

#ifdef PARALELO
      call h5dget_space_f(d3set, d3filespace, error)
      call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)
      
#endif
      physicbuff(imin:imax, jmin:jmax, kmin:kmax) =  physic(imin:imax, jmin:jmax, kmin:kmax)%tracer
      call h5dwrite_f(d3set, H5T_NATIVE_DOUBLE, physicbuff(imin:imax, jmin:jmax, kmin:kmax), d3dims, error &
#ifdef PARALELO
           ,file_space_id = d3filespace, mem_space_id = d3memspace &
           ,xfer_prp = data_transfer &
#endif
           )
      call h5dclose_f(d3set, error)
#ifdef PARALELO
      call h5sclose_f(d3filespace, error)
#endif

      if (basenm(1:1) == 'J' .and. basenm(2:2) == '3') then
!         if (c1 == '0' .and. c2 == '0' .and. c3 == '0' .and. c4 == '0' .and. c5 == '0') then
!     gravx
            call h5dcreate_f(group_id, 'gravx', H5T_NATIVE_DOUBLE, d3space, d3set, error)

#ifdef PARALELO
            call h5dget_space_f(d3set, d3filespace, error)
            call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)
      
#endif
            physicbuff(imin:imax, jmin:jmax, kmin:kmax) =  physic(imin:imax, jmin:jmax, kmin:kmax)%gravx
            call h5dwrite_f(d3set, H5T_NATIVE_DOUBLE, physicbuff(imin:imax, jmin:jmax, kmin:kmax), d3dims, error &
#ifdef PARALELO
             ,file_space_id = d3filespace, mem_space_id = d3memspace &
             ,xfer_prp = data_transfer &
#endif
             )
           call h5dclose_f(d3set, error)
#ifdef PARALELO
           call h5sclose_f(d3filespace, error)
#endif

!     gravy
           call h5dcreate_f(group_id, 'gravy', H5T_NATIVE_DOUBLE, d3space, d3set, error)

#ifdef PARALELO
           call h5dget_space_f(d3set, d3filespace, error)
           call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)
      
#endif
           physicbuff(imin:imax, jmin:jmax, kmin:kmax) =  physic(imin:imax, jmin:jmax, kmin:kmax)%gravy
           call h5dwrite_f(d3set, H5T_NATIVE_DOUBLE, physicbuff(imin:imax, jmin:jmax, kmin:kmax), d3dims, error &
#ifdef PARALELO
             ,file_space_id = d3filespace, mem_space_id = d3memspace &
             ,xfer_prp = data_transfer &
#endif
             )
           call h5dclose_f(d3set, error)
#ifdef PARALELO
           call h5sclose_f(d3filespace, error)
#endif

!     gravz
           call h5dcreate_f(group_id, 'gravz', H5T_NATIVE_DOUBLE, d3space, d3set, error)

#ifdef PARALELO
           call h5dget_space_f(d3set, d3filespace, error)
           call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error) 
#endif
           physicbuff(imin:imax, jmin:jmax, kmin:kmax) =  physic(imin:imax, jmin:jmax, kmin:kmax)%gravz
           call h5dwrite_f(d3set, H5T_NATIVE_DOUBLE, physicbuff(imin:imax, jmin:jmax, kmin:kmax), d3dims, error &
#ifdef PARALELO
               ,file_space_id = d3filespace, mem_space_id = d3memspace &
               ,xfer_prp = data_transfer &
#endif
               )
           call h5dclose_f(d3set, error)
#ifdef PARALELO
           call h5sclose_f(d3filespace, error)
#endif

!        endif
      endif


      deallocate (physicbuff)

      ELSE  !if ird=0

   
! AXIAL PLANES
      
!   denstya
      call h5dcreate_f(group_id, 'denstya', H5T_NATIVE_DOUBLE, d2spacea, d2seta, error)       

#ifdef PARALELO
      call h5dget_space_f(d2seta, d2filespacea, error)
      call h5sselect_hyperslab_f(d2filespacea, H5S_SELECT_SET_F, d2offseta, d2dimsa, error)
      
#endif
      call h5dwrite_f(d2seta, H5T_NATIVE_DOUBLE, dbufa, d2dimsa, error &
#ifdef PARALELO
           ,file_space_id = d2filespacea, mem_space_id = d2memspacea &
           ,xfer_prp = data_transfer &
#endif
           )
      call h5dclose_f(d2seta, error)
#ifdef PARALELO
      call h5sclose_f(d2filespacea, error)
#endif

            

!     epsa
      call h5dcreate_f(group_id, 'epsa', H5T_NATIVE_DOUBLE, d2spacea, d2seta, error)

#ifdef PARALELO
      call h5dget_space_f(d2seta, d2filespacea, error)
      call h5sselect_hyperslab_f(d2filespacea, H5S_SELECT_SET_F, d2offseta, d2dimsa, error)
      
#endif
      call h5dwrite_f(d2seta, H5T_NATIVE_DOUBLE, dbufa, d2dimsa, error &
#ifdef PARALELO
           ,file_space_id = d2filespacea, mem_space_id = d2memspacea &
           ,xfer_prp = data_transfer &
#endif
          )
      call h5dclose_f(d2seta, error)
#ifdef PARALELO
      call h5sclose_f(d2filespacea, error)
#endif


      
!     velxa
      call h5dcreate_f(group_id, 'velxa', H5T_NATIVE_DOUBLE, d2spacea, d2seta, error)

#ifdef PARALELO
      call h5dget_space_f(d2seta, d2filespacea, error)
      call h5sselect_hyperslab_f(d2filespacea, H5S_SELECT_SET_F, d2offseta, d2dimsa, error)
      
#endif
      call h5dwrite_f(d2seta, H5T_NATIVE_DOUBLE, dbufa, d2dimsa, error &
#ifdef PARALELO
           ,file_space_id = d2filespacea, mem_space_id = d2memspacea &
           ,xfer_prp = data_transfer &
#endif
           )
      call h5dclose_f(d2seta, error)
#ifdef PARALELO
      call h5sclose_f(d2filespacea, error)
#endif

      if (nsdim.gt.1) then

        
!     velya
         call h5dcreate_f(group_id, 'velya', H5T_NATIVE_DOUBLE, d2spacea, d2seta, error)

#ifdef PARALELO
        call h5dget_space_f(d2seta, d2filespacea, error)
        call h5sselect_hyperslab_f(d2filespacea, H5S_SELECT_SET_F, d2offseta, d2dimsa, error)
      
#endif
        call h5dwrite_f(d2seta, H5T_NATIVE_DOUBLE, dbufa, d2dimsa, error &
#ifdef PARALELO
           ,file_space_id = d2filespacea, mem_space_id = d2memspacea &
           ,xfer_prp = data_transfer &
#endif
           )
        call h5dclose_f(d2seta, error)
#ifdef PARALELO
        call h5sclose_f(d2filespacea, error)
#endif

      endif

      if (nsdim.gt.2) then

       

!     velza
        call h5dcreate_f(group_id, 'velza', H5T_NATIVE_DOUBLE, d2spacea, d2seta, error)

#ifdef PARALELO
        call h5dget_space_f(d2seta, d2filespacea, error)
        call h5sselect_hyperslab_f(d2filespacea, H5S_SELECT_SET_F, d2offseta, d2dimsa, error) 
#endif
        call h5dwrite_f(d2seta, H5T_NATIVE_DOUBLE, dbufa, d2dimsa, error &
#ifdef PARALELO
           ,file_space_id = d2filespacea, mem_space_id = d2memspacea &
           ,xfer_prp = data_transfer &
#endif
           )
        call h5dclose_f(d2seta, error)
#ifdef PARALELO
        call h5sclose_f(d2filespacea, error)
#endif
      endif

    

!     tracera
      call h5dcreate_f(group_id, 'tracera', H5T_NATIVE_DOUBLE, d2spacea, d2seta, error)

#ifdef PARALELO
      call h5dget_space_f(d2seta, d2filespacea, error)
      call h5sselect_hyperslab_f(d2filespacea, H5S_SELECT_SET_F, d2offseta, d2dimsa, error)
      
#endif
      call h5dwrite_f(d2seta, H5T_NATIVE_DOUBLE, physic(imin:imax, jmin:jmax, kmin:kmax)%tracer, d2dimsa, error &
#ifdef PARALELO
           ,file_space_id = d2filespacea, mem_space_id = d2memspacea &
           ,xfer_prp = data_transfer &
#endif
           )
      call h5dclose_f(d2seta, error)
#ifdef PARALELO
      call h5sclose_f(d2filespacea, error)
#endif

      
      


      ENDIF


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
         call h5sclose_f(dzmemspace, error)
      endif
      call h5sclose_f(dxmemspace, error)
      if (ird == 0) then
         call h5sclose_f(d3memspace, error)
      endif
#endif /* PARALELO */

!     close total memory spaces
      if (nsdim.gt.1) then
         call h5sclose_f(dyspace, error)
         call h5sclose_f(dzspace, error)
      endif
      call h5sclose_f(dxspace, error)
      if (ird == 0) then
         call h5sclose_f(d3space, error)
      endif   

!     create memory spaces for header
      call h5sclose_f(intspace, error)
      call h5sclose_f(dblspace, error)
     
! WE WILL USE ONLY ONE GROUP FOR THE RESTART. PERUCHO
      if (ird /= 0) then
         curR_R = curR_R + 1
      endif

!     close file         
      call h5fclose_f(file_id, error)

!     prepare for writing a new file
!      if (curR_R.gt.RECORDER_R) then
        if (ird == 0) then          
          print*, &
              "Closing file "//TRIM(outRoot)//rstfil//"-" &
              //c1//c2//c3//c4//c5//".h5"
          r_filer = r_filer + 1
          H_openr = 0
        else 
           if (plane == 3) then
           print*, &
              "Closing file "//TRIM(outRoot)//datfil//"-" &
              //c1//c2//c3//c4//c5//".h5"
             r_filed = r_filed + 1
             H_opend = 0          
          endif
        endif
!      endif

        if (cartrank==0) then
            print*,'---------------------------------------'
        endif

#ifdef CASERO
      !!!call flush(6)
#else    
      !call flush_(6)
#endif
      
      RETURN
      END SUBROUTINE rrhdf5_record

      SUBROUTINE rrhdf5_finalize
            
      use hdf5
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

#ifdef PARALELO
      integer ierror
#endif

      integer :: error

      call h5close_f(error)
      

      RETURN
      END SUBROUTINE rrhdf5_finalize
#endif

