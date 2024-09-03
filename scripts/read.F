!===============================================================================
!     Read HDF5 file header
!module load hdf5/1.14.1_one
!h5pfc -o read read.F
!./read
!===============================================================================

      program read_hd
          use hdf5
          implicit none

          integer(HID_T) :: file_id,group_id
          integer(HID_T) :: dblspace,intspace
          integer(HSIZE_T),dimension(1) :: dbldims = (/1/),intdims=(/1/)
          integer(HID_T) :: intset,dblset
          integer :: error
          real(8) :: timej,dt
          real(8) :: time,c,Rb,pctocm
          integer(4) :: nsdim
          character(len=256) :: file_name,a,b

          !Constants and code units
          c = 2.9979e10
          pctocm= 3.086e18
          Rb = 500*pctocm

          !File path
          a = "/storage/scratch/lv72/lv72805"
          b = "/D02/J3C0/previous/J3C0RAA-00005.h5"

          file_name = trim(a)//trim(b)

          !Create dataspace
          call h5screate_simple_f(1,intdims,intspace,error)
          call h5screate_simple_f(1,dbldims,dblspace,error)

          !Open file
          call h5open_f(error)
          call h5fopen_f(trim(file_name),H5F_ACC_RDWR_F,file_id,error)
          !Open data group within the file
          call h5gopen_f(file_id,"00001",group_id,error)

          !Open and read variables
          call h5dopen_f(group_id,"timej",dblset,error)
          call h5dread_f(dblset,H5T_NATIVE_DOUBLE,timej,dbldims,error)
          call h5dclose_f(dblset,error)

          !Close file
          call h5close_f(error)

          !Print variables
          print*,'Time in code units',timej
          time = timej*(Rb/c)/(365*24*3600)/(1e6)
          print*,'Time [Myr] = ', time

      end program read_hd
