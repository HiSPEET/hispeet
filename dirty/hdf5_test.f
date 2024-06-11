program HDF5_Test
  use Kind_Parameters
  USE ISO_C_BINDING   ! required for C_Loc
  use HDF5_Binding
  use Mesh_Element__3D
  implicit none

  

  ! HDF5 MeshElement type .................................................................

  integer(hid_t) :: H5T_MeshElement_3D     = -1

  ! auxiliary variables ........................................................

  integer          :: err, i, j, e
  integer(hid_t)   :: file_id, dataspace_id, dataset_id, dtype_id, group_id
  integer(hsize_t) :: dims(1)
  type(c_ptr)      :: f_ptr

  type(MeshElement_3D), target :: element_orig(2)
  type(MeshElement_3D), target :: element_read(2)

  ! HDF5 initialization ........................................................

  call H5open_f(err)

  call Init_HDF5_Binding()
  print '(A,I0)', 'H5T_INTEGER         =  ', H5T_INTEGER
  print '(A,I0)', 'H5T_INTEGER_IXS     =  ', H5T_INTEGER_IXS
  print '(A,I0)', 'H5T_INTEGER_IXL     =  ', H5T_INTEGER_IXL
  print '(A,I0)', 'H5T_REAL_RSP        =  ', H5T_REAL_RSP
  print '(A,I0)', 'H5T_REAL_RDP        =  ', H5T_REAL_RDP
  print '(A,I0)', 'H5T_REAL_RHP        =  ', H5T_REAL_RHP
  print '(A,I0)', 'H5T_REAL_RNP        =  ', H5T_REAL_RNP
  print '(A,I0)', 'H5T_CHARACTER       =  ', H5T_CHARACTER
  print '(A,I0)', 'H5T_FORTRAN_S1      =  ', H5T_FORTRAN_S1
  print '(A,I0)', 'H5T_LOGICAL         =  ', H5T_LOGICAL

  call Get_H5T_MeshElement_3D(H5T_MeshElement_3D)
  print '(A,I0)', 'H5T_MeshElement     =  ', H5T_MeshElement_3D

  ! create dummy data ..........................................................
  
   do e = 1,2
     element_orig(e)%id = e
     element_orig(e)%cluster_id = 0
     element_orig(e)%frozen = .true.
   
     do i = 1,8
       element_orig(e)%vertex(i)%id         = i * element_orig(e)%id
       element_orig(e)%vertex(i)%n_neighbor = 1
     end do

     do j = 1,3
     do i = 0,3
       element_orig(e)%geometry%x_c(i,j) = i*j
     end do
     end do

     allocate(element_orig(e)%neighbor(1))
     element_orig(e)%neighbor(1)%id        = 2
     element_orig(e)%neighbor(1)%component = 2
    
   end do
       


  ! save data ..................................................................
  
  ! Create a new file (or open an existing one)
  call h5fcreate_f('elements.h5', H5F_ACC_TRUNC_F, file_id, err)
  ! Create a group
  call h5gcreate_f(file_id, "mesh_data", group_id, err)
  ! Create dataspace for the dataset
  dims = size(element_orig)
  !call h5screate_simple_f(rank, dims, dataspace_id, err)
  call h5screate_simple_f(1, dims, dataspace_id, err)
  ! Create the dataset
  call h5dcreate_f(group_id, 'mesh_partition', H5T_MeshElement_3D, &
                   dataspace_id, dataset_id, err)
  ! Write the dataset
  call h5dwrite_f(dataset_id, H5T_MeshElement_3D, C_LOC(element_orig(1)), err)
  ! Close the dataset
  call h5dclose_f(dataset_id, err)
  ! Close the dataspace
  call h5sclose_f(dataspace_id, err)
  ! Close the group
  call h5gclose_f(group_id, err)
  ! Close the file
  call h5fclose_f(file_id, err)


  ! read data ..................................................................
  
  ! Open HDF5 file
  call h5fopen_f('elements.h5', H5F_ACC_RDWR_F, file_id, err)
  ! Open the group
  call h5gopen_f(file_id, 'mesh_data', group_id, err)
  ! Open the dataspace
  call h5dopen_f(group_id, 'mesh_partition', dataset_id, err)
  ! Read the dataset
  call h5dget_type_f(dataset_id, dtype_id, err)
  f_ptr = C_LOC(element_read(1))
  call h5dread_f(dataset_id, dtype_id, f_ptr, err)
  ! Close the dataset
  call h5dclose_f(dataset_id, err)
  ! Close the group
  call h5gclose_f(group_id, err)
  ! Close the file
  call h5fclose_f(file_id, err)

  ! check data .................................................................

  print '(A,99(F5.1,1X))', 'orig: e1%x_c =', element_orig(1)%geometry%x_c
  print '(A,99(F5.1,1X))', 'read: e1%x_c =', element_read(1)%geometry%x_c

  call H5close_f(err)

  !=============================================================================

end program HDF5_Test
