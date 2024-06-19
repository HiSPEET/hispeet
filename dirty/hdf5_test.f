program HDF5_Test
  USE ISO_C_Binding
  use HDF5
  implicit none

  type Element
    integer :: id = -1
    integer, allocatable :: nb(:)
  end type Element

  type(Element), target :: elem_orig(2)
  type(Element), target :: elem_read(2)

  integer(HID_T)   :: file_id, space_id, data_id, type_id, group_id
  integer(HSIZE_T) :: dims(1)
  integer(SIZE_T)  :: offset
  type(C_Ptr)      :: buf
  integer          :: err, i
  logical          :: alloc

  call H5open_f(err)

  ! initialize data, leaving dynamic component deallocated
  do i = 1, size(elem_orig)
    elem_orig(i)%id =  i
  end do
  allocate(elem_orig(2)%nb, source = [1,2])

  ! HDF5 compound to write static component
  offset = H5offsetof(C_Loc(elem_orig(1)), C_Loc(elem_orig(2)))
  call H5Tcreate_f(H5T_COMPOUND_F, offset, type_id, err)
  offset = H5offsetof(C_Loc(elem_orig(1)), C_Loc(elem_orig(1)%id))
  call H5Tinsert_f(type_id, 'id', offset, H5T_NATIVE_INTEGER, err)

  ! write data
  dims = size(elem_orig)
  call H5Fcreate_f('elements.h5', H5F_ACC_TRUNC_F, file_id, err)
  call H5Gcreate_f(file_id, 'data', group_id, err)
  call H5Screate_simple_f(1, dims, space_id, err)
  call H5Dcreate_f(group_id, 'elements', type_id, space_id, data_id, err)
  call H5Dwrite_f(data_id, type_id, C_Loc(elem_orig(1)), err)
  call H5Dclose_f(data_id, err)
  call H5Sclose_f(space_id, err)
  call H5Gclose_f(group_id, err)
  call H5Fclose_f(file_id, err)

  ! read data
  buf = C_Loc(elem_read(1))
  call H5Fopen_f('elements.h5', H5F_ACC_RDWR_F, file_id, err)
  call H5Gopen_f(file_id, 'data', group_id, err)
  call H5Dopen_f(group_id, 'elements', data_id, err)
  call H5Dget_type_f(data_id, type_id, err)
  call H5Dread_f(data_id, type_id, buf, err)
  call H5Dclose_f(data_id, err)
  call H5Gclose_f(group_id, err)
  call H5Fclose_f(file_id, err)

  ! check data .................................................................

  print '(9(G0,1X))', 'elem_orig%id =', elem_orig%id
  print '(9(G0,1X))', 'elem_read%id =', elem_read%id

  do i = 1, size(elem_read)
    alloc = allocated(elem_read(i)%nb)
    print '(9(G0,1X))', 'allocated(elem_read(',i,')%nb) =', alloc
    if (alloc) then
      print '(9(G0,1X))', 'size(elem_read(',i,')%nb) =', size(elem_read(i)%nb)
      print '(9(G0,1X))', 'elem_read(',i,')%nb =', elem_read(i)%nb
    end if
  end do

  call H5close_f(err)

end program HDF5_Test
