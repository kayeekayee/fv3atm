module FV3GFS_io_generic_mod
  use GFS_typedefs, only: kind_phys

  use fms2_io_mod,        only: FmsNetcdfDomainFile_t => FmsNetcdfDomainFile_t, unlimited => unlimited,      &
                                f_open_file => open_file, f_close_file => close_file,                 &
                                f_register_axis => register_axis, f_register_restart_field => register_restart_field, &
                                f_register_variable_attribute => register_variable_attribute, f_register_field => register_field, &
                                f_read_restart => read_restart, f_write_restart => write_restart, f_write_data => write_data,     &
                                f_get_global_io_domain_indices => get_global_io_domain_indices, f_variable_exists => variable_exists


  use FV3GFS_io_netcdf_mod,only: GFS_io_netCDF_type,      &
                                g_open_file => open_file, g_close_file => close_file,                 &
                                g_register_axis => register_axis, g_register_restart_field => register_restart_field, &
                                g_register_variable_attribute => register_variable_attribute, g_register_field => register_field, &
                                g_read_restart => read_restart, g_write_restart => write_restart, g_write_data => write_data,     &
                                g_get_global_io_domain_indices => get_global_io_domain_indices, g_variable_exists => variable_exists

  use mpp_domains_mod,    only: domain2d

  implicit none

  private

  public :: unlimited ! from fms2_io_mod

  public :: open_file, close_file, register_field, get_global_io_domain_indices
  public :: variable_exists, read_restart, write_restart, register_variable_attribute

  type, public :: GFS_io_generic_type
    type(FmsNetcdfDomainFile_t) :: fms2
    type(GFS_io_netCDF_type) :: ionet
    logical :: use_fms
  end type GFS_io_generic_type

  logical, parameter :: super_verbose = .false.
  logical :: module_initialized = .false.
  logical :: use_io_netcdf = .false.
  logical, parameter :: use_fms2_io = .true.

  public :: register_restart_field
  interface register_restart_field
     module procedure register_restart_field_int1
     module procedure register_restart_field_int2
     module procedure register_restart_field_int3
     module procedure register_restart_field_phys1
     module procedure register_restart_field_phys2
     module procedure register_restart_field_phys3
  end interface register_restart_field

  public :: register_axis
  interface register_axis
     module procedure register_axis_x_or_y
     module procedure register_axis_len
  end interface register_axis

  public :: write_data
  interface write_data
     module procedure write_data_int
     module procedure write_data_int_scalar
     module procedure write_data_phys
  end interface write_data

contains

  subroutine initialize_module
    use module_fv3_config,  only: fcst_mpi_comm
    use mpp_domains_mod, only: domain2d
    use atmosphere_mod, only: atmosphere_domain
    use mpi
    implicit none

    type(domain2d) :: fv_domain_for_read
    type(domain2d) :: fv_domain
    integer :: layout(2)
    logical :: regional
    logical :: nested
    integer :: ngrids_atmos
    integer :: mygrid_atmos
    integer, pointer :: pelist(:)
    logical :: moving_nest_parent
    logical :: is_moving_nest

    logical :: any_nested
    integer :: rank_in_fcst,ierr,iostat,unit

    namelist /FV3GFS_io/ use_io_netcdf

    if(module_initialized) then
      return
    endif

    module_initialized=.true.

    open(form="FORMATTED",file="input.nml",iostat=iostat,status='OLD',newunit=unit)
    if(iostat/=0) then
      if(super_verbose) then
        write(0,*) 'Could not open input.nml'
      endif
      return
    endif

    read(unit=unit,nml=FV3GFS_io,iostat=iostat)
    if(iostat/=0) then
      if(super_verbose) then
        write(0,*) 'could not read FV3GFS_io namelist',iostat
      endif
      close(unit)
      return
    endif

    if(super_verbose) then
8     format('Successfully read FV3GFS_io namelist with use_io_netcdf=',L1)
      print 8,use_io_netcdf
    endif
    close(unit)

    if(.not.use_io_netcdf) return

    call MPI_Comm_rank(fcst_mpi_comm,rank_in_fcst,ierr)

    call atmosphere_domain ( fv_domain=fv_domain, rd_domain=fv_domain_for_read, layout=layout, regional=regional, nested=nested, &
                                moving_nest_parent=moving_nest_parent, is_moving_nest=is_moving_nest, &
                                ngrids_atmos=ngrids_atmos, mygrid_atmos=mygrid_atmos, pelist=pelist )

38  format(I0,': is ',A)
    if(nested) then
      if(super_verbose) then
        print 38,rank_in_fcst,'nested'
      else
        print 38,rank_in_fcst,'not nested'
      endif
    endif

    call MPI_Allreduce(nested,any_nested,1,MPI_LOGICAL,MPI_LAND,fcst_mpi_comm,ierr)

    if(any_nested) then
      if(rank_in_fcst==0) then
        write(0,'(A)') 'FV3GFS_io_netcdf does not support nesting. Will use FMS2 IO for this simulation.'
      endif
      use_io_netcdf=.false.
    endif

  end subroutine initialize_module

  logical function open_file(restart,infile,mode,domain,is_restart,dont_add_res_to_filename)
    use module_fv3_config,  only: fcst_mpi_comm
    use mpi
    implicit none
    type(GFS_io_generic_type) :: restart
    character(*), intent(in) :: infile, mode
    type(domain2d), intent(in) :: domain
    logical, optional, intent(in) :: is_restart, dont_add_res_to_filename
    integer :: ierr, rank

    call MPI_Comm_rank(fcst_mpi_comm,rank,ierr)

    call initialize_module

    open_file=.false.
    restart%use_fms=.true.

    if(use_io_netcdf) then
      ! Try to use the faster parallel netcdf implementation first:
      open_file = g_open_file(restart%ionet,infile,mode,domain=domain,is_restart=is_restart, &
           dont_add_res_to_filename=dont_add_res_to_filename)
      if(open_file) then
        if(rank==0) then
83        format(A,': using NetCDF I/O in mode ',A)
          write(0,83) trim(infile),trim(mode)
        endif
        restart%use_fms=.false.
        return
      else
308   format(A)
      write(0,803) 'g_open_file failed'
      endif
    elseif(rank==0) then
803   format(A)
      write(0,803) 'IO NetCDF is disabled by namelist'
    endif

    if(use_fms2_io) then
      ! Fall back to the more general fms2 io:
      if(rank==0) then
38      format(A,': falling back to fms2 io in mode ',A)
        write(0,38) trim(infile),trim(mode)
      endif
      open_file = f_open_file(restart%fms2,infile,mode,domain=domain,is_restart=is_restart, &
           dont_add_res_to_filename=dont_add_res_to_filename)
    endif
  end function open_file

  subroutine close_file(restart)
    type(GFS_io_generic_type) :: restart
    if(restart%use_fms) then
      call f_close_file(restart%fms2)
    else
      call g_close_file(restart%ionet)
    end if
  end subroutine close_file

  subroutine write_restart(restart)
    type(GFS_io_generic_type) :: restart
    if(restart%use_fms) then
      call f_write_restart(restart%fms2)
    else
      call g_write_restart(restart%ionet)
    end if
  end subroutine write_restart

  subroutine read_restart(restart,ignore_checksum )
    type(GFS_io_generic_type) :: restart
    logical, intent(in), optional :: ignore_checksum
    if(restart%use_fms) then
      if(present(ignore_checksum)) then
        call f_read_restart(restart%fms2,ignore_checksum=ignore_checksum)
      else
        call f_read_restart(restart%fms2)
      endif
    else
      call g_read_restart(restart%ionet)
    end if
  end subroutine read_restart

  subroutine get_global_io_domain_indices(restart, name, global_start, global_end, indices)
    implicit none
    type(GFS_io_generic_type), intent(inout) :: restart
    character(len=*), intent(in) :: name
    integer, intent(out) :: global_start, global_end
    integer, allocatable :: indices(:)

    if(restart%use_fms) then
      call f_get_global_io_domain_indices(restart%fms2, name, global_start, global_end, indices=indices)
    else
      call g_get_global_io_domain_indices(restart%ionet, name, global_start, global_end, indices=indices)
    end if
  end subroutine get_global_io_domain_indices

  logical function variable_exists(restart,name)
    type(GFS_io_generic_type) :: restart
    character(len=*), intent(in) :: name
    if(restart%use_fms) then
      variable_exists=f_variable_exists(restart%fms2,name)
    else
      variable_exists=g_variable_exists(restart%ionet,name)
    end if
  end function variable_exists

  subroutine write_data_int(restart,name,buffer)
    use mpi
    use netcdf
    implicit none
    type(GFS_io_generic_type), intent(inout) :: restart
    character(*), intent(in) :: name
    integer, intent(in) :: buffer(:)
    if(restart%use_fms) then
      call f_write_data(restart%fms2,name,buffer)
    else
      call g_write_data(restart%ionet,name,buffer)
    end if
  end subroutine write_data_int

  subroutine write_data_int_scalar(restart,name,buffer)
    use mpi
    use netcdf
    implicit none
    type(GFS_io_generic_type), intent(inout) :: restart
    character(*), intent(in) :: name
    integer, intent(in) :: buffer
    if(restart%use_fms) then
      call f_write_data(restart%fms2,name,buffer)
    else
      call g_write_data(restart%ionet,name,buffer)
    end if
  end subroutine write_data_int_scalar

  subroutine write_data_phys(restart,name,buffer)
    use mpi
    use netcdf
    implicit none
    type(GFS_io_generic_type), intent(inout) :: restart
    character(*), intent(in) :: name
    real(kind_phys), intent(in) :: buffer(:)
    if(restart%use_fms) then
      call f_write_data(restart%fms2,name,buffer)
    else
      call g_write_data(restart%ionet,name,buffer)
    end if
  end subroutine write_data_phys

  subroutine register_axis_x_or_y(restart,name,axis_name)
    implicit none
    type(GFS_io_generic_type), intent(inout) :: restart
    character(len=*), intent(in) :: name, axis_name

    if(restart%use_fms) then
      call f_register_axis(restart%fms2,name,axis_name)
    else
      call g_register_axis(restart%ionet,name,axis_name)
    endif
  end subroutine register_axis_x_or_y

  subroutine register_axis_len(restart,name,dimension_length)
    implicit none
    type(GFS_io_generic_type), intent(inout) :: restart
    character(len=*), intent(in) :: name
    integer, intent(in) :: dimension_length
    if(restart%use_fms) then
      call f_register_axis(restart%fms2,name,dimension_length=dimension_length)
    else
      call g_register_axis(restart%ionet,name,dimension_length=dimension_length)
    endif
  end subroutine register_axis_len

  subroutine register_variable_attribute(restart,var_name,att_name,att_value,str_len)
    implicit none
    type(GFS_io_generic_type), intent(inout) :: restart
    character(*), intent(in) :: var_name,att_name,att_value
    integer, intent(in) :: str_len
    if(restart%use_fms) then
      call f_register_variable_attribute(restart%fms2,var_name,att_name,att_value,str_len=str_len)
    else
      call g_register_variable_attribute(restart%ionet,var_name,att_name,att_value,str_len=str_len)
    endif
  end subroutine register_variable_attribute

  subroutine register_restart_field_int1(restart, name, data, dimensions, is_optional)
    implicit none
    type(GFS_io_generic_type), intent(inout) :: restart
    integer, pointer :: data(:)
    character(*), intent(in) :: dimensions(:)
    character(*), intent(in) :: name
    logical, optional :: is_optional
    if(restart%use_fms) then
      call f_register_restart_field(restart%fms2,name,data,dimensions=dimensions,is_optional=is_optional)
    else
      call g_register_restart_field(restart%ionet,name,data,dimensions=dimensions,is_optional=is_optional)
    endif
  end subroutine register_restart_field_int1

  subroutine register_restart_field_int2(restart, name, data, dimensions, is_optional)
    implicit none
    type(GFS_io_generic_type), intent(inout) :: restart
    integer, pointer :: data(:,:)
    character(*), intent(in) :: dimensions(:)
    character(*), intent(in) :: name
    logical, optional :: is_optional
    if(restart%use_fms) then
      call f_register_restart_field(restart%fms2,name,data,dimensions=dimensions,is_optional=is_optional)
    else
      call g_register_restart_field(restart%ionet,name,data,dimensions=dimensions,is_optional=is_optional)
    endif
  end subroutine register_restart_field_int2

  subroutine register_restart_field_int3(restart, name, data, dimensions, is_optional)
    implicit none
    type(GFS_io_generic_type), intent(inout) :: restart
    integer, pointer :: data(:,:,:)
    character(*), intent(in) :: dimensions(:)
    character(*), intent(in) :: name
    logical, optional :: is_optional
    if(restart%use_fms) then
      call f_register_restart_field(restart%fms2,name,data,dimensions=dimensions,is_optional=is_optional)
    else
      call g_register_restart_field(restart%ionet,name,data,dimensions=dimensions,is_optional=is_optional)
    endif
  end subroutine register_restart_field_int3

  subroutine register_restart_field_phys1(restart, name, data, dimensions, is_optional)
    implicit none
    type(GFS_io_generic_type), intent(inout) :: restart
    real(kind_phys), pointer :: data(:)
    character(*), intent(in) :: dimensions(:)
    character(*), intent(in) :: name
    logical, optional :: is_optional
    if(restart%use_fms) then
      call f_register_restart_field(restart%fms2,name,data,dimensions=dimensions,is_optional=is_optional)
    else
      call g_register_restart_field(restart%ionet,name,data,dimensions=dimensions,is_optional=is_optional)
    endif
  end subroutine register_restart_field_phys1

  subroutine register_restart_field_phys2(restart, name, data, dimensions, is_optional)
    implicit none
    type(GFS_io_generic_type), intent(inout) :: restart
    real(kind_phys), pointer :: data(:,:)
    character(*), intent(in) :: dimensions(:)
    character(*), intent(in) :: name
    logical, optional :: is_optional
    if(restart%use_fms) then
      call f_register_restart_field(restart%fms2,name,data,dimensions=dimensions,is_optional=is_optional)
    else
      call g_register_restart_field(restart%ionet,name,data,dimensions=dimensions,is_optional=is_optional)
    endif
  end subroutine register_restart_field_phys2

  subroutine register_restart_field_phys3(restart, name, data, dimensions, is_optional)
    implicit none
    type(GFS_io_generic_type), intent(inout) :: restart
    real(kind_phys), pointer :: data(:,:,:)
    character(*), intent(in) :: dimensions(:)
    character(*), intent(in) :: name
    logical, optional :: is_optional
    if(restart%use_fms) then
      call f_register_restart_field(restart%fms2,name,data,dimensions=dimensions,is_optional=is_optional)
    else
      call g_register_restart_field(restart%ionet,name,data,dimensions=dimensions,is_optional=is_optional)
    endif
  end subroutine register_restart_field_phys3

  subroutine register_field(restart, name, type, dims)
    implicit none
    type(GFS_io_generic_type), intent(inout) :: restart
    character(*), intent(in) :: dims(:)
    character(len=*), intent(in) :: name, type
    if(restart%use_fms) then
      call f_register_field(restart%fms2,name,type,dims)
    else
      call g_register_field(restart%ionet,name,type,dims)
    end if
  end subroutine register_field

end module FV3GFS_io_generic_mod
