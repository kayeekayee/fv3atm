module FV3GFS_io_netCDF_mod
  use GFS_typedefs, only: kind_phys
  use netcdf, only: unlimited => nf90_unlimited
  implicit none

  private

  public :: unlimited
  public :: open_file, close_file, register_field, get_global_io_domain_indices
  public :: variable_exists, read_restart, write_restart, register_variable_attribute

  integer, parameter, public :: GFS_io_max_dims = 10
  integer, parameter, public :: GFS_io_max_var_dims = 4
  integer, parameter, public :: GFS_io_dimname_len = 10
  integer, parameter, public :: GFS_io_max_vars = 200
  integer, parameter, public :: GFS_io_varname_len = 40
  integer, parameter, public :: GFS_io_invalid_id = -1
  integer, parameter, public :: GFS_io_max_var_attr = 5

  logical, parameter :: super_verbose = .true.

  type, public :: GFS_io_att_type
     character(:), pointer :: name, value
  end type GFS_io_att_type

  type, public :: GFS_io_dim_type
     character(len=GFS_io_dimname_len) :: name ! dimension name
     integer :: dimid ! NetCDF dimension id
     integer :: dimlen ! Dimension length
     integer :: local_start ! start of data within this rank
     integer :: local_end ! end of data within this rank
  end type GFS_io_dim_type

  type, public :: GFS_io_var_type
     character(len=GFS_io_varname_len) :: name ! variable name
     integer :: varid ! netcdf variable id for this variable
     integer :: ndims ! number of dimensions of this variable
     integer :: xtype ! NetCDF type
     integer :: dimids(GFS_io_max_var_dims) ! netcdf dimension ids
     integer :: dimindex(GFS_io_max_var_dims) ! index of dims in GFS_io_netCDF_type%dims

     integer :: nattr = 0
     type(GFS_io_att_type) :: attr(GFS_io_max_var_attr)

     ! Data storage for read/write of this variable
     real(kind_phys), pointer :: phys1(:),   phys2(:,:), phys3(:,:,:)
     real,            pointer :: real1(:),   real2(:,:), real3(:,:,:)
     integer,         pointer ::  int1(:),    int2(:,:),  int3(:,:,:)

     logical :: is_optional, write_data_called
  end type GFS_io_var_type

  type, public :: GFS_io_netCDF_type
     character(:), pointer :: filename

     ! MPI communication information
     integer :: comm, rank

     ! mpp_domains_mod tile count
     integer :: ntiles

     ! atmosphere_mod decomposition information using atmos_mod naming
     integer :: nx, ny ! local dimensions
     integer :: cnx, cny ! tile dimensions
     integer :: tile_num, isc, iec, jsc, jec, kt

     ! nf90_create or nf90_open information
     integer :: ncid
     logical :: write

     ! For write mode, have we called nf90_enddef yet?
     logical :: is_defined

     ! Dimension storage
     integer :: ndims
     type(GFS_io_dim_type) :: dims(GFS_io_max_dims)

     ! Variable storage
     integer :: nvars
     type(GFS_io_var_type) :: vars(GFS_io_max_vars)
  end type GFS_io_netCDF_type

  public :: register_restart_field
  interface register_restart_field
     module procedure register_restart_field_int1
     module procedure register_restart_field_int2
     module procedure register_restart_field_int3
     module procedure register_restart_field_phys1
     module procedure register_restart_field_phys2
     module procedure register_restart_field_phys3
     module procedure register_restart_field_real1
     module procedure register_restart_field_real2
     module procedure register_restart_field_real3
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
     module procedure write_data_real
     module procedure write_data_phys
  end interface write_data

contains ! ------------------------------------------------------------

  logical function open_file(restart,infile,mode,domain,is_restart,dont_add_res_to_filename)
    use mpp_domains_mod,    only: domain2d, mpp_get_ntile_count
    use module_fv3_config,  only: fcst_mpi_comm
    use atmosphere_mod,     only: atmosphere_resolution, atmosphere_control_data
    use mpi
    use netcdf

    implicit none

    type(GFS_io_netCDF_type) :: restart
    character(*), intent(in) :: infile, mode
    type(domain2d), intent(in) :: domain
    logical, optional, intent(in) :: is_restart, dont_add_res_to_filename

    character(len=:), pointer :: fullname
    integer :: need, info, ierr, ncerr, ntiles, idot ! FIXME: GET COMM

    open_file=.true.

    ! Get basic communication information
    restart%comm = fcst_mpi_comm
    call MPI_Comm_rank(restart%comm,restart%rank,ierr)

    ! Sanity checks
    if(mode/='read' .and. mode/='overwrite') then
       if(restart%rank==0) then
98        format('Invalid mode "',A,'"')
          write(0,98) trim(mode)
       endif
       open_file=.false.
    endif
    if(present(is_restart)) then
       if(.not.is_restart) then
          if(restart%rank==0) then
             write(0,'(A)') "FV3GFS_io_netCDF_mod can only handle read and overwrite modes"
          endif
          open_file=.false.
       endif
    endif
    if(present(dont_add_res_to_filename)) then
       if(.not.dont_add_res_to_filename) then
          if(restart%rank==0) then
             write(0,'(A)') "FV3GFS_io_netCDF_mod dont_add_res_to_filename must be .true."
          endif
          open_file=.false.
       endif
    endif

    ! Give up if any process on this tile failed.
    if(.not.allreduce_and(open_file,restart%comm)) return

    ! Get tile count from mpp_domains_mod
    restart%ntiles = mpp_get_ntile_count(domain)

    ! Get decomposition information from atmos_mod
    call atmosphere_resolution(restart%nx, restart%ny, global=.false.)
    call atmosphere_resolution(restart%cnx, restart%cny, global=.true.)
    call atmosphere_control_data(restart%isc,restart%iec,restart%jsc,restart%jec, &
         restart%kt, tile_num=restart%tile_num)

    ! Construct the filename, adding .tileN if needed:
    idot=index(infile,'.',.true.)
    need=len_trim(infile)+len('.tile999.nc') ! longest string we may need
    allocate(character(len=need) :: fullname)
    if(ntiles>1) then
       if(idot==0) then
          ! Need to add .nc extension
91        format(A,".tile",I0,".nc")
          write(fullname,91) trim(infile),restart%tile_num
       else
          ! Has extension
92        format(A,".tile",I0,A)
       write(fullname,92) infile(1:idot-1),restart%tile_num,infile(idot:len_trim(infile))
    else
       fullname=infile
    endif

    if(super_verbose) then
104    format('Open file "',A,'"')
       if(restart%rank==0) then
          write(0,104) trim(fullname) 
    endif

    ! Open the file. We'll temporarily need an MPI_Info for this.
    call MPI_Info_create(info,ierr)
16     format(A,": ",A," tile ",I0)
    if(mode=='overwrite') then
       if(restart%rank==0) then
          write(0,) trim(fullname),'write',restart%tile_num
       endif
       ncerr=nf90_create_par(path=trim(fullname), &
            cmode=ior(NF90_CLOBBER,NF90_NETCDF4), &
            ncid=restart%ncid,comm=restart%comm,info=info)
       restart%write=.true.
    else ! mode=='read'
       if(restart%rank==0) then
          write(0,) trim(fullname),'read',restart%tile_num
       endif
       ncerr=nf90_open_par(path=trim(fullname),cmode=NF90_NOWRITE, &
            ncid=restart%ncid,comm=restart%comm,info=info)
       restart%write=.false.
    endif
    call mpi_info_free(info,ierr)

    ! Make sure the file was created or opened.
    open_file = ncerr==NF90_NOERR
    if(.not.allreduce_and(open_file,restart%comm)) then
       if(ncerr==NF90_NOERR) then
          ncerr=nf90_close(restart%ncid)
       endif
       if(restart%rank==0) then
          write(0,'(A,A,I0,A,A)') trim(fullname),": could not open (tile ",restart%tile_num,"): ",nf90_strerror(ncerr)
       endif
       deallocate(fullname)
       nullify(restart%filename)
    else
       restart%filename => fullname
    endif

    restart%is_defined = .false.

    ! Haven't read or written anything yet.
    restart%ndims=0
    restart%nvars=0

    ! Note: restart%dims and restart%vars don't need to be initialized
    ! here because they are initialized in next_dim_idx, next_var_idx,
    ! and next_attr_idx


  end function open_file

  ! --------------------------------------------------------------------

  subroutine close_file(restart)
    use netcdf
    implicit none
    type(GFS_io_netCDF_type), intent(inout) :: restart
    integer :: v, a
    call handle_ncerr(restart,'close file',nf90_close(restart%ncid))

    deallocate(restart%filename)
    nullify(restart%filename)
    restart%ncid=GFS_io_invalid_id

    if(super_verbose) then
       write(0,*) 'deallocate data (begin)'
    endif

    do v=1,restart%nvars
       if(restart%vars(v)%write_data_called) then
          if(associated(restart%vars(v)%phys1)) then
             deallocate(restart%vars(v)%phys1)
             nullify(restart%vars(v)%phys1)
          endif
          if(associated(restart%vars(v)%phys2)) then
             deallocate(restart%vars(v)%phys2)
             nullify(restart%vars(v)%phys2)
          endif
          if(associated(restart%vars(v)%phys3)) then
             deallocate(restart%vars(v)%phys3)
             nullify(restart%vars(v)%phys3)
          endif
          if(associated(restart%vars(v)%real1)) then
             deallocate(restart%vars(v)%real1)
             nullify(restart%vars(v)%real1)
          endif
          if(associated(restart%vars(v)%real2)) then
             deallocate(restart%vars(v)%real2)
             nullify(restart%vars(v)%real2)
          endif
          if(associated(restart%vars(v)%real3)) then
             deallocate(restart%vars(v)%real3)
             nullify(restart%vars(v)%real3)
          endif
          if(associated(restart%vars(v)%int1)) then
             deallocate(restart%vars(v)%int1)
             nullify(restart%vars(v)%int1)
          endif
          if(associated(restart%vars(v)%int2)) then
             deallocate(restart%vars(v)%int2)
             nullify(restart%vars(v)%int2)
          endif
          if(associated(restart%vars(v)%int3)) then
             deallocate(restart%vars(v)%int3)
             nullify(restart%vars(v)%int3)
          endif
       endif
       do a=1,restart%vars(v)%nattr
          if(associated(restart%vars(v)%attr(a)%name)) then
             deallocate(restart%vars(v)%attr(a)%name)
             nullify(restart%vars(v)%attr(a)%name)
          endif
          if(associated(restart%vars(v)%attr(a)%value)) then
             deallocate(restart%vars(v)%attr(a)%value)
             nullify(restart%vars(v)%attr(a)%value)
          endif
       enddo
    enddo

    if(super_verbose) then
       write(0,*) 'deallocate data (success)'
    endif
  end subroutine close_file

  ! --------------------------------------------------------------------

  subroutine define_restart(restart)
    use netcdf
    use mpi
    implicit none
    type(GFS_io_netCDF_type), intent(inout) :: restart
    integer ivar, idim

    if(restart%is_defined) then
       ! Ensure we don't redefine in a second write_restart() for the
       ! same restart file.
       return
    endif

    ! This is where we would sort dims, vars, and attributes, but fms2
    ! io does not do that, so FV3GFS_io_netcdf won't either.

    ! Define all dimensions:
    do idim=1,restart%ndims
       call define_dim(restart,restart%dims(idim))
    enddo

    ! Define variables and attributes:
    do ivar=1,restart%nvars
       call define_var(restart,restart%vars(ivar))
    enddo

    ! Exit definition mode. NetCDF writes faster if you do this.
    call handle_ncerr(restart,"end definition",nf90_enddef(restart%ncid))

    restart%is_defined=.true.
  end subroutine define_restart

  ! --------------------------------------------------------------------

  subroutine define_dim(restart,dim)
    use netcdf
    use mpi
    implicit none
    type(GFS_io_netCDF_type), intent(inout) :: restart
    type(GFS_io_dim_type), intent(inout) :: dim

    call handle_ncerr(restart,"define dim "//trim(dim%name),&
         nf90_def_dim(restart%ncid,trim(dim%name),len=dim%dimlen,dimid=dim%dimid))
  end subroutine define_dim

  ! --------------------------------------------------------------------

  subroutine define_var(restart,var)
    use netcdf
    implicit none
    type(GFS_io_netCDF_type), intent(inout) :: restart
    type(GFS_io_var_type), intent(inout) :: var
    integer :: i,idim
    do i=1,var%ndims
       idim = var%dimindex(i)
       var%dimids(i) = restart%dims(idim)%dimid
    enddo
    call handle_ncerr(restart,"def var "//trim(var%name), &
         nf90_def_var(ncid=restart%ncid,name=trim(var%name),xtype=var%xtype, &
         dimids=var%dimids(1:var%ndims), & ! input
         varid=var%varid))   ! output

    ! if(var%write_data_called) then
    !    ! Write whole array from first rank with write_data
    !    call handle_ncerr(restart,"turn off parallel access for "//trim(var%name), &
    !         nf90_var_par_access(restart%ncid,var%varid,NF90_INDEPENDENT))
    ! else
       ! Distributed write for register_*_field
       call handle_ncerr(restart,"turn on parallel access for "//trim(var%name), &
            nf90_var_par_access(restart%ncid,var%varid,NF90_COLLECTIVE))
    ! endif
    do i=1,var%nattr
       call define_attr(restart,var,var%attr(i))
    end do
  end subroutine define_var

  ! --------------------------------------------------------------------

  subroutine define_attr(restart,var,att)
    use netcdf
    implicit none
    type(GFS_io_netCDF_type), intent(inout) :: restart
    type(GFS_io_var_type), intent(inout) :: var
    type(GFS_io_att_type), intent(inout) :: att
    call handle_ncerr(restart,'def att '//att%name//' in var '//att%value, &
         nf90_put_att(restart%ncid,var%varid,att%name,att%value))
  end subroutine define_attr

  ! --------------------------------------------------------------------

  subroutine write_restart(restart)
    use netcdf
    use mpi
    implicit none
    type(GFS_io_netCDF_type), intent(inout) :: restart
    integer ivar

    if(.not. restart%write) then
       return
    endif

    call define_restart(restart)

    do ivar=1,restart%nvars
       call write_var(restart,restart%vars(ivar))
    enddo
  end subroutine write_restart

  ! --------------------------------------------------------------------

  subroutine write_var(restart,var)
    use netcdf
    implicit none
    type(GFS_io_netCDF_type), intent(inout) :: restart
    type(GFS_io_var_type), intent(inout) :: var
    integer :: start(var%ndims), count(var%ndims)
    integer :: idim, expected_size, ierr

    ! write(0,'(A,A)') 'enter barrier before ',trim(var%name)
    ! call MPI_Barrier(restart%comm,ierr)
    ! write(0,'(A,A)') 'exit  barrier before ',trim(var%name)

    if(var%write_data_called) then
       expected_size=1
       do idim=1,var%ndims
          start(idim) = 1
          count(idim) = restart%dims(var%dimindex(idim))%dimlen
          expected_size = expected_size*count(idim)
       enddo
    else
       expected_size=1
       do idim=1,var%ndims
          start(idim) = restart%dims(var%dimindex(idim))%local_start
          count(idim) = restart%dims(var%dimindex(idim))%local_end - start(idim) + 1
          expected_size = expected_size*count(idim)
       enddo
    endif

    if(associated(var%phys1)) then
       call check_size('phys1',size(var%phys1))
       call handle_ncerr(restart,'write '//trim(var%name), &
           nf90_put_var(restart%ncid,var%varid,var%phys1,start,count))
    else if(associated(var%phys2)) then
       call check_size('phys2',size(var%phys2))
       call handle_ncerr(restart,'write '//trim(var%name), &
           nf90_put_var(restart%ncid,var%varid,var%phys2,start,count))
    else if(associated(var%phys3)) then
       call check_size('phys3',size(var%phys3))
       call handle_ncerr(restart,'write '//trim(var%name), &
           nf90_put_var(restart%ncid,var%varid,var%phys3,start,count))
    else if(associated(var%real1)) then
       call check_size('real1',size(var%real1))
       call handle_ncerr(restart,'write '//trim(var%name), &
           nf90_put_var(restart%ncid,var%varid,var%real1,start,count))
    else if(associated(var%real2)) then
       call check_size('real2',size(var%real2))
       call handle_ncerr(restart,'write '//trim(var%name), &
           nf90_put_var(restart%ncid,var%varid,var%real2,start,count))
    else if(associated(var%real3)) then
       call check_size('real3',size(var%real3))
       call handle_ncerr(restart,'write '//trim(var%name), &
           nf90_put_var(restart%ncid,var%varid,var%real3,start,count))
    else if(associated(var%int1)) then
       call check_size('int1',size(var%int1))
       call handle_ncerr(restart,'write '//trim(var%name), &
           nf90_put_var(restart%ncid,var%varid,var%int1,start,count))
    else if(associated(var%int2)) then
       call check_size('int2',size(var%int2))
       call handle_ncerr(restart,'write '//trim(var%name), &
           nf90_put_var(restart%ncid,var%varid,var%int2,start,count))
    else if(associated(var%int3)) then
       call check_size('int3',size(var%int3))
       call handle_ncerr(restart,'write '//trim(var%name), &
           nf90_put_var(restart%ncid,var%varid,var%int3,start,count))
    else ! if(restart%rank==0) then
101    format(A,': No array provided for writing variable "',A,'"')
       write(0,101) restart%filename,trim(var%name)
    endif

    ! write(0,'(A,A)') 'enter barrier after  ',trim(var%name)
    ! call MPI_Barrier(restart%comm,ierr)
    ! write(0,'(A,A)') 'exit  barrier after  ',trim(var%name)

  contains

    subroutine check_size(what,actual_size)
      use mpi
      implicit none
      integer, intent(in) :: actual_size
      character(*), intent(in) :: what
      integer :: ierr
      if(actual_size/=expected_size) then
         if(restart%rank==0) then
            write(0,38) restart%filename,trim(var%name),trim(what),actual_size,expected_size
38          format(A,': var "',A,'" field ',A,' has size ',I0,' instead of expected size ',I0)
         endif
         call MPI_Abort(MPI_COMM_WORLD,1,ierr)
      endif
    end subroutine check_size
  end subroutine write_var

  ! --------------------------------------------------------------------

  subroutine read_restart(restart)
    implicit none
    type(GFS_io_netCDF_type), intent(inout) :: restart
    integer ivar

    if(restart%write) then
       return
    endif

    do ivar=1,restart%nvars
       call read_var(restart,restart%vars(ivar))
    enddo
  end subroutine read_restart

  ! --------------------------------------------------------------------

  subroutine read_var(restart,var)
    use netcdf
    implicit none
    type(GFS_io_netCDF_type), intent(inout) :: restart
    type(GFS_io_var_type), intent(inout) :: var
    integer :: start(var%ndims), count(var%ndims)
    integer :: idim, expected_size,ierr

    ! write(0,'(A,A)') 'enter barrier for ',trim(var%name)
    ! call MPI_Barrier(restart%comm,ierr)
    ! write(0,'(A,A)') 'exit  barrier for ',trim(var%name)

    if(var%write_data_called) then
       expected_size=1
       do idim=1,var%ndims
          start(idim) = 1
          count(idim) = restart%dims(var%dimindex(idim))%dimlen
          expected_size = expected_size*count(idim)
       enddo
    else
       expected_size=1
       do idim=1,var%ndims
          start(idim) = restart%dims(var%dimindex(idim))%local_start
          count(idim) = restart%dims(var%dimindex(idim))%local_end - start(idim) + 1
          expected_size = expected_size*count(idim)
       enddo
    endif

    if(associated(var%phys1)) then
       call check_size('phys1',size(var%phys1))
       call handle_ncerr(restart,'read '//trim(var%name), &
           nf90_get_var(restart%ncid,var%varid,var%phys1,start,count))
    else if(associated(var%phys2)) then
       call check_size('phys2',size(var%phys2))
       call handle_ncerr(restart,'read '//trim(var%name), &
           nf90_get_var(restart%ncid,var%varid,var%phys2,start,count))
    else if(associated(var%phys3)) then
       call check_size('phys3',size(var%phys3))
       call handle_ncerr(restart,'read '//trim(var%name), &
           nf90_get_var(restart%ncid,var%varid,var%phys3,start,count))
    else if(associated(var%real1)) then
       call check_size('real1',size(var%real1))
       call handle_ncerr(restart,'read '//trim(var%name), &
           nf90_get_var(restart%ncid,var%varid,var%real1,start,count))
    else if(associated(var%real2)) then
       call check_size('real2',size(var%real2))
       call handle_ncerr(restart,'read '//trim(var%name), &
           nf90_get_var(restart%ncid,var%varid,var%real2,start,count))
    else if(associated(var%real3)) then
       call check_size('real3',size(var%real3))
       call handle_ncerr(restart,'read '//trim(var%name), &
           nf90_get_var(restart%ncid,var%varid,var%real3,start,count))
    else if(associated(var%int1)) then
       call check_size('int1',size(var%int1))
       call handle_ncerr(restart,'read '//trim(var%name), &
           nf90_get_var(restart%ncid,var%varid,var%int1,start,count))
    else if(associated(var%int2)) then
       call check_size('int2',size(var%int2))
       call handle_ncerr(restart,'read '//trim(var%name), &
           nf90_get_var(restart%ncid,var%varid,var%int2,start,count))
    else if(associated(var%int3)) then
       call check_size('int3',size(var%int3))
       call handle_ncerr(restart,'read '//trim(var%name), &
           nf90_get_var(restart%ncid,var%varid,var%int3,start,count))
    else if(super_verbose) then
       ! This is expected for non-restart variables.
101    format(A,': No destination for reading variable "',A,'"')
       write(0,101) restart%filename,trim(var%name)
    endif

  contains

    subroutine check_size(what,actual_size)
      use mpi
      implicit none
      integer, intent(in) :: actual_size
      character(*), intent(in) :: what
      integer :: ierr
      if(actual_size/=expected_size) then
         if(restart%rank==0) then
            write(0,38) restart%filename,trim(var%name),trim(what),actual_size,expected_size
38          format(A,': var "',A,'" field ',A,' has size ',I0,' instead of expected size ',I0)
         endif
         call MPI_Abort(MPI_COMM_WORLD,1,ierr)
      endif
    end subroutine check_size

  end subroutine read_var

  ! --------------------------------------------------------------------

  subroutine get_global_io_domain_indices(restart, name, global_start, global_end, indices)
    implicit none
    type(GFS_io_netCDF_type), intent(inout) :: restart
    character(len=*), intent(in) :: name
    integer, intent(out) :: global_start, global_end
    integer, allocatable :: indices(:)

    integer :: idx, i, dimlen

    idx = find_dim(restart,name)
    dimlen = restart%dims(idx)%dimlen
    global_start = 1
    global_end = dimlen
    allocate(indices(dimlen))
    do i=1,dimlen
       indices(i) = i
    enddo
  end subroutine get_global_io_domain_indices

  ! --------------------------------------------------------------------

  logical function variable_exists(restart,name)
    use netcdf
    implicit none
    type(GFS_io_netCDF_type), intent(inout) :: restart
    character(*), intent(in) :: name

    integer :: varid, ierr

    if(restart%write) then
       variable_exists = find_var(restart,name,abort=.false.) > 0
    else
       variable_exists = NF90_NOERR==nf90_inq_varid(restart%ncid,trim(name),varid)
    endif
  end function variable_exists

  ! --------------------------------------------------------------------

  subroutine handle_ncerr(restart,what,err)
    use mpi
    use netcdf

    implicit none

    type(GFS_io_netCDF_type) :: restart
    character(len=*), intent(in) :: what
    integer, intent(in) :: err
    integer :: maxerr, mpierr

    if(err/=NF90_NOERR) then
38     format(A,": ",A,": NetCDF error ",I0,": ",A)
       write(0,38) trim(restart%filename),trim(what),err,trim(nf90_strerror(err))
       call MPI_Abort(MPI_COMM_WORLD,1,mpierr)
    endif

    if(super_verbose) then
11     format(A,": ",A,' (success)')
       write(0,11) trim(restart%filename),trim(what)
    endif
  end subroutine handle_ncerr

  ! --------------------------------------------------------------------

  integer function next_dim_idx(restart)
    use mpi
    implicit none
    type(GFS_io_netCDF_type), intent(inout) :: restart
    integer :: ierr

    restart%ndims=restart%ndims+1
    next_dim_idx=restart%ndims

    if(restart%ndims>GFS_io_max_dims) then
       if(restart%rank==0) then
38        format("FV3GFS_io_netcdf supports only ",I0," dimensions. Increase GFS_io_max_dims in FV3GFS_io_netcdf.F90.")
          write(0,38) GFS_io_max_dims
       endif
       call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    endif

    restart%dims(next_dim_idx)%name(:)=' '
    restart%dims(next_dim_idx)%dimid=GFS_io_invalid_id
    restart%dims(next_dim_idx)%dimlen=0
    restart%dims(next_dim_idx)%local_start=0
    restart%dims(next_dim_idx)%local_end=0
  end function next_dim_idx

  ! --------------------------------------------------------------------

  integer function next_var_idx(restart)
    use mpi
    implicit none
    type(GFS_io_netCDF_type), intent(inout) :: restart
    integer :: ierr, a

    restart%nvars=restart%nvars+1
    next_var_idx=restart%nvars

    if(restart%nvars>GFS_io_max_vars) then
       if(restart%rank==0) then
38        format("FV3GFS_io_netcdf supports only ",I0," varaibles. Increase GFS_io_max_vars in FV3GFS_io_netcdf.F90.")
          write(0,38) GFS_io_max_vars
       endif
       call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    endif

    restart%vars(next_var_idx)%name(:)=' '
    restart%vars(next_var_idx)%varid=GFS_io_invalid_id
    restart%vars(next_var_idx)%ndims=0
    restart%vars(next_var_idx)%xtype=GFS_io_invalid_id
    restart%vars(next_var_idx)%dimids=GFS_io_invalid_id
    restart%vars(next_var_idx)%dimindex=GFS_io_invalid_id
    restart%vars(next_var_idx)%nattr=0

    do a=1,GFS_io_max_var_attr
       nullify(restart%vars(next_var_idx)%attr(a)%name)
       nullify(restart%vars(next_var_idx)%attr(a)%value)
    enddo
    
    nullify(restart%vars(next_var_idx)%phys1)
    nullify(restart%vars(next_var_idx)%phys2)
    nullify(restart%vars(next_var_idx)%phys3)
    nullify(restart%vars(next_var_idx)%real1)
    nullify(restart%vars(next_var_idx)%real2)
    nullify(restart%vars(next_var_idx)%real3)
    nullify(restart%vars(next_var_idx)%int1)
    nullify(restart%vars(next_var_idx)%int2)
    nullify(restart%vars(next_var_idx)%int3)
    
    restart%vars(next_var_idx)%is_optional = .false.
    restart%vars(next_var_idx)%write_data_called = .false.
  end function next_var_idx

  ! --------------------------------------------------------------------

  subroutine write_data_int(restart,name,buffer)
    use mpi
    use netcdf
    implicit none
    type(GFS_io_netCDF_type), intent(inout) :: restart
    character(*), intent(in) :: name
    integer, intent(in) :: buffer(:)
    integer :: idx, idim, expected_size, ierr

    if(restart%write) then
       idx=find_var(restart,name)

       expected_size=1
       do idim=1,restart%vars(idx)%ndims
          expected_size = expected_size*restart%dims(restart%vars(idx)%dimindex(idim))%dimlen
       end do

       if(size(buffer)/=expected_size) then
          if(restart%rank==0) then
38           format(A,': In write_data_int, array size ',I0,' does not match expected size ',I0)
             write(0,38) restart%filename,size(buffer),expected_size
          endif
          call MPI_Abort(MPI_COMM_WORLD,1,ierr)
       endif

       restart%vars(idx)%write_data_called=.true.
       allocate(restart%vars(idx)%int1(size(buffer)))
       restart%vars(idx)%int1=buffer
    endif
  end subroutine write_data_int

  ! --------------------------------------------------------------------

  subroutine write_data_int_scalar(restart,name,buffer)
    use mpi
    use netcdf
    implicit none
    type(GFS_io_netCDF_type), intent(inout) :: restart
    character(*), intent(in) :: name
    integer, intent(in) :: buffer
    integer :: idx, idim, expected_size, ierr

    if(restart%write) then
       idx=find_var(restart,name)

       expected_size=1
       do idim=1,restart%vars(idx)%ndims
          expected_size = expected_size*restart%dims(restart%vars(idx)%dimindex(idim))%dimlen
       end do

       if(1/=expected_size) then
          if(restart%rank==0) then
38           format(A,': In write_data_int, array size ',I0,' does not match expected size ',I0)
             write(0,38) restart%filename,1,expected_size
          endif
          call MPI_Abort(MPI_COMM_WORLD,1,ierr)
       endif

       restart%vars(idx)%write_data_called=.true.
       allocate(restart%vars(idx)%int1(1))
       restart%vars(idx)%int1=buffer
    endif
  end subroutine write_data_int_scalar

  ! --------------------------------------------------------------------

  subroutine write_data_real(restart,name,buffer)
    use mpi
    use netcdf
    implicit none
    type(GFS_io_netCDF_type), intent(inout) :: restart
    character(*), intent(in) :: name
    real, intent(in) :: buffer(:)
    integer :: idx, idim, expected_size, ierr

    if(restart%write) then
       idx=find_var(restart,name)

       expected_size=1
       do idim=1,restart%vars(idx)%ndims
          expected_size = expected_size*restart%dims(restart%vars(idx)%dimindex(idim))%dimlen
       end do

       if(size(buffer)/=expected_size) then
          if(restart%rank==0) then
38           format(A,': In write_data_real, array size ',I0,' does not match expected size ',I0)
             write(0,38) restart%filename,size(buffer),expected_size
          endif
          call MPI_Abort(MPI_COMM_WORLD,1,ierr)
       endif

       restart%vars(idx)%write_data_called=.true.
       allocate(restart%vars(idx)%real1(size(buffer)))
       restart%vars(idx)%real1=buffer
    endif
  end subroutine write_data_real

  ! --------------------------------------------------------------------

  subroutine write_data_phys(restart,name,buffer)
    use mpi
    use netcdf
    implicit none
    type(GFS_io_netCDF_type), intent(inout) :: restart
    character(*), intent(in) :: name
    real(kind_phys), intent(in) :: buffer(:)
    integer :: idx, idim, expected_size, ierr

    if(restart%write) then
       idx=find_var(restart,name)

       expected_size=1
       do idim=1,restart%vars(idx)%ndims
          expected_size = expected_size*restart%dims(restart%vars(idx)%dimindex(idim))%dimlen
       end do

       if(size(buffer)/=expected_size) then
          if(restart%rank==0) then
38           format(A,': In write_data_phys, array size ',I0,' does not match expected size ',I0)
             write(0,38) restart%filename,size(buffer),expected_size
          endif
          call MPI_Abort(MPI_COMM_WORLD,1,ierr)
       endif

       restart%vars(idx)%write_data_called=.true.
       allocate(restart%vars(idx)%phys1(size(buffer)))
       restart%vars(idx)%phys1=buffer
    endif
  end subroutine write_data_phys

  ! --------------------------------------------------------------------

  subroutine register_axis_x_or_y(restart,name,axis_name)
    use mpi
    use netcdf

    implicit none
    type(GFS_io_netCDF_type), intent(inout) :: restart
    character(len=*), intent(in) :: name, axis_name
    integer :: idx, ierr, dlen

    idx=next_dim_idx(restart)

    restart%dims(idx)%name=name

    if(axis_name=='X') then
       dlen=restart%cnx
       restart%dims(idx)%local_start = restart%isc
       restart%dims(idx)%local_end = restart%iec
    else if(axis_name=='Y') then
       dlen=restart%cny
       restart%dims(idx)%local_start = restart%jsc
       restart%dims(idx)%local_end = restart%jec
    else
       if(restart%rank==0) then
          write(0,'(A,A,A)') 'ERROR: invalid axis_name "',trim(axis_name),'"'
       endif
       call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    endif

    if(restart%write) then
       restart%dims(idx)%dimlen=dlen
    else
       call handle_ncerr(restart,"inquire dim "//trim(name),&
            nf90_inq_dimid(restart%ncid,trim(name),restart%dims(idx)%dimid))
       call handle_ncerr(restart,"inquire dim len "//trim(name), &
            nf90_inquire_dimension(restart%ncid,restart%dims(idx)%dimid,&
            len=restart%dims(idx)%dimlen))
       if(dlen/=restart%dims(idx)%dimlen) then
          if(restart%rank==0) then
38           format('Dimension "',A,'" for axis "',A,'" had length ',I0,' instead of expected ',I0)
             write(0,38) trim(name),trim(axis_name),restart%dims(idx)%dimlen,dlen
          endif
          call MPI_Abort(MPI_COMM_WORLD,1,ierr)
       endif
    endif
  end subroutine register_axis_x_or_y

  ! --------------------------------------------------------------------

  subroutine register_axis_len(restart,name,dimension_length)
    use mpi
    use netcdf

    implicit none
    type(GFS_io_netCDF_type), intent(inout) :: restart
    character(len=*), intent(in) :: name
    integer, intent(in) :: dimension_length
    integer :: idx, ierr

    idx=next_dim_idx(restart)

    restart%dims(idx)%name=name

    if(restart%write) then
       if(dimension_length == nf90_unlimited) then
          restart%dims(idx)%dimlen=1
       else
          restart%dims(idx)%dimlen=dimension_length
       endif
    else
       call handle_ncerr(restart,"inquire dim id "//trim(name),&
            nf90_inq_dimid(restart%ncid,trim(name),restart%dims(idx)%dimid))
       call handle_ncerr(restart,"inquire dim len "//trim(name), &
            nf90_inquire_dimension(restart%ncid,restart%dims(idx)%dimid,&
            len=restart%dims(idx)%dimlen))
       if(dimension_length/=restart%dims(idx)%dimlen) then
          if(restart%rank==0) then
38           format('Dimension "',A,'" had length ',I0,' instead of expected ',I0)
             write(0,38) trim(name),restart%dims(idx)%dimlen,dimension_length
          endif
          call MPI_Abort(MPI_COMM_WORLD,1,ierr)
       endif
    endif

    restart%dims(idx)%local_start = 1

    if(dimension_length == nf90_unlimited) then
       restart%dims(idx)%local_end = 1
    else
       restart%dims(idx)%local_end = restart%dims(idx)%dimlen
    endif
  end subroutine register_axis_len

  ! --------------------------------------------------------------------

  function phys_size()
    use iso_c_binding, only: c_sizeof
    implicit none
    character(len=6) :: phys_size
    real(kind_phys), parameter :: test_phys = 123.456
    if(c_sizeof(test_phys)==4) then
       phys_size='float '
    else
       phys_size='double'
    endif
  end function phys_size

  ! --------------------------------------------------------------------

  function real_size()
    use iso_c_binding, only: c_sizeof
    implicit none
    character(len=6) :: real_size
    real, parameter :: test_real = 123.456
    if(c_sizeof(test_real)==4) then
       real_size='float '
    else
       real_size='double'
    endif
  end function real_size

  ! --------------------------------------------------------------------

  integer function next_attr_idx(restart,var)
    use mpi
    implicit none
    type(GFS_io_netCDF_type), intent(inout) :: restart
    type(GFS_io_var_type), intent(inout) :: var
    integer :: ierr

    var%nattr = var%nattr+1
    next_attr_idx = var%nattr
    if(var%nattr>=GFS_io_max_var_attr) then
       if(restart%rank==0) then
38        format("FV3GFS_io_netcdf supports only ",I0," attributes per variable. Increase GFS_io_max_var_attr in FV3GFS_io_netcdf.F90.")
          write(0,38) GFS_io_max_var_attr
       endif
       call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    endif
  end function next_attr_idx

  ! --------------------------------------------------------------------

  subroutine register_variable_attribute(restart,var_name,att_name,att_value,str_len)
    use netcdf
    implicit none
    type(GFS_io_netCDF_type), intent(inout) :: restart
    character(*), intent(in) :: var_name,att_name,att_value
    integer, intent(in) :: str_len

    integer :: vidx, aidx

    if(restart%write) then
       vidx = find_var(restart,var_name)
       aidx = next_attr_idx(restart,restart%vars(vidx))
       allocate(character(len_trim(att_name)) :: restart%vars(vidx)%attr(aidx)%name)
       restart%vars(vidx)%attr(aidx)%name = att_name
       allocate(character(str_len) :: restart%vars(vidx)%attr(aidx)%value)
       restart%vars(vidx)%attr(aidx)%value = att_value
    end if
  end subroutine register_variable_attribute

  ! --------------------------------------------------------------------

  subroutine register_restart_field_int1(restart, name, data, dimensions, is_optional)
    implicit none
    type(GFS_io_netCDF_type), intent(inout) :: restart
    integer, pointer :: data(:)
    character(*), intent(in) :: dimensions(:)
    character(*), intent(in) :: name
    logical, optional :: is_optional
    call register_variable(restart,name,'integer',dimensions,is_optional=is_optional)
    restart%vars(restart%nvars)%int1 => data
  end subroutine register_restart_field_int1

  ! --------------------------------------------------------------------

  subroutine register_restart_field_int2(restart, name, data, dimensions, is_optional)
    implicit none
    type(GFS_io_netCDF_type), intent(inout) :: restart
    integer, pointer :: data(:,:)
    character(*), intent(in) :: dimensions(:)
    character(*), intent(in) :: name
    logical, optional :: is_optional
    call register_variable(restart,name,'integer',dimensions,is_optional=is_optional)
    restart%vars(restart%nvars)%int2 => data
  end subroutine register_restart_field_int2

  ! --------------------------------------------------------------------

  subroutine register_restart_field_int3(restart, name, data, dimensions, is_optional)
    implicit none
    type(GFS_io_netCDF_type), intent(inout) :: restart
    integer, pointer :: data(:,:,:)
    character(*), intent(in) :: dimensions(:)
    character(*), intent(in) :: name
    logical, optional :: is_optional
    call register_variable(restart,name,'integer',dimensions,is_optional=is_optional)
    restart%vars(restart%nvars)%int3 => data
  end subroutine register_restart_field_int3

  ! --------------------------------------------------------------------

  subroutine register_restart_field_phys1(restart, name, data, dimensions, is_optional)
    implicit none
    type(GFS_io_netCDF_type), intent(inout) :: restart
    real(kind_phys), pointer :: data(:)
    character(*), intent(in) :: dimensions(:)
    character(*), intent(in) :: name
    logical, optional :: is_optional
    call register_variable(restart,name,phys_size(),dimensions,is_optional=is_optional)
    restart%vars(restart%nvars)%phys1 => data
  end subroutine register_restart_field_phys1

  ! --------------------------------------------------------------------

  subroutine register_restart_field_phys2(restart, name, data, dimensions, is_optional)
    implicit none
    type(GFS_io_netCDF_type), intent(inout) :: restart
    real(kind_phys), pointer :: data(:,:)
    character(*), intent(in) :: dimensions(:)
    character(*), intent(in) :: name
    logical, optional :: is_optional
    call register_variable(restart,name,phys_size(),dimensions,is_optional=is_optional)
    restart%vars(restart%nvars)%phys2 => data
  end subroutine register_restart_field_phys2

  ! --------------------------------------------------------------------

  subroutine register_restart_field_phys3(restart, name, data, dimensions, is_optional)
    implicit none
    type(GFS_io_netCDF_type), intent(inout) :: restart
    real(kind_phys), pointer :: data(:,:,:)
    character(*), intent(in) :: dimensions(:)
    character(*), intent(in) :: name
    logical, optional :: is_optional
    call register_variable(restart,name,phys_size(),dimensions,is_optional=is_optional)
    restart%vars(restart%nvars)%phys3 => data
  end subroutine register_restart_field_phys3

  ! --------------------------------------------------------------------

  subroutine register_restart_field_real1(restart, name, data, dimensions, is_optional)
    implicit none
    type(GFS_io_netCDF_type), intent(inout) :: restart
    real, pointer :: data(:)
    character(*), intent(in) :: dimensions(:)
    character(*), intent(in) :: name
    logical, optional :: is_optional
    call register_variable(restart,name,real_size(),dimensions,is_optional=is_optional)
    restart%vars(restart%nvars)%real1 => data
  end subroutine register_restart_field_real1

  ! --------------------------------------------------------------------

  subroutine register_restart_field_real2(restart, name, data, dimensions, is_optional)
    implicit none
    type(GFS_io_netCDF_type), intent(inout) :: restart
    real, pointer :: data(:,:)
    character(*), intent(in) :: dimensions(:)
    character(*), intent(in) :: name
    logical, optional :: is_optional
    call register_variable(restart,name,real_size(),dimensions,is_optional=is_optional)
    restart%vars(restart%nvars)%real2 => data
  end subroutine register_restart_field_real2

  ! --------------------------------------------------------------------

  subroutine register_restart_field_real3(restart, name, data, dimensions, is_optional)
    implicit none
    type(GFS_io_netCDF_type), intent(inout) :: restart
    real, pointer :: data(:,:,:)
    character(*), intent(in) :: dimensions(:)
    character(*), intent(in) :: name
    logical, optional :: is_optional
    call register_variable(restart,name,real_size(),dimensions,is_optional=is_optional)
    restart%vars(restart%nvars)%real3 => data
  end subroutine register_restart_field_real3

  ! --------------------------------------------------------------------

  subroutine register_field(restart, name, type, dims)
    implicit none
    type(GFS_io_netCDF_type), intent(inout) :: restart
    character(*), intent(in) :: dims(:)
    character(len=*), intent(in) :: name, type
    call register_variable(restart,name,type,dims,is_optional=.true.)
    ! Note: not associating any pointers means the varaible will
    ! not be automatically read.
  end subroutine register_field

  ! --------------------------------------------------------------------

  subroutine register_variable(restart, name, type, dims, is_optional)
    use mpi
    use netcdf
    implicit none

    type(GFS_io_netCDF_type), intent(inout) :: restart
    character(*), intent(in) :: name, type
    character(*), intent(in) :: dims(:)
    logical, optional :: is_optional

    integer :: idx, ierr, xtype, i, idim, didx
    integer :: ncerr
    logical :: var_missing, discard

    idx = next_var_idx(restart)

    restart%vars(idx)%name = name

    if(type=='int' .or. type=='integer') then
       xtype = NF90_INT
    else if(type=='double') then
       xtype = NF90_DOUBLE
    else if(type=='float') then
       xtype = NF90_FLOAT
    else
       if(restart%rank==0) then
31        format('Variable "',A,'" has invalid type name "',A,'". I only understand double, float, int, and integer.')
          write(0,31) trim(name), trim(type)
       endif
       call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    endif
    restart%vars(idx)%xtype=xtype

    restart%vars(idx)%ndims=size(dims)
    if(restart%vars(idx)%ndims>GFS_io_max_var_dims) then
       if(restart%rank==0) then
38        format('Variable "',A,'" has ',I0,' dimensions which is larger than the limit of ',I0,'. Increase GFS_io_max_var_dims in FV3GFS_io_netcdf.F90.')
          write(0,38) trim(name),size(dims),GFS_io_max_var_dims
       endif
       call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    endif

    restart%vars(idx)%is_optional = .false.
    if(present(is_optional)) then
       if(is_optional) then
          restart%vars(idx)%is_optional=is_optional
       endif
    endif

    if(restart%write) then
       do i=1,size(dims)
          restart%vars(idx)%dimindex(i) = find_dim(restart,dims(i))
       enddo
    else ! Read.
       ncerr = nf90_inq_varid(restart%ncid,trim(name),restart%vars(idx)%varid)
       var_missing = ncerr/=NF90_NOERR

       if(var_missing) then
          if(.not.restart%vars(idx)%is_optional) then
             if(restart%rank==0) then
83              format(A,': mandatory restart variable ',A,' is missing')
                write(0,83) trim(restart%filename), trim(name)
             endif
             call MPI_Abort(MPI_COMM_WORLD,1,ierr)
          endif
          restart%nvars = restart%nvars - 1
          return
       endif
       call handle_ncerr(restart,"inquire variable dimensions "//trim(name), &
            nf90_inquire_variable(ncid=restart%ncid, &
            varid=restart%vars(idx)%varid,  & ! input
            dimids=restart%vars(idx)%dimids, & ! output
            ndims=restart%vars(idx)%ndims))   ! output
       do i=1,size(dims)
          restart%vars(idx)%dimindex(i) = find_dimid(restart,restart%vars(idx)%dimids(i))
       enddo
    endif
  end subroutine register_variable

  ! --------------------------------------------------------------------

  integer function find_dim(restart,name,abort)
    ! Find index within restart%dims of dimension with specified name
    use mpi
    implicit none
    type(GFS_io_netCDF_type), intent(inout) :: restart
    character(len=*), intent(in) :: name
    logical, intent(in), optional :: abort
    integer :: idim,ierr
    logical :: do_abort

    do_abort=.true.
    if(present(abort)) do_abort=abort
    
    find_dim=-1

    do idim=1,restart%ndims
       if(restart%dims(idim)%name==name) then
          find_dim=idim
          return
       endif
    enddo

    if(do_abort) then
       if(restart%rank==0) then
38        format(A': No dimension defined with name "',A,'"')
          write(0,38) restart%filename,trim(name)
       endif
       call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    endif
  end function find_dim

  ! --------------------------------------------------------------------

  integer function find_var(restart,name,abort)
    ! Find index within restart%vars of variable with specified name
    use mpi
    implicit none
    type(GFS_io_netCDF_type), intent(inout) :: restart
    character(len=*), intent(in) :: name
    logical, intent(in), optional :: abort

    integer :: ivar,ierr
    logical :: do_abort

    do_abort=.true.
    if(present(abort)) do_abort=abort
    
    find_var=-1

    do ivar=1,restart%nvars
       if(restart%vars(ivar)%name==name) then
          find_var=ivar
          return
       endif
    enddo

    if(abort) then
       if(restart%rank==0) then
38        format(A': No variable defined with name "',A,'"')
          write(0,38) restart%filename,trim(name)
       endif
       call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    endif
  end function find_var

  ! --------------------------------------------------------------------

  integer function find_dimid(restart,dimid,abort)
    ! Find index within restart%dims of dimension with specified NetCDF dimension id
    use mpi
    implicit none
    type(GFS_io_netCDF_type), intent(inout) :: restart
    integer, intent(in) :: dimid
    logical, intent(in), optional :: abort
    integer :: idim,ierr
    logical :: do_abort

    do_abort=.true.
    if(present(abort)) do_abort=abort
    
    find_dimid=-1

    do idim=1,restart%ndims
       if(restart%dims(idim)%dimid==dimid) then
          find_dimid=idim
          return
       endif
    enddo

    if(do_abort) then
       if(restart%rank==0) then
38        format(A': No dimension defined with NetCDF ID "',I0,'"')
          write(0,38) restart%filename,dimid
       endif
       call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    endif
  end function find_dimid

  ! --------------------------------------------------------------------

  logical function allreduce_and(flag,comm)
    ! Replaces flag with a logical "and" of all values of that flag
    ! across the communicator, and returns that value.
    use mpi
    implicit none
    integer, intent(in) :: comm
    logical, intent(inout) :: flag
    integer :: ierr
    call MPI_Allreduce(flag,allreduce_and,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    flag=allreduce_and
  end function allreduce_and
end module FV3GFS_io_netCDF_mod
