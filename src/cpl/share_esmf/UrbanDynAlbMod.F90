module UrbanDynAlbMod
!----------------------------------------------------------------------- 
!
!DESCRIPTION:
!Transient Urban Albedo Input Stream Data
!The time-varing urban albedo is read in using this module (a stream)
!
!USES:
 use ESMF            , only : ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU, ESMF_Finalize, ESMF_END_ABORT
 use dshr_strdata_mod, only : shr_strdata_type
 use shr_kind_mod    , only : r8 => shr_kind_r8, CL => shr_kind_CL
 use shr_log_mod     , only : errMsg => shr_log_errMsg
 use abortutils      , only : endrun
 use decompMod       , only : bounds_type, subgrid_level_landunit
 use clm_varctl      , only : iulog, FL => fname_len
 use LandunitType    , only : lun
 use GridcellType    , only : grc
 use clm_varcon      , only : spval
 use landunit_varcon , only : isturb_MIN, isturb_MAX         ! isturb_MIN = 7, isturb_MAX = 9 or isturb_MAX = 16 (use_lcz=.true.)
 use clm_varpar      , only : numrad
 use UrbanParamsType , only : transient_urbanalbedo_roof
 use UrbanParamsType , only : transient_urbanalbedo_improad
 use UrbanParamsType , only : transient_urbanalbedo_wall
 ! 
 implicit none
 save
 private
 
 ! !PUBLIC TYPE
 type, public :: urbanalbtv_type
    ! urban wall albedo inputs
    real(r8), public, pointer :: dyn_alb_roof_dir        (:,:) ! dynamic lun direct  roof albedo
    real(r8), public, pointer :: dyn_alb_roof_dif        (:,:) ! dynamic lun diffuse roof albedo
    real(r8), public, pointer :: dyn_alb_improad_dir     (:,:) ! dynamic lun direct roof albedo
    real(r8), public, pointer :: dyn_alb_improad_dif     (:,:) ! dynamic lun diffuse roof albedo
    real(r8), public, pointer :: dyn_alb_wall_dir        (:,:) ! dynamic lun direct wall albedo
    real(r8), public, pointer :: dyn_alb_wall_dif        (:,:) ! dynamic lun diffuse wall albedo
    ! 
    type(shr_strdata_type)    :: sdat_urbanalbtvroof         ! urban time varying roof albedo data stream
    type(shr_strdata_type)    :: sdat_urbanalbtvimproad      ! urban time varying improad albedo data stream
    type(shr_strdata_type)    :: sdat_urbanalbtvwall         ! urban time varying wall albedo data stream
    
   contains
     ! !PUBLIC MEMBER FUNCTIONS:
     procedure, public :: dynalbinit                         ! Allocate and initialize urbanalbtv
     procedure, public :: urbanalbtvroof_init                ! Initialize urban wall albedo time varying stream
     procedure, public :: urbanalbtvroof_interp              ! Interpolate urban roof alebdo time varying stream
     procedure, public :: urbanalbtvimproad_init             ! Initialize urban improad albedo time varying stream
     procedure, public :: urbanalbtvimproad_interp           ! Interpolate urban improad alebdo time varying stream
     procedure, public :: urbanalbtvwall_init                ! Initialize urban wall albedo time varying stream
     procedure, public :: urbanalbtvwall_interp              ! Interpolate urban wall alebdo time varying stream
 end type urbanalbtv_type

  integer      , private              :: stream_varname_MIN       ! minimum index
  integer      , private              :: stream_varname_MAX       ! maximum index
  character(30), private, pointer     :: stream_var_name_roof(:)
  character(30), private, pointer     :: stream_var_name_improad(:)
  character(30), private, pointer     :: stream_var_name_wall(:)

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  !----------------------------------------------------------------------- 
 contains
  !-----------------------------------------------------------------------
  !
  subroutine dynalbinit(this, bounds, NLFilename)
  !
  ! !DESCRIPTION:
  ! Initialize data stream information for dynamic urban albedo
  !
  ! !USES:
  use shr_infnan_mod   , only : nan => shr_infnan_nan, assignment(=)
  use histFileMod      , only : hist_addfld2d
  use clm_varctl       , only : use_lcz
  !
  ! !ARGUMENTS:
  class(urbanalbtv_type)                 :: this
  type(bounds_type)      , intent(in)    :: bounds
  character(len=*)       , intent(in)    :: NLFilename   ! Namelist filename
  !
  ! !LOCAL VARIABLES:  
  integer             :: begl, endl
  !---------------------------------------------------------------------
  !
  begl = bounds%begl; endl = bounds%endl                        ! beginning and ending landunit index
  !
  ! Determine the minimum and maximum indices
  stream_varname_MIN = 1
  if (use_lcz) then
      stream_varname_MAX = 10
  else
      stream_varname_MAX = 3
  end if
  ! 
  ! Allocate urbanalbtv data structures
  ! 
  allocate(this%dyn_alb_roof_dir        (begl:endl,numrad))   ; this%dyn_alb_roof_dir        (:,:) = nan
  allocate(this%dyn_alb_roof_dif        (begl:endl,numrad))   ; this%dyn_alb_roof_dif        (:,:) = nan 
  allocate(this%dyn_alb_improad_dir     (begl:endl,numrad))   ; this%dyn_alb_improad_dir     (:,:) = nan   
  allocate(this%dyn_alb_improad_dif     (begl:endl,numrad))   ; this%dyn_alb_improad_dif     (:,:) = nan
  allocate(this%dyn_alb_wall_dir        (begl:endl,numrad))   ; this%dyn_alb_wall_dir        (:,:) = nan   
  allocate(this%dyn_alb_wall_dif        (begl:endl,numrad))   ; this%dyn_alb_wall_dif        (:,:) = nan
  allocate(stream_var_name_roof(stream_varname_MIN:stream_varname_MAX))
  allocate(stream_var_name_improad(stream_varname_MIN:stream_varname_MAX))
  allocate(stream_var_name_wall(stream_varname_MIN:stream_varname_MAX))

  if (transient_urbanalbedo_roof) then 
     call this%urbanalbtvroof_init(bounds, NLFilename)
     call this%urbanalbtvroof_interp(bounds)
     ! Add history fields
     call hist_addfld2d (fname='DYNALB_ROOF_DIR', units='',      &
            avgflag='A', long_name='time varing urban roof albedo dir',  type2d='numrad', &
            ptr_lunit=this%dyn_alb_roof_dir, default='inactive', set_nourb=spval, &
            l2g_scale_type='unity')
     call hist_addfld2d (fname='DYNALB_ROOF_DIF', units='',      &
            avgflag='A', long_name='time varing urban roof albedo dif',  type2d='numrad',  &
            ptr_lunit=this%dyn_alb_roof_dif, default='inactive', set_nourb=spval, &
            l2g_scale_type='unity')
  end if
  
  if (transient_urbanalbedo_improad) then    
     call this%urbanalbtvimproad_init(bounds, NLFilename)
     call this%urbanalbtvimproad_interp(bounds)
     ! Add history fields
     call hist_addfld2d (fname='DYNALB_IMPROAD_DIR', units='',      &
            avgflag='A', long_name='time varing urban improad albedo dir',  type2d='numrad', &
            ptr_lunit=this%dyn_alb_improad_dir, default='inactive', set_nourb=spval, &
            l2g_scale_type='unity')
     call hist_addfld2d (fname='DYNALB_IMPROAD_DIF', units='',      &
            avgflag='A', long_name='time varing urban improad albedo dif',  type2d='numrad',  &
            ptr_lunit=this%dyn_alb_improad_dif, default='inactive', set_nourb=spval, &
            l2g_scale_type='unity')
  end if
  
  if (transient_urbanalbedo_wall) then
     call this%urbanalbtvwall_init(bounds, NLFilename)
     call this%urbanalbtvwall_interp(bounds)
     ! Add history fields 
     call hist_addfld2d (fname='DYNALB_WALL_DIR', units='',      &
            avgflag='A', long_name='time varing urban wall albedo dir',  type2d='numrad', &
            ptr_lunit=this%dyn_alb_wall_dir, default='inactive', set_nourb=spval, &
            l2g_scale_type='unity')
     call hist_addfld2d (fname='DYNALB_WALL_DIF', units='',      &
            avgflag='A', long_name='time varing urban wall albedo dif',  type2d='numrad',  &
            ptr_lunit=this%dyn_alb_wall_dif, default='inactive', set_nourb=spval, &
            l2g_scale_type='unity')
  end if   
  
 end subroutine dynalbinit

 !---------------------------------------------------------------------
 subroutine urbanalbtvroof_init(this, bounds, NLFileName)
   !
   ! !DESCRIPTION:
   ! Initialize data stream information for urban time varying roof albedo
   !
   ! !USES:
   use clm_nlUtilsMod   , only : find_nlgroup_name
   use spmdMod          , only : masterproc, mpicom, iam
   use shr_mpi_mod      , only : shr_mpi_bcast
   use dshr_strdata_mod , only : shr_strdata_init_from_inline
   use lnd_comp_shr     , only : mesh, model_clock
   use clm_varctl       , only : use_lcz
   use landunit_varcon  , only : isturb_tbd, isturb_hd, isturb_md          
   use landunit_varcon  , only : isturb_lcz1, isturb_lcz2, isturb_lcz3, &
                                 isturb_lcz4, isturb_lcz5, isturb_lcz6, &
                                 isturb_lcz7, isturb_lcz8, isturb_lcz9, &
                                 isturb_lcz10
   !
   ! !ARGUMENTS:
   implicit none
   class(urbanalbtv_type)         :: this
   type(bounds_type), intent(in)  :: bounds
   character(len=*),  intent(in)  :: NLFilename   ! Namelist filename
   ! 
   ! !LOCAL VARIABLES:
   integer            :: n
   integer            :: stream_year_first_urbanalbtvroof            ! first year in urban roof albedo tv stream to use
   integer            :: stream_year_last_urbanalbtvroof             ! last year in urban roof albedo tv stream to use
   integer            :: model_year_align_urbanalbtvroof             ! align stream_year_first_urbanalbtvroof with this model year
   integer            :: nu_nml                                      ! unit for namelist file 
   integer            :: nml_error                                   ! namelist i/o error flag
   character(len=CL)  :: stream_fldFileName_urbanalbtvroof           ! urban roof albedo time-varying streams filename
   character(len=FL)  :: stream_meshfile_urbanalbtvroof              ! urban roof albedo time-varying mesh filename
   character(len=CL)  :: urbanalbtvroofmapalgo = 'nn'                ! mapping alogrithm for urban ac
   character(len=CL)  :: urbanalbtvroof_tintalgo = 'linear'          ! time interpolation alogrithm 
   integer            :: rc                                          ! error code
   character(*), parameter :: subName = "('urbanalbtvroof_init')"
   !-----------------------------------------------------------------------
   namelist /urbanalbtvroof_streams/       &
        stream_year_first_urbanalbtvroof,  &  
        stream_year_last_urbanalbtvroof,   &  
        model_year_align_urbanalbtvroof,   &  
        urbanalbtvroofmapalgo,             &  
        stream_fldFileName_urbanalbtvroof, &  
        stream_meshfile_urbanalbtvroof,    &          
        urbanalbtvroof_tintalgo  
   !-----------------------------------------------------------------------       
   !               
   ! Default values for namelist
   stream_year_first_urbanalbtvroof  = 1      ! first year in stream to use
   stream_year_last_urbanalbtvroof   = 1      ! last  year in stream to use
   model_year_align_urbanalbtvroof   = 1      ! align stream_year_first_urbanalbtvroof with this model year
   stream_fldFileName_urbanalbtvroof = ' '
   stream_meshfile_urbanalbtvroof    = ' '
   
   ! create the field list for urban albedo fields
   if (.not. use_lcz) then 
      stream_var_name_roof(isturb_tbd -6) = "dyn_alb_roof_TBD"
      stream_var_name_roof(isturb_hd -6)  = "dyn_alb_roof_HD"
      stream_var_name_roof(isturb_md -6)  = "dyn_alb_roof_MD"   
   else
      stream_var_name_roof(isturb_lcz1 -6) = "dyn_alb_roof_LCZ1"
      stream_var_name_roof(isturb_lcz2 -6) = "dyn_alb_roof_LCZ2"
      stream_var_name_roof(isturb_lcz3 -6) = "dyn_alb_roof_LCZ3"
      stream_var_name_roof(isturb_lcz4 -6) = "dyn_alb_roof_LCZ4"
      stream_var_name_roof(isturb_lcz5 -6) = "dyn_alb_roof_LCZ5"
      stream_var_name_roof(isturb_lcz6 -6) = "dyn_alb_roof_LCZ6"
      stream_var_name_roof(isturb_lcz7 -6) = "dyn_alb_roof_LCZ7"
      stream_var_name_roof(isturb_lcz8 -6) = "dyn_alb_roof_LCZ8"
      stream_var_name_roof(isturb_lcz9 -6) = "dyn_alb_roof_LCZ9"
      stream_var_name_roof(isturb_lcz10-6) = "dyn_alb_roof_LCZ10"
   end if

   ! Read urbanalbtvroof_streams namelist
   if (masterproc) then
      open( newunit=nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
      call find_nlgroup_name(nu_nml, 'urbanalbtvroof_streams', status=nml_error)
      if (nml_error == 0) then
         read(nu_nml, nml=urbanalbtvroof_streams,iostat=nml_error) 
         if (nml_error /= 0) then
            call endrun(msg='ERROR reading urbanalbtvroof_streams namelist'//errMsg(sourcefile, __LINE__))
         end if
      else
          call endrun(subname // ':: ERROR finding urbanalbtvroof_streams namelist')   
      end if
      close(nu_nml)
   endif

   call shr_mpi_bcast(stream_year_first_urbanalbtvroof  , mpicom)
   call shr_mpi_bcast(stream_year_last_urbanalbtvroof   , mpicom)
   call shr_mpi_bcast(model_year_align_urbanalbtvroof   , mpicom)
   call shr_mpi_bcast(stream_fldFileName_urbanalbtvroof , mpicom)
   call shr_mpi_bcast(stream_meshfile_urbanalbtvroof    , mpicom)
   call shr_mpi_bcast(urbanalbtvroof_tintalgo           , mpicom)
   
   if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'Attempting to read time varying urban roof albedo parameters......'
       write(iulog,'(a)') 'urbanalbtvroof_streams settings:'
       write(iulog,'(a,i8)') '  stream_year_first_urbanalbtvroof  = ',stream_year_first_urbanalbtvroof
       write(iulog,'(a,i8)') '  stream_year_last_urbanalbtvroof   = ',stream_year_last_urbanalbtvroof
       write(iulog,'(a,i8)') '  model_year_align_urbanalbtvroof   = ',model_year_align_urbanalbtvroof
       write(iulog,'(a,a)' ) '  stream_fldFileName_urbanalbtvroof = ',stream_fldFileName_urbanalbtvroof
       write(iulog,'(a,a)' ) '  stream_meshfile_urbanalbtvroof    = ',stream_meshfile_urbanalbtvroof
       write(iulog,'(a,a)' ) '  urbanalbtvroof_tintalgo           = ',urbanalbtvroof_tintalgo
       write(iulog,*) 'Read in urbanalbtvroof_streams namelist from:',trim(NLFilename)
       do n = stream_varname_MIN,stream_varname_MAX
          write(iulog,'(a,a)' ) '  stream_var_name_roof         = ',trim(stream_var_name_roof(n))
       end do
   endif
    
    call shr_strdata_init_from_inline(this%sdat_urbanalbtvroof,                  &
         my_task             = iam,                                              &
         logunit             = iulog,                                            &
         compname            = 'LND',                                            &
         model_clock         = model_clock,                                      &
         model_mesh          = mesh,                                             &
         stream_meshfile     = trim(stream_meshfile_urbanalbtvroof),             &
         stream_lev_dimname  = 'null',                                           &
         stream_mapalgo      = trim(urbanalbtvroofmapalgo),                      &
         stream_filenames    = (/trim(stream_fldfilename_urbanalbtvroof)/),      &
         stream_fldlistFile  = stream_var_name_roof(stream_varname_MIN:stream_varname_MAX), &
         stream_fldListModel = stream_var_name_roof(stream_varname_MIN:stream_varname_MAX), &
         stream_yearFirst    = stream_year_first_urbanalbtvroof,                 &
         stream_yearLast     = stream_year_last_urbanalbtvroof,                  &
         stream_yearAlign    = model_year_align_urbanalbtvroof,                  &
         stream_offset       = 0,                                                &
         stream_taxmode      = 'extend',                                         &
         stream_dtlimit      = 1.0e30_r8,                                        &
         stream_tintalgo     = urbanalbtvroof_tintalgo,                          &
         stream_name         = 'Urban time varying roof albedo data',            &
         rc                  = rc)

   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
   end if
    
  end subroutine urbanalbtvroof_init

  !==============================================================================
 
  subroutine urbanalbtvroof_interp(this, bounds)
    ! !DESCRIPTION:
    ! Interpolate data stream information for urban time varying albedo.
    ! 
    ! !USES:
    use clm_time_manager  , only : get_curr_date
    use clm_instur        , only : urban_valid
    use dshr_methods_mod  , only : dshr_fldbun_getfldptr
    use dshr_strdata_mod  , only : shr_strdata_advance
    use shr_infnan_mod    , only : nan => shr_infnan_nan, assignment(=)
    ! 
    ! !ARGUMENTS:
    ! 
    class(urbanalbtv_type)           :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    !
    logical :: found
    integer :: l, ig, g, ip, n, ib    
    integer :: year    ! year (0, ...) for nstep+1
    integer :: mon     ! month (1, ..., 12) for nstep+1
    integer :: day     ! day of month (1, ..., 31) for nstep+1
    integer :: sec     ! seconds into current date for nstep+1
    integer :: mcdate  ! Current model date (yyyymmdd)
    integer :: lindx   ! landunit index
    integer :: gindx   ! gridcell index
    integer :: lsize
    integer :: rc
    real(r8), pointer :: dataptr1d(:)
    real(r8), pointer :: dataptr2d(:,:)
    ! 
    !-----------------------------------------------------------------------
    ! 
    ! Advance sdat stream
    !
    call get_curr_date(year, mon, day, sec)
    !
    ! packing the date into an integer
    mcdate = year*10000 + mon*100 + day

    call shr_strdata_advance(this%sdat_urbanalbtvroof, ymd=mcdate, tod=sec, logunit=iulog, istr='hdmdyn', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if
    !
    ! Create 2d array for all stream variable data
    lsize = bounds%endg - bounds%begg + 1
    allocate(dataptr2d(lsize, stream_varname_MIN:stream_varname_MAX))
    do n = stream_varname_MIN,stream_varname_MAX
       call dshr_fldbun_getFldPtr(this%sdat_urbanalbtvroof%pstrm(1)%fldbun_model, trim(stream_var_name_roof(n)), &
            fldptr1=dataptr1d, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if
       ! Note that the size of dataptr1d includes ocean points so it will be around 3x larger than lsize
       ! So an explicit loop is required here
       do g = 1,lsize
          dataptr2d(g,n) = dataptr1d(g)
       end do
    end do

    ! Determine this%tbuilding_max (and this%p_ac, if applicable) for all landunits
    do l = bounds%begl,bounds%endl
       if (lun%urbpoi(l)) then
          ! Note that since l is within [begl, endl] bounds, we can assume
          ! lun%gricell(l) is within [begg, endg]
          ig = lun%gridcell(l) - bounds%begg + 1
          do ib = 1,numrad 
             do n = stream_varname_MIN,stream_varname_MAX
                if (stream_var_name_roof((lun%itype(l)-6)) == stream_var_name_roof(n)) then
                   this%dyn_alb_roof_dir(l,ib) = dataptr2d(ig,n)
                   this%dyn_alb_roof_dif(l,ib) = dataptr2d(ig,n)
                end if
             end do
          end do 
       else
           do ib = 1,numrad
              this%dyn_alb_roof_dir(l,ib) = spval
              this%dyn_alb_roof_dif(l,ib) = spval  
           end do  
       end if
    end do
    deallocate(dataptr2d)

    ! Error check
    found = .false.
    do l = bounds%begl,bounds%endl
       if (lun%urbpoi(l)) then
          do g = bounds%begg,bounds%endg
             if (g == lun%gridcell(l)) exit
          end do
          ! Check for valid urban data
          do ib = 1,numrad
             if ( .not. urban_valid(g) .or. (this%dyn_alb_roof_dir(l,ib) <= 0._r8) .or. (this%dyn_alb_roof_dif(l,ib) <= 0._r8)) then
                found = .true.
                gindx = g
                lindx = l
                exit
             end if
          end do   
       end if
    end do
    if ( found ) then
       write(iulog,*)'ERROR: no valid urban data for g= ',gindx
       write(iulog,*)'landunit type:   ',lun%itype(lindx)
       write(iulog,*)'urban_valid:     ',urban_valid(gindx)
       write(iulog,*)'dyn_alb_roof_dir:  ',this%dyn_alb_roof_dir(lindx,:)
       write(iulog,*)'dyn_alb_roof_dif:  ',this%dyn_alb_roof_dif(lindx,:)
       call endrun(subgrid_index=lindx, subgrid_level=subgrid_level_landunit, &
            msg=errmsg(sourcefile, __LINE__))
    end if

  end subroutine urbanalbtvroof_interp

 !---------------------------------------------------------------------
 subroutine urbanalbtvwall_init(this, bounds, NLFilename)
   !
   ! !DESCRIPTION:
   ! Initialize data stream information for urban time varying wall albedo
   !
   ! !USES:
   use clm_nlUtilsMod   , only : find_nlgroup_name
   use spmdMod          , only : masterproc, mpicom, iam
   use shr_mpi_mod      , only : shr_mpi_bcast
   use dshr_strdata_mod , only : shr_strdata_init_from_inline
   use lnd_comp_shr     , only : mesh, model_clock
   use clm_varctl       , only : use_lcz
   use landunit_varcon  , only : isturb_tbd, isturb_hd, isturb_md          
   use landunit_varcon  , only : isturb_lcz1, isturb_lcz2, isturb_lcz3, &
                                 isturb_lcz4, isturb_lcz5, isturb_lcz6, &
                                 isturb_lcz7, isturb_lcz8, isturb_lcz9, &
                                 isturb_lcz10
   !
   ! !ARGUMENTS:
   implicit none
   class(urbanalbtv_type)         :: this
   type(bounds_type), intent(in)  :: bounds
   character(len=*),  intent(in)  :: NLFilename   ! Namelist filename
   ! 
   ! !LOCAL VARIABLES:
   integer            :: n
   integer            :: stream_year_first_urbanalbtvwall            ! first year in urban wall albedo tv stream to use
   integer            :: stream_year_last_urbanalbtvwall             ! last year in urban wall albedo tv stream to use
   integer            :: model_year_align_urbanalbtvwall             ! align stream_year_first_urbanalbtvwall with this model year
   integer            :: nu_nml                                      ! unit for namelist file 
   integer            :: nml_error                                   ! namelist i/o error flag
   character(len=CL)  :: stream_fldFileName_urbanalbtvwall           ! urban wall albedo time-varying streams filename
   character(len=FL)  :: stream_meshfile_urbanalbtvwall              ! urban wall albedo time-varying mesh filename
   character(len=CL)  :: urbanalbtvwallmapalgo = 'nn'                ! mapping alogrithm for urban ac
   character(len=CL)  :: urbanalbtvwall_tintalgo = 'linear'          ! time interpolation alogrithm 
   integer            :: rc                                          ! error code
   character(*), parameter :: subName = "('urbanalbtvwall_init')"
   !-----------------------------------------------------------------------
   namelist /urbanalbtvwall_streams/       &
        stream_year_first_urbanalbtvwall,  &  
        stream_year_last_urbanalbtvwall,   &  
        model_year_align_urbanalbtvwall,   &  
        urbanalbtvwallmapalgo,             &  
        stream_fldFileName_urbanalbtvwall, &     
        stream_meshfile_urbanalbtvwall,    &  
        urbanalbtvwall_tintalgo  
   !-----------------------------------------------------------------------       
   !               
   ! Default values for namelist
   stream_year_first_urbanalbtvwall  = 1      ! first year in stream to use
   stream_year_last_urbanalbtvwall   = 1      ! last  year in stream to use
   model_year_align_urbanalbtvwall   = 1      ! align stream_year_first_urbanalbtvwall with this model year
   stream_fldFileName_urbanalbtvwall = ' '
   stream_meshfile_urbanalbtvwall    = ' '
   
   ! create the field list for urban albedo fields
   if (.not. use_lcz) then 
      stream_var_name_wall(isturb_tbd -6) = "dyn_alb_wall_TBD"
      stream_var_name_wall(isturb_hd -6)  = "dyn_alb_wall_HD"
      stream_var_name_wall(isturb_md -6)  = "dyn_alb_wall_MD"   
   else
      stream_var_name_wall(isturb_lcz1 -6) = "dyn_alb_wall_LCZ1"
      stream_var_name_wall(isturb_lcz2 -6) = "dyn_alb_wall_LCZ2"
      stream_var_name_wall(isturb_lcz3 -6) = "dyn_alb_wall_LCZ3"
      stream_var_name_wall(isturb_lcz4 -6) = "dyn_alb_wall_LCZ4"
      stream_var_name_wall(isturb_lcz5 -6) = "dyn_alb_wall_LCZ5"
      stream_var_name_wall(isturb_lcz6 -6) = "dyn_alb_wall_LCZ6"
      stream_var_name_wall(isturb_lcz7 -6) = "dyn_alb_wall_LCZ7"
      stream_var_name_wall(isturb_lcz8 -6) = "dyn_alb_wall_LCZ8"
      stream_var_name_wall(isturb_lcz9 -6) = "dyn_alb_wall_LCZ9"
      stream_var_name_wall(isturb_lcz10-6) = "dyn_alb_wall_LCZ10"
   end if

   ! Read urbanalbtvwall_streams namelist
   if (masterproc) then
      open( newunit=nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
      call find_nlgroup_name(nu_nml, 'urbanalbtvwall_streams', status=nml_error)
      if (nml_error == 0) then
         read(nu_nml, nml=urbanalbtvwall_streams,iostat=nml_error) 
         if (nml_error /= 0) then
            call endrun(msg='ERROR reading urbanalbtvwall_streams namelist'//errMsg(sourcefile, __LINE__))
         end if
      else
          call endrun(subname // ':: ERROR finding urbanalbtvwall_streams namelist')   
      end if
      close(nu_nml)
   endif

   call shr_mpi_bcast(stream_year_first_urbanalbtvwall  , mpicom)
   call shr_mpi_bcast(stream_year_last_urbanalbtvwall   , mpicom)
   call shr_mpi_bcast(model_year_align_urbanalbtvwall   , mpicom)
   call shr_mpi_bcast(stream_fldFileName_urbanalbtvwall , mpicom)
   call shr_mpi_bcast(stream_meshfile_urbanalbtvwall    , mpicom)
   call shr_mpi_bcast(urbanalbtvwall_tintalgo           , mpicom)
   
   if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'Attempting to read time varying urban wall albedo parameters......'
       write(iulog,'(a)') 'urbanalbtvwall_streams settings:'
       write(iulog,'(a,i8)') '  stream_year_first_urbanalbtvwall  = ',stream_year_first_urbanalbtvwall
       write(iulog,'(a,i8)') '  stream_year_last_urbanalbtvwall   = ',stream_year_last_urbanalbtvwall
       write(iulog,'(a,i8)') '  model_year_align_urbanalbtvwall   = ',model_year_align_urbanalbtvwall
       write(iulog,'(a,a)' ) '  stream_fldFileName_urbanalbtvwall = ',stream_fldFileName_urbanalbtvwall
       write(iulog,'(a,a)' ) '  stream_meshfile_urbanalbtvwall    = ',stream_meshfile_urbanalbtvwall
       write(iulog,'(a,a)' ) '  urbanalbtvwall_tintalgo           = ',urbanalbtvwall_tintalgo
       write(iulog,*) 'Read in urbanalbtvwall_streams namelist from:',trim(NLFilename)
       do n = stream_varname_MIN,stream_varname_MAX
          write(iulog,'(a,a)' ) '  stream_var_name_wall         = ',trim(stream_var_name_wall(n))
       end do
   endif
    
    call shr_strdata_init_from_inline(this%sdat_urbanalbtvwall,                  &
         my_task             = iam,                                              &
         logunit             = iulog,                                            &
         compname            = 'LND',                                            &
         model_clock         = model_clock,                                      &
         model_mesh          = mesh,                                             &
         stream_meshfile     = trim(stream_meshfile_urbanalbtvwall),             &
         stream_lev_dimname  = 'null',                                           &
         stream_mapalgo      = trim(urbanalbtvwallmapalgo),                      &
         stream_filenames    = (/trim(stream_fldfilename_urbanalbtvwall)/),      &
         stream_fldlistFile  = stream_var_name_wall(stream_varname_MIN:stream_varname_MAX), &
         stream_fldListModel = stream_var_name_wall(stream_varname_MIN:stream_varname_MAX), &
         stream_yearFirst    = stream_year_first_urbanalbtvwall,                 &
         stream_yearLast     = stream_year_last_urbanalbtvwall,                  &
         stream_yearAlign    = model_year_align_urbanalbtvwall,                  &
         stream_offset       = 0,                                                &
         stream_taxmode      = 'extend',                                         &
         stream_dtlimit      = 1.0e30_r8,                                        &
         stream_tintalgo     = urbanalbtvwall_tintalgo,                          &
         stream_name         = 'Urban time varying wall albedo data',            &
         rc                  = rc)

   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
   end if
    
  end subroutine urbanalbtvwall_init

  !==============================================================================
 
  subroutine urbanalbtvwall_interp(this, bounds)
    ! !DESCRIPTION:
    ! Interpolate data stream information for urban time varying albedo.
    ! 
    ! !USES:
    use clm_time_manager  , only : get_curr_date
    use clm_instur        , only : urban_valid
    use dshr_methods_mod  , only : dshr_fldbun_getfldptr
    use dshr_strdata_mod  , only : shr_strdata_advance
    use shr_infnan_mod    , only : nan => shr_infnan_nan, assignment(=)
    ! 
    ! !ARGUMENTS:
    ! 
    class(urbanalbtv_type)           :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    !
    logical :: found
    integer :: l, ig, g, ip, n, ib    
    integer :: year    ! year (0, ...) for nstep+1
    integer :: mon     ! month (1, ..., 12) for nstep+1
    integer :: day     ! day of month (1, ..., 31) for nstep+1
    integer :: sec     ! seconds into current date for nstep+1
    integer :: mcdate  ! Current model date (yyyymmdd)
    integer :: lindx   ! landunit index
    integer :: gindx   ! gridcell index
    integer :: lsize
    integer :: rc
    real(r8), pointer :: dataptr1d(:)
    real(r8), pointer :: dataptr2d(:,:)
    ! 
    !-----------------------------------------------------------------------
    ! 
    ! Advance sdat stream
    !
    call get_curr_date(year, mon, day, sec)
    !
    ! packing the date into an integer
    mcdate = year*10000 + mon*100 + day

    call shr_strdata_advance(this%sdat_urbanalbtvwall, ymd=mcdate, tod=sec, logunit=iulog, istr='hdmdyn', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if
    !
    ! Create 2d array for all stream variable data
    lsize = bounds%endg - bounds%begg + 1
    allocate(dataptr2d(lsize, stream_varname_MIN:stream_varname_MAX))
    do n = stream_varname_MIN,stream_varname_MAX
       call dshr_fldbun_getFldPtr(this%sdat_urbanalbtvwall%pstrm(1)%fldbun_model, trim(stream_var_name_wall(n)), &
            fldptr1=dataptr1d, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if
       ! Note that the size of dataptr1d includes ocean points so it will be around 3x larger than lsize
       ! So an explicit loop is required here
       do g = 1,lsize
          dataptr2d(g,n) = dataptr1d(g)
       end do
    end do

    ! Determine this%tbuilding_max (and this%p_ac, if applicable) for all landunits
    do l = bounds%begl,bounds%endl
       if (lun%urbpoi(l)) then
          ! Note that since l is within [begl, endl] bounds, we can assume
          ! lun%gricell(l) is within [begg, endg]
          ig = lun%gridcell(l) - bounds%begg + 1
          do ib = 1,numrad 
             do n = stream_varname_MIN,stream_varname_MAX
                if (stream_var_name_wall((lun%itype(l)-6)) == stream_var_name_wall(n)) then
                   this%dyn_alb_wall_dir(l,ib) = dataptr2d(ig,n)
                   this%dyn_alb_wall_dif(l,ib) = dataptr2d(ig,n)
                end if
             end do
          end do 
       else
           do ib = 1,numrad
              this%dyn_alb_wall_dir(l,ib) = spval
              this%dyn_alb_wall_dif(l,ib) = spval  
           end do  
       end if
    end do
    deallocate(dataptr2d)

    ! Error check
    found = .false.
    do l = bounds%begl,bounds%endl
       if (lun%urbpoi(l)) then
          do g = bounds%begg,bounds%endg
             if (g == lun%gridcell(l)) exit
          end do
          ! Check for valid urban data
          do ib = 1,numrad
             if ( .not. urban_valid(g) .or. (this%dyn_alb_wall_dir(l,ib) <= 0._r8) .or. (this%dyn_alb_wall_dif(l,ib) <= 0._r8)) then
                found = .true.
                gindx = g
                lindx = l
                exit
             end if
          end do   
       end if
    end do
    if ( found ) then
       write(iulog,*)'ERROR: no valid urban data for g= ',gindx
       write(iulog,*)'landunit type:   ',lun%itype(lindx)
       write(iulog,*)'urban_valid:     ',urban_valid(gindx)
       write(iulog,*)'dyn_alb_wall_dir:  ',this%dyn_alb_wall_dir(lindx,:)
       write(iulog,*)'dyn_alb_wall_dif:  ',this%dyn_alb_wall_dif(lindx,:)
       call endrun(subgrid_index=lindx, subgrid_level=subgrid_level_landunit, &
            msg=errmsg(sourcefile, __LINE__))
    end if

  end subroutine urbanalbtvwall_interp

 !---------------------------------------------------------------------
 subroutine urbanalbtvimproad_init(this, bounds, NLFilename)
   !
   ! !DESCRIPTION:
   ! Initialize data stream information for urban time varying impervious road albedo
   !
   ! !USES:
   use clm_nlUtilsMod   , only : find_nlgroup_name
   use spmdMod          , only : masterproc, mpicom, iam
   use shr_mpi_mod      , only : shr_mpi_bcast
   use dshr_strdata_mod , only : shr_strdata_init_from_inline
   use lnd_comp_shr     , only : mesh, model_clock
   use clm_varctl       , only : use_lcz
   use landunit_varcon  , only : isturb_tbd, isturb_hd, isturb_md          
   use landunit_varcon  , only : isturb_lcz1, isturb_lcz2, isturb_lcz3, &
                                 isturb_lcz4, isturb_lcz5, isturb_lcz6, &
                                 isturb_lcz7, isturb_lcz8, isturb_lcz9, &
                                 isturb_lcz10
   !
   ! !ARGUMENTS:
   implicit none
   class(urbanalbtv_type)         :: this
   type(bounds_type), intent(in)  :: bounds
   character(len=*),  intent(in)  :: NLFilename   ! Namelist filename
   ! 
   ! !LOCAL VARIABLES:
   integer            :: n
   integer            :: stream_year_first_urbanalbtvimproad            ! first year in urban impervious road albedo tv stream to use
   integer            :: stream_year_last_urbanalbtvimproad             ! last year in urban impervious road albedo tv stream to use
   integer            :: model_year_align_urbanalbtvimproad             ! align stream_year_first_urbanalbtvimproad with this model year
   integer            :: nu_nml                                      ! unit for namelist file 
   integer            :: nml_error                                   ! namelist i/o error flag
   character(len=CL)  :: stream_fldFileName_urbanalbtvimproad           ! urban impervious road albedo time-varying streams filename
   character(len=FL)  :: stream_meshfile_urbanalbtvimproad              ! urban impervious road albedo time-varying mesh filename
   character(len=CL)  :: urbanalbtvimproadmapalgo = 'nn'                ! mapping alogrithm for urban ac
   character(len=CL)  :: urbanalbtvimproad_tintalgo = 'linear'          ! time interpolation alogrithm 
   integer            :: rc                                          ! error code
   character(*), parameter :: subName = "('urbanalbtvimproad_init')"
   !-----------------------------------------------------------------------
   namelist /urbanalbtvimproad_streams/       &
        stream_year_first_urbanalbtvimproad,  &  
        stream_year_last_urbanalbtvimproad,   &  
        model_year_align_urbanalbtvimproad,   &  
        urbanalbtvimproadmapalgo,             &  
        stream_fldFileName_urbanalbtvimproad, &   
        stream_meshfile_urbanalbtvimproad,    &   
        urbanalbtvimproad_tintalgo  
   !-----------------------------------------------------------------------       
   !               
   ! Default values for namelist
   stream_year_first_urbanalbtvimproad  = 1      ! first year in stream to use
   stream_year_last_urbanalbtvimproad   = 1      ! last  year in stream to use
   model_year_align_urbanalbtvimproad   = 1      ! align stream_year_first_urbanalbtvimproad with this model year
   stream_fldFileName_urbanalbtvimproad = ' '
   stream_meshfile_urbanalbtvimproad    = ' '
   
   ! create the field list for urban albedo fields
   if (.not. use_lcz) then 
      stream_var_name_improad(isturb_tbd -6) = "dyn_alb_improad_TBD"
      stream_var_name_improad(isturb_hd -6)  = "dyn_alb_improad_HD"
      stream_var_name_improad(isturb_md -6)  = "dyn_alb_improad_MD"   
   else
      stream_var_name_improad(isturb_lcz1 -6) = "dyn_alb_improad_LCZ1"
      stream_var_name_improad(isturb_lcz2 -6) = "dyn_alb_improad_LCZ2"
      stream_var_name_improad(isturb_lcz3 -6) = "dyn_alb_improad_LCZ3"
      stream_var_name_improad(isturb_lcz4 -6) = "dyn_alb_improad_LCZ4"
      stream_var_name_improad(isturb_lcz5 -6) = "dyn_alb_improad_LCZ5"
      stream_var_name_improad(isturb_lcz6 -6) = "dyn_alb_improad_LCZ6"
      stream_var_name_improad(isturb_lcz7 -6) = "dyn_alb_improad_LCZ7"
      stream_var_name_improad(isturb_lcz8 -6) = "dyn_alb_improad_LCZ8"
      stream_var_name_improad(isturb_lcz9 -6) = "dyn_alb_improad_LCZ9"
      stream_var_name_improad(isturb_lcz10-6) = "dyn_alb_improad_LCZ10"
   end if

   ! Read urbanalbtvimproad_streams namelist
   if (masterproc) then
      open( newunit=nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
      call find_nlgroup_name(nu_nml, 'urbanalbtvimproad_streams', status=nml_error)
      if (nml_error == 0) then
         read(nu_nml, nml=urbanalbtvimproad_streams,iostat=nml_error) 
         if (nml_error /= 0) then
            call endrun(msg='ERROR reading urbanalbtvimproad_streams namelist'//errMsg(sourcefile, __LINE__))
         end if
      else
          call endrun(subname // ':: ERROR finding urbanalbtvimproad_streams namelist')   
      end if
      close(nu_nml)
   endif

   call shr_mpi_bcast(stream_year_first_urbanalbtvimproad  , mpicom)
   call shr_mpi_bcast(stream_year_last_urbanalbtvimproad   , mpicom)
   call shr_mpi_bcast(model_year_align_urbanalbtvimproad   , mpicom)
   call shr_mpi_bcast(stream_fldFileName_urbanalbtvimproad , mpicom)
   call shr_mpi_bcast(stream_meshfile_urbanalbtvimproad    , mpicom)
   call shr_mpi_bcast(urbanalbtvimproad_tintalgo           , mpicom)
   
   if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'Attempting to read time varying urban impervious road albedo parameters......'
       write(iulog,'(a)') 'urbanalbtvimproad_streams settings:'
       write(iulog,'(a,i8)') '  stream_year_first_urbanalbtvimproad  = ',stream_year_first_urbanalbtvimproad
       write(iulog,'(a,i8)') '  stream_year_last_urbanalbtvimproad   = ',stream_year_last_urbanalbtvimproad
       write(iulog,'(a,i8)') '  model_year_align_urbanalbtvimproad   = ',model_year_align_urbanalbtvimproad
       write(iulog,'(a,a)' ) '  stream_fldFileName_urbanalbtvimproad = ',stream_fldFileName_urbanalbtvimproad
       write(iulog,'(a,a)' ) '  stream_meshfile_urbanalbtvimproad    = ',stream_meshfile_urbanalbtvimproad
       write(iulog,'(a,a)' ) '  urbanalbtvimproad_tintalgo           = ',urbanalbtvimproad_tintalgo
       write(iulog,*) 'Read in urbanalbtvimproad_streams namelist from:',trim(NLFilename)
       do n = stream_varname_MIN,stream_varname_MAX
          write(iulog,'(a,a)' ) '  stream_var_name_improad         = ',trim(stream_var_name_improad(n))
       end do
   endif
    
    call shr_strdata_init_from_inline(this%sdat_urbanalbtvimproad,                  &
         my_task             = iam,                                                 &
         logunit             = iulog,                                               &
         compname            = 'LND',                                               &
         model_clock         = model_clock,                                         &
         model_mesh          = mesh,                                                &
         stream_meshfile     = trim(stream_meshfile_urbanalbtvimproad),             &
         stream_lev_dimname  = 'null',                                              &
         stream_mapalgo      = trim(urbanalbtvimproadmapalgo),                      &
         stream_filenames    = (/trim(stream_fldfilename_urbanalbtvimproad)/),      &
         stream_fldlistFile  = stream_var_name_improad(stream_varname_MIN:stream_varname_MAX), &
         stream_fldListModel = stream_var_name_improad(stream_varname_MIN:stream_varname_MAX), &
         stream_yearFirst    = stream_year_first_urbanalbtvimproad,                 &
         stream_yearLast     = stream_year_last_urbanalbtvimproad,                  &
         stream_yearAlign    = model_year_align_urbanalbtvimproad,                  &
         stream_offset       = 0,                                                   &
         stream_taxmode      = 'extend',                                            &
         stream_dtlimit      = 1.0e30_r8,                                           &
         stream_tintalgo     = urbanalbtvimproad_tintalgo,                          &
         stream_name         = 'Urban time varying impervious road albedo data',               &
         rc                  = rc)

   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
   end if
    
  end subroutine urbanalbtvimproad_init

  !==============================================================================
 
  subroutine urbanalbtvimproad_interp(this, bounds)
    ! !DESCRIPTION:
    ! Interpolate data stream information for urban time varying albedo.
    ! 
    ! !USES:
    use clm_time_manager  , only : get_curr_date
    use clm_instur        , only : urban_valid
    use dshr_methods_mod  , only : dshr_fldbun_getfldptr
    use dshr_strdata_mod  , only : shr_strdata_advance
    use shr_infnan_mod    , only : nan => shr_infnan_nan, assignment(=)
    ! 
    ! !ARGUMENTS:
    ! 
    class(urbanalbtv_type)           :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    !
    logical :: found
    integer :: l, ig, g, ip, n, ib    
    integer :: year    ! year (0, ...) for nstep+1
    integer :: mon     ! month (1, ..., 12) for nstep+1
    integer :: day     ! day of month (1, ..., 31) for nstep+1
    integer :: sec     ! seconds into current date for nstep+1
    integer :: mcdate  ! Current model date (yyyymmdd)
    integer :: lindx   ! landunit index
    integer :: gindx   ! gridcell index
    integer :: lsize
    integer :: rc
    real(r8), pointer :: dataptr1d(:)
    real(r8), pointer :: dataptr2d(:,:)
    ! 
    !-----------------------------------------------------------------------
    ! 
    ! Advance sdat stream
    !
    call get_curr_date(year, mon, day, sec)
    !
    ! packing the date into an integer
    mcdate = year*10000 + mon*100 + day

    call shr_strdata_advance(this%sdat_urbanalbtvimproad, ymd=mcdate, tod=sec, logunit=iulog, istr='hdmdyn', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if
    !
    ! Create 2d array for all stream variable data
    lsize = bounds%endg - bounds%begg + 1
    allocate(dataptr2d(lsize, stream_varname_MIN:stream_varname_MAX))
    do n = stream_varname_MIN,stream_varname_MAX
       call dshr_fldbun_getFldPtr(this%sdat_urbanalbtvimproad%pstrm(1)%fldbun_model, trim(stream_var_name_improad(n)), &
            fldptr1=dataptr1d, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if
       ! Note that the size of dataptr1d includes ocean points so it will be around 3x larger than lsize
       ! So an explicit loop is required here
       do g = 1,lsize
          dataptr2d(g,n) = dataptr1d(g)
       end do
    end do

    ! Determine this%tbuilding_max (and this%p_ac, if applicable) for all landunits
    do l = bounds%begl,bounds%endl
       if (lun%urbpoi(l)) then
          ! Note that since l is within [begl, endl] bounds, we can assume
          ! lun%gricell(l) is within [begg, endg]
          ig = lun%gridcell(l) - bounds%begg + 1
          do ib = 1,numrad 
             do n = stream_varname_MIN,stream_varname_MAX
                if (stream_var_name_improad((lun%itype(l)-6)) == stream_var_name_improad(n)) then
                   this%dyn_alb_improad_dir(l,ib) = dataptr2d(ig,n)
                   this%dyn_alb_improad_dif(l,ib) = dataptr2d(ig,n)
                end if
             end do
          end do 
       else
           do ib = 1,numrad
              this%dyn_alb_improad_dir(l,ib) = spval
              this%dyn_alb_improad_dif(l,ib) = spval  
           end do  
       end if
    end do
    deallocate(dataptr2d)

    ! Error check
    found = .false.
    do l = bounds%begl,bounds%endl
       if (lun%urbpoi(l)) then
          do g = bounds%begg,bounds%endg
             if (g == lun%gridcell(l)) exit
          end do
          ! Check for valid urban data
          do ib = 1,numrad
             if ( .not. urban_valid(g) .or. (this%dyn_alb_improad_dir(l,ib) <= 0._r8) .or. (this%dyn_alb_improad_dif(l,ib) <= 0._r8)) then
                found = .true.
                gindx = g
                lindx = l
                exit
             end if
          end do   
       end if
    end do
    if ( found ) then
       write(iulog,*)'ERROR: no valid urban data for g= ',gindx
       write(iulog,*)'landunit type:   ',lun%itype(lindx)
       write(iulog,*)'urban_valid:     ',urban_valid(gindx)
       write(iulog,*)'dyn_alb_improad_dir:  ',this%dyn_alb_improad_dir(lindx,:)
       write(iulog,*)'dyn_alb_improad_dif:  ',this%dyn_alb_improad_dif(lindx,:)
       call endrun(subgrid_index=lindx, subgrid_level=subgrid_level_landunit, &
            msg=errmsg(sourcefile, __LINE__))
    end if

  end subroutine urbanalbtvimproad_interp

end module UrbanDynAlbMod
