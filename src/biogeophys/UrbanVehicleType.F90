module UrbanVehicleType
!----------------------------------------------------------------------- 
! DESCRIPTION:
! Urban Traffic input data
! USES:
 use ESMF             , only : ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU, ESMF_Finalize, ESMF_END_ABORT
 use shr_kind_mod     , only : r8 => shr_kind_r8, CL => shr_kind_CL
 use shr_log_mod      , only : errMsg => shr_log_errMsg
 use abortutils       , only : endrun
 use LandunitType     , only : lun
 use dshr_strdata_mod , only : shr_strdata_type
 use clm_varctl       , only : iulog
 use decompMod        , only : bounds_type, subgrid_level_landunit
 use clm_varcon       , only : spval
 use UrbanParamsType  , only : urbanparams_type, urban_traffic
 use WaterType        , only : water_type
 implicit none
 private
 ! PUBLIC MEMBER
    public  :: traffic_flux      ! Calculate urban traffic heat flux
    public  :: vehicle_speed_env ! Calculate vehicle speed with environmental factors
    public  :: ev_heat_scale     ! Calculate electric vehicle scaler
    public  :: vehicle_flow_env  ! Calculate vehicle flow in impervious road
 ! ! PUBLIC TYPE
 type, public :: urbanvehicle_type
    !
    real(r8), pointer         :: vehicle_flow            (:)   ! traffic flow input
    real(r8), pointer         :: heat_per_vehicle        (:)   ! heat emission per vehicle (J/m)
    real(r8), pointer         :: ev_scaler_grc           (:)   ! electric vehicle heat scaler
    real(r8), pointer         :: vehicle_speed_grc       (:)   ! vehicle speed considering environmental factors
    real(r8), pointer         :: vehicle_flow_lun        (:)   ! vehicle flow in impervious road
    type(shr_strdata_type)    :: sdat_vehicletv                ! annual urban time varying traffic data stream
  contains
    !! PUBLIC MEMBER FUNCTIONS
    procedure, public  :: Init              ! Allocate and initialize urbantraffic
    procedure, public  :: vehicletv_init    ! Initialize urban time varying vehicle counts
    procedure, public  :: vehicletv_interp  ! Interpolate urban time varying vehicle counts
    procedure, private :: InitHistory       ! setup history fields
  end type urbanvehicle_type

  integer          , private              :: stream_varname_MIN       ! minimum index for stream_varnames
  integer          , private              :: stream_varname_MAX       ! maximum index for stream_varnames
  character(25)    , private, pointer     :: stream_varnames(:)       ! urban time varying variable names
  character(len=*) , parameter, private   :: sourcefile = &
       __FILE__
!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------
  subroutine Init(this, bounds, NLFilename)
    !
    ! ! USES:
    use shr_infnan_mod   , only : nan => shr_infnan_nan, assignment(=)
    !
    ! ! ARGUMENTS:
    class(urbanvehicle_type)                   :: this
    type(bounds_type)        , intent(in)      :: bounds
    character(len=*)         , intent(in)      :: NLFilename         ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer		:: begl, endl
    integer		:: begg, endg
    character(*), parameter :: subName = "('Init')"
    !---------------------------------------------------------------------
   
    begl = bounds%begl; endl = bounds%endl
    begg = bounds%begg; endg = bounds%endg
    ! allocate  
    allocate(this%vehicle_speed_grc(begg:endg))                    ; this%vehicle_speed_grc  (:)  = nan
    allocate(this%vehicle_flow_lun(begl:endl))                     ; this%vehicle_flow_lun   (:)  = nan
    allocate(this%ev_scaler_grc(begg:endg))                        ; this%ev_scaler_grc      (:)  = 1.0
    
    if (urban_traffic) then 
        ! Determine the minimum and maximum indices for stream_varnames
        stream_varname_MIN = 1
        stream_varname_MAX = 7

        allocate(stream_varnames(stream_varname_MIN:stream_varname_MAX))
        allocate(this%heat_per_vehicle(begg:endg))                 ; this%heat_per_vehicle   (:)  = 0._r8    ! should be 0._r8 rather than nan               
        allocate(this%vehicle_flow(begl:endl))                     ; this%vehicle_flow       (:)  = nan 
       
        call this%vehicletv_init(bounds, NLFilename)
        call this%vehicletv_interp(bounds)
    end if

    call this%InitHistory(bounds)    

  end subroutine Init

  !==============================================================================
  subroutine vehicletv_init(this, bounds, NLFilename)
    !
    ! ! DESCRIPTION:
    ! Initialize data stream information for urban time varying traffic data
    !
    ! ! USES:
    use clm_nlUtilsMod   , only : find_nlgroup_name
    use spmdMod          , only : masterproc, mpicom, iam
    use shr_mpi_mod      , only : shr_mpi_bcast
    use landunit_varcon  , only : isturb_tbd, isturb_hd, isturb_md
    use dshr_strdata_mod , only : shr_strdata_init_from_inline
    use lnd_comp_shr     , only : mesh, model_clock
    !
    ! ! ARGUMENTS:
    implicit none
    class(urbanvehicle_type)       :: this
    type(bounds_type) , intent(in) :: bounds
    character(len=*)  , intent(in) :: NLFilename   ! Namelist filename
    !
    ! ! LOCAL VARIABLES:
    integer            :: n
    integer            :: stream_year_first_vehicletv       ! first year in urban vehicle tv stream to use
    integer            :: stream_year_last_vehicletv        ! last year in urban vehicle tv stream to use
    integer            :: model_year_align_vehicletv        ! align stream_year_first_vehicletv with this model year  
    integer            :: nu_nml                            ! unit for namelist file
    integer            :: nml_error                         ! namelist i/o error flag
    character(len=CL)  :: stream_fldFileName_vehicletv      ! urban vehicle tv streams filename
    character(len=CL)  :: stream_meshfile_vehicletv         ! urban vehicle tv streams filename
    character(len=CL)  :: vehicletvmapalgo = 'nn'           ! mapping alogrithm for urban vehicle
    character(len=CL)  :: vehicletv_tintalgo = 'linear'     ! time interpolation alogrithm
    integer            :: rc                                ! error code
    character(*), parameter :: subName = "('vehicletv_init')"
    !-----------------------------------------------------------------------

    namelist /vehicletv_streams/       &
         stream_year_first_vehicletv,  &
         stream_year_last_vehicletv,   &
         model_year_align_vehicletv,   &
         vehicletvmapalgo,             &
         stream_fldFileName_vehicletv, &
         stream_meshfile_vehicletv,    &
         vehicletv_tintalgo

    ! Default values for namelist
    stream_year_first_vehicletv  = 1                  ! first year in stream to use
    stream_year_last_vehicletv   = 1                  ! last  year in stream to use
    model_year_align_vehicletv   = 1                  ! align stream_year_first_vehicletv with this model year
    stream_fldFileName_vehicletv = ' '
    stream_meshfile_vehicletv    = ' '
    stream_varnames(1) = "vehicle_flow_TBD"           ! vehicle flow in TBD (unit: vehicles per lane per hour)
    stream_varnames(2) = "vehicle_flow_HD"            ! vehicle flow in HD (unit: vehicles per lane per hour)
    stream_varnames(3) = "vehicle_flow_MD"            ! vehicle flow in MD (unit: vehicles per lane per hour)
    !
    ! considering electric vehicles in future
    stream_varnames(4) = "vehicle_percent_PETROL"     ! percentage of vehicles using petrol (unitless: 0-1)
    stream_varnames(5) = "vehicle_percent_DIESEL"     ! percentage of vehicles using diesel (unitless: 0-1)
    stream_varnames(6) = "vehicle_percent_ELECTRIC"   ! percentage of vehicles using electric (unitless: 0-1)
    stream_varnames(7) = "vehicle_percent_HYBRID"     ! percentage of vehicles using hybrid (unitless: 0-1)
    
    ! Read vehicletv_streams namelist
    if (masterproc) then
       open(newunit=nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
       call find_nlgroup_name(nu_nml, 'vehicletv_streams', status=nml_error)
       if (nml_error == 0) then
          read(nu_nml, nml=vehicletv_streams,iostat=nml_error)
          if (nml_error /= 0) then
             call endrun(msg='ERROR reading vehicletv_streams namelist'//errMsg(sourcefile, __LINE__))
          end if
       end if
       close(nu_nml)
    endif

    call shr_mpi_bcast(stream_year_first_vehicletv  , mpicom)
    call shr_mpi_bcast(stream_year_last_vehicletv   , mpicom)
    call shr_mpi_bcast(model_year_align_vehicletv   , mpicom)
    call shr_mpi_bcast(stream_fldFileName_vehicletv , mpicom)
    call shr_mpi_bcast(stream_meshfile_vehicletv    , mpicom)
    call shr_mpi_bcast(vehicletv_tintalgo           , mpicom)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,'(a)')    '  vehicletv_streams settings:    '
       write(iulog,'(a,i8)') '  stream_year_first_vehicletv  = ',stream_year_first_vehicletv
       write(iulog,'(a,i8)') '  stream_year_last_vehicletv   = ',stream_year_last_vehicletv
       write(iulog,'(a,i8)') '  model_year_align_vehicletv   = ',model_year_align_vehicletv
       write(iulog,'(a,a)' ) '  stream_fldFileName_vehicletv = ',stream_fldFileName_vehicletv
       write(iulog,'(a,a)' ) '  stream_meshfile_vehicletv    = ',stream_meshfile_vehicletv
       write(iulog,'(a,a)' ) '  vehicletv_tintalgo           = ',vehicletv_tintalgo
       do n = stream_varname_MIN,stream_varname_MAX
          write(iulog,'(a,a)' ) '  stream_varname            = ',trim(stream_varnames(n))
       end do
       write(iulog,*) ' '
    endif

    ! Initialize the cdeps data type this%sdat_vehicletv
    call shr_strdata_init_from_inline(this%sdat_vehicletv,             &
         my_task             = iam,                                    &
         logunit             = iulog,                                  &
         compname            = 'LND',                                  &
         model_clock         = model_clock,                            &
         model_mesh          = mesh,                                   &
         stream_meshfile     = trim(stream_meshfile_vehicletv),        &
         stream_lev_dimname  = 'null',                                 &
         stream_mapalgo      = trim(vehicletvmapalgo),                 &
         stream_filenames    = (/trim(stream_fldfilename_vehicletv)/), &
         stream_fldlistFile  = stream_varnames(stream_varname_MIN:stream_varname_MAX), &
         stream_fldListModel = stream_varnames(stream_varname_MIN:stream_varname_MAX), &
         stream_yearFirst    = stream_year_first_vehicletv,            &
         stream_yearLast     = stream_year_last_vehicletv,             &
         stream_yearAlign    = model_year_align_vehicletv,             &
         stream_offset       = 0,                                      &
         stream_taxmode      = 'extend',                               &
         stream_dtlimit      = 1.0e30_r8,                              &
         stream_tintalgo     = vehicletv_tintalgo,                     &
         stream_name         = 'Urban time varying traffic data',      &
         rc                  = rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

  end subroutine vehicletv_init

  !==============================================================================
  subroutine vehicletv_interp(this, bounds)
    !
    ! ! DESCRIPTION:
    ! Interpolate data stream information for urban time varying vehicle data.
    !
    ! ! USES:
    use clm_time_manager , only : get_curr_date
    use clm_instur       , only : urban_valid
    use dshr_methods_mod , only : dshr_fldbun_getfldptr
    use dshr_strdata_mod , only : shr_strdata_advance
    use shr_infnan_mod   , only : nan => shr_infnan_nan, assignment(=)
    use landunit_varcon  , only : isturb_MIN, isturb_MAX
    use GridcellType     , only : grc
    use clm_time_manager , only : get_local_time
    ! isturb_MIN = 7
    ! isturb_MAX = 9
    !
    ! ! ARGUMENTS:
    class(urbanvehicle_type)           :: this
    type(bounds_type), intent(in)      :: bounds
    !
    ! ! LOCAL VARIABLES:
    logical :: found
    integer :: l, ig, g, il, n
    integer :: year    ! year (0, ...) for nstep+1
    integer :: mon     ! month (1, ..., 12) for nstep+1
    integer :: day     ! day of month (1, ..., 31) for nstep+1
    integer :: sec     ! seconds into current date for nstep+1
    integer :: hour    ! hour of the day
    integer :: mcdate  ! Current model date (yyyymmdd)
    integer :: lindx   ! landunit index
    integer :: gindx   ! gridcell index
    integer :: gsize
    integer :: rc
    integer :: begl, endl
    integer :: begg, endg
    integer :: num_urbl
    real(r8), pointer                     :: scaler(:)      
    real(r8), pointer                     :: dataptr1d(:)
    real(r8), pointer                     :: dataptr2d(:,:)
    real(r8), dimension(4)                :: heat_release                       ! Watt
    real, dimension(24,3)                 :: diurnal_factors                    ! diurnal factors data (multiple dimensions)
    real(r8), pointer                     :: vehicle_percent(:)                 ! percentage of vehicles
    real(r8), pointer                     :: heat_sum(:)                        ! vehicle heat at grid-cell level
    character(*), parameter               :: subName = "('vehicletv_interp')"
    !-----------------------------------------------------------------------
    
    begl = bounds%begl; endl = bounds%endl
    begg = bounds%begg; endg = bounds%endg
    num_urbl = 1 + isturb_MAX - isturb_MIN
    !
    ! Check time zone if necessary
    ! Initialize the diurnal factors
    data diurnal_factors /  &
    ! Tall building district (TBD)
    0.6, 0.4, 0.3, 0.3, 0.5, 0.8, 1.2, 1.5, 1.8, 1.6, 1.4, 1.3, &
    1.2, 1.3, 1.4, 1.5, 1.7, 1.9, 2.0, 1.8, 1.4, 1.2, 0.9, 0.7,  &
    ! High-density area (HD)
    0.5, 0.3, 0.2, 0.2, 0.4, 0.7, 1.1, 1.3, 1.5, 1.4, 1.3, 1.2, &
    1.1, 1.2, 1.3, 1.4, 1.6, 1.8, 1.9, 1.7, 1.3, 1.1, 0.8, 0.6,  &
    ! Medium-density area (MD)
    0.04, 0.02, 0.02, 0.01, 0.01, 0.01, 0.02, 0.04, 0.08, 0.07, 0.06, 0.06, &
    0.06, 0.05, 0.06, 0.06, 0.05, 0.06, 0.05, 0.05, 0.05, 0.05, 0.04, 0.04 /    
    
    data heat_release / 18900, 19340, 5240, 670 /

    ! Advance sdat stream
    call get_curr_date(year, mon, day, sec)
    mcdate = year*10000 + mon*100 + day

    call shr_strdata_advance(this%sdat_vehicletv, ymd=mcdate, tod=sec, logunit=iulog, istr='hdmdyn', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    ! Create 2d array for all stream variable data
    ! number of grids
    gsize = endg - begg + 1 

    ! allocate 
    allocate(vehicle_percent(begg:endg)) ; vehicle_percent(:)  = nan
    allocate(heat_sum(begg:endg)) ; heat_sum(:) = 0._r8
    allocate(dataptr2d(gsize, stream_varname_MIN:stream_varname_MAX))
    allocate(scaler(begg:endg)); scaler(:) = 1.0
    !
    do n = stream_varname_MIN, stream_varname_MAX
       call dshr_fldbun_getFldPtr(this%sdat_vehicletv%pstrm(1)%fldbun_model, trim(stream_varnames(n)), &
            fldptr1=dataptr1d, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if
       do g = 1, gsize
          dataptr2d(g,n) = dataptr1d(g)
       end do
    end do

    ! Determine 2D data for all landunits
    do l = begl, endl
       ! Note that since l is within [begl, endl] bounds, we can assume
       ! lun%gricell(l) is within [begg, endg]
       if (lun%urbpoi(l)) then
          ig=lun%gridcell(l)-begg+1
          g=lun%gridcell(l)
          hour=get_local_time(grc%londeg(g)) / 3600 + 1
          hour = max(1, min(24, hour))
          do n = stream_varname_MIN, num_urbl
             ! adjust vehicle flow with diurnal variation
             ! the if condition is necessary!
             if (stream_varnames((lun%itype(l)-isturb_MIN+1)) == stream_varnames(n)) then
                this%vehicle_flow(l) = dataptr2d(ig,n) * diurnal_factors(hour, n) 
                !write(iulog,*)'g,ig:            ',g,ig
                !write(iulog,*)'n:              ',n
                !write(iulog,*)'hour:           ',hour
                !write(iulog,*)'l, lun%itype(l):',l, lun%itype(l)
                !write(iulog,*)'dataptr2d(ig,n):',dataptr2d(ig,n) 
                !write(iulog,*)'diurnal_factors:',diurnal_factors(hour,n) 
                !write(iulog,*)'lun%itype(l)-isturb_MIN+1:',lun%itype(l)-isturb_MIN+1
                !write(iulog,*)'vehicle_flow:   ',this%vehicle_flow(l)   
             end if
          end do  
        end if
    end do
    
    ! Determine 2D data at grid-cell level
    do g = begg, endg
       do n = (num_urbl+1), stream_varname_MAX
           vehicle_percent(g)=dataptr2d(g,n)
           if (n == stream_varname_MAX) then 
              scaler = this%ev_scaler_grc(g)
           else
              scaler(g) = 1.0
           end if       
           heat_sum(g)=heat_sum(g)+vehicle_percent(g)*heat_release(n-num_urbl)*scaler(g)
           !write(iulog,*)'scaler:   ',scaler(g)
           !write(iulog,*)'g:        ',g
       end do
       this%heat_per_vehicle(g) = heat_sum(g)
       !write(iulog,*)'scaler:   ',scaler          
    end do

    deallocate(dataptr2d)
    deallocate(vehicle_percent)
    deallocate(heat_sum)
    deallocate(scaler)

    ! Error check
    found = .false.
    do l = begl, endl
       if (lun%urbpoi(l)) then
          do g = begg, endg
             if (g == lun%gridcell(l)) exit
          end do
          ! Check for valid urban data
          if ((.not. urban_valid(g)) .or. (this%heat_per_vehicle(g) <= 0._r8) &
             .or. (this%vehicle_flow(l) <= 0._r8)) then   
             found = .true.
             gindx = g
             lindx = l
             exit
          end if   
       end if
    end do

    if ( found ) then
       write(iulog,*)'ERROR: no valid urban traffic flux data for g= ',gindx
       write(iulog,*)'landunit type:   ',lun%itype(lindx)
       write(iulog,*)'urban_valid:     ',urban_valid(gindx)
       write(iulog,*)'heat_per_vehicle ',this%heat_per_vehicle(gindx)
       write(iulog,*)'vehicle_flow:    ',this%vehicle_flow(lindx)
       call endrun(subgrid_index=lindx, subgrid_level=subgrid_level_landunit, &
            msg=errmsg(sourcefile, __LINE__))
    end if
  end subroutine vehicletv_interp
   
   !==============================================================================
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Setup fields that can be output to history files
    !
    ! !USES:
    use shr_infnan_mod   , only : nan => shr_infnan_nan, assignment(=)
    use histFileMod      , only : hist_addfld1d
    implicit none
    !
    ! !ARGUMENTS:
    class(urbanvehicle_type)        :: this
    type(bounds_type), intent(in)   :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer                 :: begl, endl
    integer                 :: begg, endg
    character(*), parameter :: subName = "('InitHistory')"
    !------------------------------------------------------------------------
    begl = bounds%begl; endl = bounds%endl
    begg = bounds%begg; endg = bounds%endg
    this%vehicle_flow_lun(begl:endl) = spval
    call hist_addfld1d (fname='VEHICLE_FLOW', units='vehicle per hour',  &
         avgflag='A', long_name='vehicle flow', &
         ptr_lunit=this%vehicle_flow_lun, set_nourb=0._r8, &
         default='inactive', l2g_scale_type='unity')
    
    this%vehicle_speed_grc(begg:endg) = spval
    call hist_addfld1d (fname='VEHICLE_SPEED', units='m/s',  &
         avgflag='A', long_name='vehicle speed', &
         ptr_lnd=this%vehicle_speed_grc, set_nourb=0._r8, &
         default='inactive')   

  end subroutine InitHistory  

   !==============================================================================
  subroutine vehicle_speed_env(forc_rain, forc_snow, vehicle_speed_out)
    !
    ! ! DESCRIPTION:
    ! Calculate traffic heat flux
    !
    ! ! USES:
    use shr_kind_mod , only: r8 => shr_kind_r8
    implicit none
    ! ! ARGUMENTS:
    real(r8), intent(in)        :: forc_rain                   ! unit: mm/s
    real(r8), intent(in)        :: forc_snow                   ! unit: mm/s
    real(r8), intent(out)       :: vehicle_speed_out
    !
    ! ! LOCAL:
    real(r8)                    :: rain_scale_factor
    real(r8)                    :: snow_scale_factor
    real(r8), target            :: vehicle_speed = 11.1        ! unit: m/s
    real(r8), target            :: rain_threshold = 0.00083    ! unit: mm/s
    real(r8), target            :: snow_threshold = 0.000353   ! unit: mm/s
    !-----------------------------------------------------------------------
    !
    ! secondary weather impact
    ! https://doi.org/10.1260/2046-0430.1.1.25: 
    ! 5% reduction over 0.3 cm/h = 0.00083 mm/s 
    ! 8% 
       
    if ((forc_rain > 0._r8) .and. (forc_rain <= rain_threshold)) then
        rain_scale_factor = 1.0 - 60 * forc_rain
    elseif (forc_rain > rain_threshold) then
        rain_scale_factor = 1.0 - (90 * forc_rain + 0.0425)
    elseif (forc_rain .eq. 0._r8) then
        rain_scale_factor = 1.0 
    end if 
    
    rain_scale_factor = max(rain_scale_factor, 0._r8)

    ! https://doi.org/10.1080/01441647.2017.1293188
    ! 
    if ((forc_snow > 0._r8) .and. (forc_snow <= snow_threshold)) then
        snow_scale_factor = 0.96
    elseif ((forc_snow > snow_threshold) .and. (forc_snow <= snow_threshold *2)) then
        snow_scale_factor = 0.92
    elseif ((forc_snow > snow_threshold *2) .and. (forc_snow <= snow_threshold *10)) then
        snow_scale_factor = 0.91
    elseif (forc_snow > snow_threshold*10) then
        snow_scale_factor = 0.87
    elseif (forc_snow .eq. 0._r8) then 
        snow_scale_factor = 1.0            
    end if
    snow_scale_factor = max(snow_scale_factor, 0._r8)

    vehicle_speed_out = rain_scale_factor * snow_scale_factor * vehicle_speed

  end subroutine vehicle_speed_env

   !==============================================================================
  subroutine vehicle_flow_env(vehicle_flow_in, nlane_traffic, vehicle_flow_out)
    !
    ! ! DESCRIPTION:
    ! Calculate traffic heat flux
    !
    ! ! USES:
    use shr_kind_mod , only: r8 => shr_kind_r8
    implicit none
    ! ! ARGUMENTS:
    real(r8), intent(in)        :: vehicle_flow_in
    Integer , intent(in)        :: nlane_traffic
    real(r8), intent(out)       :: vehicle_flow_out
    !
    ! ! LOCAL:
    !-----------------------------------------------------------------------
    !
    vehicle_flow_out = vehicle_flow_in * nlane_traffic

  end subroutine vehicle_flow_env

   !==============================================================================
  subroutine ev_heat_scale(forc_t, ev_scaler_out)
    !
    ! ! DESCRIPTION:
    ! Calculate electric vehicle heat scaler
    !
    ! ! USES:
    use shr_kind_mod , only: r8 => shr_kind_r8
    implicit none
    ! ! ARGUMENTS:
    real(r8), intent(in)        :: forc_t
    real(r8), intent(out)       :: ev_scaler_out
    !
    ! ! LOCAL:
    real(r8), target            :: t_threshold = 293.15           ! unit: K
    real(r8), target            :: t_threshold_zero = 273.15      ! unit: K
    !-----------------------------------------------------------------------
    !
    if ((forc_t > t_threshold_zero) .and. (forc_t < t_threshold)) then
        ev_scaler_out = 1.0 + 0.0165 * (t_threshold - forc_t)
    elseif ((forc_t <= t_threshold_zero) .and. (forc_t >263.15)) then
        ev_scaler_out =1.33
    elseif ((forc_t <= 263.15) .and. (forc_t >253.15)) then
        ev_scaler_out = 1.4
    elseif (forc_t <= 253.15) then
        ev_scaler_out = 1.55
    else
        ev_scaler_out = 1.0   
    end if    
  end subroutine ev_heat_scale

   !==============================================================================
  subroutine traffic_flux(total_vehicle_flow, heat_per_vehicle_grid, vehicle_speed_out, improad_width, eflx_traffic)
    !
    ! ! DESCRIPTION:
    ! Calculate traffic heat flux
    ! ! USES:
    use shr_kind_mod , only: r8 => shr_kind_r8

    implicit none
    ! ! ARGUMENTS:
    real(r8), intent(in)        :: total_vehicle_flow
    real(r8), intent(in)        :: heat_per_vehicle_grid
    real(r8), intent(in)        :: vehicle_speed_out
    real(r8), intent(in)        :: improad_width
    real(r8), intent(out)       :: eflx_traffic
    !
    ! ! LOCAL:
    real(r8), target            :: minimum_vehicle_speed = 0.3        ! unit: m/s
    !-----------------------------------------------------------------------
    
    if (vehicle_speed_out >= minimum_vehicle_speed) then
        !write(iulog,*) 'total_vehicle_flow', total_vehicle_flow
        !write(iulog,*) 'heat_per_vehicle_grid', heat_per_vehicle_grid
        eflx_traffic = (total_vehicle_flow / 3600._r8 * heat_per_vehicle_grid) / (improad_width * vehicle_speed_out)
    else
        eflx_traffic = 0._r8
    end if    
  end subroutine traffic_flux
end module UrbanVehicleType