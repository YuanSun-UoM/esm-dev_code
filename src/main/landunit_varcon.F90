module landunit_varcon

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing landunit indices and associated variables and routines.
  !
  ! !USES:
#include "shr_assert.h"
!YS
  use clm_varctl   , only : use_lcz
!YS  
  !
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  
  !------------------------------------------------------------------
  ! Initialize landunit type constants
  !------------------------------------------------------------------

  integer, parameter, public :: istsoil    = 1  !soil         landunit type (natural vegetation)
  integer, parameter, public :: istcrop    = 2  !crop         landunit type
  integer, parameter, public :: istocn     = 3  !ocean        landunit type
  integer, parameter, public :: istice     = 4  !land ice landunit type
  integer, parameter, public :: istdlak    = 5  !deep lake    landunit type (now used for all lakes)
  integer, parameter, public :: istwet     = 6  !wetland      landunit type (swamp, marsh, etc.)

  integer, parameter, public :: isturb_MIN = 7  !minimum urban type index
  integer, parameter, public :: isturb_tbd = 7  !urban tbd    landunit type
  integer, parameter, public :: isturb_hd  = 8  !urban hd     landunit type
  integer, parameter, public :: isturb_md  = 9  !urban md     landunit type
!YS  integer, parameter, public :: isturb_MAX = 9  !maximum urban type index
!YS
  ! 10 lCZs urban landunits  
  integer, parameter, public :: isturb_lcz1  = 7     !LCZ 1      urban landunit type
  integer, parameter, public :: isturb_lcz2  = 8     !LCZ 2      urban landunit type
  integer, parameter, public :: isturb_lcz3  = 9     !LCZ 3      urban landunit type
  integer, parameter, public :: isturb_lcz4  = 10    !LCZ 4      urban landunit type
  integer, parameter, public :: isturb_lcz5  = 11    !LCZ 5      urban landunit type
  integer, parameter, public :: isturb_lcz6  = 12    !LCZ 6      urban landunit type
  integer, parameter, public :: isturb_lcz7  = 13    !LCZ 7      urban landunit type
  integer, parameter, public :: isturb_lcz8  = 14    !LCZ 8      urban landunit type
  integer, parameter, public :: isturb_lcz9  = 15    !LCZ 9      urban landunit type
  integer, parameter, public :: isturb_lcz10 = 16    !LCZ 10     urban landunit type
  integer, public :: max_lunit   !maximum value that lun%itype can have
  integer, public :: isturb_MAX  !maximum urban type index
  integer, public :: numurbl
!YS  
!YS   integer, parameter, public :: max_lunit  = 9  !maximum value that lun%itype can have
!YS                                         !(i.e., largest value in the above list)

  integer, parameter, public                   :: landunit_name_length = 40  ! max length of landunit names
!YS  character(len=landunit_name_length), public  :: landunit_names(max_lunit)  ! name of each landunit type
  character(len=landunit_name_length), allocatable, public  :: landunit_names(:)  ! name of each landunit type 
  ! parameters that depend on the above constants

!YS  integer, parameter, public :: numurbl = isturb_MAX - isturb_MIN + 1   ! number of urban landunits

  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: landunit_varcon_init  ! initialize constants in this module
  public :: landunit_is_special   ! returns true if this is a special landunit

  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: set_landunit_names   ! set the landunit_names vector
!-----------------------------------------------------------------------

contains
  
  !-----------------------------------------------------------------------
  subroutine landunit_varcon_init()
    !
    ! !DESCRIPTION:
    ! Initialize constants in landunit_varcon
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'landunit_varcon_init'
    !-----------------------------------------------------------------------
!YS    
    ! parameters that depend on the above constants

    if (.not. use_lcz) then
       max_lunit = 9
       isturb_MAX = 9
    else if (use_lcz) then
       max_lunit = 16
       isturb_MAX = 16
    end if 

    numurbl = isturb_MAX - isturb_MIN + 1 
    allocate(landunit_names(max_lunit))
!YS
    call set_landunit_names()

  end subroutine landunit_varcon_init

  !-----------------------------------------------------------------------
  function landunit_is_special(ltype) result(is_special)
    !
    ! !DESCRIPTION:
    ! Returns true if the landunit type ltype is a special landunit; returns false otherwise
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    logical :: is_special  ! function result
    integer :: ltype       ! landunit type of interest
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'landunit_is_special'
    !-----------------------------------------------------------------------

    SHR_ASSERT((ltype >= 1 .and. ltype <= max_lunit), subname//': ltype out of bounds')

    if (ltype == istsoil .or. ltype == istcrop) then
       is_special = .false.
    else
       is_special = .true.
    end if

  end function landunit_is_special

  
  !-----------------------------------------------------------------------
  subroutine set_landunit_names
    !
    ! !DESCRIPTION:
    ! Set the landunit_names vector
    !
    ! !USES:
    use shr_sys_mod, only : shr_sys_abort
    !
    character(len=*), parameter :: not_set = 'NOT_SET'
    character(len=*), parameter :: unused  = 'UNUSED'
    character(len=*), parameter :: subname = 'set_landunit_names'
    !-----------------------------------------------------------------------
    
    landunit_names(:) = not_set

    landunit_names(istsoil) = 'vegetated_or_bare_soil'
    landunit_names(istcrop) = 'crop'
    landunit_names(istocn) = 'ocean'
    landunit_names(istice) = 'landice'
    landunit_names(istdlak) = 'deep_lake'
    landunit_names(istwet) = 'wetland'
!YS    landunit_names(isturb_tbd) = 'urban_tbd'
!YS    landunit_names(isturb_hd) = 'urban_hd'
!YS    landunit_names(isturb_md) = 'urban_md'
!YS    
    if (.not. use_lcz) then
       landunit_names(isturb_tbd) = 'urban_tbd'
       landunit_names(isturb_hd) = 'urban_hd'
       landunit_names(isturb_md) = 'urban_md'
    else if (use_lcz) then
       landunit_names(isturb_lcz1)  = 'urban_lcz1'
       landunit_names(isturb_lcz2)  = 'urban_lcz2'
       landunit_names(isturb_lcz3)  = 'urban_lcz3'
       landunit_names(isturb_lcz4)  = 'urban_lcz4'
       landunit_names(isturb_lcz5)  = 'urban_lcz5'
       landunit_names(isturb_lcz6)  = 'urban_lcz6'
       landunit_names(isturb_lcz7)  = 'urban_lcz7'
       landunit_names(isturb_lcz8)  = 'urban_lcz8'
       landunit_names(isturb_lcz9)  = 'urban_lcz9'
       landunit_names(isturb_lcz10) = 'urban_lcz10'
       !landunit_names(isturb_lcz11) = 'urban_lcz11'
    end if
!YS
    if (any(landunit_names == not_set)) then
       call shr_sys_abort(trim(subname)//': Not all landunit names set')
    end if

  end subroutine set_landunit_names

end module landunit_varcon
