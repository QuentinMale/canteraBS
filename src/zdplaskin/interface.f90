subroutine zdplaskinInit() bind(C, name='zdplaskinInit')
  use ZDPlasKin
  implicit none
  call ZDPlasKin_init()
end subroutine zdplaskinInit

subroutine zdplaskinSetDensity(cstring, DENS) bind(C, name='zdplaskinSetDensity')
  use, intrinsic :: iso_c_binding
  use ZDPlasKin
  implicit none
  type(c_ptr), target, intent(in) :: cstring
  character(kind=c_char), pointer :: fstring(:)
  character(10) :: string
  integer length, i
  real(c_double), intent(in) :: DENS

  interface
   function strlen(s) bind(C, name='strlen')
     use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
     implicit none
     type(c_ptr), intent(in), value :: s
     integer(c_size_t) :: strlen
   end function strlen
  end interface

  length = strlen(cstring)
  call c_f_pointer(cstring, fstring, [length])
  string = fstring(1)

  do i=2,size(fstring(:))
    string = trim(string)//fstring(i)
  enddo

  string = trim(string)
  call ZDPlasKin_set_density(string, DENS)
end subroutine zdplaskinSetDensity

subroutine zdplaskinSetConditions(gasTemperature, reduced_field) bind(C,name='zdplaskinSetConditions')
  use, intrinsic :: iso_c_binding
  use ZDPlasKin
  implicit none
  real(c_double), intent(in) :: gasTemperature
  real(c_double), intent(in) :: reduced_field
  call ZDPlasKin_set_conditions(GAS_TEMPERATURE=gasTemperature, REDUCED_FIELD=reduced_field)
end subroutine zdplaskinSetConditions

function zdplaskinGetElecTemp() bind(C,name='zdplaskinGetElecTemp') result(eTemperature)
  use, intrinsic :: iso_c_binding
  use ZDPlasKin
  implicit none
  real(c_double) :: eTemperature
  call ZDPlasKin_get_conditions(ELEC_TEMPERATURE=eTemperature)
end function zdplaskinGetElecTemp

function zdplaskinGetElecDiffCoeff() bind(C,name='zdplaskinGetElecDiffCoeff') result(De)
  use, intrinsic :: iso_c_binding
  use ZDPlasKin
  implicit none
  real(c_double) :: De
  call ZDPlasKin_get_conditions(ELEC_DIFF_COEFF=De)
  ! convert to SI
  De = 1e-4 * De
end function zdplaskinGetElecDiffCoeff

function zdplaskinGetElecMobility(ND) bind(C,name='zdplaskinGetElecMobility') result(mu)
  use, intrinsic :: iso_c_binding
  use ZDPlasKin
  implicit none
  real(c_double) :: mu
  real(c_double), intent(in) :: ND
  call ZDPlasKin_get_conditions(ELEC_MOBILITY_N=mu)
  ! convert to SI
  mu = 100 * mu / ND
end function zdplaskinGetElecMobility

subroutine zdplaskinGetPlasmaSource(array) bind(C,name='zdplaskinGetPlasmaSource')
  use, intrinsic :: iso_c_binding
  use ZDPlasKin
  implicit none
  type(c_ptr), intent(inout) :: array
  real(c_double), target, save :: source(0:species_max-1)

  call ZDPlasKin_get_rates(SOURCE_TERMS=source)
  ! Allocate an array and make it available in C
  array = c_loc(source)
end subroutine zdplaskinGetPlasmaSource

function zdplaskinNSpecies() bind(C,name='zdplaskinNSpecies') result(nSpecies)
  use, intrinsic :: iso_c_binding
  use ZDPlasKin
  integer(c_size_t) nSpecies
  nSpecies = species_max
end function zdplaskinNSpecies

subroutine zdplaskinGetSpeciesName(cstring, index) bind(C, name='zdplaskinGetSpeciesName')
  use iso_c_binding
  use ZDPlasKin
  implicit none
  integer, intent(in) :: index
  CHARACTER(10), TARGET :: fstring = ''
  TYPE(C_PTR) :: cstring
  fstring = species_name(index+1)
  cstring = c_loc(fstring)
end subroutine zdplaskinGetSpeciesName


