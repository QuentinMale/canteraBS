subroutine zdplaskinInit() bind(C, name='zdplaskinInit')
  use ZDPlasKin
  implicit none
  call ZDPlasKin_init()
end subroutine zdplaskinInit

subroutine zdplaskinSetDensity(cstring, DENS) bind(C, name='zdplaskinSetDensity')
  use C_interface_module
  use ZDPlasKin
  implicit none
  TYPE(C_PTR), intent(inout) :: cstring
  character(10) :: fstring
  real(c_double), intent(in) :: DENS
  call C_F_string_ptr(cstring, fstring)
  call ZDPlasKin_set_density(trim(fstring), DENS)
end subroutine zdplaskinSetDensity

subroutine zdplaskinGetDensity(cstring, DENS) bind(C, name='zdplaskinGetDensity')
  use C_interface_module
  use ZDPlasKin
  implicit none
  TYPE(C_PTR), intent(inout) :: cstring
  character(10) :: fstring
  real(c_double), intent(out) :: DENS
  call C_F_string_ptr(cstring, fstring)
  call ZDPlasKin_get_density(trim(fstring),DENS)
end subroutine zdplaskinGetDensity

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

subroutine zdplaskinGetPlasmaSource(carray) bind(C,name='zdplaskinGetPlasmaSource')
  use, intrinsic :: iso_c_binding
  use ZDPlasKin
  implicit none
  type(c_ptr), intent(inout) :: carray
  real(c_double), target :: farray(species_max)
  integer i
  call ZDPlasKin_get_rates(SOURCE_TERMS=farray)
  do i = 1, species_max
    farray(i) = farray(i) * 1e6
  enddo
  carray = c_loc(farray)
end subroutine zdplaskinGetPlasmaSource

function zdplaskinNSpecies() bind(C,name='zdplaskinNSpecies') result(nSpecies)
  use, intrinsic :: iso_c_binding
  use ZDPlasKin
  integer(c_size_t) nSpecies
  nSpecies = species_max
end function zdplaskinNSpecies

subroutine zdplaskinGetSpeciesName(cstring, index) bind(C, name='zdplaskinGetSpeciesName')
  use ZDPlasKin
  use C_interface_module
  implicit none
  integer, intent(in) :: index
  TYPE(C_PTR), intent(inout) :: cstring
  call F_C_string_ptr(trim(species_name(index+1)), cstring)
end subroutine zdplaskinGetSpeciesName
