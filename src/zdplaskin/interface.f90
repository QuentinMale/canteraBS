subroutine zdplaskinInit() bind(C, name='zdplaskinInit')
  use ZDPlasKin
  implicit none
  call ZDPlasKin_init()
end subroutine zdplaskinInit

subroutine zdplaskinSetConfig(atol,rtol) bind(C, name='zdplaskinSetConfig')
    use C_interface_module
    use ZDPlasKin
    implicit none
    real(c_double), intent(in) :: atol, rtol
    call ZDPlasKin_set_config(ATOL=atol,RTOL=rtol)
end subroutine zdplaskinSetConfig

subroutine zdplaskinSetDensity(cstring, num_density_SI) bind(C, name='zdplaskinSetDensity')
  use C_interface_module
  use ZDPlasKin
  implicit none
  character(len=1,kind=C_char), intent(in) :: cstring(*)
  character(20) :: fstring
  real(c_double), intent(in) :: num_density_SI
  real(c_double) :: DENS
  DENS = num_density_SI * 1e-6
  call C_F_string_chars(cstring, fstring)
  call ZDPlasKin_set_density(trim(fstring), DENS)
end subroutine zdplaskinSetDensity

subroutine zdplaskinGetDensity(cstring, num_density_SI) bind(C, name='zdplaskinGetDensity')
  use C_interface_module
  use ZDPlasKin
  implicit none
  character(len=1,kind=C_char), intent(in) :: cstring(*)
  character(20) :: fstring
  real(c_double), intent(out) :: num_density_SI
  real(c_double) :: DENS
  call C_F_string_chars(cstring, fstring)
  call ZDPlasKin_get_density(trim(fstring),DENS)
  num_density_SI = DENS * 1e6
end subroutine zdplaskinGetDensity

subroutine zdplaskinSetElecField(elec_field, elec_frequency, num_density) bind(C,name='zdplaskinSetElecField')
  use, intrinsic :: iso_c_binding
  use ZDPlasKin
  implicit none
  ! input variables in SI unit 
  real(c_double), intent(in) :: elec_field ! [v/m] 
  real(c_double), intent(in) :: elec_frequency ! [1/s]
  real(c_double), intent(in) :: num_density ! [m-3]
  ! zdplaskin varibles
  real(c_double) :: reduced_field
  real(c_double) :: reduced_freqency
  reduced_field = elec_field / num_density * 1e21 ! [Td]
  reduced_freqency = elec_frequency / num_density * 1e6  ! [cm3 s-1]
  call ZDPlasKin_set_conditions(REDUCED_FIELD=reduced_field, REDUCED_FREQUENCY=reduced_freqency)
end subroutine zdplaskinSetElecField

subroutine zdplaskinSetReducedField(reduced_field, reduced_frequency) bind(C,name='zdplaskinSetReducedField')
  use, intrinsic :: iso_c_binding
  use ZDPlasKin
  implicit none
  ! input variables in SI unit 
  real(c_double), intent(in) :: reduced_field ! [Td]
  real(c_double), intent(in) :: reduced_frequency ! [1/s]
  call ZDPlasKin_set_conditions(REDUCED_FIELD=reduced_field, REDUCED_FREQUENCY=reduced_frequency)
end subroutine zdplaskinSetReducedField

subroutine zdplaskinSoftReset() bind(C,name='zdplaskinSoftReset')
  use, intrinsic :: iso_c_binding
  use ZDPlasKin
  implicit none
  call ZDPlasKin_set_conditions(SOFT_RESET=.true.)
end subroutine zdplaskinSoftReset

subroutine zdplaskinSetGasTemp(Tgas) bind(C,name='zdplaskinSetGasTemp')
  use, intrinsic :: iso_c_binding
  use ZDPlasKin
  implicit none
  ! input variables in SI unit 
  real(c_double), intent(in) :: Tgas ![K]
  call ZDPlasKin_set_conditions(GAS_TEMPERATURE=Tgas)
end subroutine zdplaskinSetGasTemp

subroutine zdplaskinSetElecTemp(Te) bind(C,name='zdplaskinSetElecTemp')
  use, intrinsic :: iso_c_binding
  use ZDPlasKin
  implicit none
  ! input variables in SI unit 
  real(c_double), intent(in) :: Te ![K]
  call ZDPlasKin_set_conditions(ELEC_TEMPERATURE=Te)
end subroutine zdplaskinSetElecTemp

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

function zdplaskinGetElecPower(ND) bind(C,name='zdplaskinGetElecPower') result(ep)
  use, intrinsic :: iso_c_binding
  use ZDPlasKin
  implicit none
  real(c_double) :: ep
  real(c_double), intent(in) :: ND
  real(c_double) :: ElectronCharge = 1.602176565e-19
  call ZDPlasKin_get_conditions(ELEC_POWER_N=ep)
  ep = ep * ND * ElectronCharge * 1e-6
end function zdplaskinGetElecPower

function zdplaskinGetElecPowerElastic(ND) bind(C,name='zdplaskinGetElecPowerElastic') result(ep)
  use, intrinsic :: iso_c_binding
  use ZDPlasKin
  implicit none
  real(c_double) :: ep
  real(c_double), intent(in) :: ND
  real(c_double) :: ElectronCharge = 1.602176565e-19
  call ZDPlasKin_get_conditions(ELEC_POWER_ELASTIC_N=ep)
  ep = ep * ND * ElectronCharge * 1e-6
end function zdplaskinGetElecPowerElastic

function zdplaskinGetElecPowerInelastic(ND) bind(C,name='zdplaskinGetElecPowerInelastic') result(ep)
  use, intrinsic :: iso_c_binding
  use ZDPlasKin
  implicit none
  real(c_double) :: ep
  real(c_double), intent(in) :: ND
  real(c_double) :: ElectronCharge = 1.602176565e-19
  call ZDPlasKin_get_conditions(ELEC_POWER_INELASTIC_N=ep)
  ep = ep * ND * ElectronCharge * 1e-6
end function zdplaskinGetElecPowerInelastic

subroutine zdplaskinGetPlasmaSource(carray) bind(C,name='zdplaskinGetPlasmaSource')
  use, intrinsic :: iso_c_binding
  use ZDPlasKin
  implicit none
  type(c_ptr), intent(out) :: carray
  real(c_double), target, save :: farray(species_max)
  real(c_double) :: Avogadro = 6.02214129e26
  integer i
  call ZDPlasKin_get_rates(SOURCE_TERMS=farray)
  do i = 1, species_max
    farray(i) = farray(i) * 1e6 / Avogadro
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
  integer(c_size_t), intent(in) :: index
  character(len=1,kind=C_char), dimension(*), intent(out) :: cstring
  call F_C_string_chars(trim(species_name(index+1)), cstring)
end subroutine zdplaskinGetSpeciesName
