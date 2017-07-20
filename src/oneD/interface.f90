subroutine zdplaskinInit() bind(C, name='zdplaskinInit')
  use ZDPlasKin
  call ZDPlasKin_init()
end subroutine zdplaskinInit

subroutine zdplaskinSetDensity(cstring, DENS) bind(C, name='zdplaskinSetDensity')
  use, intrinsic :: iso_c_binding
  use ZDPlasKin
  type(c_ptr), target, intent(in) :: cstring
  character(kind=c_char), pointer :: fstring(:)
  character(10) :: string
  integer length
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

  do I=2,size(fstring(:))
    string = trim(string)//fstring(I)
  enddo

  string = trim(string)
  call ZDPlasKin_set_density(string, DENS)
end subroutine zdplaskinSetDensity

subroutine zdplaskinSetConditions(gasTemperature, reduced_field) bind(C,name='zdplaskinSetConditions')
  use, intrinsic :: iso_c_binding
  use ZDPlasKin
  real(c_double), intent(in) :: gasTemperature
  real(c_double), intent(in) :: reduced_field
  call ZDPlasKin_set_conditions(GAS_TEMPERATURE=gasTemperature, REDUCED_FIELD=reduced_field)
end subroutine zdplaskinSetConditions

function getElectronTemperature() bind(C,name='getElectronTemperature') result(eTemperature)
  use, intrinsic :: iso_c_binding
  use ZDPlasKin
  real(c_double) :: eTemperature
  call ZDPlasKin_get_conditions(ELEC_TEMPERATURE=eTemperature)
end function getElectronTemperature