subroutine setDensity(cstring, DENS) bind(C,name="setDensity")
    use, intrinsic :: iso_c_binding
    use ZDPlasKin
    type(c_ptr), target, intent(in) :: cstring
    character(kind=c_char), pointer :: fstring(:)
    character(*) :: forstr
    double precision, intent(in) :: DENS

    interface
     function strlen(s) bind(C, name='strlen')
       use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
       implicit none
       type(c_ptr), intent(in), value :: s
       integer(c_size_t) :: strlen
     end function strlen
    end interface

    call c_f_pointer(cstring, fstring, [strlen(cstring)])
    
    forstr = fstring

    !call ZDPlasKin_set_density(fstring, DENS)

end subroutine setDensity

