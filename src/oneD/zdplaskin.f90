module math
    implicit none
    contains
    function getradius(degree) result(r)
        implicit none
        real*8, intent(in) :: degree
        real*8             :: r
        r = degree / 180.0 * 3.1415
    end function getradius
end module math

