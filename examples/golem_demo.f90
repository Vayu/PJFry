
! Copied from golem95, for testing purposes

module golem_demo

  contains

  pure function packb(set) result(bits)
    integer, intent(in), dimension(:) :: set
    integer :: bits
    bits = sum( ibset(0,pos=set) )
  end function packb

  pure function scalar(p1,p2)
    integer, parameter :: ki = kind(1.0d0)
    real(ki), intent (in), dimension(4) :: p1,p2
    real(ki) :: scalar
    scalar = p1(1)*p2(1) - p1(2)*p2(2) - p1(3)*p2(3) - p1(4)*p2(4)
  end function scalar

end module golem_demo
