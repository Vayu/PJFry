

module precision_golem
    !
    integer, parameter :: ki=kind(1.0d0)
    integer, parameter :: ki_lt=kind(1.0d0)
    integer, parameter :: ki_avh=kind(1.0d0)
    integer, parameter :: ki_pjf=kind(1.0d0)
    !
end module precision_golem

module pjfry90
  implicit none

  integer, parameter :: ki_pjf=kind(1.0d0)

  interface
    function pgetmusq()
      use precision_golem, only: ki_pjf
      implicit none
      real(ki_pjf) :: pgetmusq
    end function pgetmusq
  end interface

  interface
    function psetmusq(musq)
      use precision_golem, only: ki_pjf
      implicit none
      real(ki_pjf), intent(in) :: musq
      real(ki_pjf) :: psetmusq
    end function psetmusq
  end interface

  interface
    function pe0(p1, p2, p3, p4, p5, p1p2, p2p3, p3p4, p4p5, p1p5, m1, m2, m3, m4, m5, ep)
      use precision_golem, only: ki_pjf
      implicit none
      real(ki_pjf), intent(in) :: p1, p2, p3, p4, p5
      real(ki_pjf), intent(in) :: p1p2, p2p3, p3p4, p4p5, p1p5
      real(ki_pjf), intent(in) :: m1, m2, m3, m4, m5
      integer, intent(in) :: ep
      complex(ki_pjf) :: pe0
    end function pe0
  end interface

  interface
    function pe0i(i, p1, p2, p3, p4, p5, p1p2, p2p3, p3p4, p4p5, p1p5, m1, m2, m3, m4, m5, ep)
      use precision_golem, only: ki_pjf
      implicit none
      integer, intent(in) :: i
      real(ki_pjf), intent(in) :: p1, p2, p3, p4, p5
      real(ki_pjf), intent(in) :: p1p2, p2p3, p3p4, p4p5, p1p5
      real(ki_pjf), intent(in) :: m1, m2, m3, m4, m5
      integer, intent(in) :: ep
      complex(ki_pjf) :: pe0i
    end function pe0i
  end interface

  interface
    function pe0ij(i, j, p1, p2, p3, p4, p5, p1p2, p2p3, p3p4, p4p5, p1p5, m1, m2, m3, m4, m5, ep)
      use precision_golem, only: ki_pjf
      implicit none
      integer, intent(in) :: i, j
      real(ki_pjf), intent(in) :: p1, p2, p3, p4, p5
      real(ki_pjf), intent(in) :: p1p2, p2p3, p3p4, p4p5, p1p5
      real(ki_pjf), intent(in) :: m1, m2, m3, m4, m5
      integer, intent(in) :: ep
      complex(ki_pjf) :: pe0ij
    end function pe0ij
  end interface

  interface
    function pe0ijk(i, j, k, p1, p2, p3, p4, p5, p1p2, p2p3, p3p4, p4p5, p1p5, m1, m2, m3, m4, m5, ep)
      use precision_golem, only: ki_pjf
      implicit none
      integer, intent(in) :: i, j, k
      real(ki_pjf), intent(in) :: p1, p2, p3, p4, p5
      real(ki_pjf), intent(in) :: p1p2, p2p3, p3p4, p4p5, p1p5
      real(ki_pjf), intent(in) :: m1, m2, m3, m4, m5
      integer, intent(in) :: ep
      complex(ki_pjf) :: pe0ijk
    end function pe0ijk
  end interface

  interface
    function pe0ijkl(i, j, k, l, p1, p2, p3, p4, p5, p1p2, p2p3, p3p4, p4p5, p1p5, m1, m2, m3, m4, m5, ep)
      use precision_golem, only: ki_pjf
      implicit none
      integer, intent(in) :: i, j, k, l
      real(ki_pjf), intent(in) :: p1, p2, p3, p4, p5
      real(ki_pjf), intent(in) :: p1p2, p2p3, p3p4, p4p5, p1p5
      real(ki_pjf), intent(in) :: m1, m2, m3, m4, m5
      integer, intent(in) :: ep
      complex(ki_pjf) :: pe0ijkl
    end function pe0ijkl
  end interface

  interface
    function pe0ijklm(i, j, k, l, m, p1, p2, p3, p4, p5, p1p2, p2p3, p3p4, p4p5, p1p5, m1, m2, m3, m4, m5, ep)
      use precision_golem, only: ki_pjf
      implicit none
      integer, intent(in) :: i, j, k, l, m
      real(ki_pjf), intent(in) :: p1, p2, p3, p4, p5
      real(ki_pjf), intent(in) :: p1p2, p2p3, p3p4, p4p5, p1p5
      real(ki_pjf), intent(in) :: m1, m2, m3, m4, m5
      integer, intent(in) :: ep
      complex(ki_pjf) :: pe0ijklm
    end function pe0ijklm
  end interface

  interface
    function pd0(p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4, ep)
      use precision_golem, only: ki_pjf
      implicit none
      real(ki_pjf), intent(in) :: p1, p2, p3, p4
      real(ki_pjf), intent(in) :: p1p2, p2p3
      real(ki_pjf), intent(in) :: m1, m2, m3, m4
      integer, intent(in) :: ep
      complex(ki_pjf) :: pd0
    end function pd0
  end interface

  interface
    function pd0i(i, p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4, ep)
      use precision_golem, only: ki_pjf
      implicit none
      integer, intent(in) :: i
      real(ki_pjf), intent(in) :: p1, p2, p3, p4
      real(ki_pjf), intent(in) :: p1p2, p2p3
      real(ki_pjf), intent(in) :: m1, m2, m3, m4
      integer, intent(in) :: ep
      complex(ki_pjf) :: pd0i
    end function pd0i
  end interface

  interface
    function pd0ij(i, j, p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4, ep)
      use precision_golem, only: ki_pjf
      implicit none
      integer, intent(in) :: i, j
      real(ki_pjf), intent(in) :: p1, p2, p3, p4
      real(ki_pjf), intent(in) :: p1p2, p2p3
      real(ki_pjf), intent(in) :: m1, m2, m3, m4
      integer, intent(in) :: ep
      complex(ki_pjf) :: pd0ij
    end function pd0ij
  end interface

  interface
    function pd0ijk(i, j, k, p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4, ep)
      use precision_golem, only: ki_pjf
      implicit none
      integer, intent(in) :: i, j, k
      real(ki_pjf), intent(in) :: p1, p2, p3, p4
      real(ki_pjf), intent(in) :: p1p2, p2p3
      real(ki_pjf), intent(in) :: m1, m2, m3, m4
      integer, intent(in) :: ep
      complex(ki_pjf) :: pd0ijk
    end function pd0ijk
  end interface

  interface
    function pd0ijkl(i, j, k, l, p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4, ep)
      use precision_golem, only: ki_pjf
      implicit none
      integer, intent(in) :: i, j, k, l
      real(ki_pjf), intent(in) :: p1, p2, p3, p4
      real(ki_pjf), intent(in) :: p1p2, p2p3
      real(ki_pjf), intent(in) :: m1, m2, m3, m4
      integer, intent(in) :: ep
      complex(ki_pjf) :: pd0ijkl
    end function pd0ijkl
  end interface

  interface
    function pc0(p1, p2, p3, m1, m2, m3, ep)
      use precision_golem, only: ki_pjf
      implicit none
      real(ki_pjf), intent(in) :: p1, p2, p3
      real(ki_pjf), intent(in) :: m1, m2, m3
      integer, intent(in) :: ep
      complex(ki_pjf) :: pc0
    end function pc0
  end interface

  interface
    function pc0i(i, p1, p2, p3, m1, m2, m3, ep)
      use precision_golem, only: ki_pjf
      implicit none
      integer, intent(in) :: i
      real(ki_pjf), intent(in) :: p1, p2, p3
      real(ki_pjf), intent(in) :: m1, m2, m3
      integer, intent(in) :: ep
      complex(ki_pjf) :: pc0i
    end function pc0i
  end interface

  interface
    function pc0ij(i, j, p1, p2, p3, m1, m2, m3, ep)
      use precision_golem, only: ki_pjf
      implicit none
      integer, intent(in) :: i, j
      real(ki_pjf), intent(in) :: p1, p2, p3
      real(ki_pjf), intent(in) :: m1, m2, m3
      integer, intent(in) :: ep
      complex(ki_pjf) :: pc0ij
    end function pc0ij
  end interface

  interface
    function pc0ijk(i, j, k, p1, p2, p3, m1, m2, m3, ep)
      use precision_golem, only: ki_pjf
      implicit none
      integer, intent(in) :: i, j, k
      real(ki_pjf), intent(in) :: p1, p2, p3
      real(ki_pjf), intent(in) :: m1, m2, m3
      integer, intent(in) :: ep
      complex(ki_pjf) :: pc0ijk
    end function pc0ijk
  end interface

  interface
    function pb0(p1, m1, m2, ep)
      use precision_golem, only: ki_pjf
      implicit none
      real(ki_pjf), intent(in) :: p1
      real(ki_pjf), intent(in) :: m1, m2
      integer, intent(in) :: ep
      complex(ki_pjf) :: pb0
    end function pb0
  end interface

  interface
    function pb0i(i, p1, m1, m2, ep)
      use precision_golem, only: ki_pjf
      implicit none
      integer, intent(in) :: i
      real(ki_pjf), intent(in) :: p1
      real(ki_pjf), intent(in) :: m1, m2
      integer, intent(in) :: ep
      complex(ki_pjf) :: pb0i
    end function pb0i
  end interface

  interface
    function pb0ij(i, j, p1, m1, m2, ep)
      use precision_golem, only: ki_pjf
      implicit none
      integer, intent(in) :: i, j
      real(ki_pjf), intent(in) :: p1
      real(ki_pjf), intent(in) :: m1, m2
      integer, intent(in) :: ep
      complex(ki_pjf) :: pb0ij
    end function pb0ij
  end interface

  interface
    function pa0(m1, ep)
      use precision_golem, only: ki_pjf
      implicit none
      real(ki_pjf), intent(in) :: m1
      integer, intent(in) :: ep
      complex(ki_pjf) :: pa0
    end function pa0
  end interface

  contains

end module pjfry90