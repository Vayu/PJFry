
# this file is part of PJFry library
# Copyright 2011 Valery Yundin

pg_prefix = "pg"
pj_prefix = "pj"

pjff_template = """
  function %(pjfunc)s(%(args)sb_pin)
    integer, intent (in) :: %(args)sb_pin
    type(form_factor) :: %(pjfunc)s
    complex(ki), dimension(3) :: temp
    temp(1)=%(pgfunc)s(%(args)sb_pin, 2)
    temp(2)=%(pgfunc)s(%(args)sb_pin, 1)
    temp(3)=%(pgfunc)s(%(args)sb_pin, 0)
    %(pjfunc)s=temp
  end function %(pjfunc)s
"""

pjff_template6 = """
  function %(pjfunc)s(%(args)sb_pin)
    integer, intent (in) :: %(args)sb_pin
    type(form_factor) :: %(pjfunc)s
    complex(ki), dimension(3) :: temp
    temp(:)=0._ki
    %(pjfunc)s=temp
  end function %(pjfunc)s
"""

pgff_template = """
  interface
    function %(pgfunc)s(%(args)sb_pin, ep)
!      use precision_pjfry, only: ki_pjf
      implicit none
      integer, parameter :: ki_pjf = kind(1.0d0)
      integer, intent(in) :: %(args)sb_pin, ep
      complex(ki_pjf) :: %(pgfunc)s
    end function %(pgfunc)s
  end interface
"""

initgolem95_template = """
  interface
    subroutine %(pref)sinitgolem95(n)
      implicit none
      integer, intent(in) :: n
    end subroutine %(pref)sinitgolem95
  end interface
"""

getmusq_template = """
  interface
    function %(pref)sgetmusq()
!      use precision_pjfry, only: ki_pjf
      implicit none
      integer, parameter :: ki_pjf = kind(1.0d0)
      real(ki_pjf) :: %(pref)sgetmusq
    end function %(pref)sgetmusq
  end interface
"""

setmusq_template = """
  interface
    subroutine %(pref)ssetmusq(musq)
!      use precision_pjfry, only: ki_pjf
      implicit none
      integer, parameter :: ki_pjf = kind(1.0d0)
      real(ki_pjf), intent(in) :: musq
    end subroutine %(pref)ssetmusq
  end interface
"""

setmat_template = """
  interface
    subroutine %(pref)ssetmat(i, j, val)
!      use precision_pjfry, only: ki_pjf
      implicit none
      integer, parameter :: ki_pjf = kind(1.0d0)
      integer, intent(in) :: i, j
      real(ki_pjf), intent(in) :: val
    end subroutine %(pref)ssetmat
  end interface
"""

getmat_template = """
  interface
    function %(pref)sgetmat(i, j)
!      use precision_pjfry, only: ki_pjf
      implicit none
      integer, parameter :: ki_pjf = kind(1.0d0)
      integer, intent(in) :: i, j
      real(ki_pjf) :: %(pref)sgetmat
    end function %(pref)sgetmat
  end interface
"""

preparesmatrix_template = """
  interface
    subroutine %(pref)spreparesmatrix()
      implicit none
    end subroutine %(pref)spreparesmatrix
  end interface
"""

pg_funcs = [  initgolem95_template,
              getmusq_template,
              setmusq_template,
              setmat_template,
              getmat_template,
              preparesmatrix_template ]

def ff_list(N):
  args = ["i", "j", "k", "l", "m", "n"]
  for N in range(min_legs,max_legs+1):
    for r in range(N+1):
      yield {'func' : "a%d%d" % (N, r), 'args' : args[:r], 'legs' : N }
      if N < 6 and r >= 2:
        yield { 'func' : "b%d%d" % (N, r), 'args' : args[:(r-2)], 'legs' : N }
      if N < 6 and r >= 4:
        yield { 'func' : "c%d%d" % (N, r), 'args' : args[:(r-4)], 'legs' : N }

# -------------------------------------

def write_module_golem95pg(name, min_legs, max_legs, max_rank):
  f = open("%s.f90" % name, 'w')
  f.write("module %s\n" % name)
  f.write("""
  implicit none

""")

  for func_template in pg_funcs:
    f.write(func_template % { 'pref' : pg_prefix} )

  for ff in ff_list(max_legs):
    f.write(pgff_template % {
        'pjfunc' : pj_prefix+ff['func'],
        'pgfunc' : pg_prefix+ff['func'],
        'args' : ", ".join(ff['args']+[""])
        })
  f.write("""
  contains

  end module %s
  """ % name)
  f.close()

# -------------------------------------

def write_module_golem95pj(name, min_legs, max_legs, max_rank):
  f = open("%s.f90" % name, 'w')
  f.write("module %s\n" % name)
  f.write("""
  use precision_golem, only: ki
  use pjfry95pg
  use form_factor_type
  implicit none

  contains

""")
  for ff in ff_list(max_legs):
    if ff['legs'] == 6:
      f.write(pjff_template6 % {
          'pjfunc' : pj_prefix+ff['func'],
          'pgfunc' : pg_prefix+ff['func'],
          'args' : ", ".join(ff['args']+[""])
          })
    else:
      f.write(pjff_template % {
          'pjfunc' : pj_prefix+ff['func'],
          'pgfunc' : pg_prefix+ff['func'],
          'args' : ", ".join(ff['args']+[""])
          })
  f.write("end module %s\n" % name)
  f.close()

# -------------------------------------

if __name__ == "__main__":
  min_legs = 2
  max_legs = 6
  max_rank = 6
  write_module_golem95pj("pjfry95pj", min_legs, max_legs, max_rank)
  write_module_golem95pg("pjfry95pg", min_legs, max_legs, max_rank)
