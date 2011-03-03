!
! This program computes form factors for five-point functions and related
! (by pinches) diagrams with less external legs.
!
program main
  !
  use pjfry95pg
  use golem_demo
  !
  implicit none
  integer, parameter :: ki = kind(1.0d0)
  !
  complex(ki) :: res6a
  complex(ki) :: res6b
  complex(ki) :: res6c
  real(ki) :: t1,t2
  integer :: choix, i1,i2,i3,i4,i5,s1,s2,s3,s_null,b_pin
  !
  ! Opening of the files containing the results
  !
  open(unit=17,file='test5point.txt',status='unknown')
  !
  ! These are the entries of the S matrix
  ! They are related to the cuts of the following diagram
  ! All the momenta are incoming : p1+p2+p3+p4+p5 = 0
  !
  !   p1            p5
  !    \             /
  !     \   (5)    /
  !      |---<---\ (4)
  !      |           \
  !      |            \____ p4
  ! (1) |            /
  !      |           / (3)
  !      |--->---/
  !     /   (2)    \
  !    /             \
  !   p2            p3
  !
  ! S(1,3) = (p2+p3)^2
  ! S(2,4) = (p3+p4)^2
  ! S(2,5) = (p1+p2)^2
  ! S(3,5) = (p4+p5)^2
  ! S(1,4) = (p1+p5)^2
  ! S(1,2) = p2^2
  ! S(2,3) = p3^2
  ! S(3,4) = p4^2
  ! S(4,5) = p5^2
  ! S(1,5) = p1^2
  !
  !
  call pginitgolem95(5)
  s_null=0
  !
  ! Definition of the S matrix
  !
  call pgsetmat(1,1,  0.0_ki )
  call pgsetmat(1,2 ,  0.0_ki )
  call pgsetmat(1,3 ,  -3.0_ki )
  call pgsetmat(1,4 ,  -4.0_ki )
  call pgsetmat(1,5 ,  0.0_ki )
  !
  call pgsetmat(2,1 ,  pggetmat(1,2) )
  call pgsetmat(2,2 ,  0.0_ki )
  call pgsetmat(2,3 ,  0.0_ki )
  call pgsetmat(2,4 ,  6.0_ki )
  call pgsetmat(2,5 ,  15.0_ki )
  !
  call pgsetmat(3,1 ,  pggetmat(1,3) )
  call pgsetmat(3,2 ,  pggetmat(2,3) )
  call pgsetmat(3,3 ,  0.0_ki )
  call pgsetmat(3,4 ,  0.0_ki )
  call pgsetmat(3,5 ,  2.0_ki )
  !
  call pgsetmat(4,1 ,  pggetmat(1,4) )
  call pgsetmat(4,2 ,  pggetmat(2,4) )
  call pgsetmat(4,3 ,  pggetmat(3,4) )
  call pgsetmat(4,4 ,  0.0_ki )
  call pgsetmat(4,5 ,  0.0_ki )
  !
  call pgsetmat(5,1 ,  pggetmat(1,5) )
  call pgsetmat(5,2 ,  pggetmat(2,5) )
  call pgsetmat(5,3 ,  pggetmat(3,5) )
  call pgsetmat(5,4 ,  pggetmat(4,5) )
  call pgsetmat(5,5 ,  0.0_ki )
  !
  ! This call initialize cache with kinematics above
  !
  call pgpreparesmatrix()
  !
  write (*,*) 'Choose what the program should compute:'
  write (*,*) '0) form factor for five-point function, rank 0'
  write (*,*) '1) form factor for five-point function, rank 3 (z1*z2*z4)'
  write (*,*) '2) form factor for five-point function, rank 5 (z1*z2*z3*z4*z5)'
  write (*,*) '3) form factor for diagram with propagator 3 pinched, rank 0'
  write (*,*) '4) form factor for diagram with propagators 1 and 4 pinched, rank 0'
  write (*,*) '5) Test table'
  read (*,*) choix
    !
  if (choix == 0) then
    !
    write (*,*) 'calculating form factor for 5-point function rank 0'
    !
  else if   (choix == 1) then
    !
    write (*,*) 'calculating form factor A_124 for 5-point function rank 3'
    !
  else if   (choix == 2) then
    !
    write (*,*) 'calculating form factor A_12345 for 5-point function rank 5'
    !
  else if   (choix == 3) then
    !
    write (*,*) 'calculating form factor for a box stemming from '
    write (*,*) 'the pinch of propagator 3 of a 5-point funct., rank0'
    !
  else if   (choix == 4) then
    !
    write (*,*) 'calculating form factor for a triangle stemming from '
    write (*,*) 'the pinch of propagators 1 and 4 of a 5-point funct., rank 0'
    !
  end if
  !
  write (*,*) 'The result has been written to the file test5point.txt'
  !
  call cpu_time(t1)
  !
  ! To change the value of mu^2 (in GeV) (set to 1. by default)
  ! uncomment this line
  !
  !call pgsetmusq(12._ki)
  !
  write (17,*) 'The kinematics is:'
  write (17,*) ''
  write (17,*) ' p1            p5    '
  write (17,*) '   \             /        '
  write (17,*) '    \   (5)    /         '
  write (17,*) '     |---<---\ (4)     '
  write (17,*) '     |           \        '
  write (17,*) '     |            \____ p4'
  write (17,*) ' (1)|            /       '
  write (17,*) '     |           / (3)   '
  write (17,*) '     |--->---/         '
  write (17,*) '    /   (2)    \         '
  write (17,*) '   /             \        '
  write (17,*) ' p2            p3    '
  write (17,*) ''
  write (17,*) 'p1+p2+p3+p4+p5 = 0'
  write (17,*) ''
  write (17,*) ' S(1,3) = (p2+p3)^2=',pggetmat(1,3)
  write (17,*) ' S(2,4) = (p3+p4)^2=',pggetmat(2,4)
  write (17,*) ' S(2,5) = (p1+p2)^2=',pggetmat(2,5)
  write (17,*) ' S(3,5) = (p4+p5)^2=',pggetmat(3,5)
  write (17,*) ' S(1,4) = (p1+p5)^2=',pggetmat(1,4)
  write (17,*) ' S(1,2) = p2^2=',pggetmat(1,2)
  write (17,*) ' S(2,3) = p3^2=',pggetmat(2,3)
  write (17,*) ' S(3,4) = p4^2=',pggetmat(3,4)
  write (17,*) ' S(4,5) = p5^2=',pggetmat(4,5)
  write (17,*) ' S(1,5) = p1^2=',pggetmat(1,5)
  write (17,*) '(mu)^2 =',pggetmusq()
  write (17,*) ''
  !
  if (choix == 0) then
    !
    ! form factor for five-point function, rank 0
    !
    res6c = pga50(s_null, 0)
    res6b = pga50(s_null, 1)
    res6a = pga50(s_null, 2)
    !
  else if (choix == 1) then
    !
    ! form factor for five-point function, rank 3
    !
    res6c = pga53(1,2,4, s_null, 0)
    res6b = pga53(1,2,4, s_null, 1)
    res6a = pga53(1,2,4, s_null, 2)
    !
  else if (choix == 2) then
    !
    ! form factor for five-point function, rank 5
    !
    res6c = pga55(1,2,3,4,4, s_null, 0)
    res6b = pga55(1,2,3,4,4, s_null, 1)
    res6a = pga55(1,2,3,4,4, s_null, 2)
    !
  else if (choix == 3) then
    !
    ! form factor for pinched diagram, rank 0
    ! the propagator 3 is pinched
    !
    !
    !   p1            p4
    !    \             /
    !     \   (5)    /
    !      |---<---|
    !      |          |
    ! (1) |          | (4)
    !      |          |
    !      |--->---|
    !     /   (2)     \
    !    /              \
    !  p2            p3+p4
    !
    write (17,*) 'Since the propagator 3 is pinched'
    write (17,*) 'the reduced kinematics is:'
    write (17,*) ''
    write (17,*) ' p1             p5'
    write (17,*) '   \              /  '
    write (17,*) '    \   (5)     /   '
    write (17,*) '     |---<---|    '
    write (17,*) '     |          |    '
    write (17,*) ' (1)|          |(4) '
    write (17,*) '     |          |    '
    write (17,*) '     |--->---|    '
    write (17,*) '    /   (2)     \   '
    write (17,*) '   /               \  '
    write (17,*) ' p2            p3+p4'
    write (17,*) ''
    !
    !
    b_pin=packb( (/3/) )
    res6c = pga40(b_pin, 0)
    res6b = pga40(b_pin, 1)
    res6a = pga40(b_pin, 2)
    !
  else if (choix == 4) then
    !
    ! form factor for pinched diagram, rank 0
    ! the propagators 1 and 4 are pinched
    !
    !           |
    !           | p1+p2
    !           |
    !           /\
    !          /  \
    !    (2) /    \ (5)
    !        /      \
    !       /--->--\
    !      /   (3)    \
    ! p3  /            \ p4+p5
    !
    write (17,*) 'Since the propagators 1 and 4 are pinched'
    write (17,*) 'the reduced kinematics is:'
    write (17,*) ''
    write (17,*) '            |             '
    write (17,*) '            | p1+p2       '
    write (17,*) '            |             '
    write (17,*) '            /\            '
    write (17,*) '           /  \           '
    write (17,*) '     (2) /    \ (5)      '
    write (17,*) '         /      \         '
    write (17,*) '        /--->--\        '
    write (17,*) '       /   (3)    \       '
    write (17,*) ' p3 /             \ p4+p5'
    write (17,*) ''
    !
    !
    b_pin=packb( (/1,4/) )
    res6c = pga30(b_pin, 0)
    res6b = pga30(b_pin, 1)
    res6a = pga30(b_pin, 2)
    !
! --------------------------------------------------------------
! --------------------------------------------------------------
! --------------------------------------------------------------
  else if (choix == 5) then
    !
    write(17,*) " PENTAGON "
    write(17,*) " Rank 0 "
    write(17,'("E0( 0) = (",D23.15," +I* ",D23.15,")")') pga50(s_null, 0)
    write(17,'("E0(-1) = (",D23.15," +I* ",D23.15,")")') pga50(s_null, 1)
    write(17,'("E0(-2) = (",D23.15," +I* ",D23.15,")")') pga50(s_null, 2)
    write(17,*) " Rank 1 "
    do i1=1,5
      write(17,'("E",I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1, pga51(i1, s_null, 0)
      write(17,'("E",I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1, pga51(i1, s_null, 1)
      write(17,'("E",I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1, pga51(i1, s_null, 2)
    enddo
    write(17,*) " Rank 2 "
    write(17,'("E00( 0) = (",D23.15," +I* ",D23.15,")")') pgb52( s_null, 0)
    write(17,'("E00(-1) = (",D23.15," +I* ",D23.15,")")') pgb52( s_null, 1)
    write(17,'("E00(-2) = (",D23.15," +I* ",D23.15,")")') pgb52( s_null, 2)
    do i1=1,5
      do i2=i1,5
        write(17,'("E",2I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,i2, pga52(i1,i2, s_null, 0)
        write(17,'("E",2I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,i2, pga52(i1,i2, s_null, 1)
        write(17,'("E",2I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,i2, pga52(i1,i2, s_null, 2)
      enddo
    enddo
    write(17,*) " Rank 3 "
    do i1=1,5
      write(17,'("E00",I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1, pgb53(i1, s_null, 0)
      write(17,'("E00",I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1, pgb53(i1, s_null, 1)
      write(17,'("E00",I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1, pgb53(i1, s_null, 2)
    enddo
    do i1=1,5
      if (i1.eq.s1) cycle
      do i2=i1,5
        if (i2.eq.s1) cycle
        do i3=i2,5
          if (i3.eq.s1) cycle
          write(17,'("E",3I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3, pga53(i1,i2,i3, s_null, 0)
          write(17,'("E",3I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3, pga53(i1,i2,i3, s_null, 1)
          write(17,'("E",3I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3, pga53(i1,i2,i3, s_null, 2)
        enddo
      enddo
    enddo
    write(17,*) " Rank 4 "
    write(17,'("E0000( 0) = (",D23.15," +I* ",D23.15,")")') pgc54( s_null, 0)
    write(17,'("E0000(-1) = (",D23.15," +I* ",D23.15,")")') pgc54( s_null, 1)
    write(17,'("E0000(-2) = (",D23.15," +I* ",D23.15,")")') pgc54( s_null, 2)
    do i1=1,5
      if (i1.eq.s1) cycle
      do i2=i1,5
        if (i2.eq.s1) cycle
        write(17,'("E00",2I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,i2, pgb54(i1,i2, s_null, 0)
        write(17,'("E00",2I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,i2, pgb54(i1,i2, s_null, 1)
        write(17,'("E00",2I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,i2, pgb54(i1,i2, s_null, 2)
      enddo
    enddo
    do i1=1,5
      if (i1.eq.s1) cycle
      do i2=i1,5
        if (i2.eq.s1) cycle
        do i3=i2,5
          if (i3.eq.s1) cycle
          do i4=i3,5
            if (i4.eq.s1) cycle
            write(17,'("E",4I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3,i4, pga54(i1,i2,i3,i4, s_null, 0)
            write(17,'("E",4I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3,i4, pga54(i1,i2,i3,i4, s_null, 1)
            write(17,'("E",4I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3,i4, pga54(i1,i2,i3,i4, s_null, 2)
          enddo
        enddo
      enddo
    enddo
    write(17,*) " Rank 5 "
    do i1=1,5
      if (i1.eq.s1) cycle
      write(17,'("E0000",I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,pgc55(i1, s_null, 0)
      write(17,'("E0000",I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,pgc55(i1, s_null, 1)
      write(17,'("E0000",I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,pgc55(i1, s_null, 2)
    enddo
    do i1=1,5
      if (i1.eq.s1) cycle
      do i2=i1,5
        if (i2.eq.s1) cycle
        do i3=i2,5
          if (i3.eq.s1) cycle
          write(17,'("E00",3I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3, pgb55(i1,i2,i3, s_null, 0)
          write(17,'("E00",3I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3, pgb55(i1,i2,i3, s_null, 1)
          write(17,'("E00",3I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3, pgb55(i1,i2,i3, s_null, 2)
        enddo
      enddo
    enddo
    do i1=1,5
      if (i1.eq.s1) cycle
      do i2=i1,5
        if (i2.eq.s1) cycle
        do i3=i2,5
          if (i3.eq.s1) cycle
          do i4=i3,5
            if (i4.eq.s1) cycle
            do i5=i4,5
              if (i5.eq.s1) cycle
              write(17,'("E",5I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3,i4,i5, pga55(i1,i2,i3,i4,i5, s_null, 0)
              write(17,'("E",5I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3,i4,i5, pga55(i1,i2,i3,i4,i5, s_null, 1)
              write(17,'("E",5I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3,i4,i5, pga55(i1,i2,i3,i4,i5, s_null, 2)
            enddo
          enddo
        enddo
      enddo
    enddo
    write(17,*) " Rank 5 DISOREDER "
    i1=4
    i2=3
    i3=4
    i4=2
    i5=1
    write(17,'("E",5I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3,i4,i5, pga55(i1,i2,i3,i4,i5, s_null, 0)
    write(17,'("E",5I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3,i4,i5, pga55(i1,i2,i3,i4,i5, s_null, 1)
    write(17,'("E",5I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3,i4,i5, pga55(i1,i2,i3,i4,i5, s_null, 2)
! --------------------------------------------------------------
! --------------------------------------------------------------
    do s1=1,5
      b_pin=packb( (/s1/) )
      write(17,*) " BOX s = ",s1
      write(17,*) " Rank 0 "
      write(17,'("D0( 0) = (",D23.15," +I* ",D23.15,")")') pga40(b_pin, 0)
      write(17,'("D0(-1) = (",D23.15," +I* ",D23.15,")")') pga40(b_pin, 1)
      write(17,'("D0(-2) = (",D23.15," +I* ",D23.15,")")') pga40(b_pin, 2)
      write(17,*) " Rank 1 "
      do i1=1,5
        write(17,'("D",I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1, pga41(i1, b_pin, 0)
        write(17,'("D",I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1, pga41(i1, b_pin, 1)
        write(17,'("D",I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1, pga41(i1, b_pin, 2)
      enddo
      ! --------------------------------
      write(17,*) " Rank 2 "
      write(17,'("D00( 0) = (",D23.15," +I* ",D23.15,")")') pgb42(b_pin, 0)
      write(17,'("D00(-1) = (",D23.15," +I* ",D23.15,")")') pgb42(b_pin, 1)
      write(17,'("D00(-2) = (",D23.15," +I* ",D23.15,")")') pgb42(b_pin, 2)
      do i1=1,5
        if (i1.eq.s1) cycle
        do i2=i1,5
          if (i2.eq.s1) cycle
          write(17,'("D",2I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,i2, pga42(i1,i2, b_pin, 0)
          write(17,'("D",2I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,i2, pga42(i1,i2, b_pin, 1)
          write(17,'("D",2I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,i2, pga42(i1,i2, b_pin, 2)
        enddo
      enddo
      ! --------------------------------
      write(17,*) " Rank 3 "
      do i1=1,5
        write(17,'("D00",I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1, pgb43(i1, b_pin, 0)
        write(17,'("D00",I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1, pgb43(i1, b_pin, 1)
        write(17,'("D00",I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1, pgb43(i1, b_pin, 2)
      enddo
      do i1=1,5
        if (i1.eq.s1) cycle
        do i2=i1,5
          if (i2.eq.s1) cycle
          do i3=i2,5
            if (i3.eq.s1) cycle
            write(17,'("D",3I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3, pga43(i1,i2,i3, b_pin, 0)
            write(17,'("D",3I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3, pga43(i1,i2,i3, b_pin, 1)
            write(17,'("D",3I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3, pga43(i1,i2,i3, b_pin, 2)
          enddo
        enddo
      enddo
      ! --------------------------------
      write(17,*) " Rank 4 "
      write(17,'("D0000( 0) = (",D23.15," +I* ",D23.15,")")') pgc44( b_pin, 0)
      write(17,'("D0000(-1) = (",D23.15," +I* ",D23.15,")")') pgc44( b_pin, 1)
      write(17,'("D0000(-2) = (",D23.15," +I* ",D23.15,")")') pgc44( b_pin, 2)
      do i1=1,5
        if (i1.eq.s1) cycle
        do i2=i1,5
          if (i2.eq.s1) cycle
          write(17,'("D00",2I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,i2, pgb44(i1,i2, b_pin, 0)
          write(17,'("D00",2I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,i2, pgb44(i1,i2, b_pin, 1)
          write(17,'("D00",2I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,i2, pgb44(i1,i2, b_pin, 2)
        enddo
      enddo
      do i1=1,5
        if (i1.eq.s1) cycle
        do i2=i1,5
          if (i2.eq.s1) cycle
          do i3=i2,5
            if (i3.eq.s1) cycle
            do i4=i3,5
              if (i4.eq.s1) cycle
              write(17,'("D",4I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3,i4, pga44(i1,i2,i3,i4, b_pin, 0)
              write(17,'("D",4I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3,i4, pga44(i1,i2,i3,i4, b_pin, 1)
              write(17,'("D",4I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3,i4, pga44(i1,i2,i3,i4, b_pin, 2)
            enddo
          enddo
        enddo
      enddo
    ! --------------------------------------------------------------
      do s2=(s1+1),5
        b_pin=packb( (/s1,s2/) )
        write(17,*) " TRIANGLE s,t = ",s1,s2
        write(17,*) " Rank 0 "
        write(17,'("C0( 0) = (",D23.15," +I* ",D23.15,")")') pga30(b_pin, 0)
        write(17,'("C0(-1) = (",D23.15," +I* ",D23.15,")")') pga30(b_pin, 1)
        write(17,'("C0(-2) = (",D23.15," +I* ",D23.15,")")') pga30(b_pin, 2)
        write(17,*) " Rank 1 "
        do i1=1,5
          if ((i1.eq.s1).or.(i1.eq.s2)) cycle
          write(17,'("C",I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1, pga31(i1, b_pin, 0)
          write(17,'("C",I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1, pga31(i1, b_pin, 1)
          write(17,'("C",I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1, pga31(i1, b_pin, 2)
        enddo
        ! --------------------------------
        write(17,*) " Rank 2 "
        write(17,'("C00( 0) = (",D23.15," +I* ",D23.15,")")') pgb32(b_pin, 0)
        write(17,'("C00(-1) = (",D23.15," +I* ",D23.15,")")') pgb32(b_pin, 1)
        write(17,'("C00(-2) = (",D23.15," +I* ",D23.15,")")') pgb32(b_pin, 2)
        do i1=1,5
          if ((i1.eq.s1).or.(i1.eq.s2)) cycle
          do i2=i1,5
            if ((i2.eq.s1).or.(i2.eq.s2)) cycle
            write(17,'("C",2I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,i2, pga32(i1,i2, b_pin, 0)
            write(17,'("C",2I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,i2, pga32(i1,i2, b_pin, 1)
            write(17,'("C",2I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,i2, pga32(i1,i2, b_pin, 2)
          enddo
        enddo
        ! --------------------------------
        write(17,*) " Rank 3 "
        do i1=1,5
          if ((i1.eq.s1).or.(i1.eq.s2)) cycle
          write(17,'("C00",I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1, pgb33(i1, b_pin, 0)
          write(17,'("C00",I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1, pgb33(i1, b_pin, 1)
          write(17,'("C00",I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1, pgb33(i1, b_pin, 2)
        enddo
        do i1=1,5
          if ((i1.eq.s1).or.(i1.eq.s2)) cycle
          do i2=i1,5
            if ((i2.eq.s1).or.(i2.eq.s2)) cycle
            do i3=i2,5
              if ((i3.eq.s1).or.(i3.eq.s2)) cycle
              write(17,'("C",3I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3, pga33(i1,i2,i3, b_pin, 0)
              write(17,'("C",3I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3, pga33(i1,i2,i3, b_pin, 1)
              write(17,'("C",3I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3, pga33(i1,i2,i3, b_pin, 2)
            enddo
          enddo
        enddo
      ! --------------------------------------------------------------
        do s3=(s2+1),5
          b_pin=packb( (/s1,s2,s3/) )
          write(17,*) " BUBBLES s,t,u = ",s1,s2,s3
          write(17,*) " Rank 0 "
          write(17,'("B0( 0) = (",D23.15," +I* ",D23.15,")")') pga20(b_pin, 0)
          write(17,'("B0(-1) = (",D23.15," +I* ",D23.15,")")') pga20(b_pin, 1)
          write(17,'("B0(-2) = (",D23.15," +I* ",D23.15,")")') pga20(b_pin, 2)
          write(17,*) " Rank 1 "
          do i1=1,5
            if ((i1.eq.s1).or.(i1.eq.s2).or.(i1.eq.s3)) cycle
            write(17,'("B",I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1, pga21(i1, b_pin, 0)
            write(17,'("B",I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1, pga21(i1, b_pin, 1)
            write(17,'("B",I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1, pga21(i1, b_pin, 2)
          enddo
          ! --------------------------------
          write(17,*) " Rank 2 "
          write(17,'("B00( 0) = (",D23.15," +I* ",D23.15,")")') pgb22(b_pin, 0)
          write(17,'("B00(-1) = (",D23.15," +I* ",D23.15,")")') pgb22(b_pin, 1)
          write(17,'("B00(-2) = (",D23.15," +I* ",D23.15,")")') pgb22(b_pin, 2)
          do i1=1,5
            if ((i1.eq.s1).or.(i1.eq.s2).or.(i1.eq.s3)) cycle
            do i2=i1,5
              if ((i2.eq.s1).or.(i2.eq.s2).or.(i2.eq.s3)) cycle
              write(17,'("B",2I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,i2, pga22(i1,i2, b_pin, 0)
              write(17,'("B",2I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,i2, pga22(i1,i2, b_pin, 1)
              write(17,'("B",2I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,i2, pga22(i1,i2, b_pin, 2)
            enddo
          enddo
        enddo
      ! --------------------------------------------------------------
      enddo
    ! --------------------------------------------------------------
    enddo
  ! --------------------------------------------------------------
  ! --------------------------------------------------------------
  end if
! --------------------------------------------------------------
! --------------------------------------------------------------
! --------------------------------------------------------------
  ! 
  call cpu_time(t2)
  !!
  if (choix.ne.5) then
    write (17,*) 'the program gives numbers for P2,P1,P0'
    write (17,*) ''
    write (17,'("  1/epsilon^2 * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6a,ki),aimag(res6a)
    write (17,'("+ 1/epsilon   * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6b,ki),aimag(res6b)
    write (17,'("+ 1           * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6c,ki),aimag(res6c)
  endif
  write (17,*) ''
  write (17,*) 'CPU time=',t2-t1
  !
  !
  !
  close(17)
  !
end program main
