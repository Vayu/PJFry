!
! This program computes form factors for the six-point functions and related
! (by pinches) diagrams with less external legs. 
! The normalisation is as follows:
! We define I_N^n= mu^(4-n) \int d^n k/(i*Pi^(n/2))*func(k,p_i)
! = r_Gam *(P2/eps^2+P1/eps+P0)
! n=4-2*eps
! r_Gam= Gamma(1+eps)*Gamma(1-eps)^2/Gamma(1-2eps)
! the program gives numbers for P2,P1,P0
!

program main
  !
  use pjfry95pg
  use golem_demo
  !
  implicit none
  !
  integer, parameter :: ki = kind(1.0d0)
  complex(ki) :: res6a,res6b,res6c
  real(ki) :: t1,t2
  real(ki), dimension(4) :: p1,p2,p3,p4,p5,p6
  real(ki), dimension(4) :: p12,p23,p34,p45,p56,p61,p123,p234,p345
  integer :: choix, b_pin,s1,s2,s3,s4,i1,i2,i3,i4,i5
   ! loop propagator masses
!    real(ki), dimension(6), parameter :: msq6 = (/ &
!             & 0.0_ki, &
!             & .96_ki, &
!             & 100.230_ki, &
!             & 8.0234343_ki, &
!             & 1.54987724_ki, &
!             & 0.0_ki &
!    & /)
  real(ki), dimension(6), parameter :: msq6 = (/ &
           & 0.0_ki,  &
           & 0.0_ki,  &
           & 0.0_ki,  &
           & 0.0_ki,  &
           & 0.0_ki,  &
           & 0.0_ki &
  & /)
  !
  ! Opening of the error file
  !
  open(unit=19,file='error_6point.txt',status='unknown')
  !
  ! Opening of the files containing the results
  !
  open(unit=17,file='test6point.txt',status='unknown')
  !
  !
  ! This is the entries of the S matrix
  ! defined as S(i,j) = (q_i - q_j)^2 where the q's are 
  ! the momentum flowing in the propagators
  ! It is related to the cuts of the following diagram
  ! All the momenta are incoming : p1+p2+p3+p4+p5+p6 = 0
  !
  !         p6            p5
  !          \           /
  !           \   (5)   /
  !            /---<---\ (4)
  !       (6) /         \
  !   p1 ____/           \____ p4
  !          \           /
  !       (1) \         / (3)
  !            \--->---/
  !           /   (2)   \
  !          /           \
  !         p2            p3
  !
  ! S(1,3) = (p2+p3)^2
  ! S(1,4) = (p2+p3+p4)^2
  ! S(1,5) = (p1+p6)^2
  ! S(2,4) = (p3+p4)^2
  ! S(2,5) = (p3+p4+p5)^2
  ! S(2,6) = (p1+p2)^2
  ! S(3,5) = (p4+p5)^2
  ! S(3,6) = (p4+p5+p6)^2
  ! S(4,6) = (p5+p6)^2
  ! S(1,2) = p2^2
  ! S(2,3) = p3^2
  ! S(3,4) = p4^2
  ! S(4,5) = p5^2
  ! S(5,6) = p6^2
  ! S(1,6) = p1^2
  !
  ! we define a set of four momentum whcih verify the different constraints
  p1 = (/ 0.5_ki , 0.0_ki , 0.0_ki , 0.5_ki /)
  p2 = (/ 0.5_ki , 0.0_ki , 0.0_ki , -0.5_ki /)
  p3 = -(/ 0.19178191094778038_ki , 0.12741179719516801_ki , 0.08262476614744381_ki , 0.11713105190921771_ki /)
  p4 = -(/ 0.33662712284553753_ki , -0.06648281097623857_ki , -0.3189378514746887_ki , -0.08471424069583446_ki /)
  p5 = -(/ 0.21604814388379073_ki , -0.20363139428835617_ki , 0.044157623555325_ki , 0.0571065672034082_ki /)
  p6 = -(/ 0.2555428223228916_ki , 0.1427024080694266_ki , 0.19215546177191994_ki , -0.08952337841679145_ki /)
  !
  p12 = p1 + p2
  p23 = p2 + p3
  p234 = p2+p3+p4
  p61 = p1 + p6
  p34 = p3 + p4
  p345 = p34 + p5
  p45 = p4 + p5
  p123 = p12 + p3
  p56 = p5 + p6
  !  
  !
  ! Allocates memory to store the set of initial propagators, the S matrix, 
  ! its inverse and the b coefficients.
  ! This call will allocate a derived type s_mat_p object.
  ! Includes calls to allocation_s and initializes the caching system
  !
  call pginitgolem95(6)
  !
  !
  ! Definition of the S matrix
  !
  call pgsetmat(1,1, 0.0_ki             - msq6(1) - msq6(1)  )
  call pgsetmat(1,2, scalar(p2,p2)      - msq6(1) - msq6(2)  )
  call pgsetmat(1,3, scalar(p23,p23)    - msq6(1) - msq6(3)  )
  call pgsetmat(1,4, scalar(p234,p234)  - msq6(1) - msq6(4)  )
  call pgsetmat(1,5, scalar(p61,p61)    - msq6(1) - msq6(5)  )
  call pgsetmat(1,6, scalar(p1,p1)      - msq6(1) - msq6(6)  )
  !
  call pgsetmat(2,1, pggetmat(1,2))
  call pgsetmat(2,2, 0.0_ki             - msq6(2) - msq6(2)  )
  call pgsetmat(2,3, scalar(p3,p3)      - msq6(2) - msq6(3)  )
  call pgsetmat(2,4, scalar(p34,p34)    - msq6(2) - msq6(4)  )
  call pgsetmat(2,5, scalar(p345,p345)  - msq6(2) - msq6(5)  )
  call pgsetmat(2,6, scalar(p12,p12)    - msq6(2) - msq6(6)  )
  !
  call pgsetmat(3,1, pggetmat(1,3))
  call pgsetmat(3,2, pggetmat(2,3))
  call pgsetmat(3,3, 0.0_ki             - msq6(3) - msq6(3)  )
  call pgsetmat(3,4, scalar(p4,p4)      - msq6(3) - msq6(4)  )
  call pgsetmat(3,5, scalar(p45,p45)    - msq6(3) - msq6(5)  )
  call pgsetmat(3,6, scalar(p123,p123)  - msq6(3) - msq6(6)  )
  !
  call pgsetmat(4,1, pggetmat(1,4))
  call pgsetmat(4,2, pggetmat(2,4))
  call pgsetmat(4,3, pggetmat(3,4))
  call pgsetmat(4,4, 0.0_ki             - msq6(4) - msq6(4)  )
  call pgsetmat(4,5, scalar(p5,p5)      - msq6(4) - msq6(5)  )
  call pgsetmat(4,6, scalar(p56,p56)    - msq6(4) - msq6(6)  )
  !
  call pgsetmat(5,1, pggetmat(1,5))
  call pgsetmat(5,2, pggetmat(2,5))
  call pgsetmat(5,3, pggetmat(3,5))
  call pgsetmat(5,4, pggetmat(4,5))
  call pgsetmat(5,5, 0.0_ki             - msq6(5) - msq6(5)  )
  call pgsetmat(5,6, scalar(p6,p6)      - msq6(5) - msq6(6)  )
  !
  call pgsetmat(6,1, pggetmat(1,6))
  call pgsetmat(6,2, pggetmat(2,6))
  call pgsetmat(6,3, pggetmat(3,6))
  call pgsetmat(6,4, pggetmat(4,6))
  call pgsetmat(6,5, pggetmat(5,6))
  call pgsetmat(6,6, 0.0_ki             - msq6(6) - msq6(6)  )
  !
  ! This call fills the internal array s_mat_r.
  ! It also assigns the integers in s_mat_p, which encode the positions
  ! of complex mass entries and zero mass entries. It includes call to init_invs
  !
  call pgpreparesmatrix()
  !
  write (*,*) 'Choose what the program has to compute:'
  write (*,*) '0) missing'
  write (*,*) '1) missing '
  write (*,*) '2) form factor for diagram where propagator 3 is pinched, rank 0'
  write (*,*) '3) form factor for diagram where props. 2,5 are pinched, rank 0'
  write (*,*) '4) form factor for diagram where props. 2,4,6 are pinched, rank 0'
  write (*,*) '5) test table'
  read (*,*) choix
  write (*,*) 'The result has been written to the file test6point.txt'
  !
  !
  call cpu_time(t1)
  !
  !
  ! To change the value of mu^2 (in GeV) (set to 1. by default)
  ! uncomment this line
  !
  !mu2_scale_par = 12._ki
  !
  write (17,*) 'The kinematics is:'
  write (17,*) ''
  write (17,*) '           p6          p5      '
  write (17,*) '            \             /        '
  write (17,*) '             \   (5)    /         '
  write (17,*) '              /---<---\         '
  write (17,*) '        (6) /            \ (4)       '
  write (17,*) ' p1 ____/              \____ p4 '
  write (17,*) '           \               /        '
  write (17,*) '       (1) \             / (3)     '
  write (17,*) '             \--->---/          '
  write (17,*) '            /   (2)     \         '
  write (17,*) '           /              \        '
  write (17,*) '           p2            p3      '
  write (17,*) ''
  write (17,*) ' S(1,3) = (p2+p3)^2=',pggetmat(1,3)
  write (17,*) ' S(1,4) = (p2+p3+p4)^2=',pggetmat(1,4)
  write (17,*) ' S(1,5) = (p1+p6)^2=',pggetmat(1,5)
  write (17,*) ' S(2,4) = (p3+p4)^2=',pggetmat(2,4)
  write (17,*) ' S(2,5) = (p3+p4+p5)^2=',pggetmat(2,5)
  write (17,*) ' S(2,6) = (p1+p2)^2=',pggetmat(2,6)
  write (17,*) ' S(3,5) = (p4+p5)^2=',pggetmat(3,5)
  write (17,*) ' S(3,6) = (p4+p5+p6)^2=',pggetmat(3,6)
  write (17,*) ' S(4,6) = (p5+p6)^2=',pggetmat(4,6)
  write (17,*) ' S(1,2) = p2^2=',pggetmat(1,2)
  write (17,*) ' S(2,3) = p3^2=',pggetmat(2,3)
  write (17,*) ' S(3,4) = p4^2=',pggetmat(3,4)
  write (17,*) ' S(4,5) = p5^2=',pggetmat(4,5)
  write (17,*) ' S(5,6) = p6^2=',pggetmat(5,6)
  write (17,*) ' S(1,6) = p1^2=',pggetmat(6,1)
  write (17,*) '(mu)^2 =',pggetmusq()
  write (17,*) ''
  !
  if (choix == 0) then
    !
    ! form factor for six-point function, rank 0
    !
!     res6 =  a60(s_null)
    !
  else if (choix == 1) then
    !
    ! form factor for six-point function, rank 4
    !
!       res6 = a64(1,1,2,3,s_null)
    !
  else if (choix == 2) then
    !
    ! form factor for pinched diagram, rank 0
    ! the propagator 3 is pinched
    !
    !
    ! p1            p6
    !  \           /
    !   \   (6)   /
    !    |---<---\ (5)
    !    |        \
    !    |         \____ p5
    ! (1)|         /
    !    |        / (4)
    !    |--->---/
    !   /   (2)   \
    !  /           \
    ! p2            p3+p4
    !
    !
    write (17,*) 'Since the propagator 3 is pinched'
    write (17,*) 'the reduced kinematics is:'
    write (17,*) ''
    write (17,*) ' p1            p6     '
    write (17,*) '   \              /       '
    write (17,*) '    \   (6)     /        '
    write (17,*) '     |---<---\ (5)     '
    write (17,*) '     |           \        '
    write (17,*) '     |            \____ p5'
    write (17,*) '(1) |            /       '
    write (17,*) '     |           / (4)    '
    write (17,*) '     |--->---/         '
    write (17,*) '    /   (2)     \        '
    write (17,*) '   /              \       '
    write (17,*) ' p2            p3+p4  '
    !
    !
    b_pin = packb((/3/))
    res6c = pga50(b_pin, 0)
    res6b = pga50(b_pin, 1)
    res6a = pga50(b_pin, 2)
    !
  else if (choix == 3) then
    !
    ! form factor for pinched diagram, rank 0
    ! the propagators 2 and 5 are pinched
    !
    !
    ! p1            p5+p6
    !  \           /
    !   \   (6)   /
    !    |---<---|
    !    |       |
    ! (1)|       |(4)
    !    |       |
    !    |--->---|
    !   /   (3)   \
    !  /           \
    ! p2+p3          p4
    !
    write (17,*) 'Since the propagators 2 and 5 are pinched'
    write (17,*) 'the reduced kinematics is:'
    write (17,*) ''
    write (17,*) '   p1              p5+p6'
    write (17,*) '     \                /     '
    write (17,*) '      \     (6)     /      '
    write (17,*) '       |---<-----|       '
    write (17,*) '       |            |       '
    write (17,*) ' (1)  |            | (4)   '
    write (17,*) '       |            |       '
    write (17,*) '       |--->-----|       '
    write (17,*) '      /     (3)     \      '
    write (17,*) '     /                \     '
    write (17,*) '  p2+p3          p4   '
    write (17,*) ''
    !
    !
    b_pin = packb((/2,5/))
    res6c = pga40(b_pin,0)
    res6b = pga40(b_pin,1)
    res6a = pga40(b_pin,2)
    !
  else if (choix == 4) then
    !
    ! form factor for pinched diagram, rank 0
    ! the propagator 2, 4 and 6 are pinched
    !
    !           |
    !           | p1+p6
    !           |
    !           /\
    !          /  \
    !     (1) /    \ (5)
    !        /      \
    !       /--->----\
    !      /   (3)    \
    !     /            \ 
    ! p2+p3            p4+p5
    !
    write (17,*) 'Since the propagators 2, 4 and 6 are pinched'
    write (17,*) 'the reduced kinematics is:'
    write (17,*) ''
    write (17,*) '            |             '
    write (17,*) '            | p1+p6       '
    write (17,*) '            |             '
    write (17,*) '            /\            '
    write (17,*) '           /  \           '
    write (17,*) '     (1) /     \ (5)      '
    write (17,*) '        /        \         '
    write (17,*) '       /--->----\        '
    write (17,*) '      /   (3)      \       '
    write (17,*) '     /               \      '
    write (17,*) ' p2+p3            p4+p5  '
    write (17,*) ''
    !
    !
    b_pin = packb((/2,4,6/))
    res6c = pga30(b_pin,0)
    res6b = pga30(b_pin,1)
    res6a = pga30(b_pin,2)
  else if (choix == 5) then
    do s1=1,6
      b_pin = packb((/s1/))
      write(17,*) " PENTAGON s=",s1
      write(17,*) " Rank 0 "
      write(17,'("E0( 0) = (",D23.15," +I* ",D23.15,")")') pga50(b_pin, 0)
      write(17,'("E0(-1) = (",D23.15," +I* ",D23.15,")")') pga50(b_pin, 1)
      write(17,'("E0(-2) = (",D23.15," +I* ",D23.15,")")') pga50(b_pin, 2)
      write(17,*) " Rank 1 "
      do i1=1,6
        write(17,'("E",I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1, pga51(i1, b_pin, 0)
        write(17,'("E",I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1, pga51(i1, b_pin, 1)
        write(17,'("E",I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1, pga51(i1, b_pin, 2)
      enddo
      write(17,*) " Rank 2 "
      write(17,'("E00( 0) = (",D23.15," +I* ",D23.15,")")') pgb52( b_pin, 0)
      write(17,'("E00(-1) = (",D23.15," +I* ",D23.15,")")') pgb52( b_pin, 1)
      write(17,'("E00(-2) = (",D23.15," +I* ",D23.15,")")') pgb52( b_pin, 2)
      do i1=1,6
        do i2=i1,6
          write(17,'("E",2I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,i2, pga52(i1,i2, b_pin, 0)
          write(17,'("E",2I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,i2, pga52(i1,i2, b_pin, 1)
          write(17,'("E",2I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,i2, pga52(i1,i2, b_pin, 2)
        enddo
      enddo
      write(17,*) " Rank 3 "
      do i1=1,6
        write(17,'("E00",I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1, pgb53(i1, b_pin, 0)
        write(17,'("E00",I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1, pgb53(i1, b_pin, 1)
        write(17,'("E00",I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1, pgb53(i1, b_pin, 2)
      enddo
      do i1=1,6
        if (i1.eq.s1) cycle
        do i2=i1,6
          if (i2.eq.s1) cycle
          do i3=i2,6
            if (i3.eq.s1) cycle
            write(17,'("E",3I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3, pga53(i1,i2,i3, b_pin, 0)
            write(17,'("E",3I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3, pga53(i1,i2,i3, b_pin, 1)
            write(17,'("E",3I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3, pga53(i1,i2,i3, b_pin, 2)
          enddo
        enddo
      enddo
      write(17,*) " Rank 4 "
      write(17,'("E0000( 0) = (",D23.15," +I* ",D23.15,")")') pgc54( b_pin, 0)
      write(17,'("E0000(-1) = (",D23.15," +I* ",D23.15,")")') pgc54( b_pin, 1)
      write(17,'("E0000(-2) = (",D23.15," +I* ",D23.15,")")') pgc54( b_pin, 2)
      do i1=1,6
        if (i1.eq.s1) cycle
        do i2=i1,6
          if (i2.eq.s1) cycle
          write(17,'("E00",2I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,i2, pgb54(i1,i2, b_pin, 0)
          write(17,'("E00",2I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,i2, pgb54(i1,i2, b_pin, 1)
          write(17,'("E00",2I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,i2, pgb54(i1,i2, b_pin, 2)
        enddo
      enddo
      do i1=1,6
        if (i1.eq.s1) cycle
        do i2=i1,6
          if (i2.eq.s1) cycle
          do i3=i2,6
            if (i3.eq.s1) cycle
            do i4=i3,6
              if (i4.eq.s1) cycle
              write(17,'("E",4I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3,i4, pga54(i1,i2,i3,i4, b_pin, 0)
              write(17,'("E",4I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3,i4, pga54(i1,i2,i3,i4, b_pin, 1)
              write(17,'("E",4I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3,i4, pga54(i1,i2,i3,i4, b_pin, 2)
            enddo
          enddo
        enddo
      enddo
      write(17,*) " Rank 5 "
      do i1=1,6
        if (i1.eq.s1) cycle
        write(17,'("E0000",I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,pgc55(i1, b_pin, 0)
        write(17,'("E0000",I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,pgc55(i1, b_pin, 1)
        write(17,'("E0000",I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,pgc55(i1, b_pin, 2)
      enddo
      do i1=1,6
        if (i1.eq.s1) cycle
        do i2=i1,6
          if (i2.eq.s1) cycle
          do i3=i2,6
            if (i3.eq.s1) cycle
            write(17,'("E00",3I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3, pgb55(i1,i2,i3, b_pin, 0)
            write(17,'("E00",3I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3, pgb55(i1,i2,i3, b_pin, 1)
            write(17,'("E00",3I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3, pgb55(i1,i2,i3, b_pin, 2)
          enddo
        enddo
      enddo
      do i1=1,6
        if (i1.eq.s1) cycle
        do i2=i1,6
          if (i2.eq.s1) cycle
          do i3=i2,6
            if (i3.eq.s1) cycle
            do i4=i3,6
              if (i4.eq.s1) cycle
              do i5=i4,6
                if (i5.eq.s1) cycle
                write(17,'("E",5I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3,i4,i5, pga55(i1,i2,i3,i4,i5, b_pin, 0)
                write(17,'("E",5I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3,i4,i5, pga55(i1,i2,i3,i4,i5, b_pin, 1)
                write(17,'("E",5I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3,i4,i5, pga55(i1,i2,i3,i4,i5, b_pin, 2)
              enddo
            enddo
          enddo
        enddo
      enddo
! --------------------------------------------------------------
! --------------------------------------------------------------
      do s2=s1+1,6
        b_pin=packb( (/s1,s2/) )
        write(17,*) " BOX s,t = ",s1,s2
        write(17,*) " Rank 0 "
        write(17,'("D0( 0) = (",D23.15," +I* ",D23.15,")")') pga40(b_pin, 0)
        write(17,'("D0(-1) = (",D23.15," +I* ",D23.15,")")') pga40(b_pin, 1)
        write(17,'("D0(-2) = (",D23.15," +I* ",D23.15,")")') pga40(b_pin, 2)
        write(17,*) " Rank 1 "
        do i1=1,6
          write(17,'("D",I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1, pga41(i1, b_pin, 0)
          write(17,'("D",I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1, pga41(i1, b_pin, 1)
          write(17,'("D",I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1, pga41(i1, b_pin, 2)
        enddo
        ! --------------------------------
        write(17,*) " Rank 2 "
        write(17,'("D00( 0) = (",D23.15," +I* ",D23.15,")")') pgb42(b_pin, 0)
        write(17,'("D00(-1) = (",D23.15," +I* ",D23.15,")")') pgb42(b_pin, 1)
        write(17,'("D00(-2) = (",D23.15," +I* ",D23.15,")")') pgb42(b_pin, 2)
        do i1=1,6
          if (i1.eq.s1) cycle
          do i2=i1,6
            if (i2.eq.s1) cycle
            write(17,'("D",2I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,i2, pga42(i1,i2, b_pin, 0)
            write(17,'("D",2I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,i2, pga42(i1,i2, b_pin, 1)
            write(17,'("D",2I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,i2, pga42(i1,i2, b_pin, 2)
          enddo
        enddo
        ! --------------------------------
        write(17,*) " Rank 3 "
        do i1=1,6
          write(17,'("D00",I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1, pgb43(i1, b_pin, 0)
          write(17,'("D00",I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1, pgb43(i1, b_pin, 1)
          write(17,'("D00",I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1, pgb43(i1, b_pin, 2)
        enddo
        do i1=1,6
          if (i1.eq.s1) cycle
          do i2=i1,6
            if (i2.eq.s1) cycle
            do i3=i2,6
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
        do i1=1,6
          if (i1.eq.s1) cycle
          do i2=i1,6
            if (i2.eq.s1) cycle
            write(17,'("D00",2I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,i2, pgb44(i1,i2, b_pin, 0)
            write(17,'("D00",2I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,i2, pgb44(i1,i2, b_pin, 1)
            write(17,'("D00",2I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,i2, pgb44(i1,i2, b_pin, 2)
          enddo
        enddo
        do i1=1,6
          if (i1.eq.s1) cycle
          do i2=i1,6
            if (i2.eq.s1) cycle
            do i3=i2,6
              if (i3.eq.s1) cycle
              do i4=i3,6
                if (i4.eq.s1) cycle
                write(17,'("D",4I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3,i4, pga44(i1,i2,i3,i4, b_pin, 0)
                write(17,'("D",4I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3,i4, pga44(i1,i2,i3,i4, b_pin, 1)
                write(17,'("D",4I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3,i4, pga44(i1,i2,i3,i4, b_pin, 2)
              enddo
            enddo
          enddo
        enddo
      ! --------------------------------------------------------------
        do s3=(s2+1),6
          b_pin=packb( (/s1,s2,s3/) )
          write(17,*) " TRIANGLE s,t,u = ",s1,s2,s3
          write(17,*) " Rank 0 "
          write(17,'("C0( 0) = (",D23.15," +I* ",D23.15,")")') pga30(b_pin, 0)
          write(17,'("C0(-1) = (",D23.15," +I* ",D23.15,")")') pga30(b_pin, 1)
          write(17,'("C0(-2) = (",D23.15," +I* ",D23.15,")")') pga30(b_pin, 2)
          write(17,*) " Rank 1 "
          do i1=1,6
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
          do i1=1,6
            if ((i1.eq.s1).or.(i1.eq.s2)) cycle
            do i2=i1,6
              if ((i2.eq.s1).or.(i2.eq.s2)) cycle
              write(17,'("C",2I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,i2, pga32(i1,i2, b_pin, 0)
              write(17,'("C",2I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,i2, pga32(i1,i2, b_pin, 1)
              write(17,'("C",2I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,i2, pga32(i1,i2, b_pin, 2)
            enddo
          enddo
          ! --------------------------------
          write(17,*) " Rank 3 "
          do i1=1,6
            if ((i1.eq.s1).or.(i1.eq.s2)) cycle
            write(17,'("C00",I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1, pgb33(i1, b_pin, 0)
            write(17,'("C00",I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1, pgb33(i1, b_pin, 1)
            write(17,'("C00",I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1, pgb33(i1, b_pin, 2)
          enddo
          do i1=1,6
            if ((i1.eq.s1).or.(i1.eq.s2)) cycle
            do i2=i1,6
              if ((i2.eq.s1).or.(i2.eq.s2)) cycle
              do i3=i2,6
                if ((i3.eq.s1).or.(i3.eq.s2)) cycle
                write(17,'("C",3I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3, pga33(i1,i2,i3, b_pin, 0)
                write(17,'("C",3I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3, pga33(i1,i2,i3, b_pin, 1)
                write(17,'("C",3I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3, pga33(i1,i2,i3, b_pin, 2)
              enddo
            enddo
          enddo
        ! --------------------------------------------------------------
          do s4=(s3+1),6
            b_pin=packb( (/s1,s2,s3,s4/) )
            write(17,*) " BUBBLES s,t,u,v = ",s1,s2,s3,s4
            write(17,*) " Rank 0 "
            write(17,'("B0( 0) = (",D23.15," +I* ",D23.15,")")') pga20(b_pin, 0)
            write(17,'("B0(-1) = (",D23.15," +I* ",D23.15,")")') pga20(b_pin, 1)
            write(17,'("B0(-2) = (",D23.15," +I* ",D23.15,")")') pga20(b_pin, 2)
            write(17,*) " Rank 1 "
            do i1=1,6
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
            do i1=1,6
              if ((i1.eq.s1).or.(i1.eq.s2).or.(i1.eq.s3)) cycle
              do i2=i1,6
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
    enddo
    !
  end if
  ! 
  call cpu_time(t2)
  !
  write (17,*) ''
  write (17,*) 'normalisation:'
  write (17,*) 'defining I_N^n= mu^(4-n) \int d^n k/(i*Pi^(n/2))*func(k,p_i)'
  write (17,*) '= r_Gam *(P2/eps^2+P1/eps+P0),'
  write (17,*) 'n = 4-2*eps,'
  write (17,*) 'r_Gam = Gamma(1+eps)*Gamma(1-eps)^2/Gamma(1-2eps)'
  write (17,*) 'the program gives numbers for P2,P1,P0'
  write (17,*) ''
  write (17,'("  1/epsilon^2 * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6a,ki),aimag(res6a)
  write (17,'("+ 1/epsilon   * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6b,ki),aimag(res6b)
  write (17,'("+ 1           * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6c,ki),aimag(res6c)
  write (17,*) ''
  !
  write (17,*) 'CPU time=',t2-t1
  !
  !
  ! routine to free the cache and allocated memory
  !
!   call exitgolem95()
  !
  close(17)
  close(19)
  !
end program main
