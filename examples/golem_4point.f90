!
! This program computes four-point functions
!
program main
  !
  use pjfry95pg
  use golem_demo
  !
  implicit none
  !
  integer, parameter :: ki = kind(1.0d0)
!   type(form_factor) :: res6,res6a
  complex(ki) :: res6a,res6b,res6c
  complex(ki) :: res7a,res7b,res7c
  real(ki) :: ti1,ti2,zero
  real(ki) :: p1sq,p2sq,p3sq,p4sq,m1sq,m2sq,m3sq,m4sq
  real(ki) :: s_var,t_var
  integer :: choix,choix_kinem
  integer :: i1,i2,i3,i4,s1,s2,s_null,b_pin
  real(ki) :: mu2,mu02,lmu2
  !
  zero = 0.0_ki
  !
  write (*,*) 'Choose from the following kinematics:'
  write (*,*) '1) all external legs light-like, no internal masses'
  write (*,*) '2) one external leg not light-like, no internal masses'
  write (*,*) '3) two opposite external legs not light-like, no internal masses'
  write (*,*) '4) two adjacent external legs not light-like, no internal masses'
  write (*,*) '5) three external legs not light-like, no internal masses'
  write (*,*) '6) four external legs not light-like, no internal masses'
  write (*,*) '7) two external legs not light-like, internal masses, IR divergent' 
  !(example QL Box 8)
  write (*,*) '8) three external legs not light-like, internal masses, IR divergent'
  !(example QL Box 13)
  read (*,*) choix_kinem
  !
  ! Opening of the files containing the results
  !
  open(unit=17,file='test4point.txt',status='unknown')
  !
  !
  ! These are the entries of the S matrix
  ! They are related to the cuts of the following diagram
  ! S(1,3) = (p1+p4)^2 = t_var
  ! S(2,4) = (p1+p2)^2 = s_var
  ! S(1,2) = p2^2 = p2sq
  ! S(2,3) = p3^2 = p3sq
  ! S(3,4) = p4^2 = p4sq
  ! S(4,1) = p1^2 = p1sq
  !
  !  p1            p4
  !    \           /
  !     \   (4)   /
  !      |---<---|
  !      |       |
  ! (1)  |       | (3)
  !      |       |
  !      |--->---|
  !     /   (2)   \
  !    /           \
  !  p2             p3
  !
  !
  !
  call pginitgolem95(4)
  s_null=0
  !
  m1sq = zero
  m2sq = zero
  m3sq = zero
  m4sq = zero
  !
  if (choix_kinem == 1) then
    !
    p1sq = 0.0_ki
    p2sq = 0.0_ki
    p3sq = 0.0_ki
    p4sq = 0.0_ki
    !
  else if (choix_kinem == 2) then
    !
    p1sq = 0.0_ki
    p2sq = 0.0_ki
    p3sq = 0.0_ki
    p4sq = 60.0_ki
    !
  else if (choix_kinem == 3) then
    !
    p1sq = 0.0_ki
    p2sq = 50.0_ki
    p3sq = 0.0_ki
    p4sq = 60.0_ki
    !
  else if (choix_kinem == 4) then
    !
    p1sq = 0.0_ki
    p2sq = 0.0_ki
    p3sq = 50.0_ki
    p4sq = 60.0_ki
    !
  else if (choix_kinem == 5) then
    !
    p1sq = 0.0_ki
    p2sq = 50.0_ki
    p3sq = 80.0_ki
    p4sq = 60.0_ki
    !
  else if (choix_kinem == 6) then
    !
    p1sq = 20.0_ki
    p2sq = 50.0_ki
    p3sq = 80.0_ki
    p4sq = 60.0_ki
    !
  else if (choix_kinem == 7) then
    !
    p1sq = 0.0_ki
    p2sq = 0.0_ki
    p3sq = -80.0_ki
    p4sq = -60.0_ki
    m1sq=0.0_ki
    m2sq=0.0_ki
    m3sq=20.0_ki
    m4sq=0.0_ki
    !
  else if (choix_kinem == 8) then
    !
    p1sq = 0.0_ki
    p2sq = 50.0_ki
    p3sq = -80.0_ki
    p4sq = -60.0_ki
    m1sq=0.0_ki
    m2sq=5.0_ki
    m3sq=10.0_ki
    m4sq=0.0_ki
    !
  end if
  !
  s_var = 200.0_ki
  t_var = -123.0_ki
  !
  ! Definition of the S matrix
  !
  call pgsetmat(1,1 , -2._ki*m1sq          )
  call pgsetmat(1,2 , p2sq-m1sq-m2sq       )
  call pgsetmat(1,3 , t_var-m1sq-m3sq      )
  call pgsetmat(1,4 , p1sq-m1sq-m4sq       )
  !
  call pgsetmat(2,1 , pggetmat(1,2)           )
  call pgsetmat(2,2 , -2._ki*m2sq          )
  call pgsetmat(2,3 , p3sq-m2sq-m3sq       )
  call pgsetmat(2,4 , s_var-m2sq-m4sq      )
  !
  call pgsetmat(3,1 , pggetmat(1,3)           )
  call pgsetmat(3,2 , pggetmat(2,3)           )
  call pgsetmat(3,3 , -2._ki*m3sq          )
  call pgsetmat(3,4 , p4sq-m3sq-m4sq       )
  !
  call pgsetmat(4,1 , pggetmat(1,4)           )
  call pgsetmat(4,2 , pggetmat(2,4)           )
  call pgsetmat(4,3 , pggetmat(3,4)           )
  call pgsetmat(4,4 , -2._ki*m4sq          )
  !
  ! This call initialize cache with kinematics above
  !
!   call pgsetmusq(100._ki)
  !
  call pgpreparesmatrix()
  !
  write (*,*) 'Choose what the program should compute:'
  write (*,*) '0) scalar four-point function in n dimensions'
  write (*,*) '1) four-point function in n dimensions with one Feynman parameter (z1)'
  write (*,*) '2) four-point function in n dimensions with two Feynman parameters (z1*z4)'
  write (*,*) '3) four-point function in n dimensions with three Feynman parameters (z1^2*z3)'
  write (*,*) '4) four-point fctn. in n dimensions with four Feynman parameters (z1*z2*z3*z4)'
  write (*,*) '5) missing'
  write (*,*) '6) missing'
  write (*,*) '7) missing'
  write (*,*) '8) the mu dependence'
  write (*,*) '9) Test table'
  read (*,*) choix
  !  
  if (choix == 0) then
    !
    write (*,*) 'calculating n-dim scalar box with '
    !
  else if   (choix == 1) then
    !
    write (*,*) 'calculating n-dim rank one (z1) 4-point fctn. with'
    !
  else if   (choix == 2) then
    !
    write (*,*) 'calculating n-dim rank two (z1*z4) 4-point fctn. with'
    !
  else if   (choix == 3) then
    !
    write (*,*) 'calculating n-dim rank three (z1^2*z3) 4-point fctn. with'
    !
  else if   (choix == 4) then
    !
    write (*,*) 'calculating n-dim rank four (z1*z2*z3*z4) 4-point fctn. with'
    !
  else if   (choix == 5) then
    !
    write (*,*) 'calculating (n+2)-dim scalar 4-point fctn. with'
    !
  else if   (choix == 6) then
    !
    write (*,*) 'calculating (n+2)-dim rank two (z1*z2) 4-point fctn. with'
    !
  else if (choix == 7) then
    !
    write (*,*) 'calculating (n+4)-dim scalar 4-point function with'
    !
  else if (choix == 8) then
    !
    write (*,*) 'example of a test on the renormalisation scale mu dependence'
    !
  else if (choix == 9) then
    !
    write (*,*) 'generating test table'
    !
  end if
  !
  if (choix_kinem == 1) then
    !
    write (*,*) 'no off-shell leg'
    !
  else if (choix_kinem == 2) then
    !
    write (*,*) 'one off-shell leg'
    !
  else if (choix_kinem == 3) then
    !
    write (*,*) 'two opposite off-shell legs'
    !
  else if (choix_kinem == 4) then
    !
    write (*,*) 'two adjacent off-shell legs'
    !
  else if (choix_kinem == 5) then
    !
    write (*,*) 'three off-shell legs'
    !
  else if (choix_kinem == 6) then
    !
    write (*,*) 'four off-shell legs'
    !
  else if (choix_kinem == 7) then
    !
    write (*,*) 'two off-shell legs and one internal mass'
    !
  else if (choix_kinem == 8) then
    !
    write (*,*) 'three off-shell legs and two internal masses'
    !
  end if
  !
  write (*,*) 'The result has been written to the file test4point.txt'
  !
  ! 
  call cpu_time(ti1)
  !
  ! To change the value of mu^2 (in GeV) (set to 1. by default)
  ! uncomment this line
  !
  !call pgsetmusq(12._ki)
  mu02 = pggetmusq()
  !
  if (choix == 0) then
    !
    ! Result for the scalar integral in n dimension
    !
    res6c =  pga40(s_null, 0)
    res6b =  pga40(s_null, 1)
    res6a =  pga40(s_null, 2)
    !
    ! the labels 1,2,..,4 correspond to Feynman parameter z1,z2,...,z4
  else if (choix == 1) then
    !
    ! Results for integrals in n dimension with one Feynman parameter 
    ! at the numerator: z1
    !
    res6c =  - pga41(1, s_null, 0)
    res6b =  - pga41(1, s_null, 1)
    res6a =  - pga41(1, s_null, 2)
    !
  else if (choix == 2) then
    !
    ! Results for integrals in n dimension with two Feynman parameters 
    ! at the numerator: z1*z4
    !
    res6c =  pga42(1,3, s_null, 0)
    res6b =  pga42(1,3, s_null, 1)
    res6a =  pga42(1,3, s_null, 2)
    !
  else if (choix == 3) then
    !
    ! Results for integrals in n dimension with three Feynman parameters 
    ! at the numerator: z1^2*z3
    !
    res6c =  -pga43(1,1,3, s_null, 0)
    res6b =  -pga43(1,1,3, s_null, 1)
    res6a =  -pga43(1,1,3, s_null, 2)
    !
  else if (choix == 4) then
    !
    ! Results for integrals in n dimension with four Feynman parameters 
    ! at the numerator: z1*z2*z3*z4
    !
    res6c =  pga44(1,2,3,3, s_null, 0)
    res6b =  pga44(1,2,3,3, s_null, 1)
    res6a =  pga44(1,2,3,3, s_null, 2)
    ! 
  else if (choix == 5) then
    ! 
  else if (choix == 6) then
    !
  else if (choix == 7) then
    !
  else if (choix == 8) then
    !
    ! By default, mu2_scale_par = 1._ki
    !
    b_pin = packb( (/3/) )
    res6c =  pgc44(s_null, 0)+pgb33(2, b_pin, 0)
    res6b =  pgc44(s_null, 1)+pgb33(2, b_pin, 1)
    res6a =  pgc44(s_null, 2)+pgb33(2, b_pin, 2)
    ! 
    !
    mu2 = 34._ki
    lmu2 = log(mu2/pggetmusq())
    ! we change the value of mu^2
    call pgsetmusq(mu2)
    ! note that we reset the cache with setting new value of musq
    ! we need to call preparesmatrix to refresh cache
    call pgpreparesmatrix()
    !
    res7c =  pgc44(s_null, 0)+pgb33(2, b_pin, 0)
    res7b =  pgc44(s_null, 1)+pgb33(2, b_pin, 1)
    res7a =  pgc44(s_null, 2)+pgb33(2, b_pin, 2)
  !
    !
  else if (choix == 9) then
    !
    write(17,*) " BOX "
    write(17,*) " Rank 0 "
    write(17,'("D0( 0) = (",D23.15," +I* ",D23.15,")")') pga40( s_null, 0)
    write(17,'("D0(-1) = (",D23.15," +I* ",D23.15,")")') pga40( s_null, 1)
    write(17,'("D0(-2) = (",D23.15," +I* ",D23.15,")")') pga40( s_null, 2)
    write(17,*) " Rank 1 "
    do i1=1,4
      write(17,'("D",I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1, pga41(i1, s_null, 0)
      write(17,'("D",I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1, pga41(i1, s_null, 1)
      write(17,'("D",I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1, pga41(i1, s_null, 2)
    enddo
    ! --------------------------------
    write(17,*) " Rank 2 "
    write(17,'("D00( 0) = (",D23.15," +I* ",D23.15,")")') pgb42( s_null, 0)
    write(17,'("D00(-1) = (",D23.15," +I* ",D23.15,")")') pgb42( s_null, 1)
    write(17,'("D00(-2) = (",D23.15," +I* ",D23.15,")")') pgb42( s_null, 2)
    do i1=1,4
      do i2=i1,4
        write(17,'("D",2I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,i2, pga42(i1,i2, s_null, 0)
        write(17,'("D",2I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,i2, pga42(i1,i2, s_null, 1)
        write(17,'("D",2I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,i2, pga42(i1,i2, s_null, 2)
      enddo
    enddo
    ! --------------------------------
    write(17,*) " Rank 3 "
    do i1=1,4
      write(17,'("D00",I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1, pgb43(i1, s_null, 0)
      write(17,'("D00",I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1, pgb43(i1, s_null, 1)
      write(17,'("D00",I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1, pgb43(i1, s_null, 2)
    enddo
    do i1=1,4
      do i2=i1,4
        do i3=i2,4
          write(17,'("D",3I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3, pga43(i1,i2,i3, s_null, 0)
          write(17,'("D",3I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3, pga43(i1,i2,i3, s_null, 1)
          write(17,'("D",3I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3, pga43(i1,i2,i3, s_null, 2)
        enddo
      enddo
    enddo
    ! --------------------------------
    write(17,*) " Rank 4 "
    write(17,'("D0000( 0) = (",D23.15," +I* ",D23.15,")")') pgc44( s_null, 0)
    write(17,'("D0000(-1) = (",D23.15," +I* ",D23.15,")")') pgc44( s_null, 1)
    write(17,'("D0000(-2) = (",D23.15," +I* ",D23.15,")")') pgc44( s_null, 2)
    do i1=1,4
      do i2=i1,4
        write(17,'("D00",2I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,i2, pgb44(i1,i2, s_null, 0)
        write(17,'("D00",2I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,i2, pgb44(i1,i2, s_null, 1)
        write(17,'("D00",2I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,i2, pgb44(i1,i2, s_null, 2)
      enddo
    enddo
    do i1=1,4
      do i2=i1,4
        do i3=i2,4
          do i4=i3,4
            write(17,'("D",4I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3,i4, pga44(i1,i2,i3,i4, s_null, 0)
            write(17,'("D",4I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3,i4, pga44(i1,i2,i3,i4, s_null, 1)
            write(17,'("D",4I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3,i4, pga44(i1,i2,i3,i4, s_null, 2)
          enddo
        enddo
      enddo
    enddo
  ! --------------------------------------------------------------
    do s1=1,4
      b_pin = packb( (/s1/) )
      write(17,*) " TRIANGLE s = ",s1
      write(17,*) " Rank 0 "
      write(17,'("C0( 0) = (",D23.15," +I* ",D23.15,")")') pga30(b_pin, 0)
      write(17,'("C0(-1) = (",D23.15," +I* ",D23.15,")")') pga30(b_pin, 1)
      write(17,'("C0(-2) = (",D23.15," +I* ",D23.15,")")') pga30(b_pin, 2)
      write(17,*) " Rank 1 "
      do i1=1,4
        if ((i1.eq.s1)) cycle
        write(17,'("C",I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1, pga31(i1, b_pin, 0)
        write(17,'("C",I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1, pga31(i1, b_pin, 1)
        write(17,'("C",I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1, pga31(i1, b_pin, 2)
      enddo
      ! --------------------------------
      write(17,*) " Rank 2 "
      write(17,'("C00( 0) = (",D23.15," +I* ",D23.15,")")') pgb32(b_pin, 0)
      write(17,'("C00(-1) = (",D23.15," +I* ",D23.15,")")') pgb32(b_pin, 1)
      write(17,'("C00(-2) = (",D23.15," +I* ",D23.15,")")') pgb32(b_pin, 2)
      do i1=1,4
        if ((i1.eq.s1)) cycle
        do i2=i1,4
          if ((i2.eq.s1)) cycle
          write(17,'("C",2I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,i2, pga32(i1,i2, b_pin, 0)
          write(17,'("C",2I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,i2, pga32(i1,i2, b_pin, 1)
          write(17,'("C",2I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,i2, pga32(i1,i2, b_pin, 2)
        enddo
      enddo
      ! --------------------------------
      write(17,*) " Rank 3 "
      do i1=1,4
        if ((i1.eq.s1)) cycle
        write(17,'("C00",I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1, pgb33(i1, b_pin, 0)
        write(17,'("C00",I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1, pgb33(i1, b_pin, 1)
        write(17,'("C00",I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1, pgb33(i1, b_pin, 2)
      enddo
      do i1=1,4
        if ((i1.eq.s1)) cycle
        do i2=i1,4
          if ((i2.eq.s1)) cycle
          do i3=i2,4
            if ((i3.eq.s1)) cycle
            write(17,'("C",3I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3, pga33(i1,i2,i3, b_pin, 0)
            write(17,'("C",3I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3, pga33(i1,i2,i3, b_pin, 1)
            write(17,'("C",3I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3, pga33(i1,i2,i3, b_pin, 2)
          enddo
        enddo
      enddo
    ! --------------------------------------------------------------
      do s2=(s1+1),4
        b_pin = packb( (/s1,s2/) )
        write(17,*) " BUBBLES s,t = ",s1,s2
        write(17,*) " Rank 0 "
        write(17,'("B0( 0) = (",D23.15," +I* ",D23.15,")")') pga20(b_pin, 0)
        write(17,'("B0(-1) = (",D23.15," +I* ",D23.15,")")') pga20(b_pin, 1)
        write(17,'("B0(-2) = (",D23.15," +I* ",D23.15,")")') pga20(b_pin, 2)
        write(17,*) " Rank 1 "
        do i1=1,4
          if ((i1.eq.s1).or.(i1.eq.s2)) cycle
          write(17,'("B",I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1, pga21(i1, b_pin, 0)
          write(17,'("B",I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1, pga21(i1, b_pin, 1)
          write(17,'("B",I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1, pga21(i1, b_pin, 2)
        enddo
        ! --------------------------------
        write(17,*) " Rank 2 "
        write(17,'("B00( 0) = (",D23.15," +I* ",D23.15,")")') pgb22(b_pin, 0)
        write(17,'("B00(-1) = (",D23.15," +I* ",D23.15,")")') pgb22(b_pin, 1)
        write(17,'("B00(-2) = (",D23.15," +I* ",D23.15,")")') pgb22(b_pin, 2)
        do i1=1,4
          if ((i1.eq.s1).or.(i1.eq.s2)) cycle
          do i2=i1,4
            if ((i2.eq.s1).or.(i2.eq.s2)) cycle
            write(17,'("B",2I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,i2, pga22(i1,i2, b_pin, 0)
            write(17,'("B",2I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,i2, pga22(i1,i2, b_pin, 1)
            write(17,'("B",2I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,i2, pga22(i1,i2, b_pin, 2)
          enddo
        enddo
      enddo
    ! --------------------------------------------------------------
    enddo
  end if
  !
  call cpu_time(ti2)
  !
  ! The results are written to the file with unit 17
  !
  write (17,*) 'The kinematics is:'
  write (17,*) ''
  write (17,*) '  p1              p4'
  write (17,*) '    \                /  '
  write (17,*) '     \     (4)     /   '
  write (17,*) '      |----<----|    '
  write (17,*) '      |            |    '
  write (17,*) ' (1) |            |(3) '
  write (17,*) '      |            |    '
  write (17,*) '      |---->----|    '
  write (17,*) '     /     (2)     \   '
  write (17,*) '    /                \  '
  write (17,*) '  p2              p3'
  write (17,*) ''
  write (17,*) 'p1+p2+p3+p4 = 0'
  write (17,*) ''
  write (17,*) '(p1+p2)^2 =',s_var
  write (17,*) '(p2+p3)^2 =',t_var
  write (17,*) '(p1)^2 =',p1sq
  write (17,*) '(p2)^2 =',p2sq
  write (17,*) '(p3)^2 =',p3sq
  write (17,*) '(p4)^2 =',p4sq
  write (17,*) 'm1^2 =',m1sq
  write (17,*) 'm2^2 =',m2sq
  write (17,*) 'm3^2 =',m3sq
  write (17,*) 'm4^2 =',m4sq
  write (17,*) '(mu)^2 =',pggetmusq()
  write (17,*) ''
  if (choix.ne.9) then
    write (17,'("  1/epsilon^2 * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6a,ki),aimag(res6a)
    write (17,'("+ 1/epsilon   * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6b,ki),aimag(res6b)
    write (17,'("+ 1           * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6c,ki),aimag(res6c)
    !
    write (6,*) ''
    write (6,'("  1/epsilon^2 * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6a,ki),aimag(res6a)
    write (6,'("+ 1/epsilon   * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6b,ki),aimag(res6b)
    write (6,'("+ 1           * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6c,ki),aimag(res6c)
    write (17,*) ''
  endif
  !
  if ( ((choix == 5) .or. (choix == 6) .or. (choix == 7)) .and. (choix_kinem < 7) ) then
    !
    write (17,*) 'Missing test'
    !
  end if
  if (choix == 8) then
    !
    write (17,*) 'The preceding result has been computed with mu^2=',mu02
    write (17,*) ' '
    write (17,*) 'Now setting by hand mu^2=',mu2
    write (17,*) 'and expanding (',mu2,'/',mu02,')^epsilon around epsilon=0'
    write (17,'("  1/epsilon^2 * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6a,ki),aimag(res6a)
    write (17,'("+ 1/epsilon   * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6b,ki)+ &
         &real(res6a,ki)*lmu2, aimag(res6b)+aimag(res6a)*lmu2
    write (17,'("+ 1           * (",e16.10,1x,"+ I*",1x,e16.10,")")') &
    &real(res6c,ki) + lmu2*real(res6b,ki) + lmu2**2*real(res6a,ki)/2._ki,&
    &aimag(res6c) + lmu2*aimag(res6b) + lmu2**2*aimag(res6a)/2._ki
    write (17,*) ''
    write (17,*) 'check with direct calculation using the global variable mu2_scale_par=',mu2
    write (17,'("  1/epsilon^2 * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res7a,ki),aimag(res7a)
    write (17,'("+ 1/epsilon   * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res7b,ki),aimag(res7b)
    write (17,'("+ 1           * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res7c,ki),aimag(res7c)
    ! ***********************************
    write (6,*) ' '
    write (6,*) 'The preceding result has been computed with mu^2=',mu02
    write (6,*) 'Now setting by hand mu^2=',mu2
    write (6,*) 'and expanding (',mu2,'/',mu02,')^epsilon around epsilon=0'
    write (6,'("  1/epsilon^2 * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6a,ki),aimag(res6a)
    write (6,'("+ 1/epsilon   * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6b,ki)+ &
         & real(res6a,ki)*lmu2, aimag(res6b)+aimag(res6a)*lmu2
    write (6,'("+ 1           * (",e16.10,1x,"+ I*",1x,e16.10,")")') &
    &real(res6c,ki) + lmu2*real(res6b,ki) + lmu2**2*real(res6a,ki)/2._ki,&
    &aimag(res6c) + lmu2*aimag(res6b) + lmu2**2*aimag(res6a)/2._ki
    !
    write (6,*) ''
    write (6,*) 'check with direct calculation using the global variable mu2_scale_par=',mu2
    write (6,'("  1/epsilon^2 * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res7a,ki),aimag(res7a)
    write (6,'("+ 1/epsilon   * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res7b,ki),aimag(res7b)
    write (6,'("+ 1           * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res7c,ki),aimag(res7c)
    !
    write (6,*) ''
  end if
  !
  write (17,*) 'CPU time=',ti2-ti1
  !
  !
  close(17)
  !
end program main
