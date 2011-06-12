! 
! this program computes three-point functions
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
  real(ki) :: mass_sq_1,mass_sq_2,mass_sq_3
  real(ki) :: mass_int_sq_1,mass_int_sq_2,mass_int_sq_3
  !
  zero = 0.0_ki
  !
  write (*,*) 'Choose from the following kinematics:'
  write (*,*) '1) one off-shell leg, no internal masses'
  write (*,*) '2) two off-shell legs, no internal masses'
  write (*,*) '3) three off-shell legs, no internal masses (finite)'
  write (*,*) '4) two off-shell legs, one internal mass (QL 3), '
  write (*,*) '5) one off-shell leg, one on-shell massive leg (one internal mass, QL 4)'
  write (*,*) '6) two on-shell massive legs (one internal mass, QL 5)'
  write (*,*) '7) one off-shell leg, two on-shell massive legs (two internal masses, QL 6)'
  read (*,*) choix_kinem
  !
  ! opening of the files containing the results
  !
  open(unit=17,file='test3point.txt',status='unknown')
  !
  ! These are the entries of the S matrix
  ! They are related to the cuts of the following diagram
  ! All the momenta are incoming : p1+p2+p3 = 0
  !
  !            |
  !            | p1
  !            |
  !           /\
  !          /  \
  !     (1) /    \ (3)
  !        /      \
  !       /--->----\
  !      /   (2)    \
  ! p2  /            \ p3
  !
  ! S(1,2) = p2^2
  ! S(2,3) = p3^2
  ! S(3,1) = p1^2
  !
  ! Allocates memory to store the set of initial propagators, the S matrix, 
  ! its inverse and the b coefficients.
  ! This call will allocate a derived type s_mat_p object.
  ! Includes calls to allocation_s and initializes the caching system
  !
  call pginitgolem95(3)
  s_null=0
  !
  !
  if (choix_kinem == 1) then
    !
    mass_int_sq_1 = 0.0_ki
    mass_int_sq_2 = 0._ki
    mass_int_sq_3 = 0.0_ki
    mass_sq_1 = 0.0_ki
    mass_sq_2 = 0.0_ki
    mass_sq_3 = -2.0_ki
    !
  else if (choix_kinem == 2) then
    !
    mass_int_sq_1 = 0.0_ki
    mass_int_sq_2 = 0._ki
    mass_int_sq_3 = 0.0_ki
    mass_sq_1 = 0.0_ki
    mass_sq_2 = 10.0_ki
    mass_sq_3 = -60.0_ki
    !
  else if (choix_kinem == 3) then
    !
    mass_int_sq_1 = 0.0_ki
    mass_int_sq_2 = 0._ki
    mass_int_sq_3 = 0.0_ki
    mass_sq_1 = -50.0_ki
    mass_sq_2 = 10.0_ki
    mass_sq_3 = -60.0_ki
    !
  ! case p1^2, p2^2 /= m1^2, internal line 1 massive (QL3)
  else if (choix_kinem == 4) then
    !
    mass_int_sq_1 = 20.0_ki
    mass_int_sq_2 = 0._ki
    mass_int_sq_3 = 0.0_ki
    mass_sq_1 = -123.0_ki
    mass_sq_2 = -60.0_ki
    mass_sq_3 = 0.0_ki
    !
  ! case p1^2 = m1^2, p2^2 /= m1^2, internal line 1 massive (QL4)
  else if (choix_kinem == 5) then
    !
    mass_int_sq_1 = 5.0_ki
    mass_int_sq_2 = 0._ki
    mass_int_sq_3 = 0.0_ki
    mass_sq_1 = mass_int_sq_1
    mass_sq_2 = 7.0_ki
    mass_sq_3 = 0.0_ki
    !
  ! case p1^2 = m1^2, p2^2 = m1^2, internal line 1 massive (QL5)
  else if (choix_kinem == 6) then
    !
    mass_int_sq_1 = 4.0_ki
    mass_int_sq_2 = 6.0_ki
    mass_int_sq_3 = 8.0_ki
    mass_sq_1 = 12.0_ki
    mass_sq_2 = 14.0_ki
    mass_sq_3 = 16.0_ki
!     mass_int_sq_1 = 5.0_ki
!     mass_int_sq_2 = 0.0_ki
!     mass_int_sq_3 = 0.0_ki
!     mass_sq_1 = mass_int_sq_1
!     mass_sq_2 = mass_int_sq_1
!     mass_sq_3 = 0.0_ki
    !
  ! case p2^2 = m1^2; p1^2 /= m1^2,m3^2,nonzero; p3^2=m3^2, internal lines 1,2 massive (QL6)
  else if (choix_kinem == 7) then
    !
    mass_int_sq_1 = 9.0_ki
    mass_int_sq_2 = 0.0_ki
    mass_int_sq_3 = 3.0_ki
    mass_sq_1 = 25.0_ki
    mass_sq_2 = mass_int_sq_1
    mass_sq_3 = mass_int_sq_3
    !
  ! case p2^2=m1^2; p1^2 /= m1^2,nonzero; p3^2=m1^2, internal lines 1,2 massive (QL6 with internal masses equal)
  else if (choix_kinem == 8) then
    !
    mass_int_sq_1 = 5.0_ki
    mass_int_sq_2 = 0.0_ki
    mass_int_sq_3 = mass_int_sq_1
    mass_sq_1 = 25.0_ki
    mass_sq_2 = mass_int_sq_1
    mass_sq_3 = mass_int_sq_3
    !
  end if
  !
  ! Definition of the S matrix
  !
  call pgsetmat(1,1,  -2.0_ki*mass_int_sq_1 )
  call pgsetmat(1,2,  mass_sq_2 - mass_int_sq_1 - mass_int_sq_2 )
  call pgsetmat(1,3,  mass_sq_1 - mass_int_sq_1 - mass_int_sq_3 )
  !
  call pgsetmat(2,1,  pggetmat(1,2) )
  call pgsetmat(2,2,  -2.0_ki*mass_int_sq_2 )
  call pgsetmat(2,3,  mass_sq_3 - mass_int_sq_2 - mass_int_sq_3 )
  !
  call pgsetmat(3,1,  pggetmat(1,3) )
  call pgsetmat(3,2,  pggetmat(2,3) )
  call pgsetmat(3,3,  -2.0_ki*mass_int_sq_3 )
  !
  ! This call fills the internal array s_mat_r.
  ! It also assigns the integers in s_mat_p, which encode the positions
  ! of complex mass entries and zero mass entries. It includes call to init_invs
  !
  call pgpreparesmatrix()
  !
  !
  write (*,*) 'Choose what the program should compute:'
  write (*,*) '0) scalar three-point function in n dimensions'
  write (*,*) '1) three-point function in n dimensions with one Feynman parameter'
  write (*,*) '2) three-point function in n dimensions with two Feynman parameters'
  write (*,*) '3) three-point function in n dimensions with three Feynman parameters'
  write (*,*) '4) scalar three-point function in n+2 dimensions'
  write (*,*) '5) three-point function in n+2 dimensions with one Feynman parameter'
  write (*,*) '6) test of the mu independence'
  write (*,*) '7) Test table'
  read (*,*) choix
  !
  ! info for user
  if (choix == 0) then
   !
    write (*,*) 'calculating n-dim scalar 3-point fctn. with '
    !
  else if   (choix == 1) then
    !
    write (*,*) 'calculating n-dim rank one (z1) 3-point fctn. with'
    !
  else if   (choix == 2) then
    !
    write (*,*) 'calculating n-dim rank two (z1*z2) 3-point fctn. with'
    !
  else if   (choix == 3) then
    !
    write (*,*) 'calculating n-dim rank three (z1^2*z3) 3-point fctn. with'
    !
  else if   (choix == 4) then
    !
    write (*,*) 'calculating (n+2)-dim scalar 3-point fctn. with'
    !
 else if (choix == 5) then
    !
    write (*,*) 'calculating (n+2)-dim rank one (z2) 3-point fctn. with'
    !
  end if
  !
  if (choix_kinem == 1) then
    !
    write (*,*) 'one off-shell leg'
    !
  else if (choix_kinem == 2) then
    !
    write (*,*) 'two off-shell legs'
    !
  else if (choix_kinem == 3) then
    !
    write (*,*) 'three off-shell legs'
    !
  end if
  !
  write (*,*) 'the result has been written to the file test3point.txt'
  ! 
  ! start calculation    
  !
  ! To change the value of mu^2 (in GeV) (set to 1. by default)
  ! uncomment this line
  ! mu2_scale_par = 12._ki
  !
  ! store original mu^2
  mu02 = pggetmusq()
  !
  !
  ! In the following the integrals f3p_x have to be called with argument
  ! s_mat_p. This was defined with the call to preparesmatrix.
  !
  if (choix == 0) then
    !
    ! Result for the scalar integral in n dimension
    !
    ! both are working
    !
    res6c = pga30(s_null,0)
    res6b = pga30(s_null,1)
    res6a = pga30(s_null,2)
    !
    ! the labels 1,2,3 correspond to Feynman parameters z1,z2,z3
  else if (choix == 1) then
    !
    ! Results for integrals in n dimensions with one Feynman parameter 
    ! in the numerator: z1
    !
    res6c = -pga31(2,s_null, 0)
    res6b = -pga31(2,s_null, 1)
    res6a = -pga31(2,s_null, 2)
    !~ res6 = -a31(3,s_null)
    !~ verif = f3p(s_mat,b_ref,3)
    !
  else if (choix == 2) then
    !
    ! Results for integrals in n dimensions with two Feynman parameters 
    ! in the numerator: z1*z2
    !
    res6c =  pga32(1,2,s_null, 0)
    res6b =  pga32(1,2,s_null, 1)
    res6a =  pga32(1,2,s_null, 2)
    !~ res6 =  a32(2,2,s_null)
    !~ verif = f3p(s_mat,b_ref,2,2)
    !
  else if (choix == 3) then
    !
    ! Results for integrals in n dimensions with three Feynman parameters 
    ! at the numerator: z1^2*z3
    !
    res6c =  -pga33(2,2,2,s_null, 0)
    res6b =  -pga33(2,2,2,s_null, 1)
    res6a =  -pga33(2,2,2,s_null, 2)
    !
  else if (choix == 4) then
    !
    ! Results for integrals in n+2 dimensions with no Feynman parameters 
    ! at the numerator: 
    !
    res6c =  -2.0_ki*pgb32(s_null, 0)
    res6b =  -2.0_ki*pgb32(s_null, 1)
    res6a =  -2.0_ki*pgb32(s_null, 2)
    ! 
  else if (choix == 5) then
    !
    ! Results for integrals in n+2 dimensions with one Feynman parameters 
    ! at the numerator: z2
    !
    res6c =  2.0_ki*pgb33(2,s_null, 0)
    res6b =  2.0_ki*pgb33(2,s_null, 1)
    res6a =  2.0_ki*pgb33(2,s_null, 2)
    !
  else if (choix == 6) then
    !
    ! by default, mu2_scale_par = 1._ki
    ! take scalar triangle as example
    res6c=pga30(s_null, 0)
    res6b=pga30(s_null, 1)
    res6a=pga30(s_null, 2)
    ! we have to reset the cache in order that the new value
    ! of mu2_scale_par will be effective. A call to preparesmatrix
    ! is sufficient. 
    !
    call pgpreparesmatrix()
    !
    ! we change the value of mu^2
    !
    mu2 = 34._ki
    lmu2 = log(mu2/pggetmusq())
    call pgsetmusq(mu2)
    res7c =  pga30(s_null, 0)
    res7b =  pga30(s_null, 1)
    res7a =  pga30(s_null, 2)
 !
  else if (choix == 7) then
    !
    write(17,*) " TRIANGLE "
    write(17,*) " Rank 0 "
    write(17,'("C0( 0) = (",D23.15," +I* ",D23.15,")")') pga30( s_null, 0)
    write(17,'("C0(-1) = (",D23.15," +I* ",D23.15,")")') pga30( s_null, 1)
    write(17,'("C0(-2) = (",D23.15," +I* ",D23.15,")")') pga30( s_null, 2)
    write(17,*) " Rank 1 "
    do i1=1,3
      write(17,'("C",I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1, pga31(i1, s_null, 0)
      write(17,'("C",I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1, pga31(i1, s_null, 1)
      write(17,'("C",I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1, pga31(i1, s_null, 2)
    enddo
    ! --------------------------------
    write(17,*) " Rank 2 "
    write(17,'("C00( 0) = (",D23.15," +I* ",D23.15,")")') pgb32( s_null, 0)
    write(17,'("C00(-1) = (",D23.15," +I* ",D23.15,")")') pgb32( s_null, 1)
    write(17,'("C00(-2) = (",D23.15," +I* ",D23.15,")")') pgb32( s_null, 2)
    do i1=1,3
      do i2=i1,3
        write(17,'("C",2I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,i2, pga32(i1,i2, s_null, 0)
        write(17,'("C",2I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,i2, pga32(i1,i2, s_null, 1)
        write(17,'("C",2I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,i2, pga32(i1,i2, s_null, 2)
      enddo
    enddo
    ! --------------------------------
    write(17,*) " Rank 3 "
    do i1=1,3
      write(17,'("C00",I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1, pgb33(i1, s_null, 0)
      write(17,'("C00",I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1, pgb33(i1, s_null, 1)
      write(17,'("C00",I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1, pgb33(i1, s_null, 2)
    enddo
    do i1=1,3
      do i2=i1,3
        do i3=i2,3
          write(17,'("C",3I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3, pga33(i1,i2,i3, s_null, 0)
          write(17,'("C",3I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3, pga33(i1,i2,i3, s_null, 1)
          write(17,'("C",3I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,i2,i3, pga33(i1,i2,i3, s_null, 2)
        enddo
      enddo
    enddo
    ! --------------------------------
  ! --------------------------------------------------------------
    do s1=1,3
      write(17,*) " BUBBLES s = ",s1
      b_pin = packb( (/s1/) )
      write(17,*) " Rank 0 "
      write(17,'("B0( 0) = (",D23.15," +I* ",D23.15,")")') pga20(b_pin, 0)
      write(17,'("B0(-1) = (",D23.15," +I* ",D23.15,")")') pga20(b_pin, 1)
      write(17,'("B0(-2) = (",D23.15," +I* ",D23.15,")")') pga20(b_pin, 2)
      write(17,*) " Rank 1 "
      do i1=1,3
        if ((i1.eq.s1)) cycle
        write(17,'("B",I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1, pga21(i1, b_pin, 0)
        write(17,'("B",I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1, pga21(i1, b_pin, 1)
        write(17,'("B",I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1, pga21(i1, b_pin, 2)
      enddo
      ! --------------------------------
      write(17,*) " Rank 2 "
      write(17,'("B00( 0) = (",D23.15," +I* ",D23.15,")")') pgb22(b_pin, 0)
      write(17,'("B00(-1) = (",D23.15," +I* ",D23.15,")")') pgb22(b_pin, 1)
      write(17,'("B00(-2) = (",D23.15," +I* ",D23.15,")")') pgb22(b_pin, 2)
      do i1=1,3
        if ((i1.eq.s1)) cycle
        do i2=i1,3
          if ((i2.eq.s1)) cycle
          write(17,'("B",2I1,"( 0) = (",D23.15," +I* ",D23.15,")")') i1,i2, pga22(i1,i2, b_pin, 0)
          write(17,'("B",2I1,"(-1) = (",D23.15," +I* ",D23.15,")")') i1,i2, pga22(i1,i2, b_pin, 1)
          write(17,'("B",2I1,"(-2) = (",D23.15," +I* ",D23.15,")")') i1,i2, pga22(i1,i2, b_pin, 2)
        enddo
      enddo
      ! --------------------------------
    ! --------------------------------------------------------------
      do s2=(s1+1),3
        write(17,*) " TADPOLES s,t = ",s1,s2
        b_pin = packb( (/s1,s2/) )
        write(17,*) " Rank 0 "
        write(17,'("A0( 0) = (",D23.15," +I* ",D23.15,")")') pga10(b_pin, 0)
        write(17,'("A0(-1) = (",D23.15," +I* ",D23.15,")")') pga10(b_pin, 1)
        write(17,'("A0(-2) = (",D23.15," +I* ",D23.15,")")') pga10(b_pin, 2)
      enddo
    ! --------------------------------------------------------------
    enddo
  end if
  !
  write (17,*) 'The kinematics is:'
  write (17,*) ''
  write (17,*) '             |             '
  write (17,*) '             | p1        '
  write (17,*) '             |             '
  write (17,*) '            /\            '
  write (17,*) '           /  \           '
  write (17,*) '     (1) /     \ (3)     '
  write (17,*) '        /        \         '
  write (17,*) '       /--->----\       '
  write (17,*) '      /   (2)      \      '
  write (17,*) 'p2 /               \ p3'
  write (17,*) ''
  write (17,*) 'p1+p2+p3 = 0'
  write (17,*) ''
  write (17,*) '(p1)^2 =',mass_sq_1
  write (17,*) '(p2)^2 =',mass_sq_2
  write (17,*) '(p3)^2 =',mass_sq_3
  write (17,*) 'm1^2 =',mass_int_sq_1
  write (17,*) 'm2^2 =',mass_int_sq_2
  write (17,*) 'm3^2 =',mass_int_sq_3
  write (17,*) 'mu^2 =',pggetmusq()
  write (17,*) ''
  write (17,*) 'defining I_N^n= mu^(4-n) \int d^n k/(i*Pi^(n/2))*func(k,p_i)'
  write (17,*) '= r_Gam *(P2/eps^2+P1/eps+P0),'
  write (17,*) 'n = 4-2*eps,'
  write (17,*) 'r_Gam = Gamma(1+eps)*Gamma(1-eps)^2/Gamma(1-2eps)'
  write (17,*) 'the program gives numbers for P2,P1,P0'
  write (17,*) ''
    !
  write (17,*) 'result='
  write (17,'("  1/epsilon^2 * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6a,ki),aimag(res6a)
  write (17,'("+ 1/epsilon   * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6b,ki),aimag(res6b)
  write (17,'("+ 1           * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6c,ki),aimag(res6c)
  write (17,*) ''
  write (6,*) 'result='
  write (6,'("  1/epsilon^2 * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6a,ki),aimag(res6a)
  write (6,'("+ 1/epsilon   * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6b,ki),aimag(res6b)
  write (6,'("+ 1           * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6c,ki),aimag(res6c)
  !
  if ( choix < 6 ) then
  !
  !
  else if ( choix ==6 ) then
  !
    write (17,*) 'The preceding result has been computed with mu^2=',mu02
    write (17,*) ' '
    write (17,*) 'Now setting by hand mu^2=',mu2
    write (17,*) 'and expanding (',mu2,'/',mu02,')^epsilon around epsilon=0'
    write (17,'("  1/epsilon^2 * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6a,ki),aimag(res6a)
    write (17,'("+ 1/epsilon   * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6b,ki)+ &
         &real(res6a,ki)*lmu2,aimag(res6b)+aimag(res6a)*lmu2
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
    write (6,'("+ 1/epsilon   * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6b,ki) + &
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
    write (6,*) ' '
  !
  else 
  !
  write (6,*) 'invalid choice, option number must be < 7'
  !
  endif
  ! routine to free the cache and allocated memory
  !
!   call exitgolem95()
  !
  close(17)
  close(19)
  !
end program main
!
!
