program integrator 
   use planet
   use, intrinsic :: iso_fortran_env ! Use the intrinsic kind definitions
   implicit none
   ! Arguments
   type(user_input_parameters)  :: param    ! derived type containing user-defined parameters
   ! Internals 
   integer, parameter :: P = real64
   real(P) :: e, inc, long, varpi, PI, f, mu, om, a, argr, argv
   integer :: iarg, w
   logical :: new 
   real(P) :: year, timestep, tfinal 
   CHARACTER (len = 4) :: argu
   real(P), dimension(3) :: h, r, v , ecc , n, vcrossh, z
   type(planet) :: pl
   year = 24_P*3600_P*365_P
   timestep = 1_P * year
   tfinal = 1000000_P*year
   ! Obtain inputs 
   iarg = iargc()
      do  w = 1, iarg
         call getarg( w, argu )
         SELECT CASE (argu)
         CASE("-r")
            call getarg( w+1, argu )
            Read(argu,*) argr 
            r(1) = argr
            call getarg( w+2, argu )
            Read(argu,*) argr 
            r(2) = argr
            call getarg( w+3, argu )
            Read(argu,*) argr 
            r(3) = argr
         CASE("-v")
            call getarg( w+1, argu )
            Read(argu,*) argv 
            v(1)=argv
            call getarg( w+2, argu )
            Read(argu,*) argv 
            v(2)=argv
            call getarg( w+3, argu )
            Read(argu,*) argv 
            v(3)=argv
         CASE("-new")
            new = .True.
         END SELECT
      enddo

   
   call pl_allocate(pl, pl_list)


   DO i = 1, timestep, tfinal 
      convert_dh(x,v)
   END DO 

   end program 
 
   subroutine cross(x, v, xcrossv)
      use, intrinsic :: iso_fortran_env
      implicit none 
      integer, parameter :: P = real64
      real(P), dimension(3), intent(in) :: x,v
      real(P), dimension(3), intent(out) :: xcrossv 
   
      xcrossv(1) = x(2)*v(3) - x(3)*v(2)
      xcrossv(2) = x(3)*v(1) - x(1)*v(3)
      xcrossv(3) = x(1)*v(2) - x(2)*v(1)
   
   end subroutine 


   subroutine kick(x, v, xnew,vnew)
      use, intrinsic :: iso_fortran_env
      implicit none 
      integer, parameter :: P = real64
      real(P), dimension(3), intent(in) :: x,v
      real(P), dimension(3), intent(out) :: xnew, vnew
   

   
   end subroutine 

   subroutine drift_danby(x, v, xnew,vnew)
      use, intrinsic :: iso_fortran_env
      implicit none 
      integer, parameter :: P = real64
      real(P), dimension(3), intent(in) :: x,v
      real(P), dimension(3), intent(out) :: xnew, vnew
   

   
   end subroutine 

   subroutine pl_allocate(pl, pl_list)
      implicit none 
      integer, parameter :: P = real64
      REAL(P), dimension(:,:), intent(in) :: pl_list  
      type(pl), intent(inout) :: pl
      integer :: n, i 
      if (len(pl_list) <= 1) then
         n = 1
      else
         n = len(pl_list)
      end if
      DO i = 1, n 
      pl%id = i 
      pl%xh = pl_list[i][1]
      pl%yh = pl_list[i][1]
      pl%zh = pl_list[i][1]
      pl%vxh = pl_list[i][1]
      pl%vyh = pl_list[i][1]
      pl%vz = pl_list[i][1]
      END DO 
   END SUBROUTINE 

   MODULE module_pl 
      implicit none 
      integer, parameter :: P = real64
      type planet 
      integer, dimension(:),         allocatable :: id
      real(P),     dimension(:),     allocatable :: xh
      real(P),     dimension(:),     allocatable :: yh
      real(P),     dimension(:),     allocatable :: zh
      real(P),     dimension(:),     allocatable :: vxh
      real(P),     dimension(:),     allocatable :: vyh
      real(P),     dimension(:),     allocatable :: vz
      end type planet 
   END MODULE 



