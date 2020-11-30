! Modern Fortran code that converts cartesian coordiante vectors (r,v) into orbital elements (a,e,inc, long, varpi, f)
program xv2el
   use, intrinsic :: iso_fortran_env ! Use the intrinsic kind definitions
   implicit none
   integer, parameter :: P = real64
   real(P) :: e, inc, long, varpi, PI, f, mu, om, a, argr, argv
   integer :: iarg, w
   logical :: new 
   CHARACTER (len = 4) :: argu
   real(P), dimension(3) :: h, r, v , ecc , n, vcrossh, z
   PI = 2*asin(1.0_P)
   new = .False. 
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

   WRITE(*,*) 'new = ', new 
   mu = 1._P 
   z = [0.0_P, 0.0_P,  1.0_P]
   CALL cross(r, v, h)
   CALL cross(v, h, vcrossh)
   ecc = (1._P / mu)*vcrossh - r / norm2(r)
   CALL cross(z,h,n)
   IF (dot_product(r,v) >= 0 ) THEN 
      f = acos((dot_product(ecc,r)/(norm2(r)*e)))
   ELSE 
      f = 2.0_P*PI - acos((dot_product(ecc,r)/(norm2(r)*e)))
   END IF 
   WRITE(*,*) "ecc = ", ecc
   WRITE(*,*) "r = ", r 
   WRITE(*,*) "v = ", v
   inc = acos(h(3)/norm2(h))

   e = norm2(ecc)

   IF (n(2) >=0) THEN 
      long = acos(n(1)/norm2(n))
   ELSE 
      long = 2*PI - acos(n(1)/norm2(n))
   END IF
   IF (n(2) >=0) THEN 
      om = acos(dot_product(n,ecc)/norm2(n)/e)
   ELSE 
      om = 2.0_P*PI - acos(dot_product(n,ecc)/norm2(n)/e)
   END IF 

   varpi = long + om 

   a = 1._P/((2.0_P/norm2(r))- (norm2(v)**2.0_P)/mu)

   if (new) THEN
      open(15,file='xv2el.output', action = 'write')
   else 
      open(15,file='xv2el.output', status= 'old', action = 'write', position = 'append')
   end if 

   write(15,*) a,',',e,',',inc,',',long,',',varpi,',',f
   close(15)

end program xv2el

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


