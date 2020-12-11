program integrator 
   use module_pl
   use, intrinsic :: iso_fortran_env ! Use the intrinsic kind definitions
   implicit none
   ! Arguments
   type(user_input_parameters)  :: param    ! derived type containing user-defined parameters
   ! Internals 
   integer, parameter :: P = real64
   real(P) :: e, inc, long, varpi, PI, f, mu, om, a, argr, argv
   integer :: iarg, w, npl
   logical :: new 
   real(P) :: year, timestep, tfinal 
   CHARACTER (len = 4) :: argu
   real(P), dimension(3) :: h, r, v , ecc , n, vcrossh, z
   type(planet) :: pl
   year = 24_P*3600_P*365_P
   timestep = 1_P * year
   half_timestep = 0.5_P * timestep
   tfinal = 1000000_P*year
   ! Obtain inputs 
   type planet 
   integer, dimension(:),           allocatable :: id
   real(P),     dimension(:,:),     allocatable :: xh
   real(P),     dimension(:,:),     allocatable :: vh
   real(P),     dimension(:),       allocatable :: mass
   real(P),     dimension(:,:),     allocatable :: xb
   real(P),     dimension(:,:),     allocatable :: vb
   real(P),     dimension(:,:),     allocatable :: ah
   end type planet 

   id2 = [3.3473818717761439E-01,-2.1061105379199779E-01,-4.7921461216598432E-02, 9.4572394374496608E-03, 2.5106125757836419E-02, 1.1835934147070429E-03, 1.651e-07]
   id3 = [-4.6411665443250860E-01, 5.4731602408177948E-01, 3.4285183291441222E-02,-1.5497531935705990E-02,-1.3190815340356370E-02, 7.1366857195630977E-04, 2.447e-06]
   id4 = [7.8447422290361046E-01, 6.0834664588920739E-01,-1.9969120829822830E-05,-1.0818280103689040E-02, 1.3526750837437910E-02, 2.3301627510155490E-07, 3.003e-06]
   id5 = [3.2488082974635041E-01,-1.3920413191921059E+00,-3.7142241988251279E-02, 1.4156783775670190E-02, 4.3809924986283897E-03,-2.5580138684768577E-04, 3.213e-07]
   id6 = [1.8733430375446749E+00, 4.6833225281837292E+00,-6.1370789424293443E-02,-7.1040288545694674E-03, 3.1641930734429640E-03, 1.4582112113878479E-04, 9.543e-04]
   id7 = [-8.2518654799859821E+00,-5.2250086399581823E+00, 4.1939352212419062E-01, 2.6773744239338769E-03,-4.7239614719522174E-03,-2.4586433481282280E-05, 2.857e-4]
   id8 = [1.9922348963200001E+01, 2.3426195643761392E+00,-2.4935475570954871E-01,-4.9364207001745177E-04, 3.7248681728472809E-03, 2.0274397915572331E-05, 4.365e-5]
   id9 = [2.6471589364036468E+01,-1.4096521334906500E+01,-3.1964441840438462E-01, 1.4491178510995289E-03, 2.7916302817640199E-03,-9.0630886544057097E-05,5.149e-5]
   id10= [4.9134676668880539E+00,-3.1886966347017140E+01, 1.9916040459076740E+00, 3.1560668298849150E-03,-1.4588989139232661E-04,-8.9870766554492072E-04, 6.5510e-09] 
   if (len(pl_list) <= 1) then
      npl = 1
   else
      npl = len(pl_list)
   end if

   allocate(pl%id(npl))
   allocate(pl%xh(3,npl))
   allocate(pl%vh(3,npl))
   allocate(pl%xb(3,npl))
   allocate(pl%vb(3,npl))
   allocate(pl%mass(npl))
   allocate(pl%ah(3,npl))


   DO i = 1, n 
   pl%id = i 
   pl%xh(:,i) = pl_list[i][1:3]
   pl%vh = pl_list[i][4:6]
   pl%mass = pl_list[i][7]
   END DO 
   
   open(10,file='mercury.out', action = 'write')
   open(11,file='venus.out', action = 'write')
   open(12,file='earth.out', action = 'write')
   open(13,file='mars.out', action = 'write')
   open(14,file='jupiter.out', action = 'write')
   open(15,file='saturn.out', action = 'write')
   open(16,file='uranus.out', action = 'write')
   open(17,file='neptune.out', action = 'write')
   open(18,file='pluto.out', action = 'write')

   DO i = 1, timestep, tfinal 
      !convert to barycentric coordinates 
      xbtot = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
      vbtot = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
      ptot =  (/ 0.0_DP, 0.0_DP, 0.0_DP /)
      mtot = 0_P
      DO j = 2, npl 
         mtot = mtot + pl%mass
         xbtot =  = pl%xh(:,i)*pl%mass(i) +xbtot
         vbtot = pl%vh(:,i)*pl%mass(i) +xbtot
      END DO 
      ptot = xbtot
      xbtot =  - xbtot / mtot 
      vbtot = -vbtot / mtot 
      ptot = xbtot * half_timestep / pl%mass(1)
      pl%xb(:,1) = xbtot 
      pl%vb(:,1) = vbtot 
      
      DO j = 2, npl 
         pl%xb(:,i) = pl%xh(:,i) + xbtot
         pl%vb(:,i) = pl%vh(:,i) + vbtot
      END DO 
      ! Sun momemtum drift at half_timestep 
      DO j = 2, npl 
         pl%xh(:,i) = pl%xh(:,i) + ptot
      END DO 
      !! kick due to interactions of all bodies except Sun 
      DO j = 2, npl 
         DO k = 2, npl 
            dx(:) = pl%xh(:,k) - pl%xh(:,j)
            rplpl = dot_product(dx,dx)
            pl%ah(:,j) = pl%ah(:,j) + pl%mass(k)*dx(:)/rplpl**(3_P/2_P)
            pl%ah(:,k) = pl%ah(:,k) - pl%mass(j)*faci*dx(:)/rplpl**(3_P/2_P)
         END DO 
      END DO 
   
      DO j = 2, npl 
         pl%vb(:,j) = pl%vb(:,j) + pl%ah(:,j)*half_timestep
      END DO 


      !! Danby drift of all bodies  at timestep 
      DO j = 2, npl 
         call drift_danby(pl%xh(:,i),pl%vh(:,i),xnew,vnew, half_timestep)
         pl%xh(:,i) = xnew 
         pl%vh(:,i) = vnew
      END DO 
      ! Kick at half_timestep
      DO j = 2, npl 
         DO k = 2, npl 
            dx(:) = pl%xh(:,k) - pl%xh(:,j)
            rplpl = dot_product(dx,dx)
            pl%ah(:,j) = pl%ah(:,j) + pl%mass(k)*dx(:)/rplpl**(3_P/2_P)
            pl%ah(:,k) = pl%ah(:,k) - pl%mass(j)*faci*dx(:)/rplpl**(3_P/2_P)
         END DO 
      END DO 
   
      DO j = 2, npl 
         pl%vb(:,j) = pl%vb(:,j) + pl%ah(:,j)*half_timestep
      END DO  


      ! Sun momemtum drift at half_timestep and planet drift 
      ptot =  (/ 0.0_DP, 0.0_DP, 0.0_DP /)
      DO j = 2, npl 
         ptot = pl%vb(:,i)*pl%mass(i) +xbtot
      END DO
      ptot = ptot * half_timestep / pl%mass(1)

      DO i = 2, npl 
         pl%xh(:,i) = pl%xh(:,i) + ptot
      END DO 
      !convert vb into vh 
      vbtot = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
      DO j = 2 , npl
         vbtot = vbtot * pl%mass(i)*pl%vb(:,i)
      END DO 
      vbtot = vbtot / pl%mass(1) 
      DO j = 2, npl 
         pl%vh(:,i) = pl%vb(:,i) + vbtot
      END DO 

      !! output 
      if MOD(i,cadence) == 0 then 
         DO j = 1, npl
            k = j+9
         write(k,*) i, pl%xh(1,j), pl%xh(2,j), pl%xh(3,j), pl%vh(1,j), pl%vh(2,j), pl%vh(3,j)
   END DO 

   close(10)
   close(11)
   close(12)
   close(13)
   close(14)
   close(15)
   close(16)
   close(17)
   close(18)
   close(19)




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

   subroutine drift_danby(x, v, xnew,vnew, dt)
      use, intrinsic :: iso_fortran_env
      implicit none 
      integer, parameter :: P = real64
      real(P), dimension(3), intent(in) :: x,v
      real(P), dimension(3), intent(out) :: xnew, vnew
      real(P) :: P, a, ecc, f, varpi, xperi, yperi, sE, cE
      real(P) :: E, E0, dE, M, t, telapsed
      integer :: n, nsteps  ! The number of steps to generate output
      nsteps = size(xv,2)
      ! Convert cartesian initial conditions to orbital elements 
      call xv2el(x,v,a,ecc,f,varpi)
      P = 2 * pi * sqrt(a**3) 
      ! First rotate the position of the orbit into pericenter-aligned coordinates
      xperi = x(1) * cos(-varpi) - x(2) * sin(-varpi)
      yperi = x(1) * sin(-varpi) + x(2) * cos(-varpi)
      sE = yperi / (a * sqrt(1._DP - ecc**2))
      cE = xperi / a + ecc
      E = atan2(sE, cE) 
      M = E - ecc * sin(E)
      n = 10
      DO i = 1, nsteps
      t = (n - 1) * dt 
      telapsed = mod(t, P)
      M = (telapsed / P) * 2 * pi 
      ! Update eccentric anomaly (solve Kepler's equation)
      E = danby_solver(M, ecc)
      ! Compute new cartesian position and velocity vectors
      call Etoxv(a, E, telapsed, xnew, vnew)
      end do

      return
   end subroutine drift_danby

   function danby_solver(M, ecc) result(E)
      implicit none
      integer, parameter :: P = real64
      real(P), intent(in) :: M, ecc
      real(P) :: E
      real(P), parameter :: tol = epsilon(M)
      integer, parameter :: imax = 1000
      real(P), parameter :: k = 0.85_DP
      real(P) :: Eold, dE
      integer :: i
      E = M + sign(k * ecc, sin(M))
      do i = 1, imax
         Eold = E
         E = danby_step(M, ecc, E)
         dE = abs(E - Eold)
         if (dE < tol) exit
      end do
      return
   end function danby_solver

   function danby_step(M, ecc, E) result(Enew)
      implicit none
      integer, parameter :: P = real64
      real(P), intent(in) :: M, ecc, E
      real(P) :: Enew
      real(P) :: d1, d2, d3
      real(P) :: f, fp, fpp, fppp

      f = E - ecc * sin(E) - M
      fp = 1._DP - ecc * cos(E)
      fpp = ecc * sin(E)
      fppp = ecc * cos(E)

      d1 = -f / fp
      d2 = -f / (fp + d1 * fpp / 2._DP)
      d3 = -f / (fp + d2 * fpp / 2._DP + d2**2 * fppp / 6._DP)

      Enew = E + d3
      return
   end function danby_step

   subroutine Etoxv(a,dE,dt,x,v)
      ! Uses the f and g functions to solve for the cartesian position and velocity as a function of
      implicit none
      integer, parameter :: P = real64
      real(P), intent(in) :: a  ! semimajor axis and change in eccentric anomaly
      real(P), intent(in) :: dE ! Change in eccentric anomaly in one step
      real(P), intent(in) :: dt ! Step size
      real(P), dimension(:), intent(inout)  :: x, v ! Position and velocity vectors. Input is x0, v0, output is new x and v
      real(P) :: f, g, fdot, gdot ! f and g functions, and their time derivatives
      real(P) :: r0, r, n
      real(P), dimension(:), allocatable :: x0, v0

      allocate(x0, source=x)
      allocate(v0, source=v)

      n = sqrt(a**(-3)) ! Mean motion
      r0 = norm2(x0(:))

      f = a / r0 * (cos(dE) - 1._DP) + 1._DP
      g = dt + 1._DP / n * (sin(dE) - dE)

      x(:) = f * x0(:) + g * v0(:)

      r = norm2(x(:))

      fdot = -(a**2 * n / (r * r0) ) * sin(dE)
      gdot = (a / r) * (cos(dE) - 1._DP) + 1._DP

      v(:) = fdot * x0(:) + gdot * v0(:)

      return
   end subroutine Etoxv

   subroutine xv2el(x,v,a,ecc,f,varpi, inc, capom)
      ! Arguments
      integer, parameter :: P = real64
      REAL(P), INTENT(IN)                  :: mu
      REAL(P), DIMENSION(:), INTENT(IN)    :: x, v
      REAL(P), INTENT(OUT)                 :: a, e, a,ecc,f,varpi, inc, capom
         
      ! Internals
      REAL(P)     :: r, v2, hx, hy, hz, h2, h, rdotv, energy, capm, fac, u, w, cw, sw, face, cape, tmpf, capf
         
      ! Executable code
      a = 0.0_DP
      e = 0.0_DP
      inc = 0.0_DP
      capom = 0.0_DP
      omega = 0.0_DP
      capm = 0.0_DP
      r = SQRT(DOT_PRODUCT(x(:), x(:)))
      v2 = DOT_PRODUCT(v(:), v(:))
      hx = x(2)*v(3) - x(3)*v(2)
      hy = x(3)*v(1) - x(1)*v(3)
      hz = x(1)*v(2) - x(2)*v(1)
      h2 = hx*hx + hy*hy + hz*hz
      h = SQRT(h2)
      IF (h2 == 0.0_DP) RETURN
      rdotv = DOT_PRODUCT(x(:), v(:))
      energy = 0.5_DP*v2 - mu/r
      fac = hz/h
      IF (fac < -1.0_DP) THEN
           inc = PI
      ELSE IF (fac < 1.0_DP) THEN
           inc = ACOS(fac)
      END IF
      fac = SQRT(hx*hx + hy*hy)/h
      IF (fac**2 < VSMALL) THEN
           u = ATAN2(x(2), x(1))
           IF (hz < 0.0_DP) u = -u
      ELSE
           capom = ATAN2(hx, -hy)
           u = ATAN2(x(3)/SIN(inc), x(1)*COS(capom) + x(2)*SIN(capom))
      END IF
      IF (capom < 0.0_DP) capom = capom + TWOPI
      IF (u < 0.0_DP) u = u + TWOPI

       fac = 1.0_DP - h2/(mu*a)
      IF (fac > VSMALL) THEN
         e = SQRT(fac)
         cape = 0.0_DP
         face = (a - r)/(a*e)
         IF (face < -1.0_DP) THEN
            cape = PI
         ELSE IF (face < 1.0_DP) THEN
            cape = ACOS(face)
         END IF
         IF (rdotv < 0.0_DP) cape = TWOPI - cape
         fac = 1.0_DP - e*COS(cape)
         cw = (COS(cape) - e)/fac
         sw = SQRT(1.0_DP - e*e)*SIN(cape)/fac
         w = ATAN2(sw, cw)
         IF (w < 0.0_DP) w = w + TWOPI
      ELSE
         cape = u
         w = u
      END IF
      capm = cape - e*SIN(cape)
          
      omega = u - w
      f = w 
      IF (omega < 0.0_DP) omega = omega + TWOPI
      varpi = capom + omega
      RETURN 
         
   END SUBROUTINE xv2el






