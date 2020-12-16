program integrator 
   use, intrinsic :: iso_fortran_env ! Use the intrinsic kind definitions
   implicit none
   ! Internals 
   integer, parameter :: P = real64
   real(P) :: e, inc, long, varpi, PI, f, mu, om, a, argr, argv, sdebug, convert
   integer :: iarg, w, npl,i,j,k, cadence
   logical :: new 
   real(P) :: year, timestep, tfinal, half_timestep, time, mtot, rplpl, G, time_conversion
   CHARACTER (len = 4) :: argu
   real(P), dimension(3) :: h, r, v , ecc , n, vcrossh, z, dx, xbtot, vbtot, ptot
   real(P), dimension(3) :: xnew, vnew, vdebug 
   real(P), dimension(7) :: id2, id3, id4, id5, id6, id7, id8, id9, id10
   real(P), dimension(9, 7) :: pl_list
   ! type 
   type planet 
   integer, dimension(:),           allocatable :: id
   real(P),     dimension(:,:),     allocatable :: xh
   real(P),     dimension(:,:),     allocatable :: vh
   real(P),     dimension(:),       allocatable :: mass
   real(P),     dimension(:,:),     allocatable :: xb
   real(P),     dimension(:,:),     allocatable :: vb
   real(P),     dimension(:,:),     allocatable :: ah
   end type planet 
   type(planet) :: pl
   G = 4.0_P*pi**2.0_P
   cadence = 1
   convert = 365.25_P
   year = 24_P*3600_P*365.25_P
   timestep = 0.1_P * year
   half_timestep = 0.5_P * timestep
   tfinal = 1000_P*year
   time_conversion = 1.0_P / year
   ! inputs 

   ! id2 = [3.3473818717761439E-01,-2.1061105379199779E-01,-4.7921461216598432E-02, &
   ! 9.4572394374496608E-03, 2.5106125757836419E-02, 1.1835934147070429E-03, 1.651e-07]
   ! id3 = [-4.6411665443250860E-01, 5.4731602408177948E-01, 3.4285183291441222E-02,&
   ! -1.5497531935705990E-02,-1.3190815340356370E-02, 7.1366857195630977E-04, 2.447e-06]
   ! id4 = [7.8447422290361046E-01, 6.0834664588920739E-01,-1.9969120829822830E-05,&
   ! -1.0818280103689040E-02, 1.3526750837437910E-02, 2.3301627510155490E-07, 3.003e-06]
   ! id5 = [3.2488082974635041E-01,-1.3920413191921059E+00,-3.7142241988251279E-02,&
   !  1.4156783775670190E-02, 4.3809924986283897E-03,-2.5580138684768577E-04, 3.213e-07]
   ! id6 = [1.8733430375446749E+00, 4.6833225281837292E+00,-6.1370789424293443E-02,&
   ! -7.1040288545694674E-03, 3.1641930734429640E-03, 1.4582112113878479E-04, 9.543e-04]
   ! id7 = [-8.2518654799859821E+00,-5.2250086399581823E+00, 4.1939352212419062E-01,&
   !  2.6773744239338769E-03,-4.7239614719522174E-03,-2.4586433481282280E-05, 2.857e-4]
   ! id8 = [1.9922348963200001E+01, 2.3426195643761392E+00,-2.4935475570954871E-01,&
   ! -4.9364207001745177E-04, 3.7248681728472809E-03, 2.0274397915572331E-05, 4.365e-5]
   ! id9 = [2.6471589364036468E+01,-1.4096521334906500E+01,-3.1964441840438462E-01,&
   !  1.4491178510995289E-03, 2.7916302817640199E-03,-9.0630886544057097E-05,5.149e-5]
   ! id10= [4.9134676668880539E+00,-3.1886966347017140E+01, 1.9916040459076740E+00,&
   !  3.1560668298849150E-03,-1.4588989139232661E-04,-8.9870766554492072E-04, 6.5510e-09] 
   

    pl_list(1,1:7) = [3.3473818717761439E-01,-2.1061105379199779E-01,-4.7921461216598432E-02, &
    9.4572394374496608E-03, 2.5106125757836419E-02, 1.1835934147070429E-03, 1.651e-07]
    pl_list(2,1:7) = [-4.6411665443250860E-01, 5.4731602408177948E-01, 3.4285183291441222E-02,&
    -1.5497531935705990E-02,-1.3190815340356370E-02, 7.1366857195630977E-04, 2.447e-06]
    pl_list(3,1:7) = [7.8447422290361046E-01, 6.0834664588920739E-01,-1.9969120829822830E-05,&
    -1.0818280103689040E-02, 1.3526750837437910E-02, 2.3301627510155490E-07, 3.003e-06]
    pl_list(4,1:7) = [3.2488082974635041E-01,-1.3920413191921059E+00,-3.7142241988251279E-02,&
     1.4156783775670190E-02, 4.3809924986283897E-03,-2.5580138684768577E-04, 3.213e-07]
     pl_list(5,1:7) = [1.8733430375446749E+00, 4.6833225281837292E+00,-6.1370789424293443E-02,&
    -7.1040288545694674E-03, 3.1641930734429640E-03, 1.4582112113878479E-04, 9.543e-04]
    pl_list(6,1:7) = [-8.2518654799859821E+00,-5.2250086399581823E+00, 4.1939352212419062E-01,&
     2.6773744239338769E-03,-4.7239614719522174E-03,-2.4586433481282280E-05, 2.857e-4]
     pl_list(7,1:7) = [1.9922348963200001E+01, 2.3426195643761392E+00,-2.4935475570954871E-01,&
    -4.9364207001745177E-04, 3.7248681728472809E-03, 2.0274397915572331E-05, 4.365e-5]
    pl_list(8,1:7) = [2.6471589364036468E+01,-1.4096521334906500E+01,-3.1964441840438462E-01,&
     1.4491178510995289E-03, 2.7916302817640199E-03,-9.0630886544057097E-05,5.149e-5]
     pl_list(9,1:7)= [4.9134676668880539E+00,-3.1886966347017140E+01, 1.9916040459076740E+00,&
     3.1560668298849150E-03,-1.4588989139232661E-04,-8.9870766554492072E-04, 6.5510e-09] 
    
   ! if (shape(pl_list) <= 1) then
   !    npl = 1
   ! else
   !    npl = shape(pl_list)
   ! end if

   npl = 10

   allocate(pl%id(npl))
   allocate(pl%xh(3,npl))
   allocate(pl%vh(3,npl))
   allocate(pl%xb(3,npl))
   allocate(pl%vb(3,npl))
   allocate(pl%mass(npl))
   allocate(pl%ah(3,npl))


   pl%id(1) = 1 
   pl%xh(:,1) = (/ 0.0_P, 0.0_P, 0.0_P /)
   pl%vh(:,1) = (/ 0.0_P, 0.0_P, 0.0_P /)
   pl%mass(1) = 1.0_P

   do i = 2,npl
      pl%id(i) = i

      pl%xh(1,i) = pl_list(i-1,1)
      pl%xh(2,i) = pl_list(i-1,2)
      pl%xh(3,i) = pl_list(i-1,3)
      pl%ah(:,i) = (/ 0.0_P, 0.0_P, 0.0_P /)
      pl%vh(1,i) = pl_list(i-1,4)*convert ! in AU/year 
      pl%vh(2,i) = pl_list(i-1,5)*convert
      pl%vh(3,i) = pl_list(i-1,6)*convert
      pl%mass(i) = pl_list(i-1,7)
   enddo 

   open(10,file='mercury.out', action = 'write')
   open(11,file='venus.out', action = 'write')
   open(12,file='earth.out', action = 'write')
   open(13,file='mars.out', action = 'write')
   open(14,file='jupiter.out', action = 'write')
   open(15,file='saturn.out', action = 'write')
   open(16,file='uranus.out', action = 'write')
   open(17,file='neptune.out', action = 'write')
   open(18,file='pluto.out', action = 'write')

   time = 1
   do while(time<tfinal)      !convert to barycentric coordinates 
      xbtot = (/ 0.0_P, 0.0_P, 0.0_P /)
      vbtot = (/ 0.0_P, 0.0_P, 0.0_P /)
      ptot =  (/ 0.0_P, 0.0_P, 0.0_P /)
      mtot = 0._P

      pl%ah = 0._P

      DO i = 2, npl 
         mtot = mtot + pl%mass(i)
         xbtot =  pl%xh(:,i)*pl%mass(i) + xbtot
         vbtot = pl%vh(:,i)*pl%mass(i) + vbtot
      END DO 
      xbtot =  - xbtot / mtot 
      vbtot = -vbtot / mtot 
      
      pl%xb(:,1) = xbtot * mtot
      pl%vb(:,1) = vbtot * mtot 

      DO i = 2, npl 
         pl%xb(:,i) = pl%xh(:,i) + xbtot
         pl%vb(:,i) = pl%vh(:,i) + vbtot
         ptot = ptot + pl%vb(:,i)*pl%mass(i)
      END DO 
      ptot = ptot * half_timestep * time_conversion / pl%mass(1)

      ! Sun momemtum drift at half_timestep 
      DO i = 2, npl 
         pl%xh(:,i) = pl%xh(:,i) + ptot
      END DO 

      !! kick due to interactions of all bodies except Sun 

      DO j = 2, npl 
         DO k = 2, npl 
            dx(:) = pl%xh(:,k) - pl%xh(:,j)
            rplpl = dot_product(dx,dx)
            if (rplpl > 0.0_P) then 
               pl%ah(:,j) = pl%ah(:,j) + pl%mass(j)*G*dx(:)*(1.0_P/SQRT(rplpl)/rplpl)**(3.0_P)
               pl%ah(:,k) = pl%ah(:,k) - pl%mass(k)*G*dx(:)*(1.0_P/SQRT(rplpl)/rplpl)**(3.0_P)
            end if 
         END DO 
      END DO 

      DO j = 2, npl 
         pl%vb(:,j) = pl%vb(:,j) + pl%ah(:,j)*half_timestep*time_conversion
      END DO 


      !! Danby drift of all bodies  at timestep 
      DO i = 2, npl 
         call drift_danby(pl%xh(:,i),pl%vb(:,i),xnew,vnew, timestep*time_conversion)
         pl%xh(:,i) = xnew 
         pl%vb(:,i) = vnew
      END DO  



      ! Kick at half_timestep
      !! kick due to interactions of all bodies except Sun 
      Do i = 2, npl
         pl%ah(:,i) = (/ 0.0_P, 0.0_P, 0.0_P /)
      end do 

      DO j = 2, npl 
         DO k = 2, npl 
            dx(:) = pl%xh(:,k) - pl%xh(:,j)
            rplpl = dot_product(dx,dx)
            if (rplpl > 0.0_P) then 
               pl%ah(:,j) = pl%ah(:,j) + pl%mass(j)*G*dx(:)*(1.0_P/SQRT(rplpl)/rplpl)**(3.0_P)
               pl%ah(:,k) = pl%ah(:,k) - pl%mass(k)*G*dx(:)*(1.0_P/SQRT(rplpl)/rplpl)**(3.0_P)
            end if 
         END DO 
      END DO 
      DO j = 2, npl 
         pl%vb(:,j) = pl%vb(:,j) + pl%ah(:,j)*half_timestep*time_conversion
      END DO 




      ! Sun momemtum drift at half_timestep and planet drift 
      ptot =  (/ 0.0_P, 0.0_P, 0.0_P /)
      DO i = 2, npl 
         ptot = pl%vb(:,i)*pl%mass(i) + ptot
      END DO
      ptot = ptot * half_timestep*time_conversion / pl%mass(1)

      DO i = 2, npl 
         pl%xh(:,i) = pl%xh(:,i) + ptot
      END DO 
      !convert vb into vh 
      vbtot = (/ 0.0_P, 0.0_P, 0.0_P /)
      DO i = 2 , npl
         vbtot = vbtot * pl%mass(i)*pl%vb(:,i)
      END DO 
      vbtot = vbtot / pl%mass(1) 
      DO i = 2, npl 
         pl%vh(:,i) = pl%vb(:,i) + vbtot
      END DO 
  

      !! output 
      if (MOD(i,cadence) == 0) then 
         DO j = 2, npl
            k = j+8
            write(k,*) time, pl%xh(1,j), pl%xh(2,j), pl%xh(3,j), pl%vh(1,j), pl%vh(2,j), pl%vh(3,j)
         enddo
      endif
      write(*,*)time
      time = time + timestep
      enddo

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


!    subroutine kick(x, v, xnew,vnew)
!       use, intrinsic :: iso_fortran_env
!       implicit none 
!       integer, parameter :: P = real64
!       real(P), dimension(3), intent(in) :: x,v
!       real(P), dimension(3), intent(out) :: xnew, vnew
   

   
!    end subroutine 


   subroutine Etoxv(a,dE,dt,x,v)
      use, intrinsic :: iso_fortran_env

      ! Uses the f and g functions to solve for the cartesian position and velocity as a function of
      implicit none
      integer, parameter :: P = real64
      real(P), parameter :: pi = 3.1415927
      real(P), intent(in) :: a  ! semimajor axis and change in eccentric anomaly
      real(P), intent(in) :: dE ! Change in eccentric anomaly in one step
      real(P), intent(in) :: dt ! Step size
      real(P), dimension(3), intent(inout)  :: x, v ! Position and velocity vectors. Input is x0, v0, output is new x and v
      real(P) :: f, g, fdot, gdot ! f and g functions, and their time derivatives
      real(P) :: r0, r, n
      
      n = sqrt(4.0_P*pi**2.0_P / a**(3.0_P)) ! Mean motion
      r0 = norm2(x(:))

      f = a / r0 * (cos(dE) - 1._P) + 1._P
      g = dt + 1._P / n * (sin(dE) - dE)

      x(:) = f * x(:) + g * v(:)

      r = norm2(x(:))

      fdot = -(a**2.0_P * n / (r * r0) ) * sin(dE)
      gdot = (a / r) * (cos(dE) - 1._P) + 1._P

      v(:) = fdot * x(:) + gdot * v(:)
      return
   end subroutine Etoxv

   subroutine xv2el(x,v,a,e,f,varpi, inc, capom)
      use, intrinsic :: iso_fortran_env

      ! Arguments
      integer, parameter :: P = real64
      REAL(P), PARAMETER :: pi = 3.1415927
      REAL(P)                  :: mu
      REAL(P), DIMENSION(3), INTENT(IN)    :: x, v
      REAL(P), INTENT(OUT)                 :: a, e, f,varpi, inc, capom
         
      ! Internals
      REAL(P)     :: r, v2, hx, hy, hz, h2, h, rdotv, energy, capm, fac, u, w, cw, sw, face, cape, tmpf, capf
      REAL(P)     :: ex, ey, ez
      REAL(P), dimension(3) :: ecc 
      ! Executable code
      mu = 4.0_P*(pi**2.0_P)
      a = 0.0_P
      e = 0.0_P
      inc = 0.0_P
      capom = 0.0_P
      omega = 0.0_P
      capm = 0.0_P
      r = SQRT(DOT_PRODUCT(x(:), x(:)))
      v2 = DOT_PRODUCT(v(:), v(:))
      hx = x(2)*v(3) - x(3)*v(2)
      hy = x(3)*v(1) - x(1)*v(3)
      hz = x(1)*v(2) - x(2)*v(1)
      h2 = hx*hx + hy*hy + hz*hz
      h = SQRT(h2)

      rdotv = DOT_PRODUCT(x(:), v(:))
      energy = 0.5_P*v2 - mu/r
      a = -0.5_P*mu/energy
      fac = hz/h

      IF (fac < -1.0_P) THEN
           inc = pi
      ELSE IF (fac < 1.0_P) THEN
           inc = ACOS(fac)
      END IF
      fac = SQRT(hx*hx + hy*hy)/h
      IF (fac**2 < 1.0E-30_P) THEN
           u = ATAN2(x(2), x(1))
           IF (hz < 0.0_P) u = -u
      ELSE
           capom = ATAN2(hx, -hy)
           u = ATAN2(x(3)/SIN(inc), x(1)*COS(capom) + x(2)*SIN(capom))
      END IF
      IF (capom < 0.0_P) capom = capom + TWOPI
      IF (u < 0.0_P) u = u + TWOPI

      fac = 1.0_P - h2/(mu*a)
      ex = (v(2)*hz - v(3)*hy)/mu - x(1)/r
      ey = (v(3)*hx - v(1)*hz)/mu - x(2)/r
      ez = (v(1)*hy - v(2)*hx)/mu - x(3)/r
      ecc = [ex,ey,ez]
      e = norm2(ecc)

      IF (fac > 1.0E-30_P) THEN
         !e = SQRT(fac)
         cape = 0.0_P
         face = (a - r)/(a*e)
         IF (face < -1.0_P) THEN
            cape = pi
         ELSE IF (face < 1.0_P) THEN
            cape = ACOS(face)
         END IF
       IF (rdotv < 0.0_P) cape = 2.0_p*pi - cape
       fac = 1.0_P - e*COS(cape)
       cw = (COS(cape) - e)/fac
       sw = SQRT(1.0_P - e*e)*SIN(cape)/fac
       w = ATAN2(sw, cw)
       IF (w < 0.0_P) w = w + TWOPI
      ELSE
         cape = u
         w = u
      END IF
      capm = cape - e*SIN(cape)           
          
      omega = u - w
      f = w 
      IF (omega < 0.0_P) omega = omega + TWOPI
      varpi = capom + omega
      RETURN 
         
   END SUBROUTINE xv2el

   subroutine drift_danby(x, v, xnew,vnew, dt)
      use, intrinsic :: iso_fortran_env
      implicit none 
      integer, parameter :: P = real64
      REAL(P), PARAMETER :: pi = 3.1415927
      real(P), dimension(3), intent(in) :: x,v
      real(P), dimension(3), intent(out) :: xnew, vnew
      real(P) :: Period, a, ecc, f, varpi, sE, cE, GM, n
      real(P) :: E, M,  telapsed, dt, inc, capom, alpha, dM, ec, es, esq
      integer :: i ! The number of steps to generate output
      ! Convert cartesian initial conditions to orbital elements 
      call xv2el(x,v,a,ecc,f,varpi,inc,capom)
      GM = 4.0_P*pi**2
      alpha = 2.0_P * GM / norm2(x) - dot_product(v,v)
      a = GM/alpha
      n = sqrt(GM/a**3)
      Period = 2.0_P * pi * sqrt(a**3/GM) 
      ec = 1.0_P - norm2(x)/a 
      es = dot_product(x,x)/(n*a**2)
      esq = ec*ec + es*es
      dM = dt*n - INT(dt*n/2/pi)*2*pi
      sE = SQRT(1.0_P-ecc**2.0_P)*sin(f)/(1+ecc*cos(f))
      cE = (ecc + cos(f)) /(1+ecc*cos(f))
      E = atan2(sE, cE)  
      write(*,*) "Einit =", E
      M = E - ecc * sin(E)
      !M = M + (telapsed / Period) * 2.0_P * pi 
      !M = mod(M+dM, 2*pi)
      ! Update eccentric anomaly (solve Kepler's equation)
      E = danby_solver(dM, ecc, dE )
      Enew = E+dE 
      call Etoxv(E,x,v)
      ! Compute new cartesian position and velocity vectors
      xnew = x 
      vnew = v 
      call Etoxv(a, E, telapsed, xnew, vnew)
      return

      contains 

      double precision function danby_solver(M, ecc) result(E)
         use, intrinsic :: iso_fortran_env
         implicit none
         integer, parameter :: P = real64
         real(P), intent(in) :: M, ecc
         real(P), parameter :: tol = 1e-8_P
         integer, parameter :: imax = 1000
         real(P), parameter :: k = 0.85_P
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
   
      double precision function danby_step(M, ecc, E) result(Enew)
         use, intrinsic :: iso_fortran_env
         implicit none
         integer, parameter :: P = real64
         real(P), intent(in) :: M, ecc, E
         real(P) :: d1, d2, d3
         real(P) :: f, fp, fpp, fppp
   
         f = E - ecc * sin(E) - M
         fp = 1._P - ecc * cos(E)
         fpp = ecc * sin(E)
         fppp = ecc * cos(E)
   
         d1 = -f / fp
         d2 = -f / (fp + d1 * fpp / 2._P)
         d3 = -f / (fp + d2 * fpp / 2._P + d2**2.0_P * fppp / 6._P)
   
         Enew = E + d3
         return
      end function danby_step
   
   end subroutine drift_danby






