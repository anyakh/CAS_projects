parameter(jmax=10000)
implicit real*8 (a-h,o-z)
real*8, dimension(jmax) :: r(jmax),rr(jmax),vol(jmax),mnfw(jmax), rhost(jmax),rho(jmax),checkrho(jmax),mhern(jmax),rhonfw(jmax),mdark(jmax),mdarkkkk(jmax),grvnfw(jmax),lnd(jmax),mgas(jmax),mgasstar(jmax),mtot(jmax),fb(jmax),rhostar(jmax)
real*8 :: msol,mu,mp,rmin,rmax,mvir,rvir,mbgc,ahern,b
! temperature profile
real*8, dimension(jmax) :: T(jmax), xx(jmax)
real*8 :: r500, Tmg, T_tmp
!! part II - ZFe
real*8, dimension(jmax) :: rhoFe(jmax), rhoFedot(jmax),rhoFe_obs(jmax), ZFe(jmax), ZFe_obs(jmax), ZFe_st(jmax),mFe_iniz(jmax), mFe_obs(jmax), flux(jmax), D(jmax), gradZFe(jmax), mFe(jmax)
real*8 :: ZFe_out, ZFe_sol, time, tnow, t0, fac_year, vt, lt, rhojp1, rhoj, t_targer, tar, alpha_SNIa, alpha_st,ZFe_SNIa, slope, SNu

!!! ---------- constants ----------------
msol = 1.989d33 ! m of sun
cmkpc = 3.084e21 ! cm in 1 kpc
mu=0.62
boltz=1.38066e-16
guniv=6.6720e-8 ! G
mp=1.67265e-24

!!! ---------- set the grid ----------------
rmin = 0.*cmkpc
rmax = 3000.*cmkpc
do j=1,jmax
   r(j)=rmin+(j-1)*rmax/(jmax-1)
enddo
do j=1,jmax-1
   rr(j)=r(j)+0.5*(r(j+1)-r(j))
enddo
rr(jmax)=rr(jmax-1)+(rr(jmax-1)-rr(jmax-2))

open(10,file='grid.dat',status='unknown')
do j=1,jmax
   write(10,*)real(r(j)/cmkpc),real(rr(j)/cmkpc)
enddo
close(10)

!!! ---------- calculate the mass of DM ----------------
vol(1)=4.1888*r(1)**3
do j=2,jmax
   vol(j)=4.1888*(r(j)**3-r(j-1)**3)  ! delta Vj+1/2  !! centered at r(j-1/2) =  rr(j-1) !! ??????????
enddo

! these parameters of NFW profile depend on the cluster mass
rho0nfw=7.35d-26 ! rhodm0, density of DM
rs=435.7*cmkpc
mvir=1.3e15*msol ! virial mass, slides for 1.2e15 msol
rvir=2797.*cmkpc ! virial radius, slides for 2.8 Mpc
fc=1.138799 ! normalization factor for NFW profile, 使积分质量与mvir一致

!!!!! Mstar parameters for BCG - Hernquist profile
mbgc=1.d12*msol ! total mass of BCG
ahern=12.*cmkpc/(1.+2.**(0.5)) !  12 kpc (half mass radius)



! DM density profile - NFW
do j=1,jmax
   x=rr(j)/rs ! centered at r(j+1/2) = rr(j)
   rhonfw(j)=rho0nfw/(x*(1.+x)**2)
enddo

open(20,file='mass.dat')
mnfw(1)=0. ! I.C.
do j=2,jmax
   x=r(j)/rs  ! why this use r(j) instead of rr(j)???????---> formula of mnfw at j+1/2 equal to this
   mnfw(j)=mnfw(j-1)+rhonfw(j-1)*vol(j) ! cumulative mass profile of NFW
   mdark(j)=mvir*(log(1.+x)-x/(1.+x))/fc ! analytical mass profile of NFW
   mdarkkkk(j) = 4.*3.14159*rho0nfw*rs**3* (log(1.+x)-x/(1.+x)) ! check the analytical mass of NFW (mdark)
   mhern(j)=mbgc*r(j)**2/(r(j)+ahern)**2
   write(20,1005)r(j)/cmkpc,mnfw(j)/msol,mdark(j)/msol,mhern(j)/msol,mdarkkkk(j)/msol
enddo
1005 format(5(1pe12.4))
close(20)

!!! ---------- calculate the gravitational acceleration ----------------
open(20,file='grv.dat')
grvnfw(1)=0.          !! ok per alone NFW, isotermo o beta-model
do j=2,jmax
   ! grvnfw(j)=guniv * mnfw(j)/r(j)**2 ! only mDM
   grvnfw(j)=guniv*(mnfw(j)+mhern(j))/r(j)**2 ! mDM and BCG
   write(20,1002)r(j)/cmkpc,grvnfw(j)/msol
enddo
1002 format(2(1pe12.4))
close(20)

!!! ---------- before calculating the gas density profile, estimate the temperature gradient dT/dr ----------------
r500 = 1400.*cmkpc ! r500 ~ rvir/2 for massive clusters
Tmg = 8.9e7 
do j=1,jmax
   xx(j)=rr(j)/r500 ! x, centered at r(j+1/2) = rr(j)
enddo

open(20,file='temperature.dat')
do j = 1,jmax
   T(j) = Tmg * 1.35 * ((xx(j)/0.045)**1.9 + 0.45) / ((xx(j)/0.045)**1.9 + 1.) / (1. + (xx(j)/0.6)**2)**0.45
   write(20,1002)rr(j)/cmkpc,T(j)
enddo

!!! ---------- calculate the gas density profile in hydrostatic equilibrium ----------------

! these parameters is used to calculate the gas density profile in hydrostatic equilibrium

! rho0=2.882d-26 ! (old)rho0, central density of gas, must be cahnged!!!
! rho0 = 4.0e-26 ! only isothermal gas 
! rho0=9.d-26   ! with BCG !!
rho0=1.7d-25   ! with BCG and dT/dr !! 1.6d-25 get fb = 0.15

Ticm=8.9e7 ! isothermal gas temperature in K

lnd(1)=log(rho0)          ! I.C. using rho0 at the center -- rhog(0) = rho0

! !!! isothermal gas
! do j=2,jmax
!    gg=grvnfw(j) ! GM/rr2 at rr(j)
!    lnd(j)=lnd(j-1)-gg*(mu*mp)*(rr(j)-rr(j-1))/(boltz*Ticm) ! ln rho_j+1/2
! enddo

!! non-isothermal gas with dT/dr
do j=2,jmax
   gg=grvnfw(j) ! GM/rr2 at rr(j)
   lnd(j)=lnd(j-1)-gg*(mu*mp)*(rr(j)-rr(j-1))/(boltz*T(j)) - (log(T(j)) - log(T(j-1)))  ! ln rho_j+1/2
enddo

do j=1,jmax
   rho(j)=exp(lnd(j)) ! transform back to rho
enddo

! open(20,file='density_no_BCG.dat',status='unknown')
! do j=1,jmax
!    write(20,1003)rr(j)/cmkpc,rho(j),rhonfw(j) ! compare gas density with DM density
! enddo
! close(20)

! open(20,file='density_BCG.dat',status='unknown')
! do j=1,jmax
!    write(20,1003)rr(j)/cmkpc,rho(j),rhonfw(j) ! compare gas density with DM density
! enddo
! close(20)

open(20,file='density_BCG_gradT.dat',status='unknown')
do j=1,jmax
   write(20,1003)rr(j)/cmkpc,rho(j),rhonfw(j) ! compare gas density with DM density
enddo
close(20)

1003  format(3(1pe12.4))

!!! ---------- fb-Method1: integrate the gas density profile to get the gas mass profile ----------------
mgas(1)=rho(1)*4.188*r(1)**3 ! I.C.
do j=2,jmax
   mgas(j)=mgas(j-1)+rho(j-1)*vol(j) ! cumulative mass profile of gas, but the vol(j) is centered at r(j)
   mgasstar(j)= mgas(j) + mhern(j) ! total baryon mass profile = gas mass + stellar mass
   mtot(j)= mgas(j) + mhern(j) + mnfw(j) ! total mass profile = gas mass + stellar mass + DM mass
   fb(j)= mgasstar(j)/mtot(j) ! baryon fraction profile
enddo

! open(20,file='gasmass_iso.dat',status='unknown')
! do j=1,jmax
!    write(20,1005)r(j)/cmkpc,mgas(j)/msol,mgasstar(j)/msol,mtot(j)/msol,fb(j) ! mass profile
! enddo
! close(20)   
! 1005 format(5(1pe12.4))

! open(20,file='gasmass_BCG.dat',status='unknown')
! do j=1,jmax
!    write(20,1005)r(j)/cmkpc,mgas(j)/msol,mgasstar(j)/msol,mtot(j)/msol,fb(j) ! mass profile
! enddo
! close(20)   
! 1005 format(5(1pe12.4))

open(20,file='gasmass_BCG_gradT.dat',status='unknown')
do j=1,jmax
   write(20,1005)r(j)/cmkpc,mgas(j)/msol,mgasstar(j)/msol,mtot(j)/msol,fb(j) ! mass profile
enddo
close(20)   



!!!---------- check the numerical solotion for rho(r) with the analytic formula (ONLY ISOTHERMAL CASE) ----------------
! do j = 1,jmax
! b = 8 * 3.14159 * guniv * mu * mp * rho0nfw * rs**2 / (27*boltz * Ticm)
! checkrho(j) = rho0 * exp(-27*b/2) * (1+rr(j)/rs)**(27*b/(2*rr(j)/rs)) ! notice that rho is defined on rj+1/2, using rr(j)
! enddo
! open(20,file='checkrho.dat',status='unknown')
! do j=1,jmax
!    write(20,1003)rr(j)/cmkpc,rho(j),checkrho(j) ! compare numerical rho with analytic rho 
! enddo
! close(20)


!!! ---------- Addendum ----------------


!!! ----------- part II - ZFe ----------------
!!! start by setting the initial condition fot ZFe and rhoFe.
ZFe_sol=1.8e-3  ! solar iron abundance by mass
ZFe_out=0.4*ZFe_sol   !! this is the background abundance, which is from ICM outside the cluster

! initial values
do j=1,jmax
   ZFe(j) = 0.d0
   rhoFe(j)= 0.d0
enddo

do j=1,jmax
   x=rr(j)/(80.*cmkpc)
   ZFe_obs(j)=ZFe_sol*0.3*1.4*(2.2+x**3)/(1+x**3)  !Perseus profile!
   ZFe_obs(j)=ZFe_obs(j) - ZFe_out   ! subtract ZFe_out, since the background is constant
   ZFe_obs(j)=max(ZFe_obs(j),0.d0) ! set mininum abundance to zero
   ! ZFe(j)=ZFe_obs(j)!!zfeout !!zfeobs(j)  ! select initial ZFe, or do not use this line
   rhoFe(j)=rho(j)*ZFe(j)/1.4
   rhoFe_obs(j)=rho(j)*ZFe_obs(j)/1.4
enddo

do j=1,jmax
   ZFe_st(j)=1.*ZFe_sol  ! set the stellar abundance 
enddo

!!! Calculate the initial excess of iron mass
mFe_iniz(1)=rhoFe(1)*vol(1)
mFe_obs(1)=rhoFe_obs(1)*vol(1)
do j=2,jmax
   mFe_iniz(j)=mFe_iniz(j-1)+rhoFe(j-1)*vol(j)
   mFe_obs(j)=mFe_obs(j-1)+rhoFe_obs(j-1)*vol(j)
enddo

open(20,file='ZFe_initial.dat')
do j=1,jmax
   write(20,1006)rr(j)/cmkpc, ZFe(j)/ZFe_sol, ZFe_obs(j)/ZFe_sol, r(j)/cmkpc, mFe_iniz(j)/msol, mFe_obs(j)/msol
   !!! This plot need coordinate r(j) instead of rr(j)！！！
enddo
close(20)
1006 format(6(1pe12.4))

!!! boundary conditions
ZFe(1)=ZFe(2)
ZFe(jmax)=ZFe(jmax-1)
rhoFe(1)=rho(1)*ZFe(1)/1.4
rhoFe(jmax)=rho(jmax)*ZFe(jmax)/1.4

!!! Use FTCS method to integrate time
fac_year=3.156d7 ! seconds in one year

tnow = 13.7*1.e9*fac_year
t0 = tnow-5.*1.e9*fac_year
time = t0
tend=tnow+5.*1.e9*fac_year
tar = 5 ! set your target time in Gyr, Default is 1 Gyr
t_target = t0 + tar*1.e9*fac_year

vt = 260.e5 ! Perseus cluster
lt = 15.*cmkpc ! uncertainty?  Poor theory and observations????
D = 0.11 * vt * lt ! for all D(j) 


! calculate the timestep
dt = 0.4 * (r(5)-r(4))**2 / (2*D(5)) ! stability condition for diffusion equation
! time = time + dt


n = 0 

! output diffusion at t = t0 (initial) before main loop
! set the boundary conditions (outflows)
! ZFe(1)=ZFe(2)
! ZFe(jmax)=Zfe(jmax-1)
! ! rhoFe(1)=rhoFe(2)
! ! rhoFe(jmax)=rhoFe(jmax-1)
! open(20, file="diffusion_initial.dat")
! do j=1,jmax
!    write(20,1002)rr(j)/cmkpc,ZFe(j)/ZFe_sol
! enddo
! close(20)
! print*, 'Output initial Fe abundance'

do while (time .lt. tend)
   n = n + 1
   time = time + dt

   !!! source term (SNIa (dominated) + stellar winds)---------
   slope = 1.1
   ZFe_SNIa = 0.74/1.4  ! iron yield per SNIa in solar units
   ! alpha with time dependence
   SNu = 0.15
   alpha_SNIa = 5.92e-21 * SNu* (time/tnow)**(-slope) ! this is the simplified result for alpha_SNIa
   alpha_st = 4.7e-20 * (time/tnow)**(-1.26)
   ! use constant alpha and multiply iron abundance
   const_SNIa = 4.7e-22
   const_st = 6.e-23

   ! calculate rho_st
   do j=2,jmax-2
      rhostar(j) = mbgc * ahern / ((rr(j)+ahern)**3 * 2.*3.14159 * rr(j)) ! Hernquist stellar density profile
   enddo

   ! calculate rhoFe dot without time dependence (const)
   do j=2, jmax-2
      rhoFedot(j) = rhostar(j) * (const_st + const_SNIa) 
   enddo

   ! ! calculate rhoFe dot with time dependence
   ! do j=2, jmax-2
   !    rhoFedot(j) = rhostar(j) * (alpha_st * ZFe_st(j) / 1.4 + alpha_SNIa * ZFe_SNIa) 
   ! enddo

   ! calculate the rhoFe
   do j=2,jmax-2
      rhoFe(j)=rhoFe(j) + dt*rhoFedot(j)
      ZFe(j)=rhoFe(j)/rho(j) * 1.4
   enddo

   ! set the boundary conditions (outflows)
   ZFe(1)=ZFe(2)
   ZFe(jmax)=Zfe(jmax-1)
   rhoFe(1)=rhoFe(2)
   rhoFe(jmax)=rhoFe(jmax-1)

   ! ! output source in target time
   ! if (abs(time - t_target) .lt. dt/2.) then
   !    open(20,file='source5.dat',status='unknown')
   !    do j=1,jmax
   !       write(20,1002)rr(j)/cmkpc,ZFe(j)/ZFe_sol
   !    enddo
   !    close(20)
   !    print*, 'Output Fe abundance at target time with only source term'
   ! calculate Fe mass
   !    mFe(1) = rhoFe(1)*vol(1)
   !    do j=2,jmax
   !       mFe(j) = mFe(j-1) + rhoFe(j-1)*vol(j)
   !    enddo
   ! endif


   !!! diffusion term ---------
   ! calculate the gradient of ZFe
   gradZFe(1)=0.0 ! B.C.
   gradZFe(jmax)=0.0
   do j = 2,jmax-1
      gradZFe(j) = (ZFe(j+1)-ZFe(j-1)) / (rr(j+1)-rr(j-1)) ! centered at r = r(j)
   enddo

   ! diffusion evolution with time
   do j=2,jmax-1
      rhojp1 = (rho(j+1)+rho(j))/2  ! rhoj+1, since rho centered at r(j+1/2) = rr(j)
      rhoj = (rho(j-1)+rho(j))/2  ! rhoj-1
      rhoFe(j) = rhoFe(j) + (dt/1.4) * ( r(j+1)**2 * D(j+1) * rhojp1 * gradZFe(j+1) - r(j)**2 * D(j) * rhoj * gradZFe(j) ) / ((r(j+1)**3 - r(j)**3)/3)
      ZFe(j) = 1.4 * rhoFe(j) / rho(j)  ! update ZFe with time
      ! print*, 'n, j, r(j), ZFe(j)= ', n, j, r(j)/cmkpc, ZFe(j)/ZFe_sol
   enddo

   ! set the boundary conditions (outflows)
   ZFe(1)=ZFe(2)
   ZFe(jmax)=Zfe(jmax-1)
   rhoFe(1)=rhoFe(2)
   rhoFe(jmax)=rhoFe(jmax-1)

      ! output diffusion in target time
      ! if (abs(time - t_target) .lt. dt/2.) then
      !    open(20,file='diffusion2.dat',status='unknown')
      !    do j=1,jmax
      !       write(20,1002)rr(j)/cmkpc,ZFe(j)/ZFe_sol
      !    enddo
      !    close(20)
      !    print*, 'Output Fe abundance at target time with only diffusion term'
      !    mFe(1) = rhoFe(1)*vol(1)
      !    do j=2,jmax
      !       mFe(j) = mFe(j-1) + rhoFe(j-1)*vol(j)
      !    enddo
   ! endif

   ! -------------------------------------------------------------
   !output diffusion + source(const) in target time
   if (abs(time - t_target) .lt. dt/2.) then
      ! transport tar into string
 
      open(20,file='d+s5.dat',status='unknown')
      do j=1,jmax
         write(20,1002)rr(j)/cmkpc,ZFe(j)/ZFe_sol
      enddo
      close(20)
      print*, 'Output Fe abundance at target time with diffusion + source (real situation)'
      mFe(1) = rhoFe(1)*vol(1)
      do j=2,jmax
         mFe(j) = mFe(j-1) + rhoFe(j-1)*vol(j)
      enddo
   endif

   ! !!! ----------output diffusion + source(t) in target time ----------------
   ! if (abs(time - t_target) .lt. dt/2.) then
   !    ! transport tar into string
 
   !    open(20,file='t_d+s5.dat',status='unknown')
   !    do j=1,jmax
   !       write(20,1002)rr(j)/cmkpc,ZFe(j)/ZFe_sol
   !    enddo
   !    close(20)
   !    print*, 'Output Fe abundance at target time with diffusion + source (real situation)'
   !    mFe(1) = rhoFe(1)*vol(1)
   !    do j=2,jmax
   !       mFe(j) = mFe(j-1) + rhoFe(j-1)*vol(j)
   !    enddo
   ! endif

   
   
enddo

print*, 'The theoretical Fe mass at the target time is ', 5.3e-22*1e12*(t_target-t0), ' Msol'
print*, 'The numerical Fe mass at the target time is ', mFe(jmax)/msol, ' Msol'



1001  continue
end