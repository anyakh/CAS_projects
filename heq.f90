!******************************************************************
!  this program solves the hydrostatic equilibrium equation
!  for an isothermal gas in a NFW halo
!*****************************************************************
implicit none
integer :: jmax
parameter(jmax=5000)

real*8, dimension(jmax) :: r(jmax),rr(jmax),vol(jmax),mnfw(jmax),&
        rhost(jmax),rho(jmax),mhern(jmax),rhonfw(jmax),mdark(jmax),&
        grvnfw(jmax),lnd(jmax),mgas(jmax),mgas_temp(jmax),&
        fbarr(jmax),fbarr_temp(jmax),fgasr(jmax),fgasr_temp(jmax),&
        tr(jmax),lndt(jmax),consttemp(jmax),rho_temp(jmax)
        
real*8 :: msol,mu,mp,rmin,rmax,mvir,rvir,mbcg,ahern,boltz,cmkpc,fbar,&
        fbar_temp,fc,fgas,fgas_temp,gg,guniv,j,pi,r500,rho0,rho0nfw,&
        rho0t,rs,ticm,x,zfeout,zfesn,zfesol

real*8, dimension(jmax) :: u(jmax),flux(jmax),ne(jmax),zfe(jmax),& 
        kappa(jmax),lturb,rhofedot(jmax),rhofe(jmax),zfest(jmax),&
        amfeiniz(jmax),amfe(jmax),gradzfe(jmax),zfeobs(jmax),&
        amfeobs(jmax),rhofeobs(jmax)

!  constants

pi = 3.14159265359
msol = 1.989d33
cmkpc = 3.084e21
mu=0.62
boltz=1.38066e-16
guniv=6.6720e-8
mp=1.67265e-24
zfesol=1.8e-3
zfesn=0.744/1.4

!    set the grid

rmin = 0.*cmkpc
rmax = 2800.*cmkpc ! 6-11-25: 3000 -> 2800.
!       set 1st grid
do j=1,jmax
   r(j)=rmin+(j-1)*rmax/(jmax-1)
enddo
!       set 2nd grid; r_(j+1/2) = rr_j
do j=1,jmax-1
   rr(j)=r(j)+0.5*(r(j+1)-r(j))
enddo
rr(jmax)=rr(jmax-1)+(rr(jmax-1)-rr(jmax-2))
open(10,file='grid.dat',status='unknown')
do j=1,jmax
   write(10,*)real(r(j)/cmkpc),real(rr(j)/cmkpc)
enddo
close(10)

!       volume b/t 2 shells of radii r_j and r_j+1:
!       DELTA V_j+1/2 = 4/3 pi (r^3_j+1 - r^3_j)
vol(1)=4.1888*r(1)**3
do j=2,jmax
   vol(j)=4.1888*(r(j)**3-r(j-1)**3)    !! centered on rr(j-1) !!
enddo

!  problem parameters

rho0nfw=7.35d-26        !! initial DM density !!
rs=435.7*cmkpc          !! stellar radius (???) !!
!! ~~~~~~~Need to mess around with rho0 more to get fbari~0.16 ~~~~~ !!
rho0=0.9d-25            !! central density; initial condition for gas density !! 
rho0t=1.7d-25            !! central density; initial condition for gas density !! 
!!rho0=4.d-26   !! 5000 points, with BCG, isothermal gas !!old: 2.882d-26
!!rho0=9.d-26   !! 5000 points, with BCG !!
ticm=8.9e7              !! temperature of ICM !! 

rvir=2797.*cmkpc        !! virial radius !!
r500=rvir/2.            !! r500... _______________ !!
fc=1.138799             !!
mvir=1.3e15*msol        !! virial mass !!
mbcg=1.d12*msol         !! mass of bcg !!
ahern=12.*cmkpc/(1.+2.**(0.5))  !! a from Hernquist profile 6-11-25: 10 -> 12 !!

do j=1,jmax             !! Dark matter density profile !!
   x=rr(j)/rs
   rhonfw(j)=rho0nfw/(x*(1.+x)**2)
   rhost(j)=mbcg/(2.*pi)*(ahern/rr(j))/(ahern+rr(j))**3 !! Density including BCG !!
enddo

open(20,file='masse.dat')       !! Calculate mass  !!
mnfw(1)=0.
do j=2,jmax
   x=r(j)/rs
   mnfw(j)=mnfw(j-1)+rhonfw(j-1)*vol(j)         !! Numerical !!
   mdark(j)=mvir*(log(1.+x)-x/(1.+x))/fc        !! Analytical !!
   mhern(j)=mbcg*r(j)**2/(r(j)+ahern)**2        !! Hernquist profile !!
   write(20,1001)r(j)/cmkpc,mnfw(j)/msol,mdark(j)/msol,mhern(j)/msol
enddo
1001 format(4(1pe12.4))
close(20)


open(20,file='grv.dat')
grvnfw(1)=0.          !! ok for halo NFW, isotherm or beta-model 
do j=2,jmax             !! Potential as function of radius (?) !!
   grvnfw(j)=guniv*(mnfw(j)+mhern(j))/r(j)**2
   write(20,1002)r(j)/cmkpc,grvnfw(j)/msol
enddo
1002 format(2(1pe12.4))
close(20)

!     Calculating T profile for non-constant temp -> THIS IS WHERE WE LEFT OFF 11-11-25 !!
!  NEXT STEP: PLOT TR AND COMPARE TO TICM (SEE SLIDE 15 CAS PROJECT1 SLIDES !
do j=1,jmax
   x=rr(j)/r500
   consttemp(j)=ticm
   tr(j)=1.35*((((x/0.045)**1.9)+0.45)/(((x/0.045)**1.9)+1.))*(1./(1.+((x/0.6)**2.))**0.45)*ticm
enddo
! Temperature file !
open(20,file='temperature.dat',status='unknown')
do j=1,jmax
   write(20,1000)rr(j)/cmkpc,tr(j),consttemp(j)
enddo

!     calculate the gas density, assuming ticm

lnd(1)=log(rho0)          !! puts the gas in eq. with the potential
lndt(1)=log(rho0t)         !! this doesnt do anything YET !!
do j=2,jmax
   gg=grvnfw(j)
   lnd(j)=lnd(j-1)-gg*(mu*mp)*(rr(j)-rr(j-1))/(boltz*ticm) !! original, constant temp !!
   lndt(j)=lndt(j-1)-gg*(mu*mp)*(rr(j)-rr(j-1))/(boltz*tr(j)) - log(tr(j))+log(tr(j-1)) !! update with additional terms for varying temp !!
enddo

do j=1,jmax
   rho(j)=exp(lnd(j))
enddo

do j=1,jmax
   rho_temp(j)=exp(lndt(j))
enddo

open(20,file='density.dat',status='unknown')
do j=1,jmax
   write(20,1000)rr(j)/cmkpc,rho(j),rhonfw(j),rhost(j),rho_temp(j) !! grid, gas density, dark matter density, stellar density !!
enddo
1000 format(5(1pe12.4))
close(20)

!! Calculates mgas and writes to file so we can look at it !!
open(20,file='mgas.dat',status='unknown')
mgas(1)=rho(j)*4.188*r(1)**3
do j=2,jmax
   mgas(j)=mgas(j-1)+rho(j-1)*vol(j)
   write(20,1100)r(j)/cmkpc,mgas(j)/msol,mnfw(j)/msol
enddo
close(20)

!! Calculates mgas nonisothermal and writes to file so we can look at it !!
open(20,file='mgast.dat',status='unknown')
mgas_temp(1)=rho_temp(j)*4.188*r(1)**3
do j=2,jmax
   mgas_temp(j)=mgas_temp(j-1)+rho_temp(j-1)*vol(j)
   write(20,1100)r(j)/cmkpc,mgas_temp(j)/msol,mnfw(j)/msol
enddo
close(20)


!! Baryon & Gas Fraction; NOT a function of radius (yet)
!!fbar=mgas(jmax-1)/(mnfw(jmax-1)+mgas(jmax-1)) unclear why this is here
fbar=(mhern(jmax-1)+mgas(jmax-1))/(mnfw(jmax-1)+mgas(jmax-1)+mhern(jmax-1))
fgas=(mgas(jmax-1))/(mnfw(jmax-1)+mgas(jmax-1)+mhern(jmax-1))
print*,'fbari, fgas = ',real(fbar), real(fgas)

fbar_temp=(mhern(jmax-1)+mgas_temp(jmax-1))/(mnfw(jmax-1)+mgas_temp(jmax-1)+mhern(jmax-1))
fgas_temp=(mgas_temp(jmax-1))/(mnfw(jmax-1)+mgas_temp(jmax-1)+mhern(jmax-1))
print*,'fbari, fgas_temp = ',real(fbar_temp), real(fgas_temp)


open(20,file='barfrac.dat',status='unknown')
do j=2,jmax-1
   fbarr(j)=(mhern(j)+mgas(j))/(mnfw(j)+mgas(j)+mhern(j))
   fgasr(j)=(mgas(j))/(mnfw(j)+mgas(j)+mhern(j))
   write(20,1100)r(j)/cmkpc,fbarr(j),fgasr(j)
enddo
close(20)
1100 format(4(1pe12.4))

!***********************************************************************
!! At this point we have the gas density profile and we can proceed
!! with the integration of the diffusion equation for rhofe
!***********************************************************************

!! Set the initial abundance profile

zfeout=0.4*zfesol   !! this is the background abundance !!

do j=1,jmax
   x=rr(j)/(80.*cmkpc)   !!smaller step!!
   !zfeobs(j)=zfesol*0.3*1.4*1.15*(2.2+x**3)/(1+x**3)/1.15  !! observed abundance in Perseus BUT where the 1.15 comes from
   zfeobs(j)=zfesol*0.3*(2.2+x**3)/(1+x**3) !! should be like this
   !zfeobs(j)=zfeobs(j) - zfeout   !! subtract background metallicity
   zfeobs(j)=max(zfeobs(j),0.)   !! Check if observed metallicity isnt less than 0
   !zfe(j)=0. !!zfeout !!zfeobs(j)  !! which initial zfe? !!
   zfe(j)=zfeobs(j) - zfeout !! I guess it should be this one, but not sure
   rhofe(j)=rho(j)*zfe(j)/1.4  !! Just density of the Fe 
   rhofeobs(j)=rho(j)*zfeobs(j)/1.4 !! Obs means observed
enddo

do j=1,jmax
    zfest(j)=1.*zfesol    !! set the stellar abundance !!
enddo

!! Calculate the initial excess of iron mass

amfeiniz(1)=rhofe(1)*vol(1)
amfeobs(1)=rhofeobs(1)*vol(1)
do j=2,jmax
   amfeiniz(j)=amfeiniz(j-1)+rhofe(j-1)*vol(j)
   amfeobs(j)=amfeobs(j-1)+rhofeobs(j-1)*vol(j)
enddo

open(20,file='zfe_initial.dat')
do j=1,jmax
   write(20,1500)rr(j)/cmkpc,zfe(j)/zfesol,zfeobs(j)/zfesol, &
                 r(j)/cmkpc,amfeiniz(j)/msol,amfeobs(j)/msol
enddo
close(20)
1500 format(6(1pe12.4))

open(20,file='initial.dat',status='unknown')
do j=1,jmax
   write(20,3001)rr(j)/cmkpc+0.001,zfe(j)/zfesol,zfeobs(j)/zfesol,ne(j)
enddo
close(20)
3001  format(4(1pe12.4))


!! At this point we created the no sorce only diffusion profile

stop
end
