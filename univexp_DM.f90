! f95 ../univexp_DM.f90  -o univexp_DM.e -fdefault-real-8 -fdefault-integer-8

program univexp
 implicit none
 character (len=40) :: fichier
 integer :: n,m,ifirst,ilast
 real :: tecr=0.,tsor=0.,tstop,dti,dtsor
 real :: pvit,epolar,eav,eap,gravplas
 real, dimension(:), allocatable :: x,v,mi,ma
 integer, dimension(:), allocatable :: name

 n=10000
 fichier='univexp01'

 call init

 open(unit=500,file=trim(fichier)//'a.d',form='formatted',status='replace')
 open(unit=510,file=trim(fichier)//'d.d',form='formatted',status='replace')

 call wbande

 tecr=tecr+dti
 tsor=tecr+dtsor

 do while(abs(tecr-tstop).gt.dtsor/2.)
  do while(abs(tecr-tsor).gt.dti/2.)
   call avance
   call ordonne
   tecr=tecr+dti
  enddo
  call wbande
  tsor=tsor+dtsor
 enddo

 write(*,'(/" fin du programme."/)')
 stop

contains
!**************************************************************************
 subroutine avance
 implicit none
 integer :: i
 real, parameter :: r2=-sqrt(2.),tier=1./3.
 real :: AA,BB,E
 real, dimension(1:m) :: a

! calcul de l'accÃ©lÃ©ration
 a(:)=0
 eav=epolar+0.5*n
 do i=ifirst,ilast
!  eap=eav-gravplas*ma(i)
!  a(i)=0.5*(eav+eap)/mi(i)
  eap=eav-1.
  a(i)=0.5*(eav+eap)
  eav=eap
 enddo

! les particules avancent de dti
 do i=ifirst,ilast
  E=a(i)
  AA=(x(i)+E+r2*v(i))*exp(r2*dti)
  BB=2.*(x(i)+E-v(i)/r2)*exp(-dti/r2)
  x(i)=((AA+BB)*tier)-E
  v(i)=(AA*r2-BB/r2)*tier
  if(x(i)>x(m)) then
!   epolar=epolar+gravplas*ma(i)
   epolar=epolar+1.
   x(i)=x(1)+x(i)-x(m)
 endif
  if(x(i)<x(1)) then
!   epolar=epolar-gravplas*ma(i)
   epolar=epolar-1.
   x(i)=x(m)+x(i)-x(1)
  endif
 enddo

 return
 end subroutine avance
!***************************************************************************
 subroutine ordonne
 implicit none
 integer :: i,j,k,np
 real :: xp,vp,mip,map

  do i=ifirst+1,ilast
   j=i
   xp=x(i)
   vp=v(i)
!   mip=mi(i)
!   map=ma(i)
   np=name(i)
   do while(j/=1.and.x(i)<x(j-1))
    j=j-1
   enddo
   do k=i,j+1,-1
    x(k)=x(k-1)
    v(k)=v(k-1)
!    mi(k)=mi(k-1)
!    ma(k)=ma(k-1)
    name(k)=name(k-1)
   enddo
   x(j)=xp
   v(j)=vp
!   mi(j)=mip
!   ma(j)=map
   name(j)=np
  enddo

 return
 end subroutine ordonne
!**************************************************************************
 subroutine wbande
 implicit none
 integer :: i

 write(*,'(" enregistrement a tecr=",f7.3)')tecr
 write(510,'("#",4i5,f7.3)')m,n,ifirst,ilast,real(tecr)
 do i=ifirst,ilast
  write(510,*)real(x(i)),real(v(i)),name(i)
 enddo
 write(510,'(/)')

 return
 end subroutine wbande
!************************************************************************
 subroutine init
 implicit none
 integer :: i
 real    :: vmoy

! calcul de la longueur d'un enregistrement : longrec
! longrec =2*2**karbre -1 +2 = 262145 * 3 *8 = 3145752 = arrondi 3145800 + 25 * 8 (= 200) = 3146000

 m=n+2          !! => tc(m) sera introduit dans la table de tri
 ifirst=2
 ilast=n+1
 dti=.001
 dtsor=1.
 tstop=15.
 gravplas=1.
 pvit=100.

 allocate(x(1:m),v(1:m),mi(1:m),ma(1:m),name(1:m))

 write(*,'("simulation :",a)')fichier
 write(*,'(f19.15," : pas de temps"/)')dti
 write(*,'(" n=",i8," pvit=",f7.1)')n,pvit
 write(*,'(" tstop=",f7.1," dtsor=",f7.1)')tstop,dtsor

! initialisation des "particules-mur"

 x(1)=-.5*float(n)
 v(1)=0
 mi(1)=0
 ma(1)=0
 name(1)=0
 x(m)=.5*float(n)
 v(m)=0
 mi(m)=0
 ma(m)=0
 name(m)=0

! initialisation des particules

 do i=ifirst,ilast
  name(i)=i-1
  x(i)=x(1)+.5+i-ifirst
  v(i)=2.*pvit*(0.5d0-rand( ))
  mi(i)=1.
  ma(i)=1.
 enddo

 epolar=0

! calage du barycentre a zero avec une vitesse moyenne nulle

 vmoy=sum(v(ifirst:ilast))
 vmoy=vmoy/float(n)
 v(ifirst:ilast)=v(ifirst:ilast)-vmoy
 vmoy=sum(v(ifirst:ilast))
 write(*,'(" vmoyen =",f7.3)')vmoy

 return
 end subroutine init
!************************************************************************
 end program univexp

