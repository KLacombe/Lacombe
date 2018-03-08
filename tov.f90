MODULE constants

  integer,parameter :: dp = selected_real_kind(15)
  integer, parameter :: npmax=10000
  
  real(dp) :: pi,hbc,hbc3,rho0,gn,msol,convertinertia
  integer :: nr,wdflag,npc,ntab
  real(dp) :: dr,pcmin, dpc, p_crit,a0 = 1._dp ! constant for H and beta
  real(dp)  :: rtab(npmax),ptab(npmax),epstab(npmax),rhostab(npmax)

  data pi /3.14159265359_dp/
  data hbc /197.33_dp/ 
  data gn / 1.3234923e-12_dp/ ! units are fm^3/(MeV m^2)
  data msol / 1.1156912e15_dp/ ! units are MeV m^3/fm^3 
  data convertinertia / 1.782661731e-23_dp/ ! to convert units from m^5 MeV/fm^3 to 10^38 m² kg

end MODULE constants

program tov

  USE constants
  IMPLICIT NONE

  real(dp) :: pc,emk,mkalt,pcmax
  integer :: wahl_pressure,ipc
  logical first

  call input

  write(*,*) 'Use given pressure (1) or do a sequence (2)? '
  read(*,*) wahl_pressure 

  IF(wahl_pressure.eq.1) then
     write(*,*) 'Which pressure (in MeV/fm^3)?'
     read(*,*) pc
     call integrate_tov(pc,emk)
     write(*,*) 'Resulting mass (in solar masses)',emk/msol
      
  else
     mkalt = 0.d0
     first = .true.
     pcmax = log(pcmin + (npc-1)*dpc)

     do ipc = 1,npc
        IF(wdflag.eq.0) then
           pc = exp(log(pcmin) + (ipc-1)*(pcmax-log(pcmin))/(npc-1))
        else
           pc = pcmin + (ipc-1)*dpc
        end if
        call integrate_tov(pc,emk)
     
        IF((mkalt.gt.emk).and.(emk/msol.gt.1.4).and.first) then
           write(*,*) 'maximum mass (in solar masses)',emk/msol
           first = .false.
        endif
        mkalt = emk


     end do
  end IF
  
end program tov

 ! RK2 integration of TOV equations


subroutine integrate_tov(pc,emk)

  USE constants
  IMPLICIT NONE

! units:
! pressure: MeV/fm^3
! energy density : MeV/fm^3
! distances: m
! densities: 1/fm^3

  logical first_crit
  real(dp) :: pk,epsk,rhob,mb,r,emk
  real(dp) :: mbm,pkm,epskm,rhobm,rm,emkm,lambdam,xi
  real(dp) :: plimit,mu,pc,pcmax,rhoc,x,epsc,pkalt,lambda
  integer itab,ir



  mu = 931.494013d0 ! baryon mass (arbritrary unit)  in MeV
  first_crit = .true.



!     center of the star 
  r = 0.d0
  pk = pc ! central value for pressure
  emk = 0.d0
  mb = 0.d0
  itab = ntab
  lambda = 1._dp
  call eos(pk,itab)
  if (itab.gt.0) then
     x = (pk-ptab(itab)) / (ptab(itab+1)-ptab(itab))
     epsk = epstab(itab) + (epstab(itab+1)-epstab(itab)) * x
     rhob = rtab(itab) + (rtab(itab+1) - rtab(itab)) * x
     rhoc = rhob
  else
     epsk = epstab(1)*(pk/ptab(1))**0.6
  endif
  epsc = epsk
  open(7,file='tov.out',status='unknown')
  write(7,*) '#  r [km]  p [MeV/fm^3] eps [MeV/fm^3]  M_g[M_sol] n_B [1/fm^3]  M_B [M_sol]'
  first_crit = .true.
  ! center of the star

  write(7,*) r/1000.,pk,epsk,emk/msol,rhob,mu*mb/msol

  rloop: do ir = 1, nr
     pkalt = pk
           
 ! RK scheme: y(n+1) = y(n) + dr*f(r + dr/2,xn + 1/2*dr*f(r,xn))
     rm = r + dr/2

     IF(r.eq.0.) then
        lambda = 1._dp
        emkm = emk 
        mbm = mb
        pkm = pk
     else
        lambda = 1-2*gn*emk/r
        emkm = emk + dr*2*pi*r**2*epsk
        mbm = mb + 2.*pi*dr * r**2*rhob/sqrt(lambda)
        pkm = pk - dr/2*gn*(pk+ epsk)*(emk + 4*pi*r**3*pk)/(r*r*lambda)
     end if
     
     if (pkm.lt.ptab(itab)) call eos(pk,itab)
     if (itab.gt.0) then
        x = (pkm-ptab(itab)) / (ptab(itab+1)-ptab(itab))
        epskm = epstab(itab) + (epstab(itab+1)-epstab(itab)) * x
        rhobm = rtab(itab) + (rtab(itab+1) - rtab(itab)) * x
     else
        epskm = epstab(1)*(pkm/ptab(1))**0.6
     endif

     lambdam = 1-2*gn*emkm/rm
     
     emk = emk + dr*4*pi*rm**2*epskm
     mb = mb + 4.*pi*dr * rm**2*rhobm/sqrt(lambdam)
     pk = pk - dr*(pkm+ epskm)*gn*(emkm + 4*pi*rm**3*pkm)/(rm*rm*lambdam)

     if (pk.lt.ptab(itab)) call eos(pk,itab) 
     if (itab.gt.0) then
        x = (pk-ptab(itab)) / (ptab(itab+1)-ptab(itab))
        epsk = epstab(itab) + (epstab(itab+1)-epstab(itab)) * x
        rhob = rtab(itab) + (rtab(itab+1) - rtab(itab)) * x
     else
        epsk = epstab(1)*(pk/ptab(1))**0.6
     endif
     r = r + dr
     lambda = 1-2*gn*emk/r

     IF(first_crit.and.pk.lt.p_crit) then
        first_crit = .false.
        pk = p_crit
     end if
     write(7,*) r/1000.,pk,epsk,emk/msol,mu*mb/msol
     ! test: surface of the star
     if ((pk.lt.ptab(1)).or.(pk.gt.pkalt)) then
        xi = gn*emk/r ! compactness
        write(10,*) rhoc,epsc,r/1000,emk/msol,mu*mb/msol,xi,pc
        write(7,*) '# Radius of the star in km', r/1000.
        write(7,*) 
        close(7)
        exit rloop
     endif

  end do rloop
end subroutine integrate_tov


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine input

  USE constants
  implicit NONE


  character(72) text

  character(60) infile,outfile
  integer readeof,ip
  real(dp) pressure,energy,density

  write(*,*) 'From which file read EoS?'
  read(*,*) infile
  write(*,*) 'File to write the results?'
  read(*,*) outfile
  open(9,file=infile,status='old')
  open(10,file=outfile,status='unknown')
  write(10,*) '# n_B (central) [fm^{-3}]  eps (central) [MeV/fm^3] R [km]  M_g [M_sol]  M_B[M_sol] Xi  p (central) [MeV/fm^3]'      
  
  hbc3 = hbc**3
  ! some standard values for the parameters
  pcmin = 10._dp
  dpc = 1
  npc = 1500
  dr = 1
  nr = 2500000
  !
  
  
  ip = 0
  read(9,*) ! read the comment line at the beginning of the file
  readloop: DO
     read(9,*,iostat=readeof) density,energy,pressure

     IF(readeof.eq.0) then
        ip = ip + 1
        ptab(ip) = pressure
        rtab(ip) = density
        epstab(ip) = energy
     else
        ntab = ip
        write(*,*) ntab,'points read from EoS file'
        write(*,*) 'Minimum and maximum pressure of file',ptab(1),ptab(ntab)
        exit readloop
     end if
  end do readloop
  if(ntab.gt.npmax) then
     write(7,*) 'error: ntab > npmax'
     stop
  endif
  p_crit = 0._dp

  return
end subroutine input
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine eos(pk,itab)
  
  USE constants, only : dp,ptab,ntab
  implicit NONE

  real(dp) :: pk
  integer itab

  integer it

  IF(pk.lt.ptab(1)) then
     itab = 0
     return
  else 
     testloop: do it = itab, 1, -1
        if (ptab(it).lt.pk) exit testloop
     enddo testloop
     itab = it
  end IF
  if (itab.eq.ntab) then
     write(*,*) 'too large central pressure',itab,ntab,pk,ptab(ntab)
     stop
  endif

  return
end subroutine eos





























































