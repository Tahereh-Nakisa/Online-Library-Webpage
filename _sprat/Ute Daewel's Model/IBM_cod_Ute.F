c      subroutine spr_eggs(tm,dtd,ystat,istep1,ieggday)
c
c	include 'C_model.f'	
c	include 'C_cod.f'
c	include 'C_track.f'
c	real me, be
c	parameter(me=77.011,be=-1.2639)		!Thompson @al  ICES 1981
cc	common block includin the location vectors for the Particle
c      common /locuspt/ xpfort(jstat),ypfort(jstat),zpfort(jstat), 
c     *  laca(ndrei)!alok,itoa(m,n,ilo)
cc---------------------------------------------------------------------------							
cc	entwtime includes development stage vectors 
cc	istage:	live stage, vealue depends on life stage
c	common /entwtime/ tte(jstat),istage(jstat),xlj(jstat),wla(jstat)
c     *	,ttl(jstat)
c	common /help/ iidrei(m,n,ilo),iland(m,n,ilo), iizwei(m,n),
c     *    ibotlay(m,n)
c	integer ieggday(jstat)
c	dimension ylarv(m,n)
c	real dtd,tm(m,n,ilo)
c	real lie, sp,tstep
c	integer lstat
c	tstep=24./dtd		!dtd: time interval [h]
c	print*,'sprat_eggs'
cc	pause
c	do izal=1,jstat
c		if (xpfort(izal).lt.9999.)then
c			if (istage(izal) .eq. 1) then
c			
c				ii1=ypfort(izal)+0.5
c				jj1=xpfort(izal)+0.5
c				kk1=zpfort(izal)+0.5
cc				if (izal.eq.1)then
cc			print*,xpfort(izal),istage(izal),tte(izal)
cc	print*,tm(ii1,jj1,kk1)
cc	endif
c
cc------------------------------------------------------
cc	egg development depends on temperature [h]
cc	hourly 
cc	eggs are planktonic in the pelagial
cc-------------------------------------------------------
c
c				lie=(me*tm(ii1,jj1,kk1)**be)*tstep
c				tte(izal)=tte(izal)+1./lie
c				if (tte(izal) .ge.1.) then
c				ystat=ystat+1
c			    istage(izal)=2
c				ieggday(izal)=istep1
c				ttl(izal)=0.
c				ylarv(ii1,jj1)=ylarv(ii1,jj1)+1
c				endif
c	
c		endif
c	endif
c	enddo
c	return
c	end
c
c-------------------------------------------------------------------------------
c	duration of egg and yolk_sac stage 
c	development of egg and yolk sac phase according to jordaan and Kling 2003
        subroutine cod_yolksac(tm,dtd,fstat,istep1,jolkday,lstat)
c-------------------------------------------------------------------------------

        include 'C_model.f'	
        include 'C_cod.f'
        include 'C_track.f'

	real ml, bl
        parameter(ml=1079.6,bl=-0.16)		!jordaan and KLing 2003
!	common block includin the location vectors for the Particle
        common /locuspt/ xpfort(jstat),ypfort(jstat),zpfort(jstat), 
     *  laca(ndrei)!alok,itoa(m,n,ilo)

        common /memloc1/ xjfort(jstat),yjfort(jstat),zjfort(jstat)
!---------------------------------------------------------------------------							
!	entwtime includes development stage vectors 
!	istage:	live stage, vealue depends on life stage
        common /entwtime/ tte(jstat),istage(jstat),xlj(jstat),wla(jstat)
     *  ,ttl(jstat)
        common /help/ iidrei(ilo,m,n),iland(ilo,m,n), iizwei(m,n),
     *    ibotlay(m,n)

        dimension flarv(m,n)
	real dtd,tm(0:ndrei),ttec,ger
	real lie, sp,tstep,lie1
        integer fstat,jolkday(jstat),lstat,ll

        tstep=24./dtd   !dtd: time interval [h]
!        print*,'cod_yolksac'
        do izal=1,lstat
              if (xpfort(izal).lt.9999.)then
			if (istage(izal) .eq. 2) then
			
				ii1=ypfort(izal)+0.5
				jj1=xpfort(izal)+0.5
				kk1=zpfort(izal)+0.5
			
!------------------------------------------------------
!	egg development depends on temperature [h]
!	hourly 
!	eggs are planktonic in the pelagial
!-------------------------------------------------------

!				lie1=(ml*exp(bl*tm(ii1,jj1,kk1)))						![hours]

        ll=iidrei(kk1,ii1,jj1)
        if (tm(ll).lt.-1.0)then
           ttec=-1.0
        else 
           ttec=tm(ll)
        endif
        
        lie=((10.0**(1.871-0.79*log10(ttec+2.0)))+ !egg
     ,  1.0/(exp(0.22*tm(ll)-3.6)))*tstep  !ys

        ttl(izal)=ttl(izal)+1.0/lie
        if (ttl(izal) .ge.1.0) then
        fstat=fstat+1
        istage(izal)=3
        xlj(izal)=xlpy
        jolkday(izal)=istep1
        flarv(ii1,jj1)=flarv(ii1,jj1)+1
        xjfort(izal)=xpfort(izal)
        yjfort(izal)=ypfort(izal)
        zjfort(izal)=zpfort(izal)
        endif

        endif
        endif
        enddo
        return
        end


c----------------------------------------------------------------
c	feeding larvae forage, consumption, growth
c----------------------------------------------------------------
	subroutine cod_larvae(tm,dken,dtd,ext,iturb,istep1,
     & jolkday,lstat)
c----------------------------------------------------------------
		
	include 'C_model.f'
	include 'C_paths.f'
	include 'C_cod.f'
	include 'C_track.f'

c----------------------------------------------------------------------------------------
c	common block including the location vectors for the Particle
      common /locuspt/ xpfort(jstat),ypfort(jstat),zpfort(jstat), 
     *  laca(ndrei)!alok,itoa(m,n,ilo)
c----------------------------------------------------------------------------------------							
c	entwtime includes development stage vectors 
c	istage:	live stage, vealue depends on life stage
	common /entwtime/ tte(jstat),istage(jstat),xlj(jstat),wla(jstat)
     *	,ttl(jstat)
c----------------------------------------------------------------------------------------
c	inf. prey, pl: prey length, dpr: dens. of prey, wpr: prey weight, ipr: prey type
c	iipr : max. number of prey types
	common /prey/ plx(iipr),dpr(iipr),wpr(iipr)  !dpr(ipr,m,n,ilo) ECOSMO daten
	common /preybio/ fb(ndrei,nbio),pr(0:ndrei),fpr(iipr)

c----------------------------------------------------------------------------------------
	common /help/ iidrei(ilo,m,n),iland(ilo,m,n), iizwei(m,n),
     *    ibotlay(m,n)
c-----------------------------------------------------------------------------------------
	common /spout/tms(jstat),prcon(jstat),prcons(jstat),summr(jstat)
	common /gut/ p4hr(jstat,5)
     	
c-----------------------------------------------------------------------------------
	real pi,rd(iipr),ra(iipr),ss,alp,H,sp,C,Cmax,xpr(iipr),phi(iipr)
	real rank(iipr), plm(iipr),xpr1(iipr),wpr1(iipr),partpr(iipr)		!prey ranking variables
	real partpr1(iipr),wpr2(iipr)	          !IBM sens test 2.11.05		
	integer ivecii(iipr),jolkday(jstat)												!IBM sens test 2.11.05
	integer	irind(iipr),irind1(iipr),iturb
	real tm(0:ndrei), dtd,dken(0:ndrei)   
	real mr,ma,k,Growth,ae,om1,tstep,sda,ttest,mc
	real turbs,psp1,psp(iipr),put,epr(iipr),paused,pausef
	real eterm1,eterm2,eterm3,pbiomax,gph
c---------------------foraging submodel
	real CS(iipr),HT(iipr),CS1(iipr),HT1(iipr), CSmax				!capture success, handling time
	real CSa,CSx0,CSb,CSy0,CSc
	real CF !(365)
c	dimension ext(m,n,ilo+1)
	dimension ext(0:ndrei)		!IBM sens test 2.11.05
!	print*,'cod_larvae'
        
c------Constants--------------------------------------------------------
	pi = 4.*atan(1.) 
	H=dtd*60.*60.		!deltat  !Photoperiod =1, cause deltat=1h
	tstep=24./dtd
	
	CF=1000.*(1./(3.928**(1./0.199)))

c	print*,'dken',dken(100,60,1)
c------Foraging-------------------------------------------------
	do iizal=1,lstat

C---------------------------------------------------------------
 	if (xpfort(iizal).lt.9999.)then

	if(istage(iizal) .eq.3) then
		ii1=ypfort(iizal)+0.5
		jj1=xpfort(iizal)+0.5
		kk1=zpfort(iizal)+0.5
                ll=iidrei(kk1,ii1,jj1)
c-------q10 changed 10.11.06 for test rund frontal paper
		tms(iizal)=tms(iizal)+tm(ll)

	mr=0.00114*wla(iizal)**(1.029-0.00774*log(wla(iizal)))*
     ,	exp((0.1072-0.0032*log(wla(iizal)))*tm(ll))	!routine metabolism/tstep [mikrol O2/(fish h)]
	mc=0.00114*(wla(iizal)**(1.02900-0.00774*alog(wla(iizal))))
     &  *exp((0.10720-0.00320*alog(wla(iizal)))*tm(ll))


	mr=mr*0.00463*227.0									!siehe Hinckley 1996 [mikrog/h]
		
      if (ext(ll).ge.lightt) then                          !forage depends on light condition
	if (iizal.eq.1) then
		print*,'day'
	endif
c		ss=0.776*xlj(iizal)**1.07							!swimming speed of cruise predator
	ss=0.261*xlj(iizal)**(1.552*xlj(iizal)**0.92)-5.289/xlj(iizal)  !mean ss f�r cod nach Peck et al.
						
	alp=0.0167*exp((9.14)-(2.40)*log(xlj(iizal))
     *	 +(0.229)*(log(xlj(iizal)))**2.0)						!minimum angle of acuity
	alp=alp*pi/180.0									![rad]

c---------------------------------------------------------------------
c	including capture sucess and handling time cite Letcher (nur mal testen)
c--------------------------------------------------------------------

c	put=1.33*exp(-0.0035*xlj(iizal))	!pusuit time or handling time siehe HT(Fiksen 1998,Dower et al 1997)
		
c---loop over number of prey type-----------------------------------

	do ipr=1,iipr 
	HT(ipr)=exp(0.264*10.0**(7.0151*(plx(ipr)/xlj(iizal)))) !HT nach Walton (1992)
c	HT(ipr)=exp(0.174*10.**(0.1*(plx(ipr)/xlj(iizal)))) !HT nach Walton (1992) modifiziert

	dpr(ipr)=fpr(ipr)*pr(ll)				!density of particles at a certain length

	
	rd(ipr)=plx(ipr)/(2.0*tan(alp/2.0))				!reactive distance
	ssp=3.*plx(ipr)									!swimming speed of prey 3 bl/sec
c------------turbulence influence on encounter speed-------------
	if (iturb .eq. 0) then							!no turbulence
		ecs=(ssp**2.+3.*ss**2.)/(3.*ss)			!encounter speed
		psp(ipr)=1.
	
	elseif (iturb.eq.1) then
!        print *, dpr(ipr)
		if (dpr(ipr).eq.0.)then
			psp(ipr)=1.
			ecs=ss
		else
!        print *, dken(ii1,jj1,kk1)
		if (dken(ll).eq.0.)then
			psp(ipr)=1.
			ecs=sqrt(ssp**2.+ss**2.)/(3.*ss)		!encounter speed
		else
		ppsd=0.55*dpr(ipr)**(-0.33334)				!prey separation distance (MacKenzie,2000)
!		xte=3.615*(dken(ii1,jj1,kk1)*ppsd)**0.667	
		xte=1.877*(dken(ll)*rd(ipr))**0.667
	
!        print *, 'xte=',xtei,iipr
	if (xte.lt.0.)then
	print*,'xte:',xte,'.lt. 0 --> dken:',dken(ll),
     & 'ppsd',ppsd
	psp(ipr)=1.
	elseif (xte.ge.0.) then
		turbs=sqrt(xte)
		ecs=sqrt(ssp**2.+ss**2.+2.*turbs**2.) 
	
	
! 	turbs=sqrt(3.615*(dken(ii1,jj1,kk1)*rd(ipr))**0.667)  !uncorrelated velocity (Rothschild and Osborn 1988)	
		 
	ttest=rd(ipr)-turbs*HT(ipr)
!        print *, ttest,turbs,plx(ipr),rd(ipr),xlj(iizal)
	if (abs(ttest).lt.1.*10.**(-6.))then
		psp(ipr)=0.5  
	print*,'rd.eq.HT*turbs ipr',ipr,'rd',rd(ipr),'HT',HT(ipr),
     ,	'turbs',turbs, ii1,jj1,kk1
	else
	psp1=prosp(rd(ipr),turbs,HT(ipr))		!capture success due to turbulence
c	if (psp1.eq.9999.)then
cprint*,'rd.eq.HT*turbs ipr',ipr,'rd',rd(ipr),'HT',HT(ipr),
c    ,	'turbs',turbs
cc	pause
c	endif
c-----Fallunterscheidung siehe MacKenzie@al 1994
	ttt=turbs*HT(ipr)

	if (ttt.lt.rd(ipr))then
		psp(ipr)=1./rd(ipr)*(rd(ipr)-ttt+psp1)
	elseif (ttt.ge.rd(ipr).and.ttt.lt.(2.*rd(ipr)))then
	psp(ipr)=1./rd(ipr)*((rd(ipr)**3.)/(ttt**3.)*(ttt-rd(ipr))+psp1)
	elseif (ttt.gt.(2.*rd(ipr)))then
		psp(ipr)=rd(ipr)**3./(ttt**3.)
	endif
	endif
	endif		!xte
	endif		!dken=0
	endif		!dpr=0
	endif

	ra(ipr)=pi/2.*rd(ipr)**2.					!reactive area
c	epr(ipr)=ecs*ra(ipr)*dpr(ipr)				!encounter rate [no/s]  cruise predator
c	set pausefrequenzie and pauseduration
		if (xlj(iizal).lt.5.5) then
		pausef=pausefs
		paused=pauseds
	else
		pausef=pausefl
		paused=pausedl
	endif
	eterm1=2./3.*pi*rd(ipr)**3.*dpr(ipr)*pausef
	eterm3=sqrt(ssp+2.*turbs**2.)
	eterm2=pi*rd(ipr)**2.*dpr(ipr)*pausef*paused*eterm3
	epr(ipr)=eterm1+eterm2						!encounter rate pause travel predator
	xpr(ipr)=epr(ipr)							!ingestion rate#
	enddo		
	
	xprsum=0.
	do jzal=1,iipr
	xprsum=xprsum+xpr(jzal)
	enddo
	
c-------Consumption--------------------------------------------------	
			summe=0.

c	plmax=(0.1+1.1349/(1.+exp(-(xlj(iizal)-5.256)/0.8705))) !f�r cod nach loughs model
c	plmax=(0.15+1.1349/(1.+exp(-(xlj(iizal)-5.256)/0.8705))) !f�r cod nach loughs model
	plmax=(1.2936/(1.+exp(-(xlj(iizal)-2.8579)/1.2758)))     !f�r cod nach loughs daten 22.11.07
c	plmax=(0.15+1.2936/(1.+exp(-(xlj(iizal)-2.8579)/1.2758)))!f�r cod nach loughs daten 22.11.07
			CSmax=1.1
			yvar=(CSmax*xlj(iizal))/plmax	!Munk+plmax

c		yvar=13.6						!Munk

c	ranking of prey types (Letcher)
c----capture success based upon Mirjams data parameter:
	
				
			do jzal=1,iipr

	CS(jzal)=CSmax-(yvar*(plx(jzal)/xlj(iizal)))
c			CS(jzal)=0.8*CS(jzal)

        if (CS(jzal) .lt.0.) then
	CS(jzal)=0.
	endif
	wpr2(jzal)=2.*wpr(jzal)		!umrechnung von Carbon in biomass
	rank(jzal)=(wpr2(jzal)*CS(jzal))/HT(jzal)
	enddo
        do ii=1,iipr
        irind1(ii)=1
        irind(ii)=999.
        enddo

        do ii=1,iipr-1
             iy=0
             rankk=rank(ii)
                do jj=ii+1,iipr
                        if (rankk<rank(jj))then
                        iy=iy+1
                        else
                        irind1(jj) = irind1(jj)+1
                        endif
                enddo
             irind1(ii)=irind1(ii)+iy
        enddo

!	do ii=1,iipr
!		irind(ii)=999.
!	enddo

	do ii=1,iipr
	
		if (rank(ii).gt.0.) then
			irind(irind1(ii))=ii
		endif 
	enddo
	do jzal=1,iipr
		partpr(jzal)=0.
		ivecii(jzal)=0				!senstest 2.11.05
		partpr1(jzal)=0.				!senstest 2.11.05
	enddo


		summe=0.0
		summe1=0.0
		summe2=0.0

 		do jzal=1,iipr
		if (irind(jzal).ne. 999.) then

		summealt=summe
		summe2alt=summe2

	partpr(irind(jzal))=xpr(irind(jzal))*wpr2(irind(jzal))
     ,	*CS(irind(jzal))*psp(irind(jzal))
		summe1=summe1+partpr(irind(jzal))
		summe2=summe2+xpr(irind(jzal))*HT(irind(jzal))
		summe=summe1/(1.+summe2)
	

c	print*,'xlj',xlj(iizal),'pl',plx(irind(jzal)),'CS',CS(irind(jzal)),
c     , 	'ni',xpr(irind(jzal)),'HT',
c     ,	HT(irind(jzal)),'wpr',wpr2(irind(jzal)),summe1,	summe2,summe
c      pause
		if (summe.lt.summealt) then
			summe=summealt
			summe2=summe2alt
			goto 112
		endif
		ivecii(irind(jzal))=1
			endif
			enddo
112			C=summe*H

c------prey specific Ingestion nur als output f�r sensitivit�t	2.11.05
		do jzal=1,iipr
			if (ivecii(jzal).eq.1) then
		partpr1(jzal)=partpr(jzal)/(1.0+summe2)
			endif
		enddo
		
	pbiomax=3.24+0.064*wla(iizal)
	TT=tm(ll)
		CC=C
c----------------------------------------------------------------------------------------------------------------------
c  satiation
c----------------------------------------------------------------------------------------------------------------------
c	print*,C

               call satiation(iizal,TT,pbiomax,C,ger)


			
c				prcon(iizal)=C
		prcon(iizal)=prcon(iizal)+C
		Cmax=p4hr(iizal,1)
		prcons(iizal)=prcons(iizal)+Cmax
c-------Metabolism----------------------------------------------------

				if (xlj(iizal).lt.5.5)then
					k=1.4
				else
					k=2.5
				endif

				ma=k*mr						!active metabolism
c	print*,'xlj',xlj(iizal)
c	print*,'C',C,ma,C-ma
c	pause
                ae=0.7*(1.0-0.4*exp(-0.003*(wla(iizal)-50.0)))   !assimilation efficiency
!               ae=0.8*(1.0-0.4*exp(-0.002*(wla(iizal)-50.0)))
!               ae=0.8*(1.0-0.3*exp(-0.003*(wla(iizal)-50.0)))	!test 7.8.07

!               sda=0.11+0.000000491*wla(iizal)

                if (ae*C.gt.ma) then
                sda=0.35
                else
                sda=0.0
                endif

                cgr=C*ae*(1.-sda)

!	if C=0 --> cgr=0 --> no prey consumed
                Growth=cgr-ma               !growth [mikrog]	
!               elseif(isunzal.eq.0.or.ext(ii1,jj1,kk1).lt.lightt) then
                elseif (ext(ll).lt.lightt) then


                if (iizal.eq.1) then
!                print*,'night or deep'
                endif

                pbiomax=3.24+0.064*wla(iizal)
                TT=tm(ll)
                C=0.0
                call satiation(iizal,TT,pbiomax,C,ger)
                ae=0.7*(1.-0.4*exp(-0.003*(wla(iizal)-50.0)))   !assimilation efficiency

                if (ae*C.gt.ma) then
                sda=0.35
                else
                sda=0.0
                endif

                cgr=C*ae*(1.0-sda)

                Growth=cgr-mr
                endif     !day or night	

C----------growth threshold

        wla_old=wla(iizal)
        wla(iizal)=wla(iizal)+Growth
        if (wla(iizal).gt.2.*wla_old)then
        print*,'unrealistic growth ',iizal,wla(iizal),xlj(iizal),pbiomax
     &,C,TT
        endif

        if (wla(iizal).lt.0.)then
        print*,iizal,wla(iizal),xlj(iizal),C,ae,mr,cgr
        endif
!-----------------------------------------------------------------------
!	larval death
!--------------------------------------------------------------------------------
!		print*,istep1,jolkday(iizal)
            idaynow = istep1-jolkday(iizal)
            wcodmin = wcodmin0*exp(0.0282*(real(idaynow))) ! cod	vorausgesetzt wir starten die feeding phase mit hatch      
            deadcod = wla(iizal)-wcodmin
            if(deadcod.lt.0.0.or.wla(iizal).lt.0.0)then         !larvae dies
            istage(iizal)=0
            wla(iizal)=0.0
            end if
!
!   growth (percent per hour)
!
!            gph=100.*cgr/wla(iizal)   !growth per hour


               xlja=xlj(iizal)
!-------length growth
               xlj(iizal)=1.935*wla(iizal)**(0.247) !lough etal
!               xlj(iizal)=((wla(iizal)/1000.)/exp(-9.7))**(1./3.997)  !nach Otterlei etal.1999
               if (xlj(iizal).lt.xlja)then
               xlj(iizal)=xlja
               endif 
               endif!if istage
!------------1d IBM output---------------------------------------------------------------
!	write(10,'(6(1x,f12.6))')xlj(iizal),wla(iizal),C,CC,mr,ger !,xprsum		!IBM sens test 2.11.05
!	write(11,'(150(1x,f12.6))') partpr1
!----------------------------------------------------------------------------------------

               endif   !if xpfort
               enddo   !jstat

               return 
               end

       function prosp(R,w,t)

       real R,t,w,a,ttt,t1,t2,loga
       a=R
       if (a.le.0.)then
       print*,'a=R=rd=0'
       endif
       wt=w*t
       wt2=wt*wt
       wt3=wt2*wt
       wt4=wt3*wt
        R2=R*R
        R3=R2*R
        R4=R3*R
       loga=log(a)
        a2=a*a
        a4=a2*a2
       t1=-1./64.*(-a4+12.*R2*a2+12.*a2*wt2+12.*loga*R4
     - -24.*loga*R2*wt2+12.*wt4*loga-32.*R3*a
     - -32.*a*wt3)/wt3

        ttt=w*t
        if (ttt.gt.R)then
        a=ttt-R
        elseif(ttt.le.R)then
        a=R-ttt
        endif
        if (a.le.0.)then
        print*,'a=',a, 'ttt=',ttt,'R=',R,'w',w,'t',t
        prosp=9999.
        goto 171
        endif
      loga=log(a)
        a2=a*a
        a4=a2*a2
      t2=-1./64.*(-a4+12.*R2*a2+12.*a2*wt2+12.*loga*R4
     -  -24.*loga*R2*wt2+12.*wt4*loga-32.*R3*a
     -  -32.*a*wt3)/wt3
        prosp=t1-t2

171	return 
	end


c-----------------------------------------------------------------------
c-------conv_biom-------------------------------------------------------
c	converts biomass from ecosmo into number of zooplankton particles
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
	subroutine conv_biom
c----------------------------------------------------------------------
      include 'C_model.f'
      include 'C_index.f'
      include 'C_paths.f'
	include 'C_cod.f'
	include 'C_track.f'
#ifdef MPI
        include 'C_mpi.f'
        include 'C_parallel.f'
#endif
	common /preybio/ fb(ndrei,nbio),pr(0:ndrei),fpr(iipr)
	common /prey/ plx(iipr),dpr(iipr),wpr(iipr)  !dpr(ipr,m,n,ilo) ECOSMO daten
	real fun(iipr),pl1(iipr),conc(iipr),totbiom
        common /help/ iidrei(ilo,m,n),iland(ilo,m,n), iizwei(m,n),
     *    ibotlay(m,n)

	totbiom=0.0
!        call Rzero(pr,m*n*ilo)
        pr(0) = 0.0
        do ll=1,ndrei
        pr(ll)=(fb(ll,1)+fb(ll,2))*10.**(-6.)![microgC/mm�]!sum of zooplanktonbiomass 3d rechnung
!				pr(i,j,k)=pr(i,j,k)*10.**(-6.)						!1d rechnung
        enddo
!----initialize
        do ii=1,iipr
        fun(ii)=0.0
        pl1(ii)=0.0
        conc(ii)=0.0
        enddo
!-----------preylength increments---------------
        plx(1)=plstart*10.0**(-3.0)
        pl1(1)=plstart
        do ii=2,iipr
        pl1(ii)=pl1(ii-1)+real(pincr)
        plx(ii)=plx(ii-1)+real(pincr)*10.0**(-3.0)   !length[mm]
        enddo

!	print*,plx(1:10)
!	pause
        do ii=1,iipr
!		fun(ii)=pfun1*pl1(ii)**pfun2		 !partially length distribution of biom
!		wpr(ii)=16.1216*plx(ii)**(2.6011) 
        wpr(ii)=10.0**(2.7722*log10(pl1(ii))-7.47607)  ![microg C]carbon content literature review (Myron) for Acar/tem/cal/pseudo
!----reference case ---------------------------------------------
        conc(ii)=695.73*exp(-0.0083*pl1(ii)) !no/liter macht aber nischt, da das normiert wird standart case
!-------------------------------------------------------------------

!	conc(ii)=695.73*exp(-0.00575*pl1(ii))			!testII ib4
!conc(ii)=695.73*exp(-0.0112*pl1(ii))				!testIII ib5
!	conc(ii)=5.1771*10.**(-4.)*exp(0.0083*pl1(ii))	!testIV	
!	conc(ii)=200.									!testV
!	conc(ii)=-9.2862+312.3405*exp(-0.5*((pl1(ii)-850.)/282.8503)**2.) !test VI
!	conc(ii)=0.0007+300.0009*exp(-0.5*((pl1(ii)-329.164)/254.479)**2.) !test VII

!	conc(ii)=0.1061*exp(-0.006*pl1(ii))			!tm1
!	conc(ii)=0.0403*exp(-0.003*pl1(ii))			!tm2	
!	conc(ii)=0.3754*exp(-0.012*pl1(ii))			!tm3
!	conc(ii)=1.361*10.**(-7.)*exp(0.0083*pl1(ii))			!tm4
!	conc(ii)=4.514*10.**(-4.)+0.0152*exp(-0.5*((pl1(ii)-850.)
!     ,	/282.8503)**2.) !tm5
!	conc(ii)=-4.454*10.**(-8.)+0.191*exp(-0.5*((pl 1(ii)-329.164)
!     , 	/254.479)**2.) !tm6

!		fun(ii)=(conc(ii)*wpr(ii))/692.375
        fun(ii)=(conc(ii)*wpr(ii))
        totbiom=totbiom+fun(ii)
        enddo
!alok        do ii=1,iipr
!alok        totbiom=totbiom+(conc(ii)*wpr(ii))
!alok        enddo
!	print*,'totbiom',totbiom
!	pause
        do ii=1,iipr
        fun(ii)=fun(ii)/totbiom
        fpr(ii)=fun(ii)/wpr(ii)       !factor for biom --> number of particles per length class
!		fpr(ii)=fun(ii)
        enddo
        return
        end
c-------------------------------------------------------------------------
c       calculating radiation extinction
c---------------------------------------------------------------
        subroutine lightext(einstr,flag,diat,ext,df)

        include 'C_model.f'
        dimension flag(ndrei),diat(ndrei)!alok,EXtot(m,n,ilo)
        dimension BioC(2),EXwater(ilo), EXphyto(m,n),df(ilo)
        dimension ext(0:ndrei),einstr(m,n)
        common /help/ iidrei(ilo,m,n),iland(ilo,m,n), iizwei(m,n),
     *    ibotlay(m,n)
        real EXtot1,flagg,diatt
!        do i=1,ilo
!        df(i)=0.0 
!        enddo
        BioC(1)=0.05                    !light extinction [1/m]
c       BioC(2)=2.25                    !phyto self shading [m**2/(mgC)]
        BioC(2)=0.03/(6.625*12.01)                      !phyto self shading [m**2/(mgC)]
        ext(0) = 0.0
!        call Rzero(ext,m*n*ilo)
        do j=1,n
           do i=1,m
                EXphyto(i,j)=0.0
                ll = iidrei(1,i,j)        
               if (ll .gt. 0 ) ext(ll)=einstr(i,j)*5.578
           enddo
        enddo

        EXwater(1)=0.0
        do k=2,ilo
           EXwater(k)=EXwater(k-1)+BioC(1)*df(k-1)
        enddo

!        do k=2,ilo
           do j=1,n
              do i=1,m
!               EXwater(i,j,1)=EXwater(i,j,1)+df(1)
!               EXphyto(i,j,1)=EXphyto(i,j,1)+(flag(i,j,1)+diat(i,j,1))*df(1)
!                EXtot(i,j,1)=0.0

                        do k=2,ibotlay(i,j)

!*________ light limitation_______________________________
!                        EXwater(k)=EXwater(k-1)+df(k-1)
         ll=iidrei(k,i,j)
         if (ll.gt. 0 )then
!         lld=iidrei(i,j,k-1)
         lld=ll-1
!         if(lld .gt. 0) then
         flagg = flag(lld)
         diatt = diat(lld)
!         else
!         flagg = 0.0
!         diatt = 0.0
!         endif
         EXphyto(i,j)=EXphyto(i,j)+
     *               (flagg+diatt)*df(k-1) 
!         EXtot(i,j,k)= BioC(2)*EXphyto(i,j)+EXwater(k)
         EXtot1= BioC(2)*EXphyto(i,j)+EXwater(k)
         ext(ll)=einstr(i,j)* exp(-EXtot1) *5.578
         endif
              enddo
           enddo
        enddo
!          do j=1,n
!                do i=1,m
!           if (iland(i,j,k).gt.0 )then
!          EXtot(i,j,k)= BioC(2)*EXphyto(i,j)+EXwater(k)
!          EXtot1= BioC(2)*EXphyto(i,j)+EXwater(k)
!          ext(i,j,k)=einstr(i,j)* exp(-EXtot1) *5.578
!          endif
!          enddo
!                enddo

        
!       do k=1,ilo
!       do j=1,n
!       do i=1,m
!            if (iland(i,j,k).gt.0 )
!     *          ext(i,j,k)=einstr(i,j)* exp(-EXtot(i,j,k)) *5.578   !umrechnung in mikroE/(m**2*s) lichtintensit�t !einstr aus feld
!        enddo
!        enddo
!        enddo

        return
        end


cc----------------------------------------------------------------------
cc	calculate remaining food in the gut  aus Lough et al. 2006
cc    code by Karen P. Edwards cod93.f
cc----------------------------------------------------------------------
c      subroutine satiation(l,tdegC,pbiomax,pbiomass)
cc
cc  Determine how much prey will be ingested based on consumption
cc    over previous 4 hours.  If the gut content exceeds a certain
cc    value then cap the current ingestion.  If it is below the
cc    observed max gut content, then ingest it all.
c	include 'C_model.f'
c	include 'C_paths.f'
c	include 'C_cod.f'
c	include 'C_track.f'
c
c
cc
c      real pbiomax,pbiomass,p0,p1,p2,p3,p4,sum,tmpefct,tdegc
c	common /gut/ p4hr(jstat,4)
cc
cc  gut biomass is X% of one hour ago, Y% of two hours ago, and
cc     Z% of 3 hours ago, etc... (based on Tilseth and Ellersten) 
cc
c      p0=pbiomass                     ! prey from current hour
c      p1=p4hr(l,1)		      ! prey in gut from 1 hour ago
c      p2=p4hr(l,2)		      ! prey in gut from 2 hours ago
c      p3=p4hr(l,3)		      ! prey in gut from 3 hours ago
c      p4=p4hr(l,4)		      ! prey in gut from 4 hours ago
c
cc
cc Sum the gut biomass and current meal and run check against pbiomax 
cc Reset current meal to lower value if exceed cap
cc
c      sum=p0+p1+p2+p3+p4
c      if(sum.gt.pbiomax)then
c         p0=amax1(pbiomax-(p1+p2+p3+p4),0.0)
c      end if 	 
cc
cc Return the appropriate portions of the meals from last several hours
cc  In waters warmer than 5 degC more food is digested (since tmpefct>1)
cc
c      tmpefct=1.-0.10*(5.-tdegC)          ! KPE Larry 3/26/03
c      pbiomass=(0.3*p0+0.4*p1+0.6*p2+0.8*p3)*tmpefct  !KPE, Larry 2/19/03    
cc
cc Save the gut contents for next time step.
cc  For waters warmer than 5 degC, less will remain in gut
cc    due to increased digestion (see above)
cc
c      p4hr(l,1)=amax1(p0-0.3*p0*tmpefct,0.0)		!KPE, Larry 2/19/03 
c      p4hr(l,2)=amax1(p1-0.4*p1*tmpefct,0.0)		!KPE, Larry 2/19/03 
c      p4hr(l,3)=amax1(p2-0.6*p2*tmpefct,0.0)
c      p4hr(l,4)=amax1(p3-0.8*p3*tmpefct,0.0) 
cc
c	
cc     write(90,100)p0,p1,p2,p3,pbiomass
cc 100 format(5f9.4)    
cc 100 format('p0-3, pbiomass = ',5f9.4)
c     return
c      end
cc


c----------------------------------------------------------------------
c	calculate remaining food in the gut  GER aus Peck and Daewel
c    code by Karen P. Edwards cod93.f and UTe Daewel
c----------------------------------------------------------------------
      subroutine satiation(l,tdegC,pbiomax,pbiomass,ger)
c
c  Determine how much prey will be ingested based on consumption
c    over previous 4 hours.  If the gut content exceeds a certain
c    value then cap the current ingestion.  If it is below the
c    observed max gut content, then ingest it all.
	include 'C_model.f'
	include 'C_paths.f'
	include 'C_cod.f'
	include 'C_track.f'


c
      real pbiomax,pbiomass,p0,p1,p2,p3,p4,sum,tmpefct,tdegc
	common /gut/ p4hr(jstat,5)
c----------------------------------------------------------------------------------------							
c	entwtime includes development stage vectors 
c	istage:	live stage, vealue depends on life stage
	common /entwtime/ tte(jstat),istage(jstat),xlj(jstat),wla(jstat)
     *	,ttl(jstat)
	real ger,ger0(5),ger1(5),pbin

	pbin=pbiomass
c----------------------------------------------------------------------------
c	calculating gut evacuation rate
c----------------------------------------------------------------------------
	ger=(1.79*xlj(l)**(-0.83)*3.0**((tdegC-15.0)/10.0)) !see Peck and Daewel 2007 (q10=3)
!     ,	*(xlj(l)**0.2*0.8)								!ger tuning
!     ,	*(xlj(l)**0.43*0.55)
        if (pbiomass.gt.0.0)then
        ger=ger*3.0
        endif
!	print*,'ger',ger
!	ger=1.2*xlj(l)**(-0.6)*2.**((tdegC-12.)/10.)
        ger0(1)=1.0
        ger1(1)=ger
        do i=2,5
        ger0(i)=ger0(i-1)-ger
        if (ger0(i).gt.0.0)then
                ger1(i)=amin1(ger/ger0(i),1.0)
        elseif (ger0(i).le.0.)then
        ger1(i)=1.0
        endif
        enddo
!	print*,'ger1',ger1
!	print*,'pbiomass',pbiomass
!c------------test
!	ger1(1)=0.3
!	ger1(2)=0.4
!	ger1(3)=0.6
!	ger1(4)=0.8
!	ger1(5)=1.
!
!  gut biomass is X% of one hour ago, Y% of two hours ago, and
!     Z% of 3 hours ago, etc... (based on Tilseth and Ellersten) 
!
      p0=pbiomass       ! prey from current hour
      p1=p4hr(l,1)      ! prey in gut from 1 hour ago
      p2=p4hr(l,2)      ! prey in gut from 2 hours ago
      p3=p4hr(l,3)      ! prey in gut from 3 hours ago
      p4=p4hr(l,4)      ! prey in gut from 4 hours ago
      p5=p4hr(l,5)


!
! Sum the gut biomass and current meal and run check against pbiomax 
! Reset current meal to lower value if exceed cap
!
      sum=p0+p1+p2+p3+p4+p5
      if(sum.gt.pbiomax)then
         p0=amax1(pbiomax-(p1+p2+p3+p4+p5),0.0)
      end if 
!
! Return the appropriate portions of the meals from last several hours
!
      pbiomass=(ger1(1)*p0+ger1(2)*p1+ger1(3)*p2+ger1(4)*p3+ger1(5)*p4)
!	print*,p0,p1,p2,p3,p4
!	print*,pbiomax,pbiomass
!	pause

      if (pbiomass.gt.pbin) then
      pbiomass=pbin

      p4hr(l,1)=0.0         !KPE, Larry 2/19/03 
      p4hr(l,2)=0.0         !KPE, Larry 2/19/03 
      p4hr(l,3)=0.0
      p4hr(l,4)=0.0
      p4hr(l,5)=0.0
      else
!
! Save the gut contents for next time step.
!
      p4hr(l,1)=amax1(p0-ger1(1)*p0,0.0)                !KPE, Larry 2/19/03 
      p4hr(l,2)=amax1(p1-ger1(2)*p1,0.0)                !KPE, Larry 2/19/03 
      p4hr(l,3)=amax1(p2-ger1(3)*p2,0.0)
      p4hr(l,4)=amax1(p3-ger1(4)*p3,0.0) 
      p4hr(l,5)=amax1(p4-ger1(5)*p4,0.0) 
!
      endif
!     write(90,100)p0,p1,p2,p3,pbiomass
! 100 format(5f9.4)    
! 100 format('p0-3, pbiomass = ',5f9.4)
      return
      end
c

