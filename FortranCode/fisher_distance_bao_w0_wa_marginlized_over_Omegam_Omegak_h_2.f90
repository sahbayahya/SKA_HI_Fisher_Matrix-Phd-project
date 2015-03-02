PROGRAM Fisher_Distance
  USE cosmo
  USE growth
  USE linearpk
  USE angular_distance
  ! Ref: Seo & Eisenstein, ApJ, 665, 14 (2007) .. Eq. (25)
  IMPLICIT none
  integer :: npara=2 ! # of parameters
  double precision, dimension(2)   :: der, work2
  double precision, dimension(2,2) :: cov,fis,fistot
  double precision, dimension(6,6) :: fis3x3, fis6x6
  integer, dimension(24)   :: list
  integer, dimension(2)   :: work
  double precision :: sigma0=12.4d0*0.817d0/0.9d0 ! non-linearity scale in units of h Mpc^-1, rescaled by WMAP5's sigma8 value (0.817)
  double precision :: BAO_AMP=0.5817d0 ! A_0 in Seo&Eisenstein
  double precision :: kmax_ov_h,k0_ov_h, kmin_ov_h, k0_new, k_new! h Mpc^-1
  double precision :: linear_pk,pk,p01,Pb, fistest
  integer:: status
  double precision :: k_full, fbao, pk_ref, dfbao, k_ov_h
  double precision :: wdamp_perp, wdamp_para, wdamp_silk, wdamp_factor, wdamp_sinterm
  double precision :: sigma_para,sigma_perp,sigma_silk=8.38d0 ! h^-1 Mpc
  double precision :: dlnk,mu,mu2,dmu,factor,wkmu,dummy,wdamp,Rmu
  double precision :: z,zin,beta,bias,g,dgdlna,sigma_z,delta_v,fz, w,dvdz, Dz
  double precision :: area,Vsurvey,Vtot,Ntot,Ngl, ngal,dndz,da, one_over_h, Vsurvey2, dvdz2
  character(len=200) :: filename,filein,fileout, fileFish
  character(len=24) :: sarea, filename2
  integer :: n,i,j,ibin,nbins,ifile, n1, n2, k, nline,L
  external linear_pk,dgdlna,g,dvdz,da, one_over_h
  
  if(iargc().eq.3) then
  	call getarg(1,filein)
  	call getarg(2,fileFish)
	call getarg(3,fileout)

  else
  	print *,'you have the wrong number of input arguements'  	
  	stop
  endif
 write(78,*)trim(filein)
 write(78,*)trim(fileout)
 write(78,*)trim(fileFish)
  
  !******************* PLANCK ********************************
  H_0 = 67.d0
  c   = 299792.458d0
  om0 = 0.267d0
  ode0= 0.686d0
  ok0 = 0d0
  ob0 = 0.049
  ok0 = 0.d0
  w0  =-1.d0
  w_a = 0.d0
  !==========================================
  CALL setup_growth ! tabulate the growth factor
  CALL setup_da       ! tabulate the angular diameter distance
  !==============Enter the File name ================================== 
  open(2,file=trim(filein),status='unknown')
  
  !=============  Read in linear P(k) ===================================
  filename= 'deriv_mikrom_dfbao_dk.dat'
  n  =129! no. of lines in the file
  zin=1d0 ! redshift of the input power spectrum
  CALL open_linearpk(filename,n) ! This fuction just to interpolate the values in the filename	
  
  !================================================================  
  ! loop over redshift bins
  fistot=0d0   ! total Fisher matrix integrated over all the redshift bins
  fis6x6=0d0 ! total Fisher matrix for w, Omega_k, and Omega_m
  Vtot  =0d0   ! total survey volume
  Ntot  =0d0
  nbins =0    ! initialize the number of bins
  !=========== Loop over redshift bins =======================================
  do ibin=1,100
      !read(2,*,end=10)z ,dndz,bias,kmax_ov_h,kmin_ov_h, Vsurvey , dvdz2! uncomment for BOSS
     !read(2,*,end=10)z ,dndz,bias,kmax_ov_h,kmin_ov_h,delta_v, Vsurvey , dvdz2! uncomment for Euclid ref
     read(2,*,end=10) z ,dndz,bias, kmax_ov_h, Vsurvey, dvdz2! Uncomment for SKA
     nbins=nbins+1
     write(79,*)z ,dndz,bias, kmax_ov_h, Vsurvey, dvdz2
     !--------------------------------------------------------------------
     beta  = (1d0+dgdlna(z)/g(z))/bias
     w       = w0 + w_a*(z/(1d0+z))
     dvdz2 = dvdz2*(H_0/100d0)**3 	! convert the comoving volume to Mpc^3 h^-3 Uncomment for SKA
     
    !ngal  = dndz	! Uncomment for  Euclid
    ngal = dndz/((3.14159d0/180d0)**2*dvdz2) 	! Uncomment for SKA
    
    
     Ngl   = (ngal * Vsurvey) 		!for Euclid
     Ntot  = Ntot + Ngl 		! Total number of galaxies 
     Vtot  = Vtot+Vsurvey 		! h^-3 Mpc^3.
     !==================Test the growth =================================
     !write(101,*) z, ngal!, D_plus!/(1d0+z), (1d0+dgdlna(z)/g(z))
     !write(101,*) z, Vsurvey, dvdz2, dndz, ngal, bias, beta, beta*bias, g(z)
     !=================  computing the Fisher matrix... =================================
     open(11,file=filename,status='old')
     read(11,*)k0_ov_h,dummy, dummy, dummy
     k_ov_h=k0_ov_h
     fis=0d0
     fistest= 0d0
     !================ loop over k ===================================
     do while (k_ov_h <=kmax_ov_h)
        read(11,*)k_ov_h, fbao, dfbao,pk_ref
        !dlnk=dlog(k_ov_h)-dlog(k0_ov_h) ! uncomment for widder bins
        !===========================================================
         k0_new=6.8492999999999996E-005 
         k_new =9.5589999999999998E-005
         dlnk=dlog(k_new)-dlog(k0_new)   !uncomment to use with Plank's parameters 
         factor=(k_ov_h)**3d0*dlnk/(8.d0 *3.1415926535d0**2d0) ! h^3 Mpc^-3
         !=========== loop over mu..========================
         mu=0d0
         dmu=1d-3
           do while (mu<=1d0)
             mu2= mu*mu
	     ! P(k) the galaxy power spectrum  in units of h^-3 Mpc^3  
             pk=pk_ref *(1d0 + fbao)*((g(z)/g(zin))*((1d0+zin)/(1d0+z)))**2d0 &
               *bias**2d0*(1d0+beta*mu2)**2d0      
             !==========================================     
             wkmu=Vsurvey*(ngal/(1d0+ngal*pk))**2d0
             wdamp = ((dfbao/(1d0 + fbao))*pk)**2d0
             der(1)=(mu2-1d0)*k_ov_h 	! dPkdlogDa
             der(2)=mu2*k_ov_h 			! dPkdlogH
             do i=1,npara
               fistest=factor*Vsurvey*(ngal*pk/(1d0+ngal*pk))**2d0
               fis(i,:)= fis(i,:) +factor* wkmu*wdamp*der(i)*der(:)*dmu 
           enddo
           mu=mu+dmu
        enddo
  
       !======================================================================
    k0_ov_h=k_ov_h
    enddo
    fistot=fistot+fis
    print*,'=== redshift bin#=(',ibin,') ==='
    print*, ' z = ', z
    print*,'Vsur =', Vsurvey,' h^-3 Mpc^3'
    print*,'ngal =',ngal,'h^3 Mpc^-3'
    print*,'1/ngal =',1d0/ngal,'h^-3 Mpc^3'
    print*,'Vsur * n(z)=', ngal* Vsurvey!,' h^-3 Mpc^3'
    print'(1A7,1F8.5)','Bias =',bias
    print'(1A7,1F8.5)','Beta =',beta
    print'(1A7,1F8.5)','g(z) =',g(z)
    print'(1A7,1F8.5)','gin a =',g(zin)/(1d0 +zin)
    print'(1A7,1F8.5)','f(z) =',1d0+dgdlna(z)/g(z)
    print'(1A7,1F8.5,1A9)','sigz =',sigma_z,' h^-1 Mpc'
    print'(1A7,1F8.5,1A9)','kmin =',kmin_ov_h,' h Mpc^-1'
    print'(1A7,1F8.5,1A9)','kmax =',kmax_ov_h,' h Mpc^-1'
    CALL report_result(z,bias,npara,fis,L,fileout, fileFish)
    print*,''
    close(11)
    CALL transform_fisher(z,fis,fis6x6)
 enddo
10 close(2)
CALL close_linearpk
 if(nbins>1)then
     write(*,*) ,'=== combined ==='
     write(*,*) ,'Vsur =',Vtot,' h^-3 Mpc^3'
     write(*,*) , 'N(z)=',Ntot
  endif
  print*,''
 CALL report_result3x3(fis6x6)
 print*,''
 
END PROGRAM Fisher_Distance

 !====================================================================
 !===================! SUBROUTINES  ===================================

 SUBROUTINE report_result(z,bias,npara,fis,L,fileout, fileFish)
  IMPLICIT none
  Integer, intent(IN) :: npara
  double precision, intent(IN) :: fis(npara,npara), z,bias
  Double precision, allocatable, dimension(:,:) :: cov
  Double precision, allocatable, dimension(:)   :: work
  Double precision :: r12,err_lnda,err_lnh,err_lnR,beta, linear_pk,dgdlna,g
  Integer :: i,j, L
  Character (len=24):: sarea
  Character (len=200):: fileout, fileFish
  external linear_pk,dgdlna,g
  ALLOCATE(cov(2,2),work(2))
  cov=fis
  CALL DVERT(cov,2,2,work)
  beta=(1d0+dgdlna(z)/g(z))/bias
  r12=cov(1,2)/sqrt(cov(1,1)*cov(2,2))
  err_lnda=sqrt(cov(1,1))
  err_lnh=sqrt(cov(2,2))
  err_lnR=err_lnda*sqrt((1d0-r12**2d0) &
       /(1d0+2d0*r12*err_lnda/err_lnh+(err_lnda/err_lnh)**2d0))
  write(sarea,'(I5)') L    
  !==============save desired fisher files =======================

 open(12,file=trim(fileFish), status='unknown')
 open(13,file =trim(fileout), status='unknown') 
  write(13,'(6F18.5)') z, err_lnda*1d2,err_lnh*1d2
  !========================================================
  print'(1A15,1F9.5)','Err[lnDa](%) =',err_lnda*1d2
  print'(1A15,1F9.5)','Err[lnH](%)  =',err_lnh*1d2
  print'(1A15,1F9.5)','r(lnDa,lnH)  =',r12
  print'(1A15,1F9.5)','Err[lnR](%)  =',err_lnR*1d2
  print*, 'Sarea = ', L
  DEALLOCATE(cov,work)
  return
END SUBROUTINE report_result

!===========================================================
SUBROUTINE report_result3x3(fis)
  USE cosmo
  IMPLICIT none
  double precision, intent(IN) :: fis(6,6)
  double precision :: work(6),cov(6,6), A(6,6), M55DET, DET5x5,DET2x2,A05(6,6)
  integer :: i, j 
  cov=fis
 !==============================Calculate DET 2x2 ===============
  CALL DVERT(cov,6,6,work)
  A = fis
  DET2x2 =  A(1,1)*A(2,2) - A(1,2)*A(2,1)  

  write(12,*)A(1,1) , A(1,2),  A(1,3) , A(1,4),  A(1,5), A(1,6)
  write(12,*)A(2,1) , A(2,2),  A(2,3), A(2,4),  A(2,5), A(2,6)
  write(12,*)A(3,1),  A(3,2),  A(3,3), A(3,4),  A(3,5),  A(3,6)
  write(12,*)A(4,1),  A(4,2),  A(4,3), A(4,4),  A(4,5),  A(4,6)
  write(12,*)A(5,1),  A(5,2),  A(5,3), A(5,4),  A(5,5),  A(5,6)
  write(12,*)A(6,1),  A(6,2),  A(6,3), A(6,4),  A(6,5),  A(6,6)     
      
   
 !====================write the resutls to file 12 ========================
  print*,'Figure of Merit  =' , (1d0/ dsqrt((cov(1,1) * cov(2,2) - cov(1,2)* cov(2,1))) )
  print'(1A20,1F9.5)', 'Err[w] =',dsqrt(cov(1,1))
  print'(1A20,1F9.5)', 'Err[wa] =',dsqrt(cov(2,2))
END SUBROUTINE report_result3x3  

 !================================================================
 SUBROUTINE transform_fisher(z,fisDH,fis6x6)
  USE cosmo
  IMPLICIT none
  integer :: a,b,i,j
  double precision, intent(IN)    :: z,fisDH(2,2)
  double precision, intent(INOUT) :: fis6x6(6,6)
  double precision :: dpdq(6,2)
  double precision :: chi,h2,func0,func1,func2,func3,func4,rombint,fz
  external h2,func0,func1,func2,func3,func4,rombint
  chi=rombint(func0,0d0,z,1d-7)
  fz =(1d0 + z)**(3d0*(1d0 + w0+ w_a)) * exp(-3d0 * w_a *(z/(1d0+ z)))
  dpdq(1,1)=-1.5d0*ode0*rombint(func1,0d0,z,1d-7)/chi        !dlnDa/dw
  dpdq(1,2)= 1.5d0*ode0*dlog(1d0+z)*fz/h2(z)                 !dlnH/dw
  dpdq(2,1)=-1.5d0*ode0*rombint(func4,0d0,z,1d-7)/chi        !dlnDa/dwa
  dpdq(2,2)= 1.5d0*ode0*(dlog(1d0+z) - z/(1d0+ z))*fz/h2(z)  !dlnH/dwa
  dpdq(3,1)=-0.5d0*ob0*rombint(func3,0d0,z,1d-7)/chi         !dlnDa/d(Omega_b)
  dpdq(3,2)= 0.5d0*ob0*((1d0+z)**3d0 - fz)/h2(z)             !dlnH/dOmega_b
  dpdq(4,1)=-0.5d0*rombint(func2,0d0,z,1d-7)/chi+chi**2d0/6d0!dlnDa/dOmega_k
  dpdq(4,2)= 0.5d0*((1d0+z)**2d0 - fz)/h2(z)                 !dlnH/dOmega_k
  dpdq(5,1)=-0.5d0*rombint(func3,0d0,z,1d-7)/chi             !dlnDa/d(Omega_m)
  dpdq(5,2)= 0.5d0*((1d0+z)**3d0 - fz)/h2(z)                 !dlnH/dOmega_m
  dpdq(6,1)= -100d0/H_0		                             !dlnDa/dH_0
  dpdq(6,2)= 100d0/H_0		                             !dlnH/dH_0

  do a=1,6
     do b=1,6
        do i=1,2
           do j=1,2
              ! transform and accumulate fis6x6
              fis6x6(a,b)=fis6x6(a,b)+dpdq(a,i)*dpdq(b,j)*fisDH(i,j) 
           enddo
        enddo
     enddo
  enddo
 return
 END SUBROUTINE transform_fisher
 !-------------------------------------------------------------------------------
 !==============! Functions ==================================
  DOUBLE PRECISION FUNCTION h2(redshift)
  USE cosmo
  !h2(z) = Omega_matter(1+z)^3+Omega_lambda
  IMPLICIT none
  DOUBLE PRECISION, intent(IN) :: redshift
  DOUBLE PRECISION :: fz
  fz =(1d0 + redshift)**(3d0*(1d0 + w0+ w_a)) * exp(-3d0 * w_a *(redshift/(1d0+ redshift)))
  h2 =  ( (om0+ob0)*(1d0+redshift)**3d0+ ok0 * (1d0+redshift)**2d0+ode0* fz)
  return
  END FUNCTION h2
  DOUBLE PRECISION FUNCTION func0(redshift)
  USE cosmo
  !func0(z) = 1/[h2(z)]^0.5
  IMPLICIT none
  DOUBLE PRECISION, intent(IN) :: redshift
  DOUBLE PRECISION :: h2
  external :: h2
  func0 = 1d0/dsqrt(h2(redshift))
  return
  END FUNCTION func0
  DOUBLE PRECISION FUNCTION func1(redshift)
  USE cosmo
  !func1(z) = ln(1+z)/[h2(z)]^1.5
  IMPLICIT none
  DOUBLE PRECISION, intent(IN) :: redshift
  DOUBLE PRECISION :: h2, fz
  external :: h2
  fz =(1d0 + redshift)**(3d0*(1d0 + w0+ w_a)) * exp(-3d0 * w_a *(redshift/(1d0+ redshift)))
  func1 = dlog(1d0+redshift) * fz /h2(redshift)**1.5d0
  return
  END FUNCTION func1
  DOUBLE PRECISION FUNCTION func2(redshift)
  USE cosmo
  !func2(z) = (1+z)^2/[h2(z)]^1.5
  IMPLICIT none
  DOUBLE PRECISION, intent(IN) :: redshift
  DOUBLE PRECISION :: h2,fz
  external :: h2
  fz =(1d0 + redshift)**(3d0*(1d0 + w0+ w_a)) * exp(-3d0 * w_a *(redshift/(1d0+ redshift)))
  func2 =( (1d0+redshift)**2d0 - fz)/h2(redshift)**1.5d0
  return
  END FUNCTION func2
  DOUBLE PRECISION FUNCTION func3(redshift)
  USE cosmo
  !func3(z) = (1+z)^3/[h2(z)]^1.5
  IMPLICIT none
  DOUBLE PRECISION, intent(IN) :: redshift
  DOUBLE PRECISION :: h2, fz
  external :: h2
  fz =(1d0 + redshift)**(3d0*(1d0 + w0+ w_a)) * exp(-3d0 * w_a *(redshift/(1d0+ redshift)))
  func3 =( (1d0+redshift)**3d0 - fz) /h2(redshift)**1.5d0
  return
  END FUNCTION func3
  DOUBLE PRECISION FUNCTION func4(redshift)
  USE cosmo
  !func4(z) = ln(1+z) - z/1+z /[h2(z)]^1.5
  IMPLICIT none
  DOUBLE PRECISION, intent(IN) :: redshift
  DOUBLE PRECISION :: h2, fz
  external :: h2
  fz =(1d0 + redshift)**(3d0*(1d0 + w0+ w_a)) * exp(-3d0 * w_a *(redshift/(1d0+ redshift)))
  func4 = fz * (dlog(1d0+redshift) - (redshift/(1d0+ redshift)))/h2(redshift)**1.5d0
  return
  END FUNCTION func4


 !(convert fisher matrix from 6x6 to 8x8 by adding zero columns and raws on the side)
  !  B(1,1) = 0.0d0 ; B(1,2)= 0.0d0; B(1,3)=  0.0 ;B(1,4)=  0.0  ;B(1,5)=  0.0 ; B(1,6)=  0.0  ;B(1,7)=  0.0 ;B(1,8)=0.0 
  !  B(2,1) = 0.0d0 ; B(2,2)=A(1,1); B(2,3)=A(1,2);B(2,4)= A(1,3);B(2,5)=A(1,4); B(2,6)= A(1,5);B(2,7)=A(1,6);B(2,8)=0.0 
  !  B(3,1) = 0.0d0 ; B(3,2)=A(2,1); B(3,3)=A(2,2);B(3,4)= A(2,3);B(3,5)=A(2,4); B(3,6)= A(2,5);B(3,7)=A(2,6);B(3,8)=0.0 
  !  B(4,1) = 0.0d0 ; B(4,2)=A(3,1); B(4,3)=A(3,2);B(4,4)= A(3,3);B(4,5)=A(3,4); B(4,6)= A(3,5);B(4,7)=A(3,6);B(4,8)=0.0 
  !  B(5,1) = 0.0d0 ; B(5,2)=A(4,1); B(5,3)=A(4,2);B(5,4)= A(4,3);B(5,5)=A(4,4); B(5,6)= A(4,5);B(5,7)=A(4,6);B(5,8)=0.0
  !  B(6,1) = 0.0d0 ; B(6,2)=A(5,1); B(6,3)=A(5,2);B(6,4)= A(5,3);B(6,5)=A(5,4); B(6,6)= A(5,5);B(6,7)=A(5,6);B(6,8)=0.0
  !  B(7,1) = 0.0d0 ; B(7,2)=A(6,1); B(7,3)=A(6,2);B(7,4)= A(6,3);B(7,5)=A(6,4); B(7,6)= A(6,5);B(7,7)=A(6,6);B(7,8)=0.0
  !  B(8,1) = 0.0d0 ; B(8,2)=  0.0 ; B(8,3)=  0.0 ;B(8,4)=   0.0 ;B(8,5)= 0.0  ; B(8,6)=  0.0  ;B(8,7)=  0.0 ;B(8,8)=0.0
  !==========================================================================  
  !! The planck fisher matrix for (ns, w0, wa, Ob, Ok, Om, h, sigma8)
 !P(1,1) = 1.99579245e+05 ; P(1,2) = -3.73667528e+04; P(1, 3) = -1.04936812e+04 ; P(1, 4) =4.34901751e+04 ;  P(1, 5) =5.58643962e+05; P(1,6) =  -5.80810461e+05; P(1,7) =  -7.65181989e+04; P(1, 8) = -2.23806234e+03
 !P(2,1) = -3.73667528e+04;P(2, 2) =1.83928663e+05 ; P(2, 3) =  5.16525685e+04 ; P(2, 4) = -1.13475150e+05; P(2, 5) = -3.98758357e+06; P(2, 6)=  5.46829285e+06; P(2, 7)=  1.32438370e+06  ;P(2,8) =  -4.51559188e+02
 !P(3,1) =-1.04936812e+04;P(3, 2) = 5.16525685e+04 ; P(3, 3) = 1.45055577e+04; P(3, 4) =  -3.18671540e+04; P(3, 5) = -1.11983054e+06; P(3, 6) =   1.53565718e+06; P(3, 7) =  3.71925825e+05; P(3,8)  -1.26811078e+02
 !P(4, 1)= 4.34901751e+04;P(4, 2) = -1.13475150e+05; P(4,3) = -3.18671540e+04; P(4,4) =   2.56561854e+05; P(4, 5) = 2.28693697e+06 ; P(4, 6) =  -3.18100554e+06; P(4, 7) =  -6.84459220e+05; P(4, 8) =  2.96407321e+01
 !P(5, 1) =5.58643962e+05; P(5,2) =  -3.98758357e+06; P(5,3) = -1.11983054e+06; P(5,4) =   2.28693697e+06; P(5, 5) = 8.70535526e+07; P(5, 6) = -1.19698819e+08; P(5, 7) =  -2.91740427e+07 ; P(5,8) =  1.88438127e+04
 !P(6, 1) =-5.80810461e+05; P(6,2)=  5.46829285e+06 ; P(6,3) =  1.53565718e+06; P(6, 4)  = -3.18100554e+06;P(6, 5) = -1.19698819e+08; P(6, 6) = 1.64891322e+08; P(6, 7) =   4.02520034e+07; P(6, 8) = -3.29986382e+04
 !P(7, 1) =-7.65181989e+04; P(7,2) = 1.32438370e+06; P(7, 3) = 3.71925825e+05; P(7, 4) =  -6.84459220e+05; P(7, 5) =-2.91740427e+07; P(7, 6) =   4.02520034e+07;P(7,7)=  9.88949015e+06; P(7,8) =  -1.01838183e+04
!P(8, 1) = -2.23806234e+03; P(8, 2) = -4.51559188e+02; P(8,3) = -1.26811078e+02; P(8,4) =  2.96407321e+01; P(8, 5) =1.88438127e+04 ; P(8, 6) = -3.29986382e+04; P(8,7) = -1.01838183e+04; P(8,8) =  1.51709659e+04
  !=========================================================================
 ! plus_planck = B + P
  !CT(1,1)=fis(1,1);CT(2,2)=fis(2,2);CT(1,2)=fis(1,2);CT(2,1)=fis(2,1)
  !CALL DVERT(C,2,2,work2)
  !print*, CT
  !cov = plus_planck
  !print*, plus_planck(2,2), plus_planck(3,3), plus_planck(3,2), (2,3)

