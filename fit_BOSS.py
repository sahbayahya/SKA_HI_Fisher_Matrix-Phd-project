import numpy as np
import scipy.optimize as opt
from scipy import interpolate, concatenate, reshape, sqrt, savetxt, linspace, exp, sin, log
import matplotlib.pyplot as plt
import cosmolopy.constants as cc
import cosmolopy.distance as cd
import cosmolopy.perturbation as cp
from scipy.integrate import quad, cumtrapz, simps
import sys
sys.path.append("/home/nassp/sahba/Desktop/Fisher_Euclid_numeric_derivative/cosmolopy")
#=============Functions ====================
def background_evolution_splines(zmax=10., nsamples=500):
    """
    Get interpolation functions for background functions of redshift:
      * H(z), Hubble rate in km/s/Mpc
      * r(z), comoving distance in Mpc
      * D(z), linear growth factor
      * f(z), linear growth rate
    """
    cosmo = {'omega_M_0':        0.316,
			 'omega_lambda_0':   0.684,
    			'omega_b_0':        0.049,
    			'N_eff':            3.046,
   			 'h':                0.67,
   			 'ns':               0.962,
   			 'sigma_8':          0.834,
    			'gamma':            0.55,
   			 'w0':               -1.,
    			'wa':               0.,
   			 'sigma_nl':         7.}
    _z = linspace(0., zmax, nsamples)
    a = 1. / (1. + _z)
    H0 = (100.*cosmo['h']); w0 = cosmo['w0']; wa = cosmo['wa']
    om = cosmo['omega_M_0']; ol = cosmo['omega_lambda_0']
    ok = 1. - om - ol
    C= 5e3
    # Sample Hubble rate H(z) and comoving dist. r(z) at discrete points
    omegaDE = ol * exp(3.*wa*(a - 1.)) / a**(3.*(1. + w0 + wa))
    E =sqrt( om * a**(-3.) + ok * a**(-2.) + omegaDE )
    _H = H0 * E
    
    r_c = concatenate( ([0.], cumtrapz(1./E, _z)) )
    if ok > 0.:
        _r = C/(H0*sqrt(ok)) * sinh(r_c * sqrt(ok))
    elif ok < 0.:
        _r = C/(H0*sqrt(-ok)) * sin(r_c * sqrt(-ok))
    else:
        _r = (C/H0) * r_c
    
    # Integrate linear growth rate to find linear growth factor, D(z)
    # N.B. D(z=0) = 1.
    a = 1. / (1. + _z)
    Oma = cosmo['omega_M_0'] * (1.+_z)**3. * (100.*cosmo['h']/_H)**2.
    _f = Oma**cosmo['gamma']
  #  print _f
    _D = concatenate( ([0.,], cumtrapz(_f, log(a))) )
    _D = exp(_D)
    
    # Construct interpolating functions and return
    r = interpolate.interp1d(_z, _r, kind='linear', bounds_error=False)
    H =interpolate.interp1d(_z, _H, kind='linear', bounds_error=False)
    D = interpolate.interp1d(_z, _D, kind='linear', bounds_error=False)
    f = interpolate.interp1d(_z, _f, kind='linear', bounds_error=False)
    return  _z, H, r, D, f


def dvdz(z):
	''' this function is to calculate the diff comoving volume 
	the results are given in units of Mpc^3.
	to use this function you need to install cosmolopy.
	Also note that the cosmological 
	parameters are  Planck best-fit parameters.
	'''

	cosmo = {'omega_M_0':        0.316,
			 'omega_lambda_0':   0.684,
    			'omega_b_0':        0.049,
    			'N_eff':            3.046,
   			 'h':                0.67,
   			 'ns':               0.962,
   			 'sigma_8':          0.834,
    			'gamma':            0.55,
   			 'w0':               -1.,
    			'wa':               0.,
   			 'sigma_nl':         7.}
	cosmo = cd.set_omega_k_0(cosmo)
	Vc = cd.diff_comoving_volume(z, **cosmo)
	return  Vc
def da(z):
	'''This function is to calculate the angular diameter distance
        The units are in Mpc
        The cosmological parameters are Planck best-fit parameters
        Note: you need to install cosmolopy and import cosmolopy.constants as cc
        and  import cosmolopy.distance as cd
	'''
	# 
	cosmo = {'omega_M_0':        0.316,
			 'omega_lambda_0':   0.684,
    			'omega_b_0':        0.049,
    			'N_eff':            3.046,
   			 'h':                0.67,
   			 'ns':               0.962,
   			 'sigma_8':          0.834,
    			'gamma':            0.55,
   			 'w0':               -1.,
    			'wa':               0.,
   			 'sigma_nl':         7.}
	cosmo = cd.set_omega_k_0(cosmo)
	d_a = cd.angular_diameter_distance(z, **cosmo)*(cosmo['h'])
	#print "Angular diameter distance = %.1f Mpc)" % (d_a) 
	return  d_a
def V_sur(zmin, zmax,area):
	''' This function to calculate the survey volume 
	The units will be Mpc^3 per deg^2 
	to convert to Mpc^3/h^3 you need to multiply by h^3
	'''
	vol = quad(dvdz, zmin, zmax)[0]	
	vol = area*(3.1415/180.)**2.*vol
	return vol
def func_rms(p,x):
	'''This is an exponential function with 3 parameters,
	y = a x^b  exp(-c x)
	'''
   	w=10.**p[0]*  x**p[1]  * np.exp(-p[2]*x) 
   	print w.size
   	return w

def residuals(p,x,y):
	'''This calculate the  sqrt of difference between the theortical 
	function y = a x^b  exp(-c x) and the measured one (your data at spesific x)
	'''
   	w=10.**p[0]* x**(p[1] ) *  np.exp(-p[2]*x)
   	err=w-y
   	err=err**2
   	B=sum(err)
   	return B	
def Bias(z):
	'''This is the bias function which is used by Euclid survey
	b(z) = sqrt(1+z)
	'''
	return  sqrt(1+z)
	
def find_parameters_rms(p0, x, rms , rms_value, xrange):
	plsqtot = opt.fmin(residuals, p0, args=(x,rms), maxiter=10000, maxfun=10000)
	print  'rms = %d uJy' % rms_value, 'p[0] = ' ,plsqtot[0],  'p[1] = ' , plsqtot[1]
	print '---------------------------------------------------------'
	ytot=func_rms(plsqtot,xrange)
	return ytot
	
def func_bias(p,x):
   w=p[0]*np.exp(p[1]*x)
   print w.size
   return w

def residuals_bias(p,x,y):
   w=p[0]* np.exp(p[1]*x)
   err=w-y
   err=err**2
   B=sum(err)
   return B
def growthf(z, omega_M_0):
	'''f growth calculate the 
	'''
	return cp.fgrowth(z, omega_M_0, unnormed=True)
	
'''This function call opt.fmin function to  minimize the error on
you paramters and give you the best fit parameters
'''      
def find_parameters_bias(p0, x, rms , rms_value, xrange):
	plsqtot = opt.fmin(residuals_bias, p0, args=(x,rms), maxiter=10000, maxfun=10000)
	print  'rms = %d uJy' % rms_value ,': ', 'p[0] = ' ,plsqtot[0],  'p[1] = ' , plsqtot[1]
	print '---------------------------------------------------------'
	ytot=func_bias(plsqtot,xrange)
	return ytot

def k_max(z,k_NL):
    cosmo = {'omega_M_0':        0.316,
			 'omega_lambda_0':   0.684,
    			'omega_b_0':        0.049,
    			'N_eff':            3.046,
   			 'h':                0.67,
   			 'ns':               0.962,
   			 'sigma_8':          0.834,
    			'gamma':            0.55,
   			 'w0':               -1.,
    			'wa':               0.,
   			 'sigma_nl':         7.}
    kmax = k_NL*(1.+z)**(2./(2.+cosmo['ns']))
    return kmax
	      		
#==============redshift (z) range==================================

if __name__=='__main__':
	print "_"*30

	xrange = np.array([ 0.05, 0.15 , 0.25 ,0.35 , 0.45 , 0.55, 0.65 , 0.75])
	dndz   = np.array([ 8, 50, 125, 222, 332, 447, 208, 30])
	xmin = xrange-0.05
	xmax = xrange+0.05
	print 'zmin = ' , xmin
	print 'zmax = ', xmax 
	print 'zc = ',  (xmin + xmax)/2.0

	#============k [Mpc^-1 h] range===================================
	k_NL = 0.2
	kmax = k_max(xrange,k_NL)# Output survey info
	_kmax = np.linspace(0.2, 0.2, 8)
	kmin = np.linspace(0.00334,0.00334, 8)
	print 'kmax = ' , kmax
	print '_kmax=', _kmax
	'''
	caculate the gorwth fucntion 
	'''
	x = xrange
	(x, H, r,D,f) = background_evolution_splines(zmax=2.1,nsamples=500)

	D_zin  =  np.empty(len(x))

	#========V survey [Mpc^3 h^-3]=================================
	h= 0.67
	area = 10000.0
	Vsurvey= []
	
	for i in range(len(xrange)):
          	Vsurvey.append(V_sur(xmin[i],xmax[i], area)*(h**3))
	print 'Vsurvey = ', Vsurvey 
	#======== The error on z==============================
	ErrorZ= np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
	bias =  np.array([ 2.0, 2.0 , 2.0 ,2.0 , 2.0 , 2.0, 2.0 ,2.0])
	#======= Save the resutls=============================

	print len(bias), len(kmax), len(kmin), len(Vsurvey), len(dvdz(xrange)), len(ErrorZ)
	data= concatenate((reshape( xrange,(len(xrange),1)) , reshape(dndz,(len(xrange),1)),reshape(bias,(len(xrange),1)),reshape(kmax,(len(xrange),1)),reshape(kmin,(len(xrange),1)),reshape(Vsurvey,(len(xrange),1)),reshape(dvdz(xrange),(len(xrange),1))),axis=1)
	print; print "_"*30 ; print "Output Files:"
	savetxt('FotranCodes/number_BOSS_S10000.txt' , data)
	print "'number_BOSS_S10000.txt' produced"

	print "_"*30
	# ===========Plot================================
	#fig = plt.figure()
	#ax = fig.add_subplot(1,1,1)
	#p1, = ax.plot(xrange, dndzrange_ref , 'ro')
	#ax.plot(redshift_Euclid, Vsurvey_Euclid_2)
	#ax.plot(xrange,Vsurvey)
	#plt.xlabel(r"redshift ($z$)")
	#plt.ylabel(r"$b(z)$")
	#plt.savefig('inputs/Euclid_dndz.eps')
	#plt.show()
