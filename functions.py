import numpy as np
import pylab as P
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import scipy.interpolate
from scipy.integrate import quad, cumtrapz, simps

import scipy.optimize
import scipy.optimize as opt
import sys
import cosmolopy.constants as cc
import cosmolopy.distance as cd
import cosmolopy.perturbation as cp
#from experiments import ID
from constants import *

''' -------------------------------------------------------------------------------------------------------------------
    These functions related to frequency corrected Srms case
    -------------------------------------------------------------------------------------------------------------------
'''

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
def background_evolution_splines(cosmo, zmax=10., nsamples=500):
    """
	Get interpolation functions for background functions of redshift:
	* H(z), Hubble rate in km/s/Mpc
	* r(z), comoving distance in Mpc
	* D(z), linear growth factor
	* f(z), linear growth rate
	"""
    _z = np.linspace(0., zmax, nsamples)
    a = 1. / (1. + _z)
    H0 = (100.*cosmo['h']); w0 = cosmo['w0']; wa = cosmo['wa']
    om = cosmo['omega_M_0']; ol = cosmo['omega_lambda_0']
    ok = 1. - om - ol
    
    # Sample Hubble rate H(z) and comoving dist. r(z) at discrete points
    omegaDE = ol * np.exp(3.*wa*(a - 1.)) / a**(3.*(1. + w0 + wa))
    E = np.sqrt( om * a**(-3.) + ok * a**(-2.) + omegaDE )
    _H = H0 * E
    
    r_c = np.concatenate( ([0.], scipy.integrate.cumtrapz(1./E, _z)) )
    if ok > 0.:
        _r = C/(H0*np.sqrt(ok)) * np.sinh(r_c * np.sqrt(ok))
    elif ok < 0.:
        _r = C/(H0*np.sqrt(-ok)) * np.sin(r_c * np.sqrt(-ok))
    else:
        _r = (C/H0) * r_c
    
    # Integrate linear growth rate to find linear growth factor, D(z)
    # N.B. D(z=0) = 1.
    a = 1. / (1. + _z)
    Oma = cosmo['omega_M_0'] * (1.+_z)**3. * (100.*cosmo['h']/_H)**2.
    _f = Oma**cosmo['gamma']
    _D = np.concatenate( ([0.,], scipy.integrate.cumtrapz(_f, np.log(a))) )
    _D = np.exp(_D)
    
    # Construct interpolating functions and return
    r = scipy.interpolate.interp1d(_z, _r, kind='linear', bounds_error=False)
    H = scipy.interpolate.interp1d(_z, _H, kind='linear', bounds_error=False)
    D = scipy.interpolate.interp1d(_z, _D, kind='linear', bounds_error=False)
    f = scipy.interpolate.interp1d(_z, _f, kind='linear', bounds_error=False)
    return H, r, D, f

def extend_with_linear_interp(xnew, x, y):
    """
	Extend an array using a linear interpolation from the last two points.
	"""
    dx = x[-1] - x[-2]
    dy = y[-1] - y[-2]
    ynew = y[-1] + dy * (xnew - x[-1]) / dx
    y = np.concatenate((y, [ynew,]))
    return y

def n_bin(zmin, zmax, dndz, bias=None):
    """
    Number density of galaxies in a given z bin (assumes full sky). Also
    returns volume of bin. dndz argument expects an interpolation fn. in units
    of deg^-2.
    """
            
    # Calculate cosmo. functions
    H, r, D, f = background_evolution_splines(cosmo)
    _z = np.linspace(zmin, zmax, 500)
    vol = 4.*np.pi*C * scipy.integrate.simps(r(_z)**2. / H(_z), _z)
    N_bin = FULLSKY * scipy.integrate.simps(dndz(_z), _z)
    nz = N_bin / vol
    
    # Calculate mean bias (weighted by number density)
    if bias is not None:
        b = scipy.integrate.simps(bias(_z)*dndz(_z), _z) / (N_bin / FULLSKY)
        return nz, vol, b
    return nz, vol

def redshift_bins( ID,dz=0.1, Nbins=None):
    """
	Calculate redshift bins.
	"""
    zmin = NU_LINE*1e3 / numax[ID] - 1.
    zmax = NU_LINE*1e3 / numin[ID] - 1.
	#print zmin, zmax, 100000000
    if zmin < 0.: zmin = 0.
    if Nbins is not None:
        zbins = np.linspace(zmin, zmax, Nbins+1)
    else:
        Nbins = np.floor((zmax - zmin) / dz)
        zbins = np.linspace(zmin, zmin + dz*Nbins, Nbins+1)
        if zmax - np.max(zbins) > 0.04:
            zbins = np.concatenate((zbins, [zmax,]))
    return zbins

def flux_redshift(z, ID):
    """
	Flux rms as a function of redshift.
	"""
    z = np.atleast_1d(z)
    nu = NU_LINE / (1. + z)
    if nucrit[ID] is not None:
        Sz = fluxrms[ID] * Scorr[ID] * np.ones(nu.size)
        idxs = np.where(nu*1e3 > nucrit[ID])
        Sz[idxs] *= (nu[idxs]*1e3 / nucrit[ID])
    else:
        Sz = fluxrms[ID] * Scorr[ID]
        Sz = nu * Sz if not Sconst[ID] else Sz * np.ones(nu.size)
    return Sz

def k_max(z,k_NL):
    '''
    A function to determine the maximum wave number => k_max
	you need to define too paramters
	* K_Nl
	* cosmo['ns']
    '''
    kmax = k_NL*(1.+z)**(2./(2.+cosmo['ns']))
    return kmax


def Construct_2D(Srms, Sarea, ID):

	'''Construct grid of dn/dz (deg^-2) as a function of flux rms and redshift and
	then construct 2D interpolator, you need the fitting paramters and Srms as arrays
	and you need to define the Sarea
	* c1 = dNdz parameter 
	* c2 = dNdz fitting parameter
	* c3 = dNdz fitting parameter 
	* c4 = Bias fitting parameter
	* c5 = Bias fitting parameter
	'''
	# Define fitting coefficients from Mario's note (HI_specs.pdf)
	Srms = np.array([0., 1., 3., 5., 6., 7.3, 10., 23., 40., 70., 100., 150.,])#200.,])
	c1 = [6.124, 6.556, 6.532, 6.551, 6.577, 6.555, 6.443, 6.02, 5.74, 5.62, 5.63, 5.48,]# 5.00]
	c2 = [1.717, 2.0185, 1.932, 1.932, 1.946, 1.923, 1.831, 1.43, 1.22, 1.11, 1.41, 1.33,]# 1.04]
	c3 = [0.789, 3.810, 5.225, 6.225,6.685, 7.078, 7.585, 9.03, 10.58, 13.03, 15.49, 16.62,]# 17.52]
	#c4 = [0.8695, 0.5863, 0.4780, 0.5884, 0.5908, 0.5088, 0.4489, 0.5751, 0.5125, 0.6193, 0.6212, 1., 1.]#, 1.]
	#c5 = [0.2338, 0.6410, 0.9181, 0.8076, 0.8455, 1.0222, 1.2069, 0.9993, 1.1842, 1.0179, 1.0759, 0., 0.]#, 0.]
	c4 =  [0.587, 0.497, 0.530, 0.550, 0.547, 0.562, 0.593, 0.607, 0.628, 0.609, 0.605,0.637, 1.0]
	c5 =  [0.358, 0.720, 0.781, 0.801, 0.829, 0.823, 0.807, 0.852, 0.844, 0.929, 1.086,0.965, 0.0]
	c1 = np.array(c1); c2 = np.array(c2); c3 = np.array(c3)
	c4 = np.array(c4); c5 = np.array(c5)

	Smax = np.max(Srms)
	# Extrapolate fitting functions to high flux rms
	c1 = extend_with_linear_interp(SBIG, Srms, c1)
	c2 = np.concatenate((c2, [1.,])) # Asymptote to linear fn. of redshift
	c3 = extend_with_linear_interp(SBIG, Srms, c3)
	Srms = np.concatenate((Srms, [SBIG,]))
	
	# Construct grid of dn/dz (deg^-2) as a function of flux rms and redshift and
	# then construct 2D interpolator
	z = np.linspace(0., 4., 400)
	nu = NU_LINE / (1. + z)
	_dndz = np.array([10.**c1[j] * z**c2[j] * np.exp(-c3[j]*z) for j in range(Srms.size)])
	_bias = np.array([c4[j] * np.exp(c5[j]*z) for j in range(Srms.size)])
	
	#print _dndz
	dndz = scipy.interpolate.RectBivariateSpline(Srms, z, _dndz, kx=1, ky=1)
	bias = scipy.interpolate.RectBivariateSpline(Srms, z, _bias, kx=1, ky=1)
	
	#print dndz
	# Construct dndz(z) interpolation fn. for the sensitivity of actual experiment
	fsky = Sarea[ID] / FULLSKY
	Sz = flux_redshift(z, ID)
	dndz_expt = scipy.interpolate.interp1d(z, dndz.ev(Sz, z))
	bias_expt = scipy.interpolate.interp1d(z, bias.ev(Sz, z))
	#print 'dndz_expt', dndz_expt
	
	# Fit function to dn/dz [deg^-2]
	_z = np.linspace(0., 1., 100)
	dndz_vals = dndz_expt(_z)
	bias_vals = bias_expt(_z)
	p0 = [100.*np.max(dndz_vals), 2., 10.]
	def lsq(params):
		A, c2, c3 = params
		model = A * _z**c2 * np.exp(-c3*_z)
		return model - dndz_vals
	p = scipy.optimize.leastsq(lsq, p0)[0]
	
	# Fit function to bias
	p0 = [np.max(bias_vals), 0.5]
	def lsq(params):
		c4, c5 = params
		model = c4 * np.exp(c5*_z)
		return model - bias_vals
	pb = scipy.optimize.leastsq(lsq, p0)[0]
	
	# Print best-fit coefficients
	print "-"*30
	#print name[ID]
	#print "-"*30
	print "Fitting coeffs."
	print "c1: %6.4f" % np.log10(p[0])
	print "c2: %6.4f" % p[1]
	print "c3: %6.4f" % p[2]
	print "c4: %6.4f" % pb[0]
	print "c5: %6.4f" % pb[1]
	#print name[ID], '&', " & ".join(["%6.4f" % n for n in [ np.log10(p[0]), p[1], p[2], pb[0], pb[1]]])
	

	# Calculate cosmo. functions
	H, r, D, f = background_evolution_splines(cosmo)
	
	# Calculate binned number densities
	zbins = redshift_bins(ID, dz=0.1)
	zc = np.array([0.5*(zbins[i] + zbins[i+1]) for i in range(zbins.size-1)])
	nz, vol, b = np.array( [n_bin(zbins[i], zbins[i+1], dndz_expt, bias_expt)
							for i in range(zbins.size-1)] ).T
	vol *= fsky
	k_NL = 0.2
	h= 0.67
	_kmax = k_max(zc,k_NL)# Output survey info
	print "-"*30
	print "zc      n Mpc^-3 h^3    bias     kmax    Vsur Mpc^3/h^3     dvdz      Srms"
	for i in range(zc.size):
		Szz = flux_redshift(zc[i], ID)
		print "%2.1f     %3.3e    %6.3f    %6.3f       %5.3e    %5.3e          %6.2f" % \
		(zc[i],dndz_expt(zc[i]), b[i], _kmax[i], vol[i]*(h**3), dvdz(zc[i]), Szz),
		if (nz[i]*vol[i]) < NGAL_MIN:print "*",
		if Szz > Smax: print "#",
		print ""
		data = np.concatenate((np.reshape(zc,(len(zc),1)),np.reshape(dndz_expt(zc),(len(zc),1)),\
							   np.reshape(b,(len(zc),1)),np.reshape(_kmax,(len(zc),1)),np.reshape(vol*(h**3),(len(zc),1)),\
							   np.reshape(dvdz(zc),(len(zc),1))),axis=1)
		np.savetxt("Outputs/"+name[ID]+"_Corrected_nz_Srms.txt",data)
	print "-"*30
	print "Ntot: %3.3e" % np.sum(nz * vol)
	print "fsky: %3.3f" % fsky
	_zmin = (NU_LINE*1e3 / numax[ID] - 1.)
	print "zmin: %3.3f" % (_zmin if _zmin >= 0. else 0.)
	print "zmax: %3.3f" % (NU_LINE*1e3 / numin[ID] - 1.)
	print "Srms const: %s" % Sconst[ID]
	print "Exp name:", name_rel[ID], "Srms=", name[ID]
	print "-"*30
	print "\n"
	return

#--------------------------------------------------------------------------------------
#    These functions to produce the fitting parameters for the Bias
#-----------------------------------------------------------------------------------

def make_segments(x, y):
    '''
	Create list of line segments from x and y coordinates, in the correct format for LineCollection:
	an array of the form   numlines x (points per line) x 2 (x and y) array
	'''
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments

# Interface to LineCollection:

def colorline(x, y, z=None, cmap=plt.get_cmap('hot'), norm=plt.Normalize(0.0, 1.0), linewidth=3, alpha=1.0):
    '''
	Plot a colored line with coordinates x and y
	Optionally specify colors in the array z
	Optionally specify a colormap, a norm function and a line width
	'''
    # Default colors equally spaced on [0,1]:
    if z is None:
        z = np.linspace(0.0, 1.0, len(x))
    # Special case if a single number:
    if not hasattr(z, "__iter__"):  # to check for numerical input -- this is a hack
        z = np.array([z])
    z = np.asarray(z)
    segments = make_segments(x, y)
    lc = LineCollection(segments, array=z, cmap=cmap, norm=norm, linewidth=linewidth, alpha=alpha)
    ax = plt.gca()
    ax.add_collection(lc)
    return lc

def clear_frame(ax=None):
    # Taken from a post by Tony S Yu
    if ax is None:
        ax = plt.gca()
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    for spine in ax.spines.itervalues():
        spine.set_visible(False)


    '''
    This function is the function we think it will be easy
    to fit its parameters to our data
    '''
def Fit_Bias(p,x):
    w=p[0]*np.exp(p[1]*x)
    #print w.size
    return w
    '''
    The purpose of this function is finding the
    difference between the theortical and the simulated data point
    at specific point x (or redshift).
    '''
def residuals_Bias(p,x,y):
    w=p[0]* np.exp(p[1]*x)
    err=w-y
    err=err**2
    B=sum(err)
    return B
    '''
    This function call opt.fmin function to  minimize the error on
    you paramters and give you the best fit parameters
    '''
def find_Bias_parameters(p0, x, bias_rms , rms_value, xrange):
    plsqtot = opt.fmin(residuals_Bias, p0, args=(x,bias_rms), maxiter=10000, maxfun=10000)
    print  'b(z) = %d uJy' % rms_value, 'p[0] = ' ,plsqtot[0],  'p[1] = ' , plsqtot[1]
    ytot=Fit_Bias(plsqtot,xrange)
    return ytot



#=======================================================================================
# This is the funtions to fit dNdz
#===============================================================================================


def Fit_dNdz(p,x):
	''' 
	This function is the function we think it will be easy
	to fit its parameters to our data
	'''
	w=10.**p[0]*  x**p[1]  * np.exp(-p[2]*x)
	print w.size
	return w



def residuals_dNdz(p,x,y):
	'''
	The purpose of this function is finding the
	difference between the theortical and the simulated data point
	at specific point x (or redshift).
	'''
	w=10.**p[0]* x**(p[1] ) *  np.exp(-p[2]*x)
	err=w-y
	err=err**2
	B=sum(err)
	return B
def find_parameters_(p0, x, rms , rms_value, xrange):
	'''
	This function call opt.fmin function to  minimize the error on
	you paramters and give you the best fit parameters
	'''
	plsqtot = opt.fmin(residuals_dNdz, p0, args=(x,rms), maxiter=10000, maxfun=10000)
	print  'rms = %d uJy' % rms_value, 'p[0] = ' ,plsqtot[0],  'p[1] = ' , plsqtot[1],'p[2] = ', plsqtot[2]
	print '========================================'
	ytot=Fit_dNdz(plsqtot,xrange)
	return ytot, xrange
#==============================================================================
def func(p,x):
	''' 
	This function is the function we think it will be easy
	to fit its parameters to our data
	'''
	w=10.**p[0]*  x**p[1]  * np.exp(-p[2]*x); w=np.log(w)
	print w.size
	return w



def residuals(p,x,y):
	'''
	The purpose of this function is finding the
	difference between the theortical and the simulated data point
	at specific point x (or redshift).
	'''
	w=np.log(10.**p[0]* x**(p[1] ) *  np.exp(-p[2]*x))
	err=w-y
	err=err**2
	B=sum(err)
	return B

def find_parameters(p0, x, rms , rms_value, xrange):
	'''
	This function call opt.fmin function to  minimize the error on
	you paramters and give you the best fit parameters
	'''
	plsqtot = opt.fmin(residuals, p0, args=(x,rms), maxiter=10000, maxfun=10000)
	print  'rms = %d uJy' % rms_value,'p[0] =' ,plsqtot[0],  'p[1] =' , plsqtot[1], 'p[2]=', plsqtot[2]
	print '========================================'
	ytot=func(plsqtot,xrange)
	return ytot, xrange




'''---------------------------------------------------------------------------------------
	These functions concern the FoM case
   --------------------------------------------------------------------------------------
'''



from scipy import linalg
def FoM_no_prior(matrix_plus_prior, area):
	'''
		:param matrix_plus_prior: The fisher matrix that you would like
		find FoM of w and w0 for it,
		:return:FoM
	'''
	matrix_plus_prior =  linalg.inv(matrix_plus_prior)
	#============parameters========================
	w0 = np.sqrt(matrix_plus_prior[0,0])#/1e-1 #;   print 'sigma_w0 = ', w0
	wa = np.sqrt(matrix_plus_prior[1,1])   # ; print 'sigma_wa =', wa
	w0a = (matrix_plus_prior[1,0])  #;  print 'sigma w0a = ', w0a
	wa0 =((matrix_plus_prior[0,1]))   #; print' wa0 = ', wa0
	ob0 =np.sqrt(matrix_plus_prior[2,2])*0.049 #/1e-4   #print'sigma_ob0 =', ob0
	ok0 = np.sqrt((matrix_plus_prior[3,3])) #/1e-2  #print 'sigma_ok0=', ok0
	om0 = np.sqrt(matrix_plus_prior[4,4])#/1e-2  #print 'sigma_om0 = ', om0
	h = np.sqrt(matrix_plus_prior[5,5])#/1e-2    #print 'sigma_h = ',  h
	FoM2 =  1.0/np.sqrt(matrix_plus_prior[1,1] * matrix_plus_prior[0,0]
					 - matrix_plus_prior[1,0]* matrix_plus_prior[0,1])#/(pi*(sqrt(2.31)))
	return area,  FoM2

def add_cmb(M, prior_fish, matrix, n):
	'''
	This function add the prior matrix to the Matrix

	:param M: convert the parameters of the prior matrix to the parameters we desire
	:param prior_fish: The prior matrix with the initial parameters
	:param matrix: The matrix we would like to add the prior matrix to it
	:param n: the nxn number of the coloumn and raws of the matrices you want to add
	:return: the inverse of Full matrix = inverse of (prior matrix + matrix)
	'''
	#======== convert the parameters ============
	MT = M.T
	M11 = np.dot(M, prior_fish)
	Final_prior_Fisher = np.dot(M11 , MT)
	#===== Stack new columns and raws===========
	newraw = np.linspace(0.,0.,n)
	matrix = np.vstack((matrix,newraw))
	matrix = np.vstack((newraw,matrix))
	newcolumn = np.linspace(0., 0., n+2)
	matrix = np.column_stack(( matrix, newcolumn))
	matrix = np.column_stack((newcolumn, matrix))
	matrix_plus_prior = matrix +  Final_prior_Fisher
	matrix_plus_prior =  linalg.inv(matrix_plus_prior)
	return matrix_plus_prior

def FoM_with_prior(matrix_plus_prior, Sarea):
	
	'''
		:param matrix_plus_prior: The fisher matrix that you would like
		find FoM of w and w0 for it,
		:return:FoM
	'''
	#============parameters========================
	w0 = np.sqrt(matrix_plus_prior[1,1])
	wa = np.sqrt(matrix_plus_prior[2,2])
	w0a = (matrix_plus_prior[2,1])
	wa0 =((matrix_plus_prior[1,2]))
	ob0 =np.sqrt(matrix_plus_prior[2,2])*0.049 #/1e-4   #print'sigma_ob0 =', ob0
	ok0 = np.sqrt((matrix_plus_prior[3,3])) #/1e-2  #print 'sigma_ok0=', ok0
	okw =(matrix_plus_prior[1,3])
	om0 = np.sqrt(matrix_plus_prior[4,4])#/1e-2  #print 'sigma_om0 = ', om0
	h = np.sqrt(matrix_plus_prior[5,5])#/1e-2    #
	FoM2 =  1.0/np.sqrt(matrix_plus_prior[1,1] * matrix_plus_prior[2,2]
					 - matrix_plus_prior[1,2]* matrix_plus_prior[2,1])
	print Sarea,   w0, wa,  w0a, okw, om0, ob0,h, FoM2
	return Sarea, FoM2
