import numpy as np
from scipy import interpolate, concatenate, reshape, sqrt, savetxt, linspace, exp, sin, log
from scipy.integrate import quad, cumtrapz, simps
import sys
import os
import subprocess
from constants import *
from functions import dvdz, k_max, V_sur

#================================
# fuctions
#================================

def Bias_Euclid(z):
	'''This is the bias function which is used by Euclid survey
		b(z) = sqrt(1+z)
		'''
	return  sqrt(1+z)


if __name__=='__main__':
	xrange = np.array([ 0.7, 0.8,  0.9 ,  1.0, 1.1 , 1.2 ,1.3  , 1.4 ,1.5 , 1.6 ,1.7 , 1.8 ,1.9, 2.0])
	dndzrange_ref = np.array([ 1.25, 1.92, 1.83, 1.68, 1.51, 1.35, 1.20, 1.00, 0.80, 0.58, 0.38, 0.35, 0.21, 0.11])
	xmin = xrange-0.05
	xmax = xrange+0.05
	#===================================
	#k [Mpc^-1 h] range
	#===================================
	kmax = np.linspace(0.16004, 0.2, 14)
	kmin = np.linspace(0.00435,0.00334, 14)
	x = xrange
	D_zin  =  np.empty(len(x))#; kmax.fill(D(1.))

	#================================
	#V survey [Mpc^3 h^-3]
	#=================================
	area_Euclid = 15000.0
	Vsurvey= []
	for i in range(len(xrange)):
          	Vsurvey.append(V_sur(xmin[i],xmax[i], area_Euclid)*(h**3))
	#==============================
	# The error on z
	#==============================

	ErrorZ= 0.001*(1. + xrange)
	dndz = dndzrange_ref*((3.14159/180.)**2*dvdz(xrange)*h**3)*10**(-3)
	
	#==============================
	#		Save the resutls
	#=============================
	
	data= concatenate((reshape( xrange,(len(xrange),1)) , reshape(dndz,(len(xrange),1)),reshape(Bias_Euclid(xrange),(len(xrange),1)),reshape(kmax,(len(xrange),1))	,reshape(Vsurvey,(len(xrange),1)),reshape(dvdz(xrange),(len(xrange),1))),axis=1)
	print; print "_"*30 ; print "Output Files:"
	savetxt('inputs/number_EuclidmJy_ref.txt' , data)
	print "'number_EuclidmJy_ref.txt' produced"
	print "_"*30
		
	filein="/Users/sahba/Dropbox/SKA_Survey/SKA_HI_Fisher_Matrix/inputs/number_EuclidmJy_ref.txt"
	fileout="/Users/sahba/Dropbox/SKA_Survey/SKA_HI_Fisher_Matrix/Outputs/output_EuclidmJy_ref.txt"
	fileFish="/Users/sahba/Dropbox/SKA_Survey/SKA_HI_Fisher_Matrix/Outputs/Fisher_EuclidmJy_ref.txt"
	print "we are running ",filein
	print "outputting to files ",fileout
	print "="*30
		
	#===================================================================================
	#    Now lets run the Fortran code to get the resutls for the Fisher Matrix
	#===================================================================================
	cwd = os.getcwd()
	os.chdir("/Users/sahba/Dropbox/SKA_survey/SKA_HI_Fisher_Matrix/FortranCode")
	print "got to this directory:"
	output = subprocess.check_output(["./fisher_distance_bao_w0_wa_marginlized_over_Omegam_Omegak_h_2",filein,fileFish,fileout])
	os.chdir(cwd)
	#====================================================================================
	
	print "Done with Fortran Code.. Back to python"
	print "==================success!=============="

