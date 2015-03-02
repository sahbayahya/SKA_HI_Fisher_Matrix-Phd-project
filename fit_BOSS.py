import numpy as np
import scipy.optimize as opt
from scipy import interpolate, concatenate, reshape, sqrt, savetxt, linspace, exp, sin, log
import matplotlib.pyplot as plt
import cosmolopy.constants as cc
import cosmolopy.distance as cd
import cosmolopy.perturbation as cp
import sys
import os
import subprocess
from constants import *
from functions import k_max, background_evolution_splines, dvdz, da, V_sur
#sys.path.append("/home/nassp/sahba/Desktop/Fisher_Euclid_numeric_derivative/cosmolopy")


if __name__=='__main__':
	print "_"*30

	xrange = np.array([ 0.05, 0.15 , 0.25 ,0.35 , 0.45 , 0.55, 0.65 , 0.75])
	dndz   = np.array([ 8, 50, 125, 222, 332, 447, 208, 30])
	xmin = xrange-0.05
	xmax = xrange+0.05
	print 'zmin = ' , xmin
	print 'zmax = ', xmax 
	print 'zc = ',  (xmin + xmax)/2.0

	#=====================================================
	#    k [Mpc^-1 h] range
	#=====================================================
	kmax = k_max(xrange,k_NL)# Output survey info
	_kmax = np.linspace(0.2, 0.2, 8)
	kmin = np.linspace(0.00334,0.00334, 8)
	print 'kmax = ' , kmax
	print '_kmax=', _kmax
	x = xrange
	D_zin  =  np.empty(len(x))
	h= 0.67
	boss_sarea = np.array([5000, 10000,1000, ])
	ErrorZ= np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
	bias =  np.array([ 2.0, 2.0 , 2.0 ,2.0 , 2.0 , 2.0, 2.0 ,2.0])
	
	for ii in range(len(boss_sarea)):
		#===============================================================
		#    V survey [Mpc^3 h^-3]
		#===============================================================
		
		Vsurvey= []
		for i in range(len(xrange)):
			Vsurvey.append(V_sur(xmin[i],xmax[i], boss_sarea[ii])*(h**3))
		print 'Vsurvey = ', Vsurvey

		#==================================================================
		#   Save the resutls
		#==================================================================

		print len(bias), len(kmax), len(kmin), len(Vsurvey), len(dvdz(xrange)), len(ErrorZ)
		data= concatenate((reshape( xrange,(len(xrange),1)) , reshape(dndz,(len(xrange),1)),reshape(bias,(len(xrange),1)),reshape(kmax,(len(xrange),1)),reshape(Vsurvey,(len(xrange),1)),reshape(dvdz(xrange),(len(xrange),1))),axis=1)
		print; print "_"*30 ; print "Output Files:"
		
		savetxt("inputs/number_BOSS_S%d.txt"%(boss_sarea[ii],) , data)
		print " 'inputs/number_BOSS_S%d.txt' produced"%(boss_sarea[ii],)
		print "_"*30

		filein="/Users/sahba/Dropbox/SKA_Survey/SKA_HI_Fisher_Matrix/inputs/number_BOSS_S%d.txt"%(boss_sarea[ii],)
		fileout="/Users/sahba/Dropbox/SKA_Survey/SKA_HI_Fisher_Matrix/Outputs/output_BOSS_S%d.txt"%(boss_sarea[ii],)
		fileFish="/Users/sahba/Dropbox/SKA_Survey/SKA_HI_Fisher_Matrix/Outputs/Fisher_BOSS_S%d.txt"%(boss_sarea[ii],)
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
	print "=====================Success============================="
