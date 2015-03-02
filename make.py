#!/usr/bin/python
"""
    Calculate dn/dz for SKA HI galaxy redshift surveys, using dn/dz curves for
    given flux thresholds (from Mario Santos) and flux scalings with redshift for
    specific arrays.
    """
from numpy import linspace
from constants import *
from functions import *
import os
import subprocess





if __name__== "__main__" :
	#--------------------------------------------------------------------
	#This part handle dNdz
	#Read your file where dN/dz [deg^-2 per unit z] are stored
	# Fit dNdz
	# and plot dNdz vs redshift
	#-------------------------------------------------------------------
	
	(x2, to, rm00muJy,rm01muJy, rm03muJy, rm05muJy, rm06muJy, rm073muJy, rm010muJy, rm023muJy, rm040muJy, rm070muJy, rm100muJy, rm150muJy, rm200muJy) = np.loadtxt('inputs/HIdndzb_modified_high.txt', unpack=True)
	(x1, to, rm0muJy,rm1muJy, rm3muJy, rm5muJy, rm6muJy, rm73muJy, rm10muJy, rm23muJy, rm40muJy, rm70muJy, rm0100muJy, rm0150muJy, rm0200muJy) = np.loadtxt('inputs/HIdndzb_modified.txt', unpack=True)
	(x, total, rms0muJy,rms1muJy, rms3muJy, rms5muJy, rms6muJy, rms73muJy, rms10muJy, rms23muJy, rms40muJy, rms70muJy, rms100muJy, rms150muJy, rms200muJy) = np.loadtxt('inputs/HIdndzb3_corrected.txt', unpack=True)
	#=========================================================================================
	
	''' p0 and p04 are the intial guess for your parameters
		In this case its 3 parameters.
	'''
	p0=[5.52,  0.6, 4.6]
	p04=[5.74, 1.14, 3.95]
	
	'''Define x axis range (or redshift range)
	'''
	xrange = np.linspace(0, 3.0, 200)
	
	# fit high srmses in the log space
	(yt,xrange)=  find_parameters(p0,  x,np.log(total),0,xrange)
	(y0t,xrange)= find_parameters(p0,  x,np.log(rms0muJy), 0,xrange)
	(y1t,xrange)= find_parameters(p0,  x,np.log(rms1muJy), 1,xrange)
	(y3t,xrange)= find_parameters(p0,  x,np.log(rms3muJy), 3, xrange)
	(y5t,xrange)= find_parameters(p0,  x,np.log(rms5muJy) ,5, xrange)
	(y6t,xrange)= find_parameters(p0,  x,np.log(rms6muJy), 6, xrange)
	(y7t,xrange)= find_parameters(p0,  x,np.log(rms73muJy), 7, xrange)
	(y10t,xrange)=find_parameters(p0,  x,np.log(rms10muJy), 10, xrange)
	# fit the low Srmses in normal space
	(y23t,xrange)= find_parameters_(p0, x1,rm23muJy, 23, xrange)
	(y40t,xrange)= find_parameters_(p0, x1,rm40muJy, 40, xrange)
	(y70t,xrange)= find_parameters_(p0, x1,rm70muJy, 70, xrange)
	(y100t,xrange)=find_parameters_(p0, x2,rm100muJy,100, xrange)
	(y150t,xrange)=find_parameters_(p0, x2,rm150muJy,150, xrange)
	(y200t,xrange)=find_parameters_(p0, x2,rm200muJy,200, xrange)
	print '============ Program excuted successfully ==========='
	'''
	Plot the results from this program
	'''
	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)
	ax.set_yscale('log')
	yt = ([np.exp(y0t), np.exp(y1t), np.exp(y3t), np.exp(y5t), np.exp(y10t), y23t, y100t, y200t])
	tic = [r'  $  \ 0 \mu$Jy', r'  $ \ 1 \mu$Jy', r'  $ \ 3 \mu$Jy', r'  $ \ 5 \mu$Jy', r'  $  \ 10\mu$Jy' , r'  $ \ 23\mu$Jy' ,r'  $  \ 100\mu$Jy'  ,  r'  $ \ 200\mu$Jy']
	'''plot
	'''
	nt = len(yt)
	for i in range(nt):
		color = i / float(nt)
		heatmap2=colorline(xrange, yt[i], color, cmap="RdBu")
	colors = ['#67001f', '#cc595d' ,'#f2a789' , '#fce4d6', '#ddebf3','#93c4df', '#2f79b6', '#4083bc']
	Srms= [0, 1, 3, 5, 10, 23, 100, 200]
	rms = [rms0muJy,rms1muJy, rms3muJy, rms5muJy, rms10muJy, rms23muJy,
		   rms100muJy, rms200muJy]
	for i in range (len(Srms)):
		ax.scatter(x, (rms[i]),s= 35, marker= 'o', edgecolor =colors[i], facecolor=colors[i])
	'''sidebar
	'''
	cbar = plt.colorbar(heatmap2, ticks=())
	cbar.ax.set_yticklabels(tic)
	for j, lab in enumerate(tic):
		cbar.ax.text(.6, (1* j) /float(nt), lab, va='center')
	plt.xlim(0.1,2.2 ,0.2)
	plt.ylim(1, 5e6)
	'''x axis
	'''
	xticks = np.arange(min(x), max(x)+1, 0.3)
	plt.xlabel(r"${ \rm redshift} (z)$", fontsize=15)
	'''y axis
	'''
	plt.ylabel(r'$\frac{d{\rm N}}{dz}(z) \ [ {\rm deg}^{-2} \ {\rm per} \ {\rm unit} \ z ]$', fontsize= 15)
	'''save figure
	'''
	plt.savefig("plots/fittingMario_dNOverdz_using_ObreschkowFunc.pdf")
	#plt.show()
	print "==================The Number Density part computed successfully  ! ======================"

	#----------------------------------------------------------------
	# This part handle the Bias
	# Input data, note that the bais^2 in the input files
	# Fit the Bias
	# Plot Bias vs redshift
	#----------------------------------------------------------------
	def read_input(v):
			z = [] ; bz = []
			for i in range(n):
				if Srms_[i] == v:
					z.append(redshift[i]); bz.append(bias[i])
			return z, bz
	(redshift, bias, Srms_) = np.loadtxt('inputs/bias_mill.out', unpack= True)
	n = len(redshift)
	y = ['y0', 'y1', 'y3', 'y5', 'y10', 'y23', 'y100', 'y200']
	z = ['z0', 'z1', 'z3', 'z5', 'z10', 'z23', 'z100', 'z200']
	b =['b0', 'b1', 'b3', 'b5', 'b10', 'b23',  'b100', 'b200']
	Srms= [0, 1, 3, 5, 10, 23, 100, 200]
	z_ = [] ; b_ = []
	for i in range (len(Srms)):
		z[i], b[i] = read_input(Srms[i])
		z_.append(z[i]) ; b_.append(b[i])
	'''
	The initial guess
	'''
	p0=[6.3,2.]    ;    p04=[5.74, 1.14]
	''' x range
	'''
	xrange = [linspace(0, 2.5, 200),  linspace(0, 2.5, 200),  linspace(0, 2.5, 200),  linspace(0, 2.5, 200)\
	,  linspace(0, 2.5, 200),  linspace(0, 2., 200),  linspace(0, 0.8, 200), linspace(0, 0.6, 200), ]
	'''
	Fit the bias to a fucntion
	'''
	y_ = []
	for i in range(len(Srms)):
		 y_.append(find_Bias_parameters(p0, np.array(z_[i]), np.array(b_[i]) , Srms[i], xrange[i]))

	'''plot using cmap
	'''
	y = y_ ; b= b_
	tic = [ r'   $0 \mu$Jy',r'   $1 \mu$Jy', r'   $3 \mu$Jy', r'   $5 \mu$Jy', \
				   r'   $10 \mu$Jy' ,  r'   $23\mu$Jy' ,r'   $100\mu$Jy',r'   $200\mu$Jy',  ]
	p_tic =[0, 1, 3, 5,10, 23, 100, 200]

	'''plot the Bias vs redshift
	'''
	fig, ax= plt.subplots()
	n = len(y)
	for i in range(n):
		color = i / float(n)
		heatmap=colorline(xrange[i], y[i], color, cmap="RdBu" , linewidth=1.5)
	'''
	sidebar
	'''
	cbar = plt.colorbar(heatmap, ticks=())
	cbar.ax.set_yticklabels(tic)
	for j, lab in enumerate(tic):
		cbar.ax.text(.6, (1* j) /float(n), lab, va='center')
	colors = ["#67001f", "#cc595d" ,"#f2a789" , "#fce4d6", "#ddebf3","#93c4df", "#2f79b6", "#4083bc"]
	'''loop to plot the data
	'''
	for i in range (len(z_)):
		ax.scatter(z_[i], b_[i],s= 35, marker= "o", edgecolor =colors[i], facecolor=colors[i])
	'''
	y axis
	'''
	plt.ylabel(r"$b(z)$", fontsize=15)
	plt.ylim(0., 5)
	'''
	x axis
	'''
	plt.xlabel(r"$ {\rm redshift} (z)$", fontsize=15)
	plt.xlim(0.01,2.2)

	plt.savefig("plots/fitted_bias.pdf")
	#plt.show()
	print "==================The Bias part computed successfully!======================"

	import fit_Euclid
	import fit_BOSS
	#---------------------------------------------------------
	# Flux rms limits at 1GHz for various configurations
	# The output used as input for the fortran code
	#--------------------------------------------------------
	ID =0
	while ID < 18:
		Construct_2D(Srms, Sarea, ID)
		ID = ID + 1



	for ii in range(len(nameint)):
		
		filein="/Users/sahba/Dropbox/SKA_Survey/SKA_HI_Fisher_Matrix/Outputs/%d_Corrected_nz_Srms.txt"%(nameint[ii],)
		fileout="/Users/sahba/Dropbox/SKA_Survey/SKA_HI_Fisher_Matrix/Outputs/output_Srms_ska1_%d.txt"%(nameint[ii],)
		fileFish="/Users/sahba/Dropbox/SKA_Survey/SKA_HI_Fisher_Matrix/Outputs/Fisher_Srms_ska1_%d.txt"%(nameint[ii],)
		print "we are running ",filein
		print "outputting to files ",fileout
		print ""

		#===================================================================================
		#Now lets run the Fortran code to get the resutls for the Fisher Matrix
		#===================================================================================
		cwd = os.getcwd()
		os.chdir("/Users/sahba/Dropbox/SKA_Survey/Fishers/fisher_distance_bao/Fisher_bao_SKA_Euclid/diff_survey_area_SKA1")
		output = subprocess.check_output(["./fisher_distance_bao_w0_wa_marginlized_over_Omegam_Omegak_h_2",filein,fileFish,fileout])
		
		print "Done with Fortran Code.. Back to python"
		os.chdir(cwd)
		
		

	import plot_Err_Da
	import plot_Err_H
	#==========================================================================
	# Now lets calculate the FoM
	#==========================================================================


	print  '  w0', '   wa',  '   w0a', '   wa0',  '   FoM2'
	import FoM_all_with_prior
	import FoM_all_no_prior

	#====================================================================
	# plot the contours
	#====================================================================

	import plot_w_vs_wa
	import plot_Ok_vs_w
	print '---------------------- Program excuted successfully--------------------'
	print '-----The plot is saved in plots folder and the txt files in Ouputs-----'
	print '--------------------------Thanks! -------------------------------------'

	