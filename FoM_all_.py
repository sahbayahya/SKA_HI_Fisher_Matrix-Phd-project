from scipy import *
from numpy import *
import pylab as P
from scipy import linalg
from functions import *

#============================Main progarm ===========================================================================
if __name__== "__main__" :

    	ob = 0.049
    	M = array([[1., 0 , 0., 0., 0.,  0.,  0.,  0.],
               [0., 1., 0., 0., 0.,  0.,  0.,  0.],
               [0., 0., 1., 0., 0.,  0.,  0.,  0.],
               [0., 0., 0., ob,-ob, -ob,  0.,  0.],
               [0., 0., 0., 0., 1.,  0.,  0.,  0.],
               [0., 0., 0.,-ob,-1., -1.,  0.,  0.],
               [0., 0., 0., 0., 0., 0. ,  1.,  0.],
               [0., 0., 0., 0., 0., 0. ,  0.,  1.]])
	sarea_print =[ 3, 5, 23, 70, 150, 200]
	save_param= []
	for i in range(len(sarea_print)):
		A =  loadtxt('Outputs/Fisher_Srms_ska1_'+str(sarea_print[i])+'.txt', unpack=True)
		ska_refuJy  = add_cmb(M, Planks_prior, A, 6)
		save_param.append(FoM_with_prior(ska_refuJy, sarea_print[i]))
		#Sarea, w0, wa,  w0a, wa0, om0, ob0, ok0, h, FoM2
	print 
