from scipy import *
from numpy import *
import pylab as P
from scipy import linalg
#from store_Fisher_matrix_with_prior2 import *

#==================Plank's prior matrix  from DETF=====================================================================
Planks_prior = array([[1.99579245e+05,  -3.73667528e+04, -1.04936812e+04 ,  1.39977603e+06  ,  5.58643962e+05 , -4.64225267e+04 , -7.65181989e+04 , -2.23806234e+03],
                    [-3.73667528e+04,   1.83928663e+05,   5.16525685e+04 , -7.42050738e+06  , -3.98758357e+06 , -1.11710442e+06 ,  1.32438370e+06 , -4.51559188e+02],
                    [-1.04936812e+04,   5.16525685e+04,   1.45055577e+04 , -2.08389634e+06  , -1.11983054e+06 , -3.13715719e+05 ,  3.71925825e+05 , -1.26811078e+02],
                    [1.39977603e+06,  -7.42050738e+06 ,  -2.08389634e+06 ,  3.64943809e+08  ,  1.58599621e+08 ,  4.25932543e+07 , -5.16878541e+07 ,  3.20338905e+04],
                    [5.58643962e+05,  -3.98758357e+06 ,  -1.11983054e+06 ,  1.58599621e+08  ,  8.70535526e+07 ,  2.48738854e+07 , -2.91740427e+07 ,  1.88438127e+04],
                    [-4.64225267e+04,  -1.11710442e+06,  -3.13715719e+05 ,  4.25932543e+07  ,  2.48738854e+07 ,  7.49686718e+06 , -8.54525588e+06 ,  1.25851649e+04],
                    [-7.65181989e+04,   1.32438370e+06,   3.71925825e+05 , -5.16878541e+07  , -2.91740427e+07 , -8.54525588e+06 ,  9.88949015e+06 , -1.01838183e+04],
                    [-2.23806234e+03,  -4.51559188e+02,  -1.26811078e+02 ,  3.20338905e+04  ,  1.88438127e+04 ,  1.25851649e+04 , -1.01838183e+04 ,  1.51709659e+04]])

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
    M11 = dot(M, prior_fish)
    Final_prior_Fisher = dot(M11 , MT)
    #===== Stack new columns and raws===========
    newraw = linspace(0.,0.,n)
    matrix = vstack((matrix,newraw))
    matrix = vstack((newraw,matrix))
    newcolumn = linspace(0., 0., n+2)
    matrix = column_stack(( matrix, newcolumn))
    matrix = column_stack((newcolumn, matrix))
    matrix_plus_prior = matrix +  Final_prior_Fisher
    matrix_plus_prior =  linalg.inv(matrix_plus_prior)
    return matrix_plus_prior

def FoM(matrix_plus_prior, sarea):

    '''
    :param matrix_plus_prior: The fisher matrix that you would like
    find FoM of w and w0 for it,
    :return:FoM
    '''
    #============parameters========================
    w0 = sqrt(matrix_plus_prior[1,1])
    wa = sqrt(matrix_plus_prior[2,2]) 
    w0a = (matrix_plus_prior[2,1])  
    wa0 =((matrix_plus_prior[1,2]))  
    FoM2 =  1.0/sqrt(matrix_plus_prior[1,1] * matrix_plus_prior[2,2]
                     - matrix_plus_prior[1,2]* matrix_plus_prior[2,1])
    print sarea , FoM2
    return  sarea, FoM2




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
	sarea_print = [100, 250, 500, 750, 1000, 2000, 2500, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 12000, 15000, 17000, 20000, 22000, 25000,27000, 30000]
	save_area= []	; save_area_800_1300 = []
	save_fom = []	; save_fom_800_1300= []	
	for i in range(len(sarea_print)):
		A =  loadtxt('Fisher_diff_area_ska1_'+str(sarea_print[i])+'.txt', unpack=True)	
		ska_refuJy  = add_cmb(M, Planks_prior, A, 6) 
		ska_refuJy =(FoM(ska_refuJy, sarea_print[i]))
		save_area.append(ska_refuJy[0]) ; save_fom.append(ska_refuJy[1])
		save_area.append(ska_refuJy[0]) ; save_fom.append(ska_refuJy[1] )
		B =  loadtxt('Fisher_diff_area_ska1_800_1300_'+str(sarea_print[i])+'.txt', unpack=True)	
		ska_refuJy_800_1300 = add_cmb(M, Planks_prior, B, 6) 
		ska_refuJy_800_1300 =FoM(ska_refuJy_800_1300, sarea_print[i])
		save_area_800_1300.append(ska_refuJy_800_1300[0]) ; save_fom_800_1300.append(ska_refuJy_800_1300[1] )	
	#P.plot(save_area, save_fom, 'red', linewidth=1.5, label = r'$\nu=950-1670  \ MHz$')
	P.plot(save_area, save_fom,  linewidth=2.0, marker= 'o', color='purple', markeredgecolor='purple', label = r'$\nu=950-1670  \ MHz$')
	P.plot(save_area_800_1300, save_fom_800_1300, marker= 'o',color ="orange" ,markeredgecolor='orange', linewidth=2.0, label=r'$\nu =800-1300 \  MHz$')
	P.legend(loc="lower right")
	P.ylabel("FoM", fontsize= 18)
	P.xlabel(r"$S_{{\rm area} }  \ [{\rm deg}^2]$", fontsize= 18)
	#P.title(r"$\rm{with \ Planck's \ prior}$")
	P.savefig("area_vs_FoM_with_prior_2feq.pdf")
	P.show()
