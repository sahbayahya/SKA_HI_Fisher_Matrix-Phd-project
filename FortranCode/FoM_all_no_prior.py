from scipy import *
from numpy import *
import pylab as P
from scipy import linalg


def FoM_no_prior(matrix_plus_prior, area):
    '''
    :param matrix_plus_prior: The fisher matrix that you would like
    find FoM of w and w0 for it,
    :return:FoM
    '''
    matrix_plus_prior =  linalg.inv(matrix_plus_prior)
    #============parameters========================
    w0 = sqrt(matrix_plus_prior[0,0])#/1e-1 #;   print 'sigma_w0 = ', w0
    wa = sqrt(matrix_plus_prior[1,1])   # ; print 'sigma_wa =', wa
    w0a = (matrix_plus_prior[1,0])  #;  print 'sigma w0a = ', w0a
    wa0 =((matrix_plus_prior[0,1]))   #; print' wa0 = ', wa0
    ob0 =sqrt(matrix_plus_prior[2,2])*0.049 #/1e-4   #print'sigma_ob0 =', ob0
    ok0 = sqrt((matrix_plus_prior[3,3])) #/1e-2  #print 'sigma_ok0=', ok0
    om0 = sqrt(matrix_plus_prior[4,4])#/1e-2  #print 'sigma_om0 = ', om0
    h = sqrt(matrix_plus_prior[5,5])#/1e-2    #print 'sigma_h = ',  h
    FoM2 =  1.0/sqrt(matrix_plus_prior[1,1] * matrix_plus_prior[0,0]
                     - matrix_plus_prior[1,0]* matrix_plus_prior[0,1])#/(pi*(sqrt(2.31)))
    print  area , FoM2  
    return area,  FoM2




#============================Main progarm ===========================================================================
if __name__== "__main__" :
	sarea_print = [100, 250, 500, 750, 1000, 2000, 2500, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 12000, 15000, 17000, 20000, 22000,25000, 27000 ,30000]
	save_area= []	; save_area_800_1300 = []
	save_fom = []	; save_fom_800_1300= []							
	print ' Sarea  ' , '  FoM  '
	for i in range(len(sarea_print)):
		A =  loadtxt('Fisher_diff_area_ska1_'+str(sarea_print[i])+'.txt', unpack=True)	
		ska_refuJy = FoM_no_prior(A, sarea_print[i])
		save_area.append(ska_refuJy[0]) ; save_fom.append(ska_refuJy[1] )
	print len(save_area), len(save_fom)	
	for i in range(len(sarea_print)):	
		B =  loadtxt('Fisher_diff_area_ska1_800_1300_'+str(sarea_print[i])+'.txt', unpack=True)	
		ska_refuJy_800_1300 = FoM_no_prior(B, sarea_print[i])
		save_area_800_1300.append(ska_refuJy_800_1300[0]) ; save_fom_800_1300.append(ska_refuJy_800_1300[1] )
	print len(save_area), len(save_fom)	
	P.plot(save_area, save_fom, 'red', linewidth=1.5, label = r'$\nu=950-1670  \ MHz$')
	P.plot(save_area, save_fom, 'ro', markeredgecolor='red')
	P.plot(save_area_800_1300, save_fom_800_1300, 'blue', linewidth=1.5,  label=r'$\nu =800-1300 \  MHz$')
	P.plot(save_area_800_1300, save_fom_800_1300, 'bo', markeredgecolor='blue')
	P.ylim(0.0, 0.2)
	P.legend()
	P.ylabel("FoM")
	P.xlabel(r"$S_{area}$")
	P.title("No prior")
	P.savefig("area_vs_FoM_no_prior_2freq.pdf")
	P.show()

