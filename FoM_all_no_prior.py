from scipy import *
from numpy import *
import pylab as P
from scipy import linalg
from functions import *

#============================Main progarm ===========================================================================
sarea_print = [100, 250, 500, 750, 1000, 2000, 2500, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 12000, 15000, 17000, 20000, 22000,25000, 27000 ,30000]
save_area= []	; save_area_800_1300 = []
save_fom = []	; save_fom_800_1300= []							
for i in range(len(sarea_print)):
	A =  loadtxt('diff_survey_area_SKA1_outputs/Fisher_diff_area_ska1_'+str(sarea_print[i])+'.txt', unpack=True)
	ska_refuJy = FoM_no_prior(A, sarea_print[i])
	save_area.append(ska_refuJy[0]) ; save_fom.append(ska_refuJy[1])
for i in range(len(sarea_print)):	
	B =  loadtxt('diff_survey_area_SKA1_outputs/Fisher_diff_area_ska1_800_1300_'+str(sarea_print[i])+'.txt', unpack=True)
	ska_refuJy_800_1300 = FoM_no_prior(B, sarea_print[i])
	save_area_800_1300.append(ska_refuJy_800_1300[0]) ; save_fom_800_1300.append(ska_refuJy_800_1300[1] )
#print len(save_area), len(save_fom)
P.plot(save_area, save_fom,  linewidth=2.0, marker= 'o', color='purple', markeredgecolor='purple', label = r'$\nu=950-1670  \ MHz$')
P.plot(save_area_800_1300, save_fom_800_1300, marker= 'o',color ="orange" ,markeredgecolor='orange', linewidth=2.0, label=r'$\nu =800-1300 \  MHz$')

P.ylim(0.0, 0.2)
P.legend()
P.ylabel("FoM")
P.xlabel(r"$S_{area}$")
P.title("No prior")
P.savefig("plots/area_vs_FoM_no_prior_2freq.pdf")
P.show()

