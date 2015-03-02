from CLtools import *
import Fisher
from numpy import *
from scipy import linalg
from matplotlib import *
from pylab import *
reload(Fisher)

#==============# Functions===========================


def expplot(exp1, dx, dy, p, alpha, fill, lw):
    F = Fisher.Fisher(exp1+`stage`)  # Stage III / IV
    F.rename()  # parameters renamed to nicknames (w0, wa, etc.)
    print dx, dy , p
    plotellsp(xo, yo, dx, dy, p, colors=expcolors[exp1], alpha=alpha, fill=fill, lw=lw)
    
    
def FoM(dx, dy, dxy, Delta_x):
	part1 = (dx**2 + dy**2)/ 2.0
	part2 =sqrt( (((dx**2 - dy**2 )**2)/ 4.0) + dxy**2)
	a = abs(part1 + part2)
	b = abs(part1 -  part2)
	FoM =pi/( pi * Delta_x * sqrt(a)* sqrt(b) )
	return FoM

#=============# Color selection========================
expcolors = {}
expcolors['WL'] = lightyellows
expcolors['SN'] =reds #greens
expcolors['BAO'] =lightblues
stage = 4
#===========#  Variables =============================

xvar, yvar = 'w0 wa'.split(' ')
xo, yo = -1, 0  # Best fit w0, wa

#==================================================
# The Values of w_0 and w_a for a ll the sensitivites, 7.3 mJy For the SKA 
#==================================================

Delta_x1 = 1.0
Delta_x2= 2.48
dx_150 =0.414272922179 ; dy_150 =1.75683373278; dxdy_150 =  -0.719095072654; p_150 = dxdy_150/(dx_150 * dy_150)  #dxdy_100/dx_ * dy_7.3 #  100 mJy_diff_bins
dx_5 =0.0396498    ; dy_5 = 0.1337945 ; dxdy_5 = -0.0047904059692 ; p_5 =  dxdy_5/(dx_5 * dy_5) # Values of w_0 , w_a and w_0w_a for  7.3 mJy_diff_bins
dx_Euclid= 0.1135 ; dy_Euclid = 0.2991 ; dxdy_Euclid =  -0.032599781664 ; p_Euclid= dxdy_Euclid/(dx_Euclid* dy_Euclid)  # Euclid 14_diff_bins

#================================================
# Calculate the Figure of Merit as in Coe. arXiv: 0906.4123v1
#================================================
print
print "=========================================="
print "Figure of Merit using Coe DETF"
print "==========================================="
print 
fom_1mJy = FoM(dx_150, dy_150, dxdy_150, Delta_x1) ; print 'FoM  Coe DETF + SKA (100mJy) = ', fom_1mJy
fom_5 = FoM(dx_5, dy_5, dxdy_5, Delta_x1); print 'FoM coe DETF  + (7.3 mJySKA) = ', fom_5
fom_Euclid = FoM(dx_Euclid, dy_Euclid, dxdy_Euclid,  Delta_x1) ; print 'FoM  Coe DETF +Euclid = ', fom_Euclid
print "=================Thanks!=========================="
print 
#====================# Plotting ====================
from matplotlib.pyplot import *
ax = subplot(1,1,1)
exp1= 'BAO'
p02 = expplot(exp1, dx_150, dy_150, p_150, alpha= 1, fill = 1, lw = 1)
exp3 = 'WL'
p03=  expplot(exp3, dx_Euclid, dy_Euclid, p_Euclid, alpha= 1, fill = 1, lw = 1)
exp2 = 'SN'
p01 =expplot(exp2, dx_5, dy_5, p_5, alpha= 1, fill = 1, lw = 1)
#xlim(-1.7,-0.3)
ylabel('$w_a$', fontsize= 18)
xlabel('$w_0$', fontsize= 18)
#ymax = 4.00
#ylim(-ymax,ymax)
ymax =0.5
ylim(-ymax,ymax)
xlim(-1.2, -0.7)
plot(xo, yo, 'kx')
plot(xo, yo, color = 'k', linestyle= '-',linewidth=1.5)
plot(xo, yo, color = 'k', linestyle= '-',linewidth=1.5)
savefig("output_ellipse_w0_wa.pdf")
print 'You should see a plot which with ellipses of different colors.'
print '(Looks like a big X at this scale.)'
print 'KILL PLOT WINDOW TO TERMINATE'
show()
