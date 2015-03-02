from CLtools import *
import Fisher
from matplotlib import *
from pylab import *
reload(Fisher)

expcolors = {}
expcolors['CL'] = reds
expcolors['SN'] = greys
expcolors['BAO'] = blues
expcolors['WL'] = greens
stage = 4

xvar, yvar = 'w0 Om'.split(' ')
xo, yo =  0.277, -1  # Best fit w0, wa
dx = 0.0131914 ; dy = 0.02953; dxdy = -0.00261 ; p =   0.880401 #dxdy/dx * dy
#print p_test# 0 mJy_diff_bins
dx_cmb = 0.0099762; dy_cmb =  0.01983 ; dxdy_cmb = -0.00202 ; p_cmb =  0.880401 #dxdy/dx * dy
#print p_test# 0 mJy_diff_bins
dx_1 = 0.02152; dy_1 = 0.25421 ; dxdy_1 = -0.00261 ; p_1 = -0.880401#dxdy_1/dx_1 * dy_1 #  1 mJy_diff_bins
dx_2 = 0.02169; dy_2 = 0.25571 ; dxdy_2 = -0.00261 ; p_2 = -0.880401#dxdy_2/dx_2 * dy_2 #  7.3 mJy_diff_bins
dx_3 = 0.02441; dy_3 = 0.27755 ; dxdy_3 = -0.00215 ; p_3 = -0.880401#dxdy_3/dx_3 * dy_3 #  23 mJy_diff_bins
def expplot(exp1, dx, dy, p, alpha, fill, lw, zorder, colors):
    F = Fisher.Fisher(exp1+`stage`)  # Stage III / IV
    F.rename()  # parameters renamed to nicknames (w0, wa, etc.)
    #dx, dy, p= F.dxdyp(xvar, yvar)
    #dx = 0.02150; dy = 0.25400
    #dxdy = -0.00261
    #p  = -0.0244462310528 # -0.00261 / (dx * dy)
    print dx, dy , p
    plotellsp(xo, yo, dx, dy, p, colors=expcolors[exp1], alpha=alpha, fill=fill, lw=lw, zorder=zorder)

#clf()
#for exp1 in 'SN BAO CL WL'.split(' '):
#exp4= 'SN'
#p04 = expplot(exp4, dx_3, dy_3, p_3, alpha= 1, fill= 1)
#exp3= 'WL'
#p03 = expplot(exp3, dx_2, dy_2, p_2, alpha = 1, fill = 1, lw= 1)
#exp2  = 'CL'
#p02 = expplot(exp2, dx_1, dy_1, p_1, alpha = 1, fill = 1, lw= 1)
exp2 = 'CL'
exp1 = 'SN'
p02 = expplot(exp2, dx, dy, p, colors=expcolors[exp2] ,alpha= 1, fill = 1, lw = 1, zorder= 1)
p01 = expplot(exp1, dx_cmb, dy_cmb, p_cmb, colors= expcolors[exp1],alpha= 1, fill = 0, lw = 1,zorder= 1)
#legend([p01, p02, p03,p04],['$ 0\mu$Jy', '$1\mu$Jy', '$7.3\mu$Jy',' $23\mu$Jy'], loc='best')
xlim(0.23,0.32)
ylim(-1.1,-0.9)
plot(xo, yo, 'ko')
xlabel('$\Omega_m$', fontsize= 18)
ylabel('$w$', fontsize= 18)
savefig("output_ellipse_w0_om_2.eps")
#finishup(xo, yo, xvar, yvar, c='k', dc='w', sh=0)
print 'You should see a plot which with ellipses of different colors.'
print '(Looks like a big X at this scale.)'
print 'KILL PLOT WINDOW TO TERMINATE'
show()
#pause()
