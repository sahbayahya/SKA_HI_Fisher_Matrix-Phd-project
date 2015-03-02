# Compare "optimistic Stage IV" constraints from different experiments
# TD = Time Delays
# SN = Supernovae
# WL = Weak Lensing
# BAO = Baryon Acoustic Oscillations
# CL = Cluster counts
# prior = Planck + Stage II WL+SN+CL

# TD constraints calculated by CoeMoustakas09
# All others calculated by Dark Energy Task Force (Albrecht06)

# Ellipses plotted in order of increasing area

# RUN FROM COMMAND LINE (EXAMPLES):
"""
python plotexps.py h w0
python plotexps.py h Ok
python plotexps.py w0 wa
python plotexps.py w0 Ok
python plotexps.py w0 Ode

python plotexps.py h w0 -save
python plotexps.py h Ok -save
python plotexps.py w0 wa -save
python plotexps.py w0 Ok -save

python plotexps.py w0 wa -H
"""

# ~/cosmo/DETFast/plots/8/
# plotexps.py

# ~/cosmo/DETFast/plots/7/
#
# DETFtest.py
# CLtest.py
# plotexpsw.py
# pivot.py
# plotexps.py

# plotexps.py
# ~/cosmo/DETFast/plots/
# plotexps.py
# plothw0.py
# ~/glens/h0limits/cosmo/3/ensemblecosmo.py
# plothw.py
# plotell.py
# transform.py


from coetools import *
from CLtools import *
import Fisher
reload(Fisher)

verbose = 0

#http://www.scipy.org/Cookbook/Matplotlib/LaTeX_Examples
# ~/p/pylab.rcParams.txt
import pylab
pparams = {'axes.labelsize': 26,
           'font.size': 20,
           'figure.subplot.left': 0.2,
           'figure.subplot.bottom': 0.125,
           'figure.subplot.top': 0.95,
           'figure.subplot.right': 0.95,
           }
pylab.rcParams.update(pparams)

#for param in pylab.rcParams:
#   print param, pylab.rcParams[param]

concordance = (0.3,  0.7,  0.7)  # Omega_m, Omega_L, h
cosmoWMAP5  = (0.28, 0.72, 0.7)  # Omega_m, Omega_L, h
#cosmoDETF1   = (0.278, 0.722, 0.725)  # Omega_m, Omega_L, h
# USED BY DETF TO CALCULATE hstprior.fisher:
cosmoDETF   = (0.27, 0.73, 0.72)  # Omega_m, Omega_L, h

#Om, OL, h = cosmoDETF
#Om, OL, h = concordance
# CHANGED BELOW FOR PLOTTING!
# BUT SHOULD HAVE THESE VALUES FOR TRANSFORMATIONS
#Om, OL, h = 0.27, 0.73, 0.72  # USED BY DETF TO CALCULATE hstprior.fisher
Om, OL, h = cosmoDETF  # USED BY DETF TO CALCULATE hstprior.fisher
Ode = OL
Ok = 0
om = Om * h**2
w0, wa = -1, 0

#################################

sh = 1  # show, else save
pcl = params_cl()
if 'save' in pcl:
    sh = 0

exp0 = 'TD'
if 'H' in pcl:
    exp0 = 'H'

dh = 0.009
#dh = 1e-6
#TD dh = 0.00815
#H  dh = 0.00794

"""
# exp1  zp  FOM
prior 0.263935455466 1.0
H 0.285872615922 1.24335811694
H 0.293604054969 1.35308463322  (dh = 1e-6)
TD 0.308266404742 1.66887692437
CL 0.360709056915 2.54967733899
SN 0.266895885608 3.54468888502
BAO 0.443408733715 3.57875198533
WL 0.474327712049 8.10912014493
"""

legroom = 1  # make room for legend
left = 1  # legend placement: left, else right

vars = 'h Ok'
vars = 'h w0'
vars = 'w0 Ok'
vars = 'w0 wa'
#vars = 'Ode w0'
#vars = 'h wa'
#vars = 'h Ode'
#vars = 'Ode Ok'
#vars = 'Ok wa'

#print sys.argv, len(sys.argv)

if len(sys.argv) > 1:
    prog, xvar, yvar = sys.argv[:3]
else:
    xvar, yvar = strspl(vars)

exec('xo = %s' % xvar)
exec('yo = %s' % yvar)

vars = '%s %s' % (xvar, yvar)

if vars in ['w0 wa']:  # 'h Ok'
    legroom = 0  # don't make room for legend

if vars in ['w0 wa']:
    left = 0  # legend placement: right

#dTC = 0.022  # nominal, adjust using sigfac
#sigfac = 1
N = 1000
#N = 2000
#N = 2400
#N = 1000000
sigfac = sqrt(N / 100.)  # MULTIPLY F
# w0-wa: 10,000 lenses does well

stage = 4

Planck = False
Planck = True

HSTKey = False
#HSTKey = True
#dhkey = 0.08

Stage2 = False
Stage2 = True
Stage2fac = 1

SN4prior = 0

p1sig = 0

fixrest = False

fixes = None
#fixes = 'Ok wa'
#fixes = 'Ok'
#fixes = 'wa'

priors = {}
#priors['Ode'] = 0.3
#priors['Ok'] = 0.3
#priors['w0'] = 0.5
#priors['wa'] = 1

#priors['Ode'] = 0.1
#priors['Ok'] = 0.1
#priors['w0'] = 0.2
#priors['wa'] = 0.8

#priors['Ok'] = 0.001
#priors['wa'] = 0.001

alpha = 0.9

#################################
# Transformation
# om = Om * h**2
# Om + Ode + Ok = 1
# om = (1 - Ode - Ok) * h**2

# d om / d h  = 2 Om h
# d om / d Ode = -h**2
# d om / d Ok = -h**2

oldparams = strspl('om Ode Ok w0 wa ob ns lnP')
newparams = strspl('h  Ode Ok w0 wa ob ns lnP')
npar = len(oldparams)

M = identity(npar)
M[0,0] = 2 * Om * h
M[0,1] = -h**2
M[0,2] = -h**2

TDparams = strspl('h  Ode Ok w0 wa')

#################################
# Priors

Fisher.silentall = True

FP = Fisher.Fisher('Planck')
FH = Fisher.Fisher('hstprior')

# Stage II WL+SN+CL
WL2 = Fisher.Fisher('WL2')
SN2 = Fisher.Fisher('SN2')
CL2 = Fisher.Fisher('CL2')
F2 = WL2 + SN2 + CL2
#print F2.params

F2 = F2 * Stage2fac
#print F2.params

SN4 = Fisher.Fisher('SN4')

###

Fprior = FP
if not Planck:
    Fprior.data = Fprior.data * 0
if HSTKey:
    Fprior = Fprior + FH
if Stage2:
    Fprior = Fprior + F2
if SN4prior:
    Fprior = Fprior + SN4

Fprior.rename()
Fprior.reorder(oldparams)
Fprior.transform(newparams, M)

#################################
# LEGEND

def addlab(lab, colors=[(0,1,1), (0,0,1)], alpha=alpha, sx=0.05, sy=0.05, zorder=10, sh=0):
    global yl
    addell(xl, yl, colors, alpha, sxl, syl, zorder=zorder, sh=0)
    text(xl+dxl, yl, lab, va='center')
    yl += dyl
    if sh: show()

zorder = 1

#################################

def fishload(exp1):
    if exp1 == 'TD':
        indir = fishdir  # defined in Fisher.py
        Fcat = loadcat('TD.fisher', dir=indir, silent=(not verbose))
        F = Fisher.Fisher()
        F.data = Fcat.data
        F.params = Fcat.labels
    elif 0:
        indir = '~/glens/h0limits/cosmo/4/'
        FF = loadfits('allcosmo_ensemble_N100', dir=indir)  # Fisherensemble.py: dTC = 0.022
        FF = FF + 0  # doesn't like being loaded in from fits! This fixes it. (Fthis.py)
        FF = FF * sigfac**2
        F = Fisher.Fisher()
        F.data = FF
        F.params = TDparams
    elif exp1 == 'prior':
        F = Fisher.Fisher()
        F.data = zeros((npar,npar))
        F.params = newparams
    elif exp1 == 'H':
        F = Fisher.Fisher()
        F.data = zeros((npar,npar))
        F.params = newparams
        F.prior('h', dh)
        #F.prior('h', 0.022)
        #F.prior('h', 0.012)
        #print F.data
    else:
        F = Fisher.Fisher(exp1+`stage`)  # Stage III / IV
        F.rename()
        F.reorder(oldparams)
        F.transform(newparams, M)
    #print exp1, F.dx('h')
    #print exp1, F.dxdyp(xvar, yvar)
    return F

def expplot(exp1, alpha=alpha, leg=1):
    global yl, zorder
    colors = expcolors[exp1]
    F = Fdict[exp1]
    dx, dy, p = F.dxdyp(xvar, yvar)
    #print '---------------------------------'
    #print exp1, dx, dy, p
    A = dx * dy * sqrt(1 - p**2) # * pi
    #print 'FOM = ', 1/A
    if p1sig:
        colors = array(colors)
        color = (colors[0] + colors[1]) / 2.
        plotell1p(xo, yo, dx, dy, p, color=color, alpha=alpha, lw=5)
    else:
        plotellsp(xo, yo, dx, dy, p, colors=colors, alpha=alpha)
    exp1 = string.split(exp1, '3')[0]
    if leg: addlab(exp1, colors, alpha)
    zorder += 1

def Fcalc(exp1):
    F = fishload(exp1)
    F = F + Fprior
    for param in priors.keys():
        sig = priors[param]
        F.prior(param, sig)
    F.reorder(newparams)
    F.fix(fixes)
    return F

exps = string.split('WL BAO SN CL %s prior' % exp0)
#exps = string.split('TD H prior')
#exps = string.split('WL BAO SN CL')
#exps = string.split('WL BAO SN CL TD prior')
#exps = string.split('WL BAO SN CL TD H prior')
#exps = string.split('WL BAO SN CL prior')

if SN4prior and 'SN' in exps:
    exps.remove('SN')

#print Fprior.params
#pint(Fprior.data)

Fdict = {}
dds = []
dxmax = 0
dymax = 0
for exp1 in exps:
    F = Fdict[exp1] = Fcalc(exp1)
    #print exp1, xvar, yvar
    #print F.params
    #pint(F.data)
    dx, dy, p = F.dxdyp(xvar, yvar)
    if not isnan(dx) and not isnan(dy):
        dxmax = max([dx, dxmax])
        dymax = max([dy, dymax])

# RESTORE CONCORDANCE COSMO AFTER DETF VARIABLE TRANSFORMATION!
Om, Ode, h = concordance
exec('xo = %s' % xvar)
exec('yo = %s' % yvar)

# AUTOMATICALLY DETERMINE PLOT RANGES
ranges = {}
rfac = 1.1
if p1sig:
    rdx = dxmax * rfac * nsig1
    rdy = dymax * rfac * nsig1
else:
    rdx = dxmax * rfac * nsig2
    rdy = dymax * rfac * nsig2

ranges[xvar] = xo - rdx, xo + rdx
ranges[yvar] = yo - rdy, yo + rdy

# LEGEND
sxl = 0.06
syl = 0.05
dxl = sxl * 0.8
dyl = syl * 1.2
if left:
    xl, yl = 0.07, 0.95 - 5 * dyl  # Upper left
else:
    xl, yl = 0.8, 0.95 - 5 * dyl  # Upper right

if legroom:  # make room for legend
    lo, hi = ranges[xvar]
    d = hi - lo
    d = 1.2 * d
    lo = hi - d
    ranges[xvar] = lo, hi

xr = ranges[xvar]
yr = ranges[yvar]

#xl, yl = 0.75, 0.65
#xl, yl = 0.8, 0.95 - 4 * dyl
#
xl = interp([xl], array([0, 1]), array(xr))[0]
yl = interp([yl], array([0, 1]), array(yr))[0]
dxl = dxl * p2p(xr)
dyl = dyl * p2p(yr)
sxl = sxl * p2p(xr)
syl = syl * p2p(yr)

for exp1 in exps:
    F = Fdict[exp1]
    dx, dy, p = F.dxdyp(xvar, yvar)
    #print exp1, dx, dy
    if verbose >= 2:
        print exp1
        for var in newparams:
            print '%.5f' % F.dx(var), var
            # print '%3s' % var, F.dx(var)

        print

    #dx = dx / p2p(ranges[xvar])
    #dy = dy / p2p(ranges[yvar])

    dd = dx * dy
    dd = dd * sqrt(1 - abs(p))

    dd = dx * dy * sqrt(1 - p**2)
    
    A, B, ang = setell(dx, dy, p)
    #dd = B

    dds.append(dd)

SI = argsort(dds)[::-1]
exps = take(exps, SI)
#exps = string.split('WL BAO SN CL TD prior')
#exps = exps[::-1]

FOMdict = {}
zpdict = {}

if vars == 'w0 wa':
    if verbose: print '---------------------------------'
    for exp1 in exps:
        F = Fdict[exp1]
        dx, dy, p = F.dxdyp(xvar, yvar)
        zp = -1 / (1 + dy / (p * dx))
        if verbose:
            print '%6s' % exp1,
            print 'Pivot redshift = ', zp
        zpdict[exp1] = zp

    if verbose: print '---------------------------------'
    for exp1 in exps:
        F = Fdict[exp1]
        dx, dy, p = F.dxdyp(xvar, yvar)
        A = dx * dy * sqrt(1 - p**2) # * pi
        FOM = 1 / A
        FOMdict[exp1] = FOM
        if verbose:
            print '%6s' % exp1,
            print 'FOM = ', FOM / FOMdict['prior']
    if verbose: print '---------------------------------'
    
    #print '# exp1  zp  FOM'
    #for exp1 in exps:
    #    print exp1, zpdict[exp1], FOMdict[exp1] / FOMdict['prior']


expcolors = {}
expcolors['TD'] = blues
expcolors['H'] = lightblues
expcolors['CL'] = reds
expcolors['SN'] = yellows
expcolors['BAO'] = purples
expcolors['WL'] = greens
expcolors['prior'] = lightgreys

clf()
ioff()
for exp1 in exps:
    expplot(exp1)

xlim(ranges[xvar])
ylim(ranges[yvar])

#xlim(-1.5, -0.65)
#ylim(-0.03, 0.03)

if 0:
    xlim(-2, 0)
    ylim(-0.1, 0.1)
    outroot += '_zoomout'

finishup(xo, yo, xvar, yvar, c='k', dc='w', sh=0)

if fixrest:
    print "*** FIXING omega_b n_s lnP (in exp1, not Planck) ***"

outroot = 'Stage' + roman(stage)

if Planck:
    outroot += '_Planck'

if Stage2:
    outroot += '_Stage2'
    if Stage2fac <> 1:
        outroot += 'x%g' % Stage2fac

# wasn't always here:
if HSTKey:
    outroot += '_HST'

if SN4prior:
    outroot += '_SN4'

outroot += '_%s-%s' % (xvar, yvar)

if exp0 == 'H':
    outroot += '_H'

if fixes:
    if 'wa' in strspl(fixes):
        outroot += '_wconst'
    if 'Ok' in strspl(fixes):
        outroot += '_flat'

#subplots_adjust(left=0.15)  # Adjust plot margin (3/h.py hcosmo.png)

if sh:
    print 'Would have saved', outroot
    show()
else:
    print 'Saving', outroot
    savepngpdf(outroot)


#savepngpdf('StageIV_h-w_var')
#savepngpdf('StageIV_h-w_const')
#savepngpdf('StageIV_h-w_flat')
#savepngpdf('StageIV_h-w_const_flat')

#outroot = 'Stage%s_%s-%s' % (roman(stage), xvar, yvar)

