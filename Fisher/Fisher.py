# MANIPULATING FISHER MATRICES
# SEE Dark Energy Task Force (Albrecht06) Technical Appendix

"""
Dan Coe (2009)
http://arxiv.org/abs/0906.4123
Fisher Matrices and Confidence Ellipses: A Quick Start Guide and Software

used in Coe & Moustakas (2009):
http://arxiv.org/abs/0906.4108
Cosmological Constraints from Gravitational Lens Time Delays
"""

#from pylab import *  # inv
#from cosmocoe import *
#from useful import *
#from os.path import join, exists
#from coeplot2 import *
from coetools import *
from numpy.linalg.linalg import *  # inv (Matrix inverse)
from os.path import join, exists
import os

home = os.environ.get('HOME', '')
#fishdir = join(home, '/Dropbox/SKA Survey/Fishers/Fisher/cosmo/DETFast/data')  # PATH FOR FISHER MATRICES
fishdir = join(home, '/Dropbox/SKA Survey/Fishers/Fisher/data')  # PATH FOR FISHER MATRICES
fishpath = os.environ.get('FISHERPATH', '')  # setenv FISHERPATH /Users/coe/Fisher
fishdir = join(fishpath, 'data')  # /Users/coe/Fisher/data/
Fishdict = loaddict('Fisher.dict', dir=fishdir, silent=True)  # File nickanmes (exp1  inroot)
pdict = loaddict('params.dict', dir=fishdir, silent=True)  # Parameter nicknames (Omega_Q  Ode)

#DETFdir = join(home, 'cosmo/DETFast')
#fishdir = join(DETFdir, 'data')
#Fishdict = loaddict('Fisher.dict', dir=join(DETFdir, 'plots'), silent=True)  # File nickanmes (exp1  inroot)
#pdict = loaddict('params.dict', dir=join(DETFdir, 'plots'), silent=True)  # Parameter nicknames (Omega_Q  Ode)

silentall = False

def matrix_multiply(MM):
    """Multiplies a list of matrices: M[0] * M[1] * M[2]..."""
    P = MM[0]
    for M in MM[1:]:
        P = dot(P, M)
    return P

class Fisher:
    def __init__(self, inroot='', xvar='', yvar='', fixes=[], margs=[],
                 data=[], params=[], silent=False):
        self.inroot = Fishdict.get(inroot, inroot)
        self.xvar = xvar
        self.yvar = yvar
        self.fixes = fixes
        self.margs = margs
        self.data = data
        self.params = params
        self.silent = silent or silentall
        if self.inroot: self.load()

    def load(self):  # DETFast format
        txt = loadfile(self.inroot+'.fisher', dir=fishdir, silent=self.silent)
        nparam = string.atoi(string.split(txt[0])[1])
        self.params = []
        for ivar in range(nparam):
            param = string.split(txt[ivar+4])[0]
            self.params.append(param)
                
        self.data = loaddata(self.inroot+'.fisher+', dir=fishdir, headlines=4+nparam, silent=1)

    def ii(self):
        self.ix = None
        self.iy = None
        if self.xvar:
            self.ix = self.params.index(self.xvar)
        if self.yvar:
            self.iy = self.params.index(self.yvar)
        
        self.ifixes = []
        for i, param in enumerate(self.params):
            if param in self.fixes:
                self.ifixes.append(i)
        
        self.imargs = []
        for i, param in enumerate(self.params):
            if param in self.margs:
                self.imargs.append(i)

    def pindex(self, param):
        return self.params.index(param)

    def take(self, iparams):
        self.data = self.data.take(iparams, 0)
        self.data = self.data.take(iparams, 1)
        self.params = list(take(self.params, iparams))
        
    def reorder(self, params):
        """Matrix with params in new order"""
        self.repar(params)
        iparams = map(self.pindex, params)
        self.take(iparams)
        
    def rename(self, pdict1=None):
        """Rename parameters given a dictionary of names & nicknames"""
        pdict1 = pdict1 or pdict
        for i, param in enumerate(self.params):
            self.params[i] = pdict1.get(param, param)

    def fix(self, fixes=[]):
        """Fix parameters constant <==> Remove them from the Fisher matrix"""
        self.fixes = fixes or self.fixes
        self.fixes = strspl(self.fixes)
        
        self.ii()
        
        iall = arange(len(self.params))
        ikeep = set(iall) - set(self.ifixes)  # Sorts result
        ikeep = list(ikeep)
        
        self.take(ikeep)

    def marg(self, margs=[]):
        """Marginalize over variables: Remove them from the covariance matrix"""
        self.margs = margs or self.margs
        self.ii()
        
        #ikeep = invertselection(arange(len(self.params)), self.fixes)
        iall = arange(len(self.params))
        ikeep = set(iall) - set(self.imargs)  # Sorts result
        ikeep = list(ikeep)
        
        C = self.cov()
        C = C.take(ikeep, 0)
        C = C.take(ikeep, 1)
        self.data = inv(C)
        
        self.params = list(take(self.params, ikeep))

    def transform(self, params, M):
        """Transform to new set of parameters using matrix provided"""
        self.data = matrix_multiply([transpose(M), self.data, M])
        self.params = params

    def prior(self, param, sig):
        """Set a prior of sig on param"""
        ivar = self.pindex(param)
        self.data[ivar,ivar] = self.data[ivar,ivar] + 1 / sig**2

    def release(self, param):
        """Release constraints on a parameter"""
        ivar = self.pindex(param)
        print "I don't think release works in Fisher.py"
        print
        xxx[9] = 3
        #self.data[ivar,ivar] = 1e-6
        #self.data[ivar,:] = 0
        #self.data[:,ivar] = 0

    def cov(self):
        """Covariance matrix"""
        return inv(self.data)

    def dxdyp(self, xvar='', yvar=''):  # , fixes=None
        """Return uncertainty in two parameters and their correlation"""
        self.xvar = xvar or self.xvar
        self.yvar = yvar or self.yvar
        #self.fixes = strspl(fixes or self.fixes)
        self.ii()
        
        C = self.cov()
        C = C.take((self.ix,self.iy),0)
        C = C.take((self.ix,self.iy),1)
        dx = sqrt(C[0,0])
        dy = sqrt(C[1,1])
        dxy = C[0,1]
        p = dxy / (dx * dy)
        self.C = C
        return dx, dy, p

    def dx(self, xvar=''): # , fixes=None
        """Return uncertainty in parameter (if marginalizing over others)"""
        self.xvar = xvar or self.xvar
        #self.fixes = strspl(fixes or self.fixes)
        self.ii()
        
        #dx = 1 / sqrt(self.data[self.ix,self.ix])
        self.C = C = self.cov()
        dx = sqrt(C[self.ix,self.ix])
        return dx

    def addpar(self, param):
        npar = len(self.params)
        data = zeros((npar+1, npar+1))
        data[:npar,:npar] = self.data
        self.data = data
        self.params.append(param)

    def repar(self, params):
        for param in params:
            if param not in self.params:
                self.addpar(param)

    def pr(self):
        """Print contents"""
        print self.params
        pint(self.data)

    def __add__(self, sel2):
        """Add Fisher matrices"""
        pl1 = self.params
        pl2 = sel2.params
        n1 = len(pl1)
        n2 = len(pl2)
        F1 = self.data
        F2 = sel2.data
        if pl1 == pl2:
            pl = pl1
        else:
            params1 = set(pl1)
            params2 = set(pl2)
            params = params1 | params2
            pl = list(params)
        #pl = list(sort(pl))
        n = len(pl)
        
        # Put F1 in merged F format
        ii = []
        for p in pl1:
            i = pl.index(p) 
            ii.append(i)
        
        FF1 = zeros((n,n))
        for i1, i in enumerate(ii):
            FF1[i].put(ii, F1[i1])
        
        # Put F2 in merged F format
        ii = []
        for p in pl2:
            i = pl.index(p) 
            ii.append(i)
        
        FF2 = zeros((n,n))
        for i2, i in enumerate(ii):
            FF2[i].put(ii, F2[i2])

        # Add
        new = Fisher()
        new.data = FF1 + FF2
        new.params = pl
        new.xvar = self.xvar
        new.yvar = self.yvar
        new.fixes = self.fixes
        return new

    def __mul__(self, fac):
        """Multiply Fisher matrix by some factor"""
        new = Fisher()
        new.data = self.data * fac
        new.params = self.params
        new.xvar = self.xvar
        new.yvar = self.yvar
        new.fixes = self.fixes
        return new
