import numpy as np
DEBUG_PLOT = True #True # Whether to plot fitting functions or not
NU_LINE = 1.420 # HI emxission line freq. in GHz
FULLSKY = (4.*np.pi * (180./np.pi)**2.) # deg^2 in the full sky
NGAL_MIN = 1e3 # Min. no. of galaxies to tolerate in a redshift bin
CBM = 1. #np.sqrt(1.57) # Correction factor due to effective beam for MID/MK (OBSOLETE)
CTH = 0.5 # Correction factor due to taking 5 sigma (not 10 sigma) cuts for SKA1
SBIG = 500. # Flux rms to extrapolate dn/dz out to (constrains behaviour at large Srms)
C = 3e5 # Speed of light, km/s
k_NL = 0.2
h= 0.67
    
name_rel = ['0n', 'SKA1MID-B1','SKA1MID_B2','SKA1MID_B2_Opt', 'SKA1MID_B2_Real', 'SKA1MID_B2_Pess', 'MEERKAT-B1', 'MEERKAT-B2', 'MID+MK-B1','MID_MK_B2_Real', 'SKA1SUR-B1', 'SKA1SUR-B2', 'ASKAP', 'SUR+ASKAP', 'SKA2_Opt', 'SKA2_Real','SKA2_Pess' ,'SKA2']
name = ['0', '315', '187', '70',  '150', '200', '696','750', '247', '152','174','192', '645',  '179',  '3', '7',  '23',  '5' ]
nameint=[0,   315  , 187,   70 ,   150,   200,   696,  750,   247,   152,  174 , 192,   645,    179,    3,   7,    23,    5  ]
fluxrms=[0.,  315., 187.,   70 ,   150.,  200.,  696., 750.,  247.,  152., 174., 192.,  645.,   179.,   3.,  7.3,  23.,   5.4]
numin = [500, 350., 950.,  800.,   800.  ,800.,  580., 900.,  580.,  950., 350., 650.,  700.,   700.,  500., 500., 500.,  500.]
numax = [1200,1050.,1760., 1300., 1300., 1300., 1015., 1670., 1015., 1670.,900., 1670., 1800., 1670., 1200., 1200.,1200., 1200.]
nucrit= [None, None, None, None, None, None, None, None, None, None, 710., 1300., 1250., 1300., None, None, None, None]
Sarea = [30e3, 5e3, 5e3, 5e3,5e3 ,5e3 , 5e3, 5e3, 5e3, 5e3, 5e3, 5e3, 5e3, 5e3, 30e3, 30e3, 30e3, 30e3]
Sconst= [True, False, False, False, False, False, False, False, False, False, False, False, False, False,True,True ,True,True]
Scorr = [1, CBM*CTH, CBM*CTH, CBM*CTH, CBM*CTH, CBM*CTH, CBM*CTH, CBM*CTH, CBM*CTH, CBM*CTH, CTH, CTH, CTH, CTH, 1.,1., 1. ,1.]

# Planck-only best-fit parameters, from Table 2 of Planck 2013 XVI.
cosmo = {
    'omega_M_0':        0.316,
    'omega_lambda_0':   0.684,
    'omega_b_0':        0.049,
    'omega_HI_0':       6.50e-4,
    'N_eff':            3.046,
    'h':                0.67,
    'ns':               0.962,
    'sigma_8':          0.834,
    'gamma':            0.55,
    'w0':               -1.,
    'wa':               0.,
    'fNL':              0.,
    'mnu':              0.,
    'k_piv':            0.05, # n_s
    'aperp':            1.,
    'apar':             1.,
    'bHI0':             0.702,
    'A':                1.,
    'sigma_nl':         7.,
    'beta_1':           0.,         # Scale-dependent bias (k^1 term coeff. [Mpc])
    'beta_2':           0.          # Scale-dependent bias (k^2 term coeff. [Mpc^2])
}
# Plancks prior matrix 
Planks_prior = np.array([[1.99579245e+05,  -3.73667528e+04, -1.04936812e+04 ,  1.39977603e+06  ,  5.58643962e+05 , -4.64225267e+04 , -7.65181989e+04 , -2.23806234e+03],
					  [-3.73667528e+04,   1.83928663e+05,   5.16525685e+04 , -7.42050738e+06  , -3.98758357e+06 , -1.11710442e+06 ,  1.32438370e+06 , -4.51559188e+02],
					  [-1.04936812e+04,   5.16525685e+04,   1.45055577e+04 , -2.08389634e+06  , -1.11983054e+06 , -3.13715719e+05 ,  3.71925825e+05 , -1.26811078e+02],
					  [1.39977603e+06,  -7.42050738e+06 ,  -2.08389634e+06 ,  3.64943809e+08  ,  1.58599621e+08 ,  4.25932543e+07 , -5.16878541e+07 ,  3.20338905e+04],
					  [5.58643962e+05,  -3.98758357e+06 ,  -1.11983054e+06 ,  1.58599621e+08  ,  8.70535526e+07 ,  2.48738854e+07 , -2.91740427e+07 ,  1.88438127e+04],
					  [-4.64225267e+04,  -1.11710442e+06,  -3.13715719e+05 ,  4.25932543e+07  ,  2.48738854e+07 ,  7.49686718e+06 , -8.54525588e+06 ,  1.25851649e+04],
					  [-7.65181989e+04,   1.32438370e+06,   3.71925825e+05 , -5.16878541e+07  , -2.91740427e+07 , -8.54525588e+06 ,  9.88949015e+06 , -1.01838183e+04],
					  [-2.23806234e+03,  -4.51559188e+02,  -1.26811078e+02 ,  3.20338905e+04  ,  1.88438127e+04 ,  1.25851649e+04 , -1.01838183e+04 ,  1.51709659e+04]])
ob = 0.049
M = np.array([[1., 0 , 0., 0., 0.,  0.,  0.,  0.],
		   [0., 1., 0., 0., 0.,  0.,  0.,  0.],
		   [0., 0., 1., 0., 0.,  0.,  0.,  0.],
		   [0., 0., 0., ob,-ob, -ob,  0.,  0.],
		   [0., 0., 0., 0., 1.,  0.,  0.,  0.],
		   [0., 0., 0.,-ob,-1., -1.,  0.,  0.],
		   [0., 0., 0., 0., 0., 0. ,  1.,  0.],
		   [0., 0., 0., 0., 0., 0. ,  0.,  1.]])