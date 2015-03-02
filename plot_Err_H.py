from pylab import *
import matplotlib.ticker as tic
plot_these_Srms = [ 3, 5, 23, 70, 150, 200]
z = [] ; Da =[] ; hz =[]
(z_E, Da_E, H_E,R_E,B_E) = loadtxt('Outputs/output_Fisher_bao_Euclid_s15k_3.txt', unpack=True)
colors_err = [ '#a40000', '#cc0000', '#ef2929', '#204a87', '#3465a4', '#729fcf']
label_err = [r"${\rm SKA2} \ {\rm Opt}$", r"${\rm SKA2}$"
    , r"${\rm SKA2} \ {\rm Pess}$", r"${\rm SKA1}\ {\rm Opt}$", r"${\rm SKA1}$"
    ,  r"${\rm SKA1} \ {\rm Pess}$" ]
(z_BOSS, Da_BOSS, hz_BOSS) = loadtxt("/Users/sahba/Dropbox/SKA_Survey/SKA_HI_Fisher_Matrix/Outputs/output_BOSS_S10000.txt", unpack=True)
p_ = []
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
for i in range(len(plot_these_Srms)):
	(z_, Da_, hz_) = (loadtxt("/Users/sahba/Dropbox/SKA_Survey/SKA_HI_Fisher_Matrix/Outputs/output_Srms_ska1_"+str(plot_these_Srms[i])+".txt", unpack=True))
	z.append(z_); Da.append(Da_) ; hz.append(hz_)
	ax.plot(z[i], hz[i], color = colors_err[i],linewidth=1.5, linestyle="-", label=label_err[i])
	ax.scatter(z[i],hz[i], s= 50, marker='o',  edgecolor = colors_err[i], facecolor= colors_err[i])
p_E, = ax.plot(z_E,  H_E , color = '#edd400',linewidth=1.5, linestyle="-")
p_BOSS, = ax.plot(z_BOSS, hz_BOSS, color= "#6E6E6E",linewidth =1.5,linestyle="-")
ax.scatter(z_E, H_E, s= 50, marker='o',  edgecolor = '#edd400', facecolor= '#edd400')
ax.scatter(z_BOSS, hz_BOSS, s=50, marker='o', edgecolor= "#6E6E6E", facecolor="#6E6E6E")



l2 =legend([p_E],[r"${\rm Euclid} \ {\rm }$", r"${\rm BOSS}$"], loc=4, frameon = False)
legend(p_,label_err,loc=1, frameon= False)
ax.legend(loc='upper right', ncol=2, frameon=False)
ax.set_yscale('log')
plt.gca().add_artist(l2)


#======= x axis
xlim(0., 2.1)
ax.set_xlabel(r"$ { \rm redshift} (z)$", fontsize=20)
xticks = arange(0, 3.1, 0.3)
#======= y axis
ax.get_yaxis().set_major_formatter(tic.ScalarFormatter())
ax.yaxis.set_major_formatter(tic.FormatStrFormatter('%0.1f'))
ax.set_yscale('log')
ax.set_ylabel(r"$\sigma_{H}/H\%$", fontsize=20)
yticks = [0.1, 0.2, 0.5, 1 ,2, 5, 10,20,25]
plt.yticks(yticks,[r'$0.1$',r'$0.2$', r'$0.5$',r'$1$',r'$2$',  r'$5$',r'$10$',r'$20$'])
plt.ylim(0.1,  25)
plt.savefig('plots/output_lnH_mario_bias_corrected_nz.pdf')
show()
