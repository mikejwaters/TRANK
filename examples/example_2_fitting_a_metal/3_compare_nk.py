
import sys
noshow = False
if len(sys.argv)>1:
	#print (sys.argv[1])
	if sys.argv[1] == 'noshow':
		noshow = True
		from matplotlib import use
		use('Agg')
		from matplotlib.pylab import ioff
		ioff()



data_directory = 'TRANK_nk_fit/'
from TRANK import functionize_nk_file


from matplotlib.pylab import *
from numpy import loadtxt, array
lamda_smooth = (loadtxt(data_directory+'fit_nk_fine.txt').T)[0]
lamda_points = (loadtxt(data_directory+'fit_nk.txt').T)[0]
fit_nk_f = functionize_nk_file(data_directory+'fit_nk.txt', skiprows = 0, kind = 'cubic')

## you can just load the file directly
#actual_nk_f = functionize_nk_file('Au-glass_10nm_30p_effective_nk.txt', skiprows = 0, kind = 'cubic')
# or use the setup file
from basic_setup import nk_f_list, layer_index_of_fit
actual_nk_f = nk_f_list[layer_index_of_fit]

from matplotlib.pylab import figure, gca, subplots_adjust, tight_layout, show, savefig
n_color = 'teal'
k_color = 'orange'


figure(figsize =(3.2,2.5), dpi = 220*2.0/3.0 )
nax = gca()
kax = nax.twinx()

nk_points =  fit_nk_f(lamda_points)
nk_smooth =  fit_nk_f(lamda_smooth)


nax.plot(lamda_points, nk_points.real, linewidth = 1.0, color = 'k',   linestyle = '', marker = 'o', markersize = .8, zorder = 1)
nax.plot(lamda_smooth, nk_smooth.real,   linewidth = 1.0, color = n_color, label ='n', linestyle = '-', zorder = 0.9)
nax.plot(lamda_smooth,actual_nk_f(lamda_smooth).real,  linewidth = 2.0, color = n_color, label ='n', linestyle = '--', zorder = 1.1)

kax.plot(lamda_points, nk_points.imag, linewidth = 1.0, color = 'k',   linestyle = '',   marker = 'o', markersize = .8, zorder = 1)
kax.plot(lamda_smooth, nk_smooth.imag,   linewidth = 1.0, color = k_color, label ='k', linestyle = '-', zorder = 0.9)
kax.plot(lamda_smooth, actual_nk_f(lamda_smooth).imag,  linewidth = 2.0, color = k_color, label ='n', linestyle = '--', zorder = 1.1)


nax.minorticks_on()


kax.minorticks_on()
nax.tick_params(axis='y',which = 'both', colors=n_color)
kax.tick_params(axis='y',which = 'both', colors=k_color)
nax.tick_params(axis='x',which = 'both',labelsize=10)

nax.tick_params(axis='y',which = 'both',labelsize=10)
kax.tick_params(axis='y',which = 'both',labelsize=10)

nax.set_ylabel('n')
kax.set_ylabel('k')

nax.yaxis.label.set_color(n_color)
kax.yaxis.label.set_color(k_color)

nax.set_xlabel('Wavelength (nm)')

nax.spines['left'].set_color(n_color)
kax.spines['right'].set_color(k_color)
kax.spines['left'].set_visible(False)
nax.spines['right'].set_visible(False)


#nax.set_title(title_string , fontsize = 8)

subplots_adjust(top = 0.99,
		bottom = 0.18,
		left = 0.15,
		right = 0.84)




### fix limits
ylim =nax.get_ylim()
nax.set_ylim(0,ylim[1])
ylim =kax.get_ylim()
kax.set_ylim(0,ylim[1])


nax.set_xlim(min(lamda_smooth),max(lamda_smooth))
tight_layout(pad = 0.1)
savefig(data_directory +'nk_compared.pdf' , dpi = 600,transparent = True)






if noshow==False:
	show()
