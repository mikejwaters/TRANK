


def compute_outer_log_levels(error_map, base10delta = 0.2):
	from numpy import floor, ceil, log10, logspace
	def bottom_log(x):
		return base10delta * floor( log10(x)/base10delta )

	def top_log(x):
		return base10delta * ceil( log10(x)/base10delta )

	low =  bottom_log( error_map.min() )
	high = top_log(    error_map.max() )

	nlevels = (high-low)/base10delta +1
	levels =  logspace(low, high ,nlevels  )

	return levels


def compute_inner_log_levels(error_map, base10delta = 1.0):
	from numpy import floor, ceil, log10, logspace
	def bottom_log(x):
		return base10delta * ceil( log10(x)/base10delta )

	def top_log(x):
		return base10delta * floor( log10(x)/base10delta )

	low =  bottom_log( error_map.min() )
	high = top_log(    error_map.max() )

	nlevels = (high-low)/base10delta +1
	levels =  logspace(low, high ,nlevels  )

	return levels



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




from TRANK import single_lamda_TRA_error_map, functionize_nk_file, try_mkdir, find_min_indices_2d_array
from numpy import  loadtxt, linspace

data_directory = 'TRANK_nk_fit/'
map_direct = 'TRANK_error_maps/'

from basic_setup  import  TR_pair_list_generator,   parameter_list_generator
fit_nk_f =  functionize_nk_file(data_directory+'fit_nk_fine.txt', skiprows = 0)
lamda_list = loadtxt(data_directory+'fit_nk.txt' , unpack = True, usecols = [0])
lamda_fine = loadtxt(data_directory+'fit_nk_fine.txt' , unpack = True, usecols = [0])



try_mkdir(map_direct )


from matplotlib.pylab import *
from matplotlib import ticker
from matplotlib.colors import LogNorm

fit_n_max = fit_nk_f(lamda_fine).real.max()
fit_k_max = fit_nk_f(lamda_fine).imag.max()

num_n_points = 50
num_k_points = 100

nmin, nmax = 0.01,  fit_n_max *2.0
kmin, kmax = 0.0,  fit_k_max *1.1

#dn = 0.01
#dk = 0.01
use_logscale = True

nlist = linspace(nmin, nmax, num_n_points)
klist = linspace(kmin, kmax, num_k_points)
#nlist = arange(nmin+dn, nmax+dn/2.0, dn)
#klist = arange(kmin, kmax+dk/2.0, dk)


coarse_lamda_list = lamda_list # uses the spacing of mesh points
#coarse_lamda_list = arange(min(lamda_fine),max(lamda_fine)+.0001, 50) # fixed spacing
#coarse_lamda_list = [600] # single point
for lamda in coarse_lamda_list:
	print('lambda:',lamda)



	#print(len(nlist),len(klist))
	error_map = single_lamda_TRA_error_map(lamda = lamda,
			nlist = nlist,
			klist = klist,
			TR_pair_list_generator = TR_pair_list_generator,
			parameter_list_generator = parameter_list_generator) *100.0

	#print(error_map.size)




	f = figure(figsize =(3.2,2.5), dpi = 220*2.0/3.0 * 2.0)
	N, K = meshgrid(nlist, klist, indexing = 'ij')

	ylabel('k', fontsize = 10)
	xlabel('n',fontsize = 10)
	title('$\lambda$ = %.3f nm' % lamda)
	minorticks_on()
	gca().tick_params(axis='both', which='major', labelsize=10)


	if use_logscale:

		cf_levels  = compute_outer_log_levels(error_map, base10delta = 0.1)
		CF = contourf(N,K, error_map,  norm = LogNorm(), levels = cf_levels , locator=ticker.LogLocator() )
		cbar = colorbar(CF, ticks = compute_inner_log_levels(error_map, base10delta = 0.2))
		cbar.ax.set_ylabel('RMS Error (%)')
		#cbar.ax.yaxis.set_ticks(log10(cf_levels) , minor=True)
	else:
		cf_level_delta = 0.5
		error_min = error_map.min()
		error_max = error_map.max()
		cfmax = (error_max - error_min) *0.15 +error_min
		cf_levels = arange( floor(error_min), ceil(cfmax) ,cf_level_delta  )

		CF = contourf(N,K, error_map, levels = cf_levels )
		cbar = colorbar(CF)
		cbar.ax.set_ylabel('RMS Error (%)')

	min_indices = find_min_indices_2d_array(error_map)
	plot(nlist[min_indices[0]], klist[min_indices[1]], marker = 'o', markerfacecolor = 'white',markeredgecolor = 'white', markersize = 3)

	plot(fit_nk_f(lamda).real, fit_nk_f(lamda).imag, marker = '+', markerfacecolor = 'grey',markeredgecolor = 'grey')



	gcf().tight_layout(pad=0.15)
	savefig(map_direct+'error_map_lamda_%f.pdf'%lamda,dpi = 600, transparent = True)
	close(f)


#show()
