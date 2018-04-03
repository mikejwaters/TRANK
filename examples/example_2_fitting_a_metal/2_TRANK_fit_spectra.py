

from TRANK import *

show_plots = False
data_directory = 'TRANK_nk_fit/'

try_mkdir(data_directory)
if show_plots:
	from matplotlib.pylab import show

###################### structure parameters
from basic_setup  import fit_nk_f, TR_pair_list_generator,   parameter_list_generator



###########
from os import getcwd, walk, listdir
from os.path import isfile

from numpy import arange, loadtxt, sqrt, mean, array

dlamda_min = 1
dlamda_max = 50
delta_weight = 0.02
lamda_min = 300
lamda_max = 1200
lamda_fine = arange(lamda_min, lamda_max + dlamda_min/2.0 , dlamda_min)



max_rms_cutoff = 5.0 #percentage points
net_rms_cutoff = 1.0
use_old_nk = False
has_old_nk = False
old_lamda = []
if isfile(data_directory+'fit_nk_fine.txt') and isfile(data_directory+'fit_nk.txt'): # fine has the complete set
	print('Found local data.')

	old_data = loadtxt( data_directory+'fit_nk.txt').T
	fit_nk_f =  functionize_nk_file(data_directory+'fit_nk.txt', skiprows = 0, kind = 'cubic')
	old_lamda = old_data[0]
	has_old_nk = True

	if has_old_nk:

		rms_spectrum = rms_TRA_error_spectrum(lamda_list = lamda_fine,
			nk_f = fit_nk_f,
			TR_pair_list_generator = TR_pair_list_generator,
			parameter_list_generator = parameter_list_generator)

		net_rms = sqrt( mean( array(rms_spectrum)**2 ) ) * 100.0
		max_rms = 	max(rms_spectrum) * 100.0

		print('nk found! RMS (max): %.2f (%.2f)'%(net_rms, max_rms))
		ylim = max_rms_cutoff - (max_rms_cutoff/net_rms_cutoff)*net_rms
		if max_rms  < ylim:
			use_old_nk = True
			passes = 2


#use_old_nk = False
if use_old_nk == False:
	from numpy.random import rand
	min_n, max_n  = 0.0, 2.0
	min_k, max_k  = 0.0, 0.1
	rand_n = rand(lamda_fine.size)*(max_n - min_n) + min_n
	rand_k = rand(lamda_fine.size)*(max_k - min_k) + min_k
	nf = extrap(lamda_fine, rand_n)
	kf = extrap(lamda_fine, rand_k)
	old_lamda = lamda_fine
	def fit_nk_f(lamda):
		return nf(lamda)+1.0j*kf(lamda)
	def fit_nk_f(lamda):
		return 1.0+0.0*lamda



nk_plot(fit_nk_f, lamda_fine = lamda_fine, lamda_list = old_lamda, file_name = data_directory+'initial_nk.pdf',title_string='Initial nk',show_nodes = True, show_plots = show_plots)
if show_plots:
	show()



error_adaptive_iterative_fit(
			nk_f_guess = fit_nk_f,
			TR_pair_list_generator = TR_pair_list_generator,
			parameter_list_generator = parameter_list_generator,
			lamda_min = lamda_min,
			lamda_max = lamda_max,
			dlamda_min = dlamda_min,
			dlamda_max = dlamda_max,
			delta_weight = delta_weight, 
			tolerance = 1e-5, interpolation_type = 'cubic',
			adaptation_threshold_max = 0.01, adaptation_threshold_min = 0.002,
			use_reducible_error = True,
			method='least_squares',
			KK_compliant = False,
			reuse_mode = use_old_nk, lamda_list = old_lamda,
			zero_weight_extra_pass = False,
			verbose = True, make_plots = True, show_plots = show_plots )
