

from TRANK import rms_error_spectrum, reducible_rms_error_spectrum, nk_plot, try_mkdir, functionize_nk_file, extrap,  error_adaptive_iterative_fit_spectra

if __name__=='__main__':

	show_plots = False
	data_directory = 'TRANK_nk_fit/'

	try_mkdir(data_directory)
	if show_plots:
		from matplotlib.pylab import show

	###################### structure parameters
	from basic_setup  import  spectrum_list_generator,   parameter_list_generator, lamda_min, lamda_max



	###########
	from os import getcwd, walk, listdir
	from os.path import isfile

	from numpy import arange, loadtxt, sqrt, mean, array

	dlamda_min = 4
	dlamda_max = 120
	lamda_max = 1300
	scaled_weight_crossover = 21.612856
	delta_weight = 2* scaled_weight_crossover/dlamda_max
	print ('delta_weight:', delta_weight)
	use_reducible_error = True




	def fit_nk_f(lamda):
		return 2.0+0.05j+0.0*lamda
	### just for plotting
	lamda_plot = arange(lamda_min, lamda_max + dlamda_min/2.0 , dlamda_min)
	nk_plot(fit_nk_f,
			lamda_fine = lamda_plot,
			file_name = data_directory+'initial_nk.pdf',
			title_string='Initial nk',show_nodes = False, show_plots = show_plots)

	if show_plots: show()


	#This is a little bit different! we are doing performing the fit with KK compliance on for most of the steps. It's much faster with this much data.
	# we then relax the function without KK compliance in the next step but using only one pass so there is no Lambda-node refinment. This ensures a good fit at the bandwidth edges.
	fit_nk_f, lamda_list = error_adaptive_iterative_fit_spectra(
				nk_f_guess = fit_nk_f,
				spectrum_list_generator = spectrum_list_generator,
				parameter_list_generator = parameter_list_generator,
				lamda_min = lamda_min,
				lamda_max = lamda_max,
				dlamda_min = dlamda_min,
				dlamda_max = dlamda_max,
				delta_weight = delta_weight, k_weight_fraction = 1.0,
				tolerance = 1e-5, interpolation_type = 'linear',
				adaptation_threshold_max = 0.05, adaptation_threshold_min = 0.001,
				use_reducible_error = use_reducible_error,
				method='least_squares',
				KK_compliant = True,
				max_passes = 1,
				reuse_mode = False,
				zero_weight_extra_pass = False,
				interpolate_to_fine_grid_at_end = False,
				verbose = True, make_plots = True, show_plots = show_plots,
				nk_spectrum_file_format = 'TRANK_nk_KK_pass_%i.pdf', rms_spectrum_file_format = 'rms_spectrum_KK_pass_%i.pdf' )


	error_adaptive_iterative_fit_spectra(
				nk_f_guess = fit_nk_f,
				spectrum_list_generator = spectrum_list_generator,
				parameter_list_generator = parameter_list_generator,
				lamda_min = lamda_min,
				lamda_max = lamda_max,
				dlamda_min = dlamda_min,
				dlamda_max = dlamda_max,
				delta_weight = delta_weight, k_weight_fraction = 0.25,
				tolerance = 1e-5, interpolation_type = 'cubic',
				adaptation_threshold_max = 0.05, adaptation_threshold_min = 0.001,
				use_reducible_error = use_reducible_error,
				method='least_squares',
				KK_compliant = False,
				max_passes = 1,
				reuse_mode = True, lamda_list = lamda_list,
				zero_weight_extra_pass = False,
				verbose = True, make_plots = True, show_plots = show_plots,
				nk_spectrum_file_format = 'TRANK_nk_pass_%i.pdf', rms_spectrum_file_format = 'rms_spectrum_pass_%i.pdf' )
