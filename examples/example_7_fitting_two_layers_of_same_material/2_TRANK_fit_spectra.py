

from TRANK import rms_error_spectrum, nk_plot, try_mkdir, functionize_nk_file, extrap_c,  error_adaptive_iterative_fit_spectra, scan_for_scaled_weight_crossover

if __name__=='__main__':

	show_plots = False
	data_directory = 'TRANK_nk_fit/'

	try_mkdir(data_directory)
	if show_plots:
		from matplotlib.pylab import show

	###################### structure parameters
	from basic_setup  import  spectrum_list_generator,   parameter_list_generator



	###########
	from os import getcwd, walk, listdir
	from os.path import isfile

	from numpy import arange, loadtxt, sqrt, mean, array

	dlamda_min = 5
	dlamda_max = 80

	lamda_min = 300
	lamda_max = 1200

	scaled_weight_multiplier = 2.0
	auto_find_weight = True
	adaptation_threshold_min = 0.2/100 # good guess is noise/(2 sqrt(n_spectra))
	adaptation_threshold_max = 5.0/100

	# this is your guess
	def fit_nk_f(lamda):
		return 1.0+3.0j+0.0*lamda

	if auto_find_weight:
		scaled_weight_crossover = scan_for_scaled_weight_crossover(
									nk_f_guess = fit_nk_f,
									spectrum_list_generator = spectrum_list_generator,
									parameter_list_generator = parameter_list_generator,
									lamda_min = lamda_min,
									lamda_max = lamda_max,
									use_free_drude = False, # this isn't guaranteed to everfind the crossover with KK or free-druid on.
									show_plots = show_plots)
	else:
		scaled_weight_crossover = 1.510919


	lamda_fine = arange(lamda_min, lamda_max + dlamda_min/2.0 , dlamda_min)
	nk_plot(fit_nk_f, lamda_fine = lamda_fine, lamda_list = [], file_name = data_directory+'initial_nk.pdf',
			title_string='Initial nk',show_nodes = False, show_plots = show_plots)

	if show_plots: show()

	delta_weight = scaled_weight_multiplier * scaled_weight_crossover/dlamda_max

	if True: # turn off the Drude model (use_free_drude = False) to see the horrors of bandwidth edge effects on metals
		fit_nk_f, lamda_list = error_adaptive_iterative_fit_spectra(
					nk_f_guess = fit_nk_f,
					spectrum_list_generator = spectrum_list_generator,
					parameter_list_generator = parameter_list_generator,
					lamda_min = lamda_min,
					lamda_max = lamda_max,
					dlamda_min = dlamda_min,
					dlamda_max = dlamda_max,
					delta_weight = delta_weight, tolerance = 1e-5, interpolation_type = 'cubic',
					adaptation_threshold_max = adaptation_threshold_max, adaptation_threshold_min = adaptation_threshold_min,
					use_reducible_error = True,
					max_passes = 3,
					method='least_squares',
					KK_compliant = True,
					use_free_drude = True,
					reuse_mode = False, lamda_list = [],
					zero_weight_extra_pass = False,
					verbose = True, make_plots = True, show_plots = show_plots,
					nk_spectrum_file_format = 'TRANK_nk_free_drude_pass_%i.pdf', rms_spectrum_file_format = 'rms_spectrum_free_drude_pass_%i.pdf' )



	error_adaptive_iterative_fit_spectra(
				nk_f_guess = fit_nk_f,
				spectrum_list_generator = spectrum_list_generator,
				parameter_list_generator = parameter_list_generator,
				lamda_min = lamda_min,
				lamda_max = lamda_max,
				dlamda_min = dlamda_min,
				dlamda_max = dlamda_max,
				delta_weight = delta_weight, tolerance = 1e-5, interpolation_type = 'cubic',
				adaptation_threshold_max = adaptation_threshold_max, adaptation_threshold_min = adaptation_threshold_min,
				use_reducible_error = True,
				method='least_squares',
				KK_compliant = False,
				max_passes = 3,
				reuse_mode = True, lamda_list = lamda_list,
				zero_weight_extra_pass = False,
				verbose = True, make_plots = True, show_plots = show_plots,
				nk_spectrum_file_format = 'TRANK_nk_pass_%i.pdf', rms_spectrum_file_format = 'rms_spectrum_pass_%i.pdf' )
