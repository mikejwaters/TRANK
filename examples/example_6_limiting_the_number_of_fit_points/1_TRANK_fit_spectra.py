

from TRANK import rms_error_spectrum, nk_plot, try_mkdir, functionize_nk_file, extrap,  error_adaptive_iterative_fit_spectra

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

	dlamda_min = 5 # 5 is the default
	dlamda_max = 80 #
	delta_weight = 0.2/dlamda_min # 0.1 for scan
	film_thickness = 24.0

	max_points = 20 # 20
	max_passes = 10

	def fit_nk_f(lamda):
		return 4.0+1.0j +0.0*lamda
	method = 'least_squares'

	lamda_fine = arange(lamda_min, lamda_max + dlamda_min/2.0 , dlamda_min)
	nk_plot(fit_nk_f, lamda_fine = lamda_fine, lamda_list = [], file_name = data_directory+'initial_nk.pdf',
			title_string='Initial nk',show_nodes = False, show_plots = show_plots)

	if show_plots: show()


	parameter_list_generator.thickness = film_thickness
	fit_nk_f, lamda_list = error_adaptive_iterative_fit_spectra(
				nk_f_guess = fit_nk_f,
				spectrum_list_generator = spectrum_list_generator,
				parameter_list_generator = parameter_list_generator,
				lamda_min = lamda_min,
				lamda_max = lamda_max,
				dlamda_min = dlamda_min,
				dlamda_max = dlamda_max,
				delta_weight = delta_weight, tolerance = 1e-5, interpolation_type = 'linear',
				adaptation_threshold_max = 0.05, adaptation_threshold_min = 0.002,
				use_reducible_error = True,
				method = method,
				KK_compliant = True,
				max_passes = 10, # was 5
				delete_low_error_points = True,
				max_points = max_points,
				reuse_mode = False, lamda_list = [],
				zero_weight_extra_pass = False,
				verbose = True, make_plots = True, show_plots = show_plots, interpolate_to_fine_grid_at_end = False,
				nk_spectrum_file_format = 'KK_precondition_nk_pass_%i.pdf',
				rms_spectrum_file_format = 'KK_precondition_rms_spectrum_pass_%i.pdf' )


	fit_nk_f = error_adaptive_iterative_fit_spectra(
				nk_f_guess = fit_nk_f,
				spectrum_list_generator = spectrum_list_generator,
				parameter_list_generator = parameter_list_generator,
				lamda_min = lamda_min,
				lamda_max = lamda_max,
				dlamda_min = dlamda_min,
				dlamda_max = dlamda_max,
				delta_weight = delta_weight, tolerance = 1e-5, interpolation_type = 'linear',
				adaptation_threshold_max = 0.05, adaptation_threshold_min = 0.002,
				use_reducible_error = True,
				method = method,
				KK_compliant = False,
				max_passes = max_passes,
				delete_low_error_points = True,
				max_points = max_points,
				reuse_mode = True, lamda_list = lamda_list,
				zero_weight_extra_pass = False,
				verbose = True, make_plots = True, show_plots = show_plots,
				nk_spectrum_file_format = 'TRANK_nk_pass_%i.pdf',
				rms_spectrum_file_format = 'rms_spectrum_pass_%i.pdf' )
