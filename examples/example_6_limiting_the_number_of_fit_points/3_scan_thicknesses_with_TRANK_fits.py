

from TRANK import  try_mkdir,  error_adaptive_iterative_fit_spectra, extrap, nk_plot
from numpy import  inf, arange, loadtxt, pi, logspace, log10
from basic_setup import  spectrum_list_generator,   parameter_list_generator, lamda_min, lamda_max

from numpy.random import rand

#lamda_min, lamda_max = 400, 1000

min_thickness = 15
max_thickness = 31

#####uniform spacing
film_thickness_list = range(min_thickness, max_thickness , 1)

print(film_thickness_list)
####log spacing
#number_of_points = 10
#film_thickness_list = logspace(log10(min_thickness),log10(max_thickness), number_of_points)

### single point thickness
#film_thickness_list = [1.0, 2.0 ]



dlamda_min = 5
dlamda_max = 100
delta_weight = 0.2/dlamda_min
adaptation_threshold_min = 0.002
lamda_fine = arange(lamda_min, lamda_max + dlamda_min/2.0 , dlamda_min)

show_plots = False
tolerance = 1e-5
max_points = 20
max_passes = 5 # limit the passes for a quicker scan, 0 means it will guess the amount needed.

if __name__=='__main__':
	if show_plots:
		from matplotlib.pylab import show

	film_thickness = film_thickness_list[0]
	print ('------ Thickness:',film_thickness ,'------')



	data_directory = 'TRANK_nk_fit_%f_nm/'%film_thickness

	try_mkdir(data_directory)

	parameter_list_generator.thickness = film_thickness


	def fit_nk_f(lamda):
		return 4.0+1.0j +0.0*lamda

	nk_plot(fit_nk_f, lamda_fine = lamda_fine, lamda_list = lamda_fine, file_name = data_directory+'initial_nk.pdf',
				title_string='Initial nk',show_nodes = True, show_plots = show_plots)

	fit_nk_f,lamda_list = error_adaptive_iterative_fit_spectra(
				nk_f_guess = fit_nk_f,
				spectrum_list_generator = spectrum_list_generator,
				parameter_list_generator = parameter_list_generator,
				lamda_min = lamda_min,
				lamda_max = lamda_max,
				dlamda_min = dlamda_min,
				dlamda_max = dlamda_max,
				delta_weight = delta_weight, tolerance = 1e-4, interpolation_type = 'cubic',
				adaptation_threshold_max = 0.05, adaptation_threshold_min = adaptation_threshold_min,
				use_reducible_error = True,
				method='least_squares',
				KK_compliant = True, max_passes = 1, # <------ one pass with KK on to precondition
				reuse_mode = False,
				zero_weight_extra_pass = False,
				data_directory = data_directory,
				verbose = True, make_plots = True, show_plots = show_plots, interpolate_to_fine_grid_at_end = False,
				nk_spectrum_file_format =  'KK_precondition_nk_pass_%i.pdf',
				rms_spectrum_file_format = 'KK_precondition_rms_spectrum_pass_%i.pdf' )




	for film_thickness in film_thickness_list:

		print ('\n------ Thickness:',film_thickness ,'------\n')
		data_directory = 'TRANK_nk_fit_%f_nm/'%film_thickness

		try_mkdir(data_directory)

		parameter_list_generator.thickness = film_thickness


		fit_nk_f, lamda_list = error_adaptive_iterative_fit_spectra(
					nk_f_guess = fit_nk_f,
					spectrum_list_generator = spectrum_list_generator,
					parameter_list_generator = parameter_list_generator,
					lamda_min = lamda_min,
					lamda_max = lamda_max,
					dlamda_min = dlamda_min,
					dlamda_max = dlamda_max,
					max_passes = max_passes,
					delta_weight = delta_weight, tolerance = tolerance, interpolation_type = 'linear',
					adaptation_threshold_max = 0.05, adaptation_threshold_min = adaptation_threshold_min,
					use_reducible_error = True,
					method='least_squares',
					KK_compliant = False,
					reuse_mode = True,
					lamda_list = lamda_list,
					zero_weight_extra_pass = False,
					data_directory = data_directory,
					delete_low_error_points = True,
					max_points = max_points,
					verbose = True, make_plots = True, show_plots = show_plots,
					nk_spectrum_file_format = 'TRANK_nk_pass_%i.pdf',
					rms_spectrum_file_format = 'rms_spectrum_pass_%i.pdf' )
