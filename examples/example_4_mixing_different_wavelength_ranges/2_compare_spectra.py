###########################################################




show_plots = True



########### get stucture
from basic_setup  import  spectrum_list_generator,   parameter_list_generator, spectrum_name_list, spectrum_function_list


data_directory = 'TRANK_nk_fit/'
from TRANK import functionize_nk_file, TMM_spectra, error_plot, reducible_rms_error_spectrum, rms_error_spectrum


if __name__=='__main__':


	from numpy import loadtxt, array, arange
	lamda_fine = loadtxt(data_directory+'fit_nk_fine.txt', usecols = [0], unpack = True)
	lamda_list = loadtxt(data_directory+'fit_nk.txt', usecols = [0], unpack = True)
	fit_nk_f = functionize_nk_file(data_directory+'fit_nk.txt', skiprows = 0, kind = 'linear')

	#lamda_fine = arange(lamda_list.min(), lamda_list.max(), 100)

	spectra_from_fit = TMM_spectra(lamda_list = lamda_fine,
										nk_f = fit_nk_f,
										parameter_list_generator = parameter_list_generator)
	#spectra_from_fit[lamda][ param spectrum] # computational reasons why this is ordered is this way

	#now i need to build independent lamda lists and their spectra
	list_of_lamda_lists = []
	list_of_fit_spectra = []
	for spectrum_function in spectrum_function_list: # create the sub lists for these to all go into
		list_of_lamda_lists.append( [])
		list_of_fit_spectra.append( [])

	for lamda_index in range(len(lamda_fine )):
		lamda = lamda_fine[lamda_index]
		param_index = 0
		for spectrum_index in range(len(spectrum_function_list)):
			if spectrum_function_list[spectrum_index].is_in_bounds(lamda):
				(list_of_fit_spectra[spectrum_index]).append(  spectra_from_fit[lamda_index][param_index])
				(list_of_lamda_lists[spectrum_index]).append(lamda)
				param_index += 1





	from matplotlib.pylab import figure, plot, legend, minorticks_on, gca, title, tight_layout, ylabel, xlabel, savefig, show

	figure(figsize =(3.2 * 2, 2.5 * 2), dpi = 220*2.0/3.0 )

	for spectrum_index in range(len(spectrum_name_list)):


		name = spectrum_name_list[spectrum_index]

		fit_data = list_of_fit_spectra[spectrum_index]
		#print(len(fit_data))
		fit_lamda = list_of_lamda_lists[spectrum_index]
		exp_data = spectrum_function_list[spectrum_index](array(fit_lamda))

		p = plot( fit_lamda, array(fit_data)*100.0, label = 'Fit '+name, linestyle = '--',         linewidth = 1.0 )
		plot(     fit_lamda, exp_data*100.0, label = 'Exp. '+name, color = p[0].get_color(), linewidth = 1.0)


	legend(loc= 'center right', labelspacing = 0.09, handletextpad = 0.3, handlelength = 0.8, fontsize = 8 )
	minorticks_on()
	gca().tick_params(axis='both', which='major', labelsize=10)

	ylabel('(%)',fontsize =10)
	xlabel('Wavelength (nm)',fontsize =10)
	#angle = angle_list[angle_index]
	#title('%i$^{\circ}$'%angle, fontsize =12)

	gca().set_ylim(0.0,100.0)
	gca().set_xlim(min(lamda_fine),max(lamda_fine))

	tight_layout(pad=0.1)
	savefig(data_directory+'spectra_comparison.pdf' ,dpi = 600, transparent = True)



	###########################################
	rms_spectrum = rms_error_spectrum(lamda_list = lamda_list, nk_f = fit_nk_f,
					spectrum_list_generator = spectrum_list_generator,
					parameter_list_generator = parameter_list_generator)

	reducible_error_spectrum, irreducible_error_spectrum = reducible_rms_error_spectrum(lamda_list= lamda_list, nk_f = fit_nk_f,
								spectrum_list_generator = spectrum_list_generator,
								parameter_list_generator = parameter_list_generator)

	rms_spectrum_fine = rms_error_spectrum(lamda_list = lamda_fine, nk_f = fit_nk_f,
						spectrum_list_generator = spectrum_list_generator,
						parameter_list_generator = parameter_list_generator)

	reducible_error_spectrum_fine, irreducible_error_spectrum_fine = reducible_rms_error_spectrum(lamda_list= lamda_fine, nk_f = fit_nk_f,
									spectrum_list_generator = spectrum_list_generator,
									parameter_list_generator = parameter_list_generator)

	error_plot(lamda_list = lamda_list, rms_spectrum = rms_spectrum,
					adaptation_threshold = -1.0, adaptation_threshold_min = -1.0, adaptation_threshold_max = -1.0,
					reducible_error_spectrum = reducible_error_spectrum,
					lamda_fine = lamda_fine, rms_spectrum_fine = rms_spectrum_fine, reducible_error_spectrum_fine = reducible_error_spectrum_fine,
					title_string = '',
					file_name = data_directory+'final_error_spectrum_fine.pdf', zoom_window = [], show_plots = True )

	if show_plots: show()
