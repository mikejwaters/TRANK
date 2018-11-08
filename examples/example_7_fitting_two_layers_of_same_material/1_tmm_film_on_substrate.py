from basic_setup  import nk_f_list,  parameter_list_generator, layer_index_of_fit, angle_list
from TRANK import TMM_spectra, try_mkdir
from numpy import arange, savetxt, array

if __name__=='__main__':

	make_plots = True
	tmm_spectra_dir = 'tmm_predicted_spectra/'

	noise = 1.0/100.0

	dlamda_min = 1
	lamda_min = 300
	lamda_max = 1200

	lamda_list = arange(lamda_min, lamda_max + dlamda_min/2.0 , dlamda_min)

	# this little monster calculates the spectra in parallel
	spectra = array( TMM_spectra(lamda_list = lamda_list,
						nk_f = nk_f_list[layer_index_of_fit[0]],  # we may have buit our parameter_list_generator to have our layer of interest in it, but these functions are all built around analyze variations of our layer nk
						parameter_list_generator = parameter_list_generator) )

	#####spectra[lamda index][ T/R/A spectrum index] for reference!


	## for adding noise
	from numpy.random import rand

	### for plotting
	if make_plots:
		from matplotlib.pyplot import figure, plot, minorticks_on, show, xlabel, ylabel, close
		from matplotlib.pyplot import gca, subplots_adjust,legend,savefig, title, tight_layout

	### for saving data
	try_mkdir(tmm_spectra_dir)

	for angle_index in range(len(angle_list)):

		angle = angle_list[angle_index]


		T_list = spectra[:, angle_index*2]   + noise* (rand(len(lamda_list))-0.5)
		R_list = spectra[:, angle_index*2+1] + noise* (rand(len(lamda_list))-0.5)
		A_list = 1.0 - T_list - R_list


		##### just plotting now!
		file_name =  tmm_spectra_dir  + '%i_deg_spectra.txt'%angle
		savetxt(file_name, array([lamda_list,T_list, R_list, A_list]).T , header =  '# lamda T R A')



		if make_plots:
			figure(figsize=(3.2,2.5),dpi = 220.0*2/3) # put your screen dpi in here to get in-print size

			plot(lamda_list, A_list*100.0, linewidth = 1.0, color = 'r', label = 'Absorptance')
			plot(lamda_list, T_list*100.0, linewidth = 1.0, color = 'g', label = 'Transmittance')
			plot(lamda_list, R_list*100.0, linewidth = 1.0, color = 'b', label = 'Reflectance')

			legend( labelspacing = 0.09, handletextpad = 0.5, handlelength = 1.2,
				fontsize = 9 , ncol = 1, columnspacing = 0.4)# loc= 'center right' ,title = '$Y\%$, $a^{*}$, $b^{*}$')
			minorticks_on()
			gca().tick_params(axis='both', which='major', labelsize=10)

			ylabel('(%)',fontsize =12)
			xlabel('Wavelength (nm)',fontsize =12)
			title('%i$^{\circ}$'%angle, fontsize =12)

			gca().set_ylim(0.0,100.0)
			gca().set_xlim(lamda_min,lamda_max)
			tight_layout(pad=0.1)

			im_name = tmm_spectra_dir  + '%i_deg_spectra.pdf'%angle
			savefig(im_name,dpi = 600, transparent = True)


	if make_plots: show()
