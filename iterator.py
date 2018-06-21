


def error_adaptive_iterative_fit_spectra(
			nk_f_guess,
			spectrum_list_generator,
			parameter_list_generator,
			lamda_min,
			lamda_max,
			dlamda_min,
			dlamda_max,
			delta_weight = 0.1, tolerance = 1e-5, k_weight_fraction = 1.0,
			adaptation_threshold_max = 0.05, adaptation_threshold_min = 0.0005, adaptation_percentile = 85,
			max_passes = 0,
			lamda_list = [],
			use_reducible_error = True,
			reuse_mode = False,
			KK_compliant = False,
			interpolation_type = 'cubic',
			no_negative = True,
			interpolate_to_fine_grid_at_end = True,
			threads = 0,
			delete_low_error_points = False,
			max_points = 30,
			zero_weight_extra_pass = False, data_directory ='TRANK_nk_fit/', method = 'least_squares', verbose = True, write_nk_files = True,
			make_plots = True, show_plots = True, nk_spectrum_file_format = 'TRANK_nk_pass_%i.pdf', rms_spectrum_file_format = 'rms_spectrum_pass_%i.pdf' ):


	from TRANK import (fit_spectra_nk_sqr, fit_spectra_nk_sqr_KK_compliant,
						rms_error_spectrum, reducible_rms_error_spectrum, nk_plot, error_plot, try_mkdir)
	from time import time
	from numpy import floor, log2, ceil, linspace, diff, sqrt, mean, array, savetxt, percentile, argsort
	from copy import deepcopy
	if write_nk_files or make_plots: try_mkdir(data_directory)

	if show_plots:
		from matplotlib.pylab import show

	if reuse_mode == False: #picks lambda points accordingly
		from TRANK import compute_coarse_and_fine_grid
		lamda_list, lamda_fine = compute_coarse_and_fine_grid(dlamda_max, dlamda_min, lamda_max, lamda_min)
		dlamda_min = lamda_fine[1]-lamda_fine[0]
		dlamda_max = lamda_list[1]-lamda_list[0]

		power_of_2 = int(round( log2(dlamda_max/dlamda_min) ))
		if max_passes == 0:
			passes = int(power_of_2) + 1
		else:
			passes = int(max_passes)



	else:

		# Determines the fine grid based on the smallest delta lambda you have
		dlamda_min_found = min(diff(lamda_list))
		power_of_2 = int(round( log2(dlamda_min_found/dlamda_min) ))
		#print( log2(dlamda_min_found/dlamda_min)  )
		dlamda_min = dlamda_min_found/(2**power_of_2) # finds power of two spacing to match your dlamda_min

		### use what we find
		dlamda_max = dlamda_max_found = max(diff(lamda_list))

		nfine = ceil((lamda_max - lamda_min)/dlamda_min)
		lamda_fine = linspace(lamda_min,  lamda_max, nfine+1)
		#print ('lamda_fine', lamda_fine)

		if max_passes == 0:
			# this part guesses how many passes are required to reach the finest grid level
			passes = max( int(power_of_2) + 1, 2) # this makes sure that it runs on restart!
		else:
			passes = int(max_passes)


	if zero_weight_extra_pass: # this will fail if the num new points conidtion is met
		passes+=1


	fit_nk_f = nk_f_guess

	print ('dlamda_max:',dlamda_max )
	print ('dlamda_min:',dlamda_min )

	# literally jury rigging the conidtion so it starts the loop, ugly, but cleaner than the alternatives
	num_new_points = len(lamda_list)
	total_iteration_time = 0.0
	pass_number = 1
	new_lamda_list = [] # we add no new lamda points for the first pass
	indicies_to_delete = []
	while pass_number <= passes and num_new_points > 0:

		## delete pointless low error points
		if delete_low_error_points and len(indicies_to_delete)>0:
			print('Deleted Points:', array(lamda_list)[indicies_to_delete])
			for index in indicies_to_delete:
				lamda_list.pop(index)


		## add new lamda points from last pass, does nothing if it is the first pass
		lamda_list = sorted(new_lamda_list+list(lamda_list))
		if pass_number > 1:
			print('New Points:', new_lamda_list)
			print('--> Points Added: ', num_new_points)



		print('-----------> Pass %i/%i' % (pass_number,passes))
		print('--> Fitting %i Points' % len(lamda_list))


		# here we build the inputs for the fitter
		inputs = dict(lamda_list = lamda_list,
					spectrum_list_generator = spectrum_list_generator,
					parameter_list_generator = parameter_list_generator,
					nk_f_guess = fit_nk_f,
					delta_weight = delta_weight,
					tolerance = tolerance,
					no_negative = no_negative,
					k_weight_fraction = k_weight_fraction,
					interpolation_type = interpolation_type, method = method, threads = threads)

		t0 = time()
		if KK_compliant:
			inputs.update(dict(lamda_fine = lamda_fine))
			fit_nk_f = fit_spectra_nk_sqr_KK_compliant(**inputs ) # <-----
		else:
			fit_nk_f = fit_spectra_nk_sqr(**inputs)  # <-----
		pass_time = time()-t0

		total_iteration_time += pass_time
		print('Pass Time: %.1f seconds'%pass_time)


		rms_spectrum = rms_error_spectrum(lamda_list = lamda_list,
							nk_f = fit_nk_f,
							spectrum_list_generator = spectrum_list_generator,
							parameter_list_generator = parameter_list_generator, threads = threads)
		net_rms = sqrt( mean( array(rms_spectrum)**2 ) )
		max_rms = max(rms_spectrum)

		rms_spectrum_fine = rms_error_spectrum(lamda_list = lamda_fine,
						nk_f = fit_nk_f,
						spectrum_list_generator = spectrum_list_generator,
						parameter_list_generator = parameter_list_generator, threads = threads)
		net_rms_fine = sqrt( mean( array(rms_spectrum_fine)**2 ) )


		nk = fit_nk_f(lamda_list)
		if use_reducible_error == False:
			reducible_error_spectrum = []
			adaptation_threshold = max( min(percentile(rms_spectrum,adaptation_percentile),adaptation_threshold_max) , adaptation_threshold_min)
			if write_nk_files: savetxt(data_directory+'fit_nk.txt',array([lamda_list, nk.real, nk.imag, array(rms_spectrum)*100.0]).T)
		else:
			reducible_error_spectrum, irreducible_error_spectrum = reducible_rms_error_spectrum(
		 						lamda_list = lamda_list,
								nk_f = fit_nk_f,
								spectrum_list_generator = spectrum_list_generator,
								parameter_list_generator = parameter_list_generator, threads = threads)
			adaptation_threshold = max( min(percentile(reducible_error_spectrum, adaptation_percentile),adaptation_threshold_max) , adaptation_threshold_min)
			if write_nk_files: savetxt(data_directory+'fit_nk.txt',array([lamda_list, nk.real, nk.imag, array(rms_spectrum)*100.0, array(reducible_error_spectrum)*100]).T)


		print('Fine Grid Net RMS Error: %f %%' % (net_rms_fine*100))
		print('--> Net RMS Error: %f %%' % (net_rms*100))
		print('--> Adaptation Threshold: %f %%' % (adaptation_threshold* 100))


		if make_plots:
			err_fig = error_plot(lamda_list = lamda_list, rms_spectrum = rms_spectrum,
							adaptation_threshold = adaptation_threshold,
							adaptation_threshold_min = adaptation_threshold_min,
							adaptation_threshold_max = adaptation_threshold_max,
							reducible_error_spectrum = reducible_error_spectrum,
							file_name = data_directory + rms_spectrum_file_format % pass_number,
							title_string = 'Pass %i: Net RMS Error = %.3f %%' %( pass_number, net_rms*100),
							show_plots = show_plots )

			nk_fig = nk_plot(lamda_list = lamda_list, lamda_fine = lamda_fine, nkf = fit_nk_f,
				file_name = data_directory + nk_spectrum_file_format % pass_number ,title_string='TRANK Pass %i' % pass_number, show_nodes = True, show_plots = show_plots)

			if show_plots:
				show()


		############ adaptation
		if use_reducible_error:
			adaptation_spectrum = reducible_error_spectrum
		else:
			adaptation_spectrum = rms_spectrum


		refinement_method = 'interpolate_and_check_all'
		#### with our adaptation selection method set, we find new points
		if refinement_method == 'near_worst':
			new_lamda_list = []
			for i in range(len(lamda_list)-1):
				if (adaptation_spectrum[i] > adaptation_threshold) or (adaptation_spectrum[i+1] > adaptation_threshold): # should we refine?
					if (lamda_list[i+1] - lamda_list[i]) > dlamda_min: # if the gap is bigger than the minimum, then it is allowed to refine
						new_lamda = (lamda_list[i]+lamda_list[i+1])/2.0
						new_lamda_list.append( new_lamda)

		elif refinement_method == 'interpolate_and_check_all':
			test_lamda_list = []
			for i in range(len(lamda_list)-1):
				if (lamda_list[i+1] - lamda_list[i]) > dlamda_min: # if the gap is bigger than the minimum, then it is allowed to refine
					test_lamda_list.append(  (lamda_list[i]+lamda_list[i+1])/2.0 )

			if use_reducible_error :
				test_error_spectrum, test_irreducible_error_spectrum = reducible_rms_error_spectrum(
			 						lamda_list = test_lamda_list,
									nk_f = fit_nk_f,
									spectrum_list_generator = spectrum_list_generator,
									parameter_list_generator = parameter_list_generator, threads = threads)
			else:
				test_error_spectrum = rms_error_spectrum(lamda_list = test_lamda_list,
									nk_f = fit_nk_f,
									spectrum_list_generator = spectrum_list_generator,
									parameter_list_generator = parameter_list_generator, threads = threads)
			#sorted_indices =  argsort(test_error_spectrum)
			new_lamda_list = []
			for i in range(len(test_lamda_list)):
				if (test_error_spectrum[i] > adaptation_threshold) :
					new_lamda_list.append( test_lamda_list[i])

		#### we combine the new points with the old at the start of the next pass
		### this is important to have here for the termination condition
		num_new_points = len(new_lamda_list)

		############ adaptation
		if delete_low_error_points:
			if ( (num_new_points + len(lamda_list)) > max_points):
				n_delete = num_new_points+len(lamda_list) - max_points
				sorted_indices =  argsort(adaptation_spectrum)
				sorted_indices_without_edge_values = list(sorted_indices)
				sorted_indices_without_edge_values.remove(0)
				sorted_indices_without_edge_values.remove(len(adaptation_spectrum)-1)
				indicies_to_delete  = sorted_indices_without_edge_values[0:n_delete]
				indicies_to_delete.sort(reverse = True)



		#### doing the stuff for the last extra pass if there is one
		if zero_weight_extra_pass:
			if (pass_number +1) == passes: # normal zero_weight_extra_pass , just finished second to last pass
				delta_weight = 0.0
				tolerance = 1e-8
				num_new_points = 1 # jury rig it so it continues regardless of state of convergence
				pass_number += 1

			elif num_new_points == 0 and pass_number < passes: # test if terminates early, but still needs that extra pass
				delta_weight = 0.0
				tolerance = 1e-8
				num_new_points = 1 # jury rig it so it continues regardless of state of convergence
				pass_number = passes # skip to last passes
				print('--> Skipping to extra pass due to early conidtion statisfaction')
		else:
			pass_number += 1

	print('Total Iterating Time: %.1f seconds'%total_iteration_time)


	if interpolate_to_fine_grid_at_end and write_nk_files:
		print('Interpolating to fine grid and saving...')
		nk = fit_nk_f(lamda_fine)

		if use_reducible_error == False:
			savetxt(data_directory+'fit_nk_fine.txt',array([lamda_fine, nk.real, nk.imag, array(rms_spectrum_fine)*100.0]).T)
		else:
			reducible_error_spectrum_fine, irreducible_error_spectrum = reducible_rms_error_spectrum(
		 						lamda_list = lamda_fine,
								nk_f = fit_nk_f,
								spectrum_list_generator = spectrum_list_generator,
								parameter_list_generator = parameter_list_generator, threads = threads)
			savetxt(data_directory+'fit_nk_fine.txt',array([lamda_fine, nk.real, nk.imag, array(rms_spectrum_fine)*100.0,  array(reducible_error_spectrum_fine)*100.0 ]).T)


	return fit_nk_f, lamda_list
















def error_adaptive_iterative_fit(
			nk_f_guess,
			TR_pair_list_generator,
			parameter_list_generator,
			lamda_min,
			lamda_max,
			dlamda_min,
			dlamda_max,
			delta_weight = 0.1, tolerance = 1e-5,
			adaptation_threshold_max = 0.01, adaptation_threshold_min = 0.0005,
			max_passes = 0,
			extra_passes = 0,
			lamda_list = [],
			use_reducible_error = True,
			reuse_mode = False,
			KK_compliant = False,
			interpolation_type = 'cubic',
			zero_weight_extra_pass = False, data_directory ='TRANK_nk_fit/', method = 'least_squares', verbose = True,
			make_plots = True, show_plots = True, nk_spectrum_file_format = 'TRANK_nk_pass_%i.pdf', rms_spectrum_file_format = 'rms_spectrum_pass_%i.pdf'  ):


	from TRANK import (fit_TRA_nk_sqr, fit_TRA_nk_sqr_KK_compliant,
						rms_TRA_error_spectrum, reducible_rms_TRA_error_spectrum, nk_plot, error_plot, try_mkdir)
	from time import time
	from numpy import floor, log2, ceil, linspace, diff, sqrt, mean, array, savetxt, percentile
	try_mkdir(data_directory)

	if show_plots:
		from matplotlib.pylab import show

	if reuse_mode == False: #picks lambda points accordingly
		ncoarse = ceil((lamda_max - lamda_min)/dlamda_max)
		dlamda_max =  (lamda_max - lamda_min)/ncoarse
		lamda_list = linspace(lamda_min,  lamda_max, ncoarse+1)
		#print(lamda_list)
		power_of_2 = int(round( log2(dlamda_max/dlamda_min) ))
		#print(log2(dlamda_max/dlamda_min), power_of_2)
		dlamda_min = dlamda_max/(2**power_of_2)
		lamda_fine = linspace(lamda_min,  lamda_max, ncoarse*(2**power_of_2)+1)
		if max_passes == 0:
			passes = int(power_of_2) + 1
		else:
			passes = int(max_passes)

	else:
		dlamda_min_found = min(diff(lamda_list))
		power_of_2 = int(round( log2(dlamda_min_found/dlamda_min) ))
		#print( log2(dlamda_min_found/dlamda_min)  )
		dlamda_min = dlamda_min_found/(2**power_of_2)
		#print ('dlamda_min', dlamda_min)
		nfine = ceil((lamda_max - lamda_min)/dlamda_min)
		#print ('nfine', nfine)
		lamda_fine = linspace(lamda_min,  lamda_max, nfine+1)
		#print ('lamda_fine', lamda_fine)

		if max_passes == 0:
			passes = max( int(power_of_2) + 1, 2) # this makes sure that it runs on restart!
		else:
			passes = int(max_passes)



	if zero_weight_extra_pass: # this will fail if the num new points conidtion is met
		passes+=1


	fit_nk_f = nk_f_guess

	print ('dlamda_max:',dlamda_max )
	print ('dlamda_min:',dlamda_min )

	num_new_points = len(lamda_list)
	total_iteration_time = 0.0
	pass_number = 1
	while pass_number <= passes and num_new_points > 0:

		print('-----------> Pass %i/%i' % (pass_number,passes))
		print('--> Fitting %i Points' % len(lamda_list))


		# here we build the inputs for the fitter
		inputs = dict(lamda_list = lamda_list,
					TR_pair_list_generator = TR_pair_list_generator,
					parameter_list_generator = parameter_list_generator,
					nk_f_guess = fit_nk_f,
					delta_weight = delta_weight,
					tolerance = tolerance,
					interpolation_type = interpolation_type, method = method)

		t0 = time()
		if KK_compliant:
			inputs.update(dict(lamda_fine = lamda_fine))
			fit_nk_f = fit_TRA_nk_sqr_KK_compliant(**inputs ) # <-----
		else:
			fit_nk_f = fit_TRA_nk_sqr(**inputs)  # <-----
		pass_time = time()-t0

		total_iteration_time += pass_time
		print('Pass Time: %.1f seconds'%pass_time)


		rms_spectrum = rms_TRA_error_spectrum(lamda_list = lamda_list,
							nk_f = fit_nk_f,
							TR_pair_list_generator = TR_pair_list_generator,
							parameter_list_generator = parameter_list_generator)
		net_rms = sqrt( mean( array(rms_spectrum)**2 ) )
		max_rms = max(rms_spectrum)

		rms_spectrum_fine = rms_TRA_error_spectrum(lamda_list = lamda_fine,
						nk_f = fit_nk_f,
						TR_pair_list_generator = TR_pair_list_generator,
						parameter_list_generator = parameter_list_generator)
		net_rms_fine = sqrt( mean( array(rms_spectrum_fine)**2 ) )

		### saving the pass data
		nk = fit_nk_f(lamda_list)
		savetxt(data_directory+'fit_nk.txt',array([lamda_list, nk.real, nk.imag, array(rms_spectrum)*100.0]).T)

		if use_reducible_error:
			reducible_error_spectrum, irreducible_error_spectrum = reducible_rms_TRA_error_spectrum(
		 						lamda_list = lamda_list,
								nk_f = fit_nk_f,
								TR_pair_list_generator = TR_pair_list_generator,
								parameter_list_generator = parameter_list_generator)
			adaptation_threshold = max( min(percentile(reducible_error_spectrum,85),adaptation_threshold_max) , adaptation_threshold_min)
		else:
			reducible_error_spectrum = []
			adaptation_threshold = max( min(percentile(rms_spectrum,85),adaptation_threshold_max) , adaptation_threshold_min)


		print('Fine Grid Net RMS Error: %f %%' % (net_rms_fine*100))
		print('--> Net RMS Error: %f %%' % (net_rms*100))
		print('--> Adaptation Threshold: %f %%' % (adaptation_threshold* 100))


		if make_plots:
			err_fig = error_plot(lamda_list = lamda_list, rms_spectrum = rms_spectrum,
							adaptation_threshold = adaptation_threshold,
							adaptation_threshold_min = adaptation_threshold_min,
							adaptation_threshold_max = adaptation_threshold_max,
							reducible_error_spectrum = reducible_error_spectrum,
							file_name = data_directory+ rms_spectrum_file_format % pass_number,
							title_string = 'Pass %i: Net RMS Error = %.3f %%' %( pass_number, net_rms*100),
							show_plots = show_plots )

			nk_fig = nk_plot(lamda_list = lamda_list, lamda_fine = lamda_fine, nkf = fit_nk_f,
				file_name = data_directory + nk_spectrum_file_format % pass_number ,title_string='TRANK Pass %i' % pass_number, show_nodes = True, show_plots = show_plots)

			if show_plots:
				show()


		if use_reducible_error:
			adaptation_spectrum = reducible_error_spectrum
		else:
			adaptation_spectrum = rms_spectrum
		############ adaptation
		new_lamda_list = []
		#adaptation_threshold = max(rms_spectrum )/2.0
		for i in range(len(lamda_list)-1):
			if (adaptation_spectrum[i] > adaptation_threshold) or (adaptation_spectrum[i+1] > adaptation_threshold): # should we refine?
				if (lamda_list[i+1] - lamda_list[i]) > dlamda_min: # if the gap is bigger than the minimum, then it is allowed to refine
					new_lamda = (lamda_list[i]+lamda_list[i+1])/2.0
					new_lamda_list.append( new_lamda)

		#### now we combine the new points with the old
		num_new_points = len(new_lamda_list)
		print('New Points:', new_lamda_list)
		print('--> Points Added: ', num_new_points)

		lamda_list = sorted(new_lamda_list+list(lamda_list))

		#### doing the stuff for the last extra pass if there is one
		if zero_weight_extra_pass:
			if (pass_number +1) == passes: # normal zero_weight_extra_pass , just finished second to last pass
				delta_weight = 0.0
				tolerance = 1e-8
				num_new_points = 1 # jury rig it so it continues regardless of state of convergence
				pass_number += 1

			elif num_new_points == 0 and pass_number < passes: # test if terminates early, but still needs that extra pass
				delta_weight = 0.0
				tolerance = 1e-8
				num_new_points = 1 # jury rig it so it continues regardless of state of convergence
				pass_number = passes # skip to last passes
				print('--> Skipping to zero weight extra pass due to early conidtion statisfaction')
		else:
			pass_number += 1




	print('Total Iterating Time: %.1f seconds'%total_iteration_time)
	nk = fit_nk_f(lamda_fine)
	savetxt(data_directory+'fit_nk_fine.txt',array([lamda_fine, nk.real, nk.imag, array(rms_spectrum_fine)*100.0]).T)


	return fit_nk_f
