###########################################################



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




########### get stucture
from basic_setup  import  TR_pair_list_generator,   parameter_list_generator, angle_list


data_directory = 'TRANK_nk_fit/'
from TRANK import functionize_nk_file, TR, TRA_spectra


from matplotlib.pylab import *
from numpy import loadtxt, array
lamda_fine = loadtxt(data_directory+'fit_nk_fine.txt', usecols = [0], unpack = True)
fit_nk_f = functionize_nk_file(data_directory+'fit_nk_fine.txt', skiprows = 0, kind = 'cubic')



spectra_from_fit = array(TRA_spectra(lamda_list = lamda_fine,
									nk_f = fit_nk_f,
									parameter_list_generator = parameter_list_generator))
#spectra_from_fit[lamda][ param][ spectrum] # computational reasons why this is ordered is this way


original_spectra = TR_pair_list_generator(lamda_fine)
#original_spectra[parameter][spectrum][lamda] # again, computational reasons why this is ordered is this way


from multiprocessing import Pool, cpu_count

# you could do this if you had a mixed bag of conifurations
#num_parameters = len(parameter_list_generator(lamda = lamda_fine[0]))
#for parameter_number in range(num_parameters):

for parameter_number in range(len(angle_list)):

	figure(figsize =(3.2,2.5), dpi = 220*2.0/3.0 )

	T_fit = spectra_from_fit[:,parameter_number,0]
	R_fit = spectra_from_fit[:,parameter_number,1]
	A_fit = spectra_from_fit[:,parameter_number,2]


	T_data = original_spectra[parameter_number][0]
	R_data = original_spectra[parameter_number][1]
	A_data = 1.0 - T_data - R_data

	### T
	plot(lamda_fine,T_fit *100.0, label = 'T fit',  color = 'g', linewidth = 1.0, linestyle = '--')
	plot(lamda_fine,T_data*100.0, label = 'T data', color = 'g', linewidth = 1.0)

	### R
	plot(lamda_fine,R_fit *100.0, label = 'R fit',  color = 'b', linewidth = 1.0, linestyle = '--')
	plot(lamda_fine,R_data*100.0, label = 'R data', color = 'b', linewidth = 1.0)

	### A
	plot(lamda_fine,A_fit *100.0, label = 'A fit',  color = 'r', linewidth = 1.0, linestyle = '--')
	plot(lamda_fine,A_data*100.0, label = 'A data', color = 'r', linewidth = 1.0)


	legend(loc= 'center right', labelspacing = 0.09, handletextpad = 0.3, handlelength = 0.8, fontsize = 10 )
	minorticks_on()
	gca().tick_params(axis='both', which='major', labelsize=10)

	ylabel('(%)',fontsize =10)
	xlabel('Wavelength (nm)',fontsize =10)
	angle = angle_list[parameter_number]
	title('%i$^{\circ}$'%angle, fontsize =12)

	gca().set_ylim(0.0,100.0)
	gca().set_xlim(min(lamda_fine),max(lamda_fine))

	tight_layout(pad=0.1)
	savefig(data_directory+'spectrum_comparison_%i.pdf'%parameter_number ,dpi = 600, transparent = True)






if noshow==False:
	show()
