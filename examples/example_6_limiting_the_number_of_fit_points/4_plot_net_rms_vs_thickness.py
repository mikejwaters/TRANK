





from numpy import loadtxt, sqrt, mean, arange, argsort
from matplotlib.pylab import *
from os.path import isfile
from glob import glob
directory_list = glob('TRANK_nk_fit_*_nm/')

#print (directory_list)
#film_thickness_list = arange(1,200, 0.5)

use_reducible_error = True
plot_reducible_error = False

net_rms_list = []
net_rms_fine_list = []

reducible_net_rms_list = []
reducible_net_rms_fine_list = []

irreducible_net_rms_list = []
irreducible_net_rms_fine_list = []

found_thicknesses = []
for data_directory in directory_list:

	if isfile(data_directory+'fit_nk_fine.txt'):
		film_thickness = float(data_directory.split('_')[3])
		found_thicknesses.append(film_thickness)
		#######################
		data = loadtxt(data_directory+'fit_nk.txt').T

		net_rms_list.append(sqrt( mean(data[3]**2)))
		if use_reducible_error:
			reducible_net_rms_list.append(sqrt( mean(data[4]**2)))
			irreducible_net_rms_list.append(sqrt( mean(data[3]**2 - data[4]**2 )))

		#######################
		data = loadtxt(data_directory+'fit_nk_fine.txt').T

		net_rms_fine_list.append(sqrt( mean(data[3]**2)))
		if use_reducible_error:
			reducible_net_rms_fine_list.append(sqrt( mean(data[4]**2)))
			irreducible_net_rms_fine_list.append(sqrt( mean(data[3]**2 - data[4]**2)))



order = argsort(found_thicknesses)
figure(figsize = (3.2,2.5), dpi= 220*2/3)
plot(array(found_thicknesses)[order], array(net_rms_list)[order], label = 'Non-Uniform $\lambda$ Points', marker = 'o', markersize = 2)
plot(array(found_thicknesses)[order], array(net_rms_fine_list)[order], label =  'Uniform Fine $\lambda$ Grid',  marker = 'o', markersize = 2)
if use_reducible_error:
	plot(array(found_thicknesses)[order], array(irreducible_net_rms_list)[order], label = 'Irreducible Non-Uniform $\lambda$ Points', marker = 'o', markersize = 2)
	plot(array(found_thicknesses)[order], array(irreducible_net_rms_fine_list)[order], label =  'Irreducible Uniform Fine $\lambda$ Grid',  marker = 'o', markersize = 2)
	if plot_reducible_error:
		plot(array(found_thicknesses)[order], array(reducible_net_rms_list)[order], label = 'Reducible Non-Uniform $\lambda$ Points', marker = 'o', markersize = 2)
		plot(array(found_thicknesses)[order], array(reducible_net_rms_fine_list)[order], label =  'Reducible Uniform Fine $\lambda$ Grid',  marker = 'o', markersize = 2)



xlabel('Film Thickness Used for Fit (nm)')
ylabel('Net RMS Error (%)')
legend(loc = 'upper left', fontsize = 8)
minorticks_on()
tight_layout(pad = 0.1)
subplots_adjust(left = 0.19)
savefig('RMS_error_vs_thickness.pdf',transparent = True)
show()
