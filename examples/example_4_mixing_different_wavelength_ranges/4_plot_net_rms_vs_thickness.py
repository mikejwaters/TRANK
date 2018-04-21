




from numpy import loadtxt, sqrt, mean, arange, argsort
from matplotlib.pylab import *
from os.path import isfile
from glob import glob
directory_list = glob('TRANK_nk_fit_*_nm/')

#print (directory_list)
#film_thickness_list = arange(1,200, 0.5)

net_rms_list = []
net_rms_fine_list = []
found_thicknesses = []
for data_directory in directory_list:

	#data_directory = 'TRANK_nk_fit_%f_nm/'%film_thickness
	#print (film_thickness)

	if isfile(data_directory+'fit_nk_fine.txt'):
		film_thickness = float(data_directory.split('_')[3])

		###########
		data = loadtxt(data_directory+'fit_nk_fine.txt').T

		rms = data[3]
		net_rms = sqrt( mean(rms**2))

		net_rms_fine_list.append(net_rms)
		found_thicknesses.append(film_thickness)

		############
		data = loadtxt(data_directory+'fit_nk.txt').T

		rms = data[3]
		net_rms = sqrt( mean(rms**2))
		net_rms_list.append(net_rms)

order = argsort(found_thicknesses)
figure(figsize = (3.2,2.5), dpi= 220*2/3)
plot(array(found_thicknesses)[order], array(net_rms_list)[order], label = 'Coarse nk grid', marker = 'o', markersize = 2)
plot(array(found_thicknesses)[order], array(net_rms_fine_list)[order], label =  'Fine nk grid',  marker = 'o', markersize = 2)
xlabel('Film Thickness Used for Fit (nm)')
ylabel('Net RMS Error (%)')
legend()
minorticks_on()
tight_layout(pad = 0.1)
savefig('RMS_error_vs_thickness.pdf',transparent = True)
show()
