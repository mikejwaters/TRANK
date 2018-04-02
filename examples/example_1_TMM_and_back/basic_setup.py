
'''This file defines inputs by pulling nk from files'''


from TRANK import functionize_nk_file, extrap
from numpy import  inf, arange, loadtxt, pi

tm_polarization_fraction = 0.5  # this is unpolarized light



layer_index_of_fit = 1



kind = 'cubic'

#### air layer
def nk_f_air(lamda):
	return 1.0+0.0j*lamda

nk_f_list = [nk_f_air] #top working down
thickness_list = [inf]
coherency_list = ['i']

##### film layers  ######
film_thickness = 50 # nm
#nk_f_list.append(functionize_nk_file('Au-glass_10nm_30p_effective_nk.txt',skiprows=1, kind = kind))
nk_f_list.append(functionize_nk_file('Mie_Au-glass_10nm_30p_effective_nk.txt',skiprows=1, kind = kind))
thickness_list.append(film_thickness)
coherency_list.append('c')

###### substrate layer
substrate_thickness = 2e-3 *(1e9) # 2 mm in nm
nk_f_list.append(functionize_nk_file('Luke_Si3N4.txt',skiprows=1, kind = kind))
thickness_list.append(substrate_thickness)
coherency_list.append('i')

##### back air layer ####
nk_f_list.append(nk_f_air)
thickness_list.append(inf)
coherency_list.append('i')

###########
fit_nk_f = nk_f_list[layer_index_of_fit]


#################  illumination paramters


angle_list = [0,80]

tmm_spectra_dir = 'tmm_predicted_spectra/'


######
# this ugly hack only defines the TR_pair_list_generator if the data actually exists
from os.path import isfile
if isfile(tmm_spectra_dir+ '0_deg_spectra.txt'):

	Tf_Rf_pair_list = []

	for angle in angle_list:
			file_name = tmm_spectra_dir+ '%i_deg_spectra.txt'%angle
			lamda, T_list, R_list = loadtxt(file_name, usecols = [0,1,2], unpack = True)
			# lnear interpoation prevents interpoation of TRA values outside 0-100%
			Tf_Rf_pair_list.append([extrap(lamda, T_list, kind = 'linear' ), extrap(lamda, R_list, kind = 'linear' )])


	### this function is what we are creating as aa
	def TR_pair_list_generator(lamda): # for measured or simulated spectra
		TR_pair_list = []
		for Tf_Rf_pair in Tf_Rf_pair_list: #effectivly loops over angles, evaluated the interpolating functions
			TR_pair_list.append( ( Tf_Rf_pair[0](lamda), Tf_Rf_pair[1](lamda) ) )

		return TR_pair_list #must return list of [T,R] pairs matching the params


def parameter_list_generator(lamda): # layer geometry and neighbor properties
	param_list = []

	# order must match the TR_pair_list_generator
	for angle in angle_list:

		params = {'lamda' : lamda,
				'snell_angle_front' : angle* pi/180,
				'layer_index_of_fit' : layer_index_of_fit,
				'nk_f_list' : nk_f_list,
				'thickness_list' : thickness_list,
				'coherency_list' : coherency_list,
				'tm_polarization_fraction' : tm_polarization_fraction}


		param_list.append(params)

	return param_list

#############################
