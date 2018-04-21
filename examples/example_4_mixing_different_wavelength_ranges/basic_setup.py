
'''This file defines inputs by pulling nk from files'''


from TRANK import functionize_nk_file, extrap
from numpy import pi, inf, loadtxt, arange





layer_index_of_fit = 1




#### air layer
def nk_f_air(lamda):
	return 1.0+0.0j*lamda

def nk_f_silica(lamda):
	return 1.5+0.0j*lamda



nk_f_list = [nk_f_air] #top working down
thickness_list = [inf]
coherency_list = ['i']

##### film layers  ######
film_thickness = 40.0 # nm
#nk_f_list.append(functionize_nk_file('Au-glass_10nm_30p_effective_nk.txt',skiprows=1, kind = kind))
nk_f_list.append(nk_f_air) # this is guess and it doesn't matter since it will be overwritten on the fly
thickness_list.append(film_thickness)
coherency_list.append('c')

###### substrate layer
substrate_thickness = 0.5e-3 *(1e9) # 2 mm in nm
nk_f_list.append(nk_f_silica)
thickness_list.append(substrate_thickness)
coherency_list.append('i')

##### back air layer ####
nk_f_list.append(nk_f_air)
thickness_list.append(inf)
coherency_list.append('i')

###########
#fit_nk_f = nk_f_list[layer_index_of_fit]


#################  experimental data
# lnear interpoation prevents interpoation of TRA values outside 0-100%

spectrum_function_list = []
spectrum_name_list = [] # this is for me to label things later
lamda_ranges = [] # [first index, last index, lamda min, lamda_max] # groups data by lamda ranges


#######
# In MS Excel, I get the best results when saving as windows .txt files and then removing any poorly converted characters in those files
T_data = loadtxt('transmittance.txt', skiprows = 2).T
lamda  = T_data[0]
spectrum_function_list.append( extrap(lamda, T_data[1] ) )
spectrum_name_list.append('0 deg Transmittance')

#######
R_data = loadtxt('normal_reflectance.txt', skiprows = 2).T
lamda = R_data[0]
spectrum_function_list.append( extrap(lamda, R_data[1] ) )
spectrum_name_list.append('0 deg unpolarized Reflectance')


#######
R_data = loadtxt('reflectance.txt', skiprows = 2).T
lamda = R_data[0]
# S => TE
# P => TM
spectrum_function_list.append( extrap(lamda, R_data[1] ) )
spectrum_name_list.append('55 deg unpolarized Reflectance')
spectrum_function_list.append( extrap(lamda, R_data[2] ) )
spectrum_name_list.append('65 deg unpolarized Reflectance')
spectrum_function_list.append( extrap(lamda, R_data[3] ) )
spectrum_name_list.append('75 deg unpolarized Reflectance')

spectrum_function_list.append( extrap(lamda, R_data[5] ) )
spectrum_name_list.append('55 deg P-Polarized Reflectance')
spectrum_function_list.append( extrap(lamda, R_data[6] ) )
spectrum_name_list.append('55 deg S-Polarized Reflectance')

spectrum_function_list.append( extrap(lamda, R_data[7] ) )
spectrum_name_list.append('65 deg P-Polarized Reflectance')
spectrum_function_list.append( extrap(lamda, R_data[8] ) )
spectrum_name_list.append('65 deg S-Polarized Reflectance')

spectrum_function_list.append( extrap(lamda, R_data[9] ) )
spectrum_name_list.append('75 deg P-Polarized Reflectance')
spectrum_function_list.append( extrap(lamda, R_data[10] ) )
spectrum_name_list.append('75 deg S-Polarized Reflectance')


lamda_min = spectrum_function_list[0].lower_bound
lamda_max = spectrum_function_list[0].upper_bound
for spectrum_function in spectrum_function_list[1:]:
	lamda_min = min(lamda_min, spectrum_function.lower_bound) #
	lamda_max = max(lamda_max, spectrum_function.upper_bound) #

### this function is what we are creating as
def spectrum_list_generator(lamda): # for measured or simulated spectra

	spectrum_list = []

	for spectrum_function in spectrum_function_list:
		if spectrum_function.is_in_bounds(lamda):
			spectrum_list.append(spectrum_function(lamda))

	return spectrum_list #must be a list of spectra  matching the params




#################################### illumination conditions


#we will attach this to
def generator_function(self,lamda): # layer geometry and neighbor properties

	thickness_list[layer_index_of_fit] = self.thickness
	# order must match the spectrum_list_generator
	param_list = []
	param_list.append( {
			'lamda' : lamda,
			'snell_angle_front' : 0 * pi/180,
			'layer_index_of_fit' : layer_index_of_fit,
			'nk_f_list' : nk_f_list,
			'thickness_list' : thickness_list,
			'coherency_list' : coherency_list,
			'tm_polarization_fraction' : 0.5,
			'spectrum' : 'T'} )

	param_list.append( {
			'lamda' : lamda,
			'snell_angle_front' : 0 * pi/180,
			'layer_index_of_fit' : layer_index_of_fit,
			'nk_f_list' : nk_f_list,
			'thickness_list' : thickness_list,
			'coherency_list' : coherency_list,
			'tm_polarization_fraction' : 0.5,
			'spectrum' : 'R'} )

	param_list.append( {
			'lamda' : lamda,
			'snell_angle_front' : 55 * pi/180,
			'layer_index_of_fit' : layer_index_of_fit,
			'nk_f_list' : nk_f_list,
			'thickness_list' : thickness_list,
			'coherency_list' : coherency_list,
			'tm_polarization_fraction' : 0.5,
			'spectrum' : 'R'} )

	param_list.append( {
			'lamda' : lamda,
			'snell_angle_front' : 65 * pi/180,
			'layer_index_of_fit' : layer_index_of_fit,
			'nk_f_list' : nk_f_list,
			'thickness_list' : thickness_list,
			'coherency_list' : coherency_list,
			'tm_polarization_fraction' : 0.5,
			'spectrum' : 'R'} )

	param_list.append( {
			'lamda' : lamda,
			'snell_angle_front' : 75 * pi/180,
			'layer_index_of_fit' : layer_index_of_fit,
			'nk_f_list' : nk_f_list,
			'thickness_list' : thickness_list,
			'coherency_list' : coherency_list,
			'tm_polarization_fraction' : 0.5,
			'spectrum' : 'R'} )

	param_list.append( {
			'lamda' : lamda,
			'snell_angle_front' : 55 * pi/180,
			'layer_index_of_fit' : layer_index_of_fit,
			'nk_f_list' : nk_f_list,
			'thickness_list' : thickness_list,
			'coherency_list' : coherency_list,
			'tm_polarization_fraction' : 1.0,
			'spectrum' : 'R'} )

	param_list.append( {
			'lamda' : lamda,
			'snell_angle_front' : 55 * pi/180,
			'layer_index_of_fit' : layer_index_of_fit,
			'nk_f_list' : nk_f_list,
			'thickness_list' : thickness_list,
			'coherency_list' : coherency_list,
			'tm_polarization_fraction' : 0.0,
			'spectrum' : 'R'} )

	param_list.append( {
			'lamda' : lamda,
			'snell_angle_front' : 65 * pi/180,
			'layer_index_of_fit' : layer_index_of_fit,
			'nk_f_list' : nk_f_list,
			'thickness_list' : thickness_list,
			'coherency_list' : coherency_list,
			'tm_polarization_fraction' : 1.0,
			'spectrum' : 'R'} )

	param_list.append( {
			'lamda' : lamda,
			'snell_angle_front' : 65 * pi/180,
			'layer_index_of_fit' : layer_index_of_fit,
			'nk_f_list' : nk_f_list,
			'thickness_list' : thickness_list,
			'coherency_list' : coherency_list,
			'tm_polarization_fraction' : 0.0,
			'spectrum' : 'R'} )

	param_list.append( {
			'lamda' : lamda,
			'snell_angle_front' : 75 * pi/180,
			'layer_index_of_fit' : layer_index_of_fit,
			'nk_f_list' : nk_f_list,
			'thickness_list' : thickness_list,
			'coherency_list' : coherency_list,
			'tm_polarization_fraction' : 1.0,
			'spectrum' : 'R'} )

	param_list.append( {
			'lamda' : lamda,
			'snell_angle_front' : 75 * pi/180,
			'layer_index_of_fit' : layer_index_of_fit,
			'nk_f_list' : nk_f_list,
			'thickness_list' : thickness_list,
			'coherency_list' : coherency_list,
			'tm_polarization_fraction' : 0.0,
			'spectrum' : 'R'} )

	good_param_list = []
	for spectrum_index in range(len(param_list)):
		if spectrum_function_list[spectrum_index].is_in_bounds(lamda):
			good_param_list.append(param_list[ spectrum_index ])


	return good_param_list

###### okay so this goofy stuff should allow us to update the thickness of our layer from other places
class updatable_parameter_list_generator:
	def __init__(self):
		self.thickness = 0.0
	def __call__(self,lamda):
		return self.generator_function(self,lamda)


parameter_list_generator = updatable_parameter_list_generator()
parameter_list_generator.thickness = thickness_list[layer_index_of_fit]
parameter_list_generator.generator_function = generator_function

#print (parameter_list_generator.thickness)
#print (parameter_list_generator(500.0))


#############################


if __name__ == '__main__':
	print ('Spectrum Ranges:')
	for i in range(len(spectrum_function_list)):
		print(spectrum_name_list[i] , spectrum_function_list[i].lower_bound, spectrum_function_list[i].upper_bound)
