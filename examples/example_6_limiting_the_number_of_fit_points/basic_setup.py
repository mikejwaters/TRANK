
'''This file defines inputs by pulling nk from files'''


from TRANK import functionize_nk_file, extrap
from numpy import  inf, arange, loadtxt, pi



layer_index_of_fit = 1



kind = 'linear'

#### air layer
def nk_f_air(lamda):
	return 1.0+0.0j*lamda





############## Structure Setup #############
nk_f_list = [nk_f_air] #top working down
thickness_list = [inf]
coherency_list = ['i']

##### film layers  ######
film_thickness = 25 # nm

nk_f_list.append(nk_f_air)
thickness_list.append(film_thickness)
coherency_list.append('c')

###### substrate layer
substrate_thickness = 520e-6 *(1e9) # 520 micons in nm
nk_f_silica = functionize_nk_file('fused_silica_substrate.txt', skiprows = 2, kind = kind)
nk_f_list.append(nk_f_silica)
thickness_list.append(substrate_thickness)
coherency_list.append('i')

##### back air layer ####
nk_f_list.append(nk_f_air)
thickness_list.append(inf)
coherency_list.append('i')

#### create the reversed lists
coherency_list_r = list(reversed(coherency_list))
nk_f_list_r = list(reversed(nk_f_list))
thickness_list_r = list(reversed(thickness_list))
layer_index_of_fit_r = len(thickness_list)-layer_index_of_fit - 1


###########
#fit_nk_f = nk_f_list[layer_index_of_fit]


#########################  experimental data

spectrum_function_list = []
spectrum_name_list = [] # this is for me to label things later


#### front side

R_data = loadtxt('Normal_Reflectance.txt', skiprows = 2).T
lamda = R_data[0]
# lnear interpoation prevents interpoation of TRA values outside 0-100%
spectrum_function_list.append( extrap(lamda, R_data[1], kind = 'linear' ) )
spectrum_name_list.append('0 deg Reflectance')

R_data = loadtxt('Reflectance.txt', skiprows = 2).T
lamda = R_data[0]


spectrum_function_list.append( extrap(lamda, R_data[1], kind = 'linear' ) )
spectrum_name_list.append('55 deg Unpolarized Reflectance')

if False:
	spectrum_function_list.append( extrap(lamda, R_data[2], kind = 'linear' ) )
	spectrum_name_list.append('65 deg Unpolarization Reflection')
	spectrum_function_list.append( extrap(lamda, R_data[3], kind = 'linear' ) )
	spectrum_name_list.append('75 deg Unpolarization Reflection')

if False:
	spectrum_function_list.append( extrap(lamda, R_data[5], kind = 'linear' ) )
	spectrum_name_list.append('55 deg P-polarized Reflectance')
	spectrum_function_list.append( extrap(lamda, R_data[6], kind = 'linear' ) )
	spectrum_name_list.append('55 deg S-polarized Reflectance')

if False:
	spectrum_function_list.append( extrap(lamda, R_data[7], kind = 'linear' ) )
	spectrum_name_list.append('65 deg P-polarized Reflectance')
	spectrum_function_list.append( extrap(lamda, R_data[8], kind = 'linear' ) )
	spectrum_name_list.append('65 deg S-polarized Reflectance')

	spectrum_function_list.append( extrap(lamda, R_data[9], kind = 'linear' ) )
	spectrum_name_list.append('75 deg P-polarized Reflectance')
	spectrum_function_list.append( extrap(lamda, R_data[10], kind = 'linear' ) )
	spectrum_name_list.append('75 deg S-polarized Reflectance')


T_data = loadtxt('Normal_Transmittance.txt', skiprows = 2).T
spectrum_function_list.append( extrap(T_data[0], T_data[1], kind = 'linear' ) )
spectrum_name_list.append('0 deg Transmittance')




########### back side ######
if True:
	R_data = loadtxt('Backside_Normal_Reflectance.txt', skiprows = 2).T
	lamda = R_data[0]
	# lnear interpoation prevents interpoation of TRA values outside 0-100%
	spectrum_function_list.append( extrap(lamda, R_data[1], kind = 'linear' ) )
	spectrum_name_list.append('0 deg Reflectance')

	R_data = loadtxt('Backside_Reflectance.txt', skiprows = 2).T
	lamda = R_data[0]

	if True:
		spectrum_function_list.append( extrap(lamda, R_data[1], kind = 'linear' ) )
		spectrum_name_list.append('55 deg Unpolarized Reflectance')

	if False:
		spectrum_function_list.append( extrap(lamda, R_data[2], kind = 'linear' ) )
		spectrum_name_list.append('65 deg Unpolarization Reflection')

		spectrum_function_list.append( extrap(lamda, R_data[3], kind = 'linear' ) )
		spectrum_name_list.append('75 deg Unpolarization Reflection')

	if False:
		spectrum_function_list.append( extrap(lamda, R_data[5], kind = 'linear' ) )
		spectrum_name_list.append('55 deg P-polarized Reflectance')
		spectrum_function_list.append( extrap(lamda, R_data[6], kind = 'linear' ) )
		spectrum_name_list.append('55 deg S-polarized Reflectance')

		spectrum_function_list.append( extrap(lamda, R_data[7], kind = 'linear' ) )
		spectrum_name_list.append('65 deg P-polarized Reflectance')
		spectrum_function_list.append( extrap(lamda, R_data[8], kind = 'linear' ) )
		spectrum_name_list.append('65 deg S-polarized Reflectance')

	if False:
		spectrum_function_list.append( extrap(lamda, R_data[9], kind = 'linear' ) )
		spectrum_name_list.append('75 deg P-polarized Reflectance')
		spectrum_function_list.append( extrap(lamda, R_data[10], kind = 'linear' ) )
		spectrum_name_list.append('75 deg S-polarized Reflectance')


	T_data = loadtxt('Backside_Normal_Transmittance.txt', skiprows = 2).T
	spectrum_function_list.append( extrap(T_data[0], T_data[1], kind = 'linear' ) )
	spectrum_name_list.append('0 deg Transmittance')





### this function is what we are creating as
def spectrum_list_generator(lamda): # for measured or simulated spectra

	spectrum_list = []

	for spectrum_function in spectrum_function_list:
		if spectrum_function.is_in_bounds(lamda):
			spectrum_list.append(spectrum_function(lamda))

	return spectrum_list #must be a list of spectra  matching the params

#################################### find Wavelength ranges

lamda_min = spectrum_function_list[0].lower_bound
lamda_max = spectrum_function_list[0].upper_bound
for spectrum_function in spectrum_function_list[1:]:
	lamda_min = min(lamda_min, spectrum_function.lower_bound) #
	lamda_max = max(lamda_max, spectrum_function.upper_bound) #





#################################### illumination conditions
# P => TM
# S => TE
def generator_function(self,lamda): # layer geometry and neighbor properties

	thickness_list[layer_index_of_fit] = self.thickness
	thickness_list_r[layer_index_of_fit_r] = self.thickness

	# order must match the spectrum_list_generator
	param_list = []

	#### frontside ####
	param_list.append( {
			'lamda' : lamda,
			'snell_angle_front' : 0 * pi/180,
			'layer_index_of_fit' : layer_index_of_fit,
			'nk_f_list' : nk_f_list,
			'thickness_list' : thickness_list,
			'coherency_list' : coherency_list,
			'tm_polarization_fraction' : 0.5,
			'spectrum' : 'R'} )

	if True:
		param_list.append( {
				'lamda' : lamda,
				'snell_angle_front' : 55 * pi/180,
				'layer_index_of_fit' : layer_index_of_fit,
				'nk_f_list' : nk_f_list,
				'thickness_list' : thickness_list,
				'coherency_list' : coherency_list,
				'tm_polarization_fraction' : 0.5,
				'spectrum' : 'R'} )

	if False:
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


	if False:
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

	if False:
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


	param_list.append( {
			'lamda' : lamda,
			'snell_angle_front' : 0 * pi/180,
			'layer_index_of_fit' : layer_index_of_fit,
			'nk_f_list' : nk_f_list,
			'thickness_list' : thickness_list,
			'coherency_list' : coherency_list,
			'tm_polarization_fraction' : 0.5,
			'spectrum' : 'T'} )


	#### backside ####
	if True:
		param_list.append( {
				'lamda' : lamda,
				'snell_angle_front' : 0 * pi/180,
				'layer_index_of_fit' : layer_index_of_fit_r,
				'nk_f_list' : nk_f_list_r,
				'thickness_list' : thickness_list_r,
				'coherency_list' : coherency_list_r,
				'tm_polarization_fraction' : 0.5,
				'spectrum' : 'R'} )


		param_list.append( {
					'lamda' : lamda,
					'snell_angle_front' : 55 * pi/180,
					'layer_index_of_fit' : layer_index_of_fit_r,
					'nk_f_list' : nk_f_list_r,
					'thickness_list' : thickness_list_r,
					'coherency_list' : coherency_list_r,
					'tm_polarization_fraction' : 0.5,
					'spectrum' : 'R'} )

		if False:
			param_list.append( {
					'lamda' : lamda,
					'snell_angle_front' : 65 * pi/180,
					'layer_index_of_fit' : layer_index_of_fit_r,
					'nk_f_list' : nk_f_list_r,
					'thickness_list' : thickness_list_r,
					'coherency_list' : coherency_list_r,
					'tm_polarization_fraction' : 0.5,
					'spectrum' : 'R'} )

			param_list.append( {
				'lamda' : lamda,
				'snell_angle_front' : 75 * pi/180,
				'layer_index_of_fit' : layer_index_of_fit_r,
				'nk_f_list' : nk_f_list_r,
				'thickness_list' : thickness_list_r,
				'coherency_list' : coherency_list_r,
				'tm_polarization_fraction' : 0.5,
				'spectrum' : 'R'} )


		if False:
			param_list.append( {
					'lamda' : lamda,
					'snell_angle_front' : 55 * pi/180,
					'layer_index_of_fit' : layer_index_of_fit_r,
					'nk_f_list' : nk_f_list_r,
					'thickness_list' : thickness_list_r,
					'coherency_list' : coherency_list_r,
					'tm_polarization_fraction' : 1.0,
					'spectrum' : 'R'} )

			param_list.append( {
					'lamda' : lamda,
					'snell_angle_front' : 55 * pi/180,
					'layer_index_of_fit' : layer_index_of_fit_r,
					'nk_f_list' : nk_f_list_r,
					'thickness_list' : thickness_list_r,
					'coherency_list' : coherency_list_r,
					'tm_polarization_fraction' : 0.0,
					'spectrum' : 'R'} )


			param_list.append( {
					'lamda' : lamda,
					'snell_angle_front' : 65 * pi/180,
					'layer_index_of_fit' : layer_index_of_fit_r,
					'nk_f_list' : nk_f_list_r,
					'thickness_list' : thickness_list_r,
					'coherency_list' : coherency_list_r,
					'tm_polarization_fraction' : 1.0,
					'spectrum' : 'R'} )

			param_list.append( {
					'lamda' : lamda,
					'snell_angle_front' : 65 * pi/180,
					'layer_index_of_fit' : layer_index_of_fit_r,
					'nk_f_list' : nk_f_list_r,
					'thickness_list' : thickness_list_r,
					'coherency_list' : coherency_list_r,
					'tm_polarization_fraction' : 0.0,
					'spectrum' : 'R'} )

		if False:
			param_list.append( {
					'lamda' : lamda,
					'snell_angle_front' : 75 * pi/180,
					'layer_index_of_fit' : layer_index_of_fit_r,
					'nk_f_list' : nk_f_list_r,
					'thickness_list' : thickness_list_r,
					'coherency_list' : coherency_list_r,
					'tm_polarization_fraction' : 1.0,
					'spectrum' : 'R'} )


			param_list.append( {
					'lamda' : lamda,
					'snell_angle_front' : 75 * pi/180,
					'layer_index_of_fit' : layer_index_of_fit_r,
					'nk_f_list' : nk_f_list_r,
					'thickness_list' : thickness_list_r,
					'coherency_list' : coherency_list_r,
					'tm_polarization_fraction' : 0.0,
					'spectrum' : 'R'} )


		param_list.append( {
				'lamda' : lamda,
				'snell_angle_front' : 0 * pi/180,
				'layer_index_of_fit' : layer_index_of_fit_r,
				'nk_f_list' : nk_f_list_r,
				'thickness_list' : thickness_list_r,
				'coherency_list' : coherency_list_r,
				'tm_polarization_fraction' : 0.5,
				'spectrum' : 'T'} )


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




############################# a quick diagnostic


if __name__ == '__main__':
	print ('Lambda Range:', lamda_min, lamda_max)
	print ('%i Spectra with Ranges:' % len(spectrum_function_list))
	for i in range(len(spectrum_function_list)):
		print(spectrum_name_list[i] , spectrum_function_list[i].lower_bound, spectrum_function_list[i].upper_bound)
