



def dual_grid_direct_KK_n_from_lamda_k(lamda_list, lamda_fine,  k,  cshift = 1e-4) :
	#KK transform in Wavelength
	# cshift is basically scaled so its universal i think
	# testing shows that cubic interpoltion and trapezoind rule excede the accurcy of just using the fine grid due to error cancelation
	from scipy.integrate import simps, trapz
	from scipy.interpolate import griddata

	from numpy import array, zeros, pi
	#print (lamda_list, k)
	#rint (len(lamda_list), len(k))
	k_fine = griddata(array(lamda_list), array(k), (array(lamda_fine)), method='cubic', fill_value = 0.0)
	#print (k_fine)
	n_out = zeros(len(lamda_list)) ### we use the fine lamda grid for the KK integral but we only need to evalute it at the coarse lamba grid points!
	k_over_lamda = k_fine/lamda_fine
	for i in range(len(lamda_list)): #parralelize this in the future!
		lamda_out = lamda_list[i]

		#n_out[i] = 1.0 + 2.0/pi * simps( k_over_lamda * (1.0/ ( (lamda_out/lamda_fine)**2 - 1.0 + cshift*1.0j)), lamda_fine).real
		n_out[i] = 1.0 + 2.0/pi * trapz( k_over_lamda * (1.0/ ( (lamda_out/lamda_fine)**2 - 1.0 + cshift*1.0j)), lamda_fine).real

	return n_out





######### DKKT
from numpy import abs, pi
from numpy import log as ln
def g(x,y):
	cshift = 0.000001j
	return (x+y)*ln(abs(x+y)) + (x-y)*ln(abs(x-y)+cshift)


def n_jj_omega(jj, omega, omega_list): # the index kernel
	term1 = g(omega,omega_list[jj-1])/(omega_list[jj] - omega_list[jj-1])
	term2 = - ( (omega_list[jj+1] - omega_list[jj-1]) * g(omega,omega_list[jj]) )/( (omega_list[jj] - omega_list[jj-1]) * (omega_list[jj+1] - omega_list[jj]))
	term3 = g(omega,omega_list[jj+1])/(omega_list[jj+1] - omega_list[jj])
	return (term1 + term2 + term3).real/pi



def DKKT_n_from_lamda_k(lamda_list, k):
	n_list = zeros(len(lamda_list))
	omega_list = 1.0/array(lamda_list)
	for i in range(len(lamda_list)):
		lamda = lamda_list[i]

		n = 0
		for jj in range(1,len(lamda_list)-1):
			n += n_jj_omega(jj = jj, omega = 1.0/lamda_list[i], omega_list = omega_list  ) * k[jj]

		n_list[i] = n

	return n_list



def single_point_DKKT_from_omega_k(inputs):
	point_index = inputs['point_index']
	omega_list = inputs['omega_list']
	k = inputs['k']

	n = 0
	for jj in range(1,len(omega_list)-1):
		n += n_jj_omega(jj = jj, omega = omega_list[point_index], omega_list = omega_list  ) * k[jj]
	return n

from numpy import array
def parallel_DKKT_n_from_lamda_k(lamda_list, k, compute_pool):
	# we could build the KK transform matrix out of the lambda points and
	# not have to recompute it later. but this is very clear for now

	omega_list = 1.0/array(lamda_list)
	input_list = []
	for i in range(len(lamda_list)):
		input_list.append({'point_index': i, 'omega_list': omega_list, 'k': k})

	n_list = compute_pool.map(single_point_DKKT_from_omega_k, input_list )
	return n_list
