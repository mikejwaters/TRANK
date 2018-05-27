force_mix = True

from time import time
from matplotlib.pylab import *
from numpy import *
c = 299792458


lamda_list = arange(300,2500.1,10)

n_0 = 1.5
epsilon_inf = n_0**2

def lorrentz_occilator_of_omega(omega,
					omega_p = 0.002, #height
					omega_0 = 1.0/800, # location
					gamma = 1.0/3000): #width

	return ( omega_p**2 )/(omega_0**2 - omega**2 + 1.0j*omega*gamma)


def drude_occilator_of_omega(omega,
					omega_p = 0.002, #height
					gamma = 1.0/3000): #location

	return  - ( omega_p**2 )/( omega**2 - 1.0j*omega*gamma)

def epsilon_of_omega_test_0(omega): #  simple lorrentzian, a sanity test that shouldn't be affected by edges
	epsilon =  epsilon_inf  + 1.0 * lorrentz_occilator_of_omega(omega, omega_0 = 1.0/1000, gamma = 1.0/10000, omega_p = 0.001) \
							+ 0.0 * lorrentz_occilator_of_omega(omega, omega_0 = 1.0/2100, gamma = 1.0/4000, omega_p = 0.0005) \
							+ 0.0 * drude_occilator_of_omega(omega, omega_p = 0.001, gamma = 1.0/2000)
	return epsilon

def epsilon_of_omega_test_1(omega): #  sharp peak and atruncated peak with negative slope, easiest to correct for
	epsilon =  epsilon_inf  + 1.0 * lorrentz_occilator_of_omega(omega, omega_0 = 1.0/800, gamma = 1.0/10000, omega_p = 0.001) \
							+ 1.0 * lorrentz_occilator_of_omega(omega, omega_0 = 1.0/2100, gamma = 1.0/4000, omega_p = 0.0005) \
							+ 0.0 * drude_occilator_of_omega(omega, omega_p = 0.001, gamma = 1.0/2000)
	return epsilon

def epsilon_of_omega_test_2(omega): # metallic conductivivty is hardest
	epsilon =  epsilon_inf  + 1.0 * lorrentz_occilator_of_omega(omega, omega_0 = 1.0/800, gamma = 1.0/10000, omega_p = 0.001) \
							+ 0.0 * lorrentz_occilator_of_omega(omega, omega_0 = 1.0/2100, gamma = 1.0/4000, omega_p = 0.0005) \
							+ 1.0 * drude_occilator_of_omega(omega, omega_p = 0.001, gamma = 1.0/2000)

	return epsilon


def epsilon_of_omega_test_3(omega): # broad peak touching the inner band width edge, no correction yet, but not as bad.
	epsilon =  epsilon_inf  + 0.0 * lorrentz_occilator_of_omega(omega, omega_0 = 1.0/800, gamma = 1.0/10000, omega_p = 0.001) \
							+ 1.0 * lorrentz_occilator_of_omega(omega, omega_0 = 1.0/400, gamma = 1.0/700, omega_p = 0.004) \
							+ 0.0 * drude_occilator_of_omega(omega, omega_p = 0.001, gamma = 1.0/2000)

	return epsilon

def epsilon_of_omega_test_4(omega): # broad peak just outside the inner band width edge, no correction yet, really hard
	epsilon =  epsilon_inf  + 0.0 * lorrentz_occilator_of_omega(omega, omega_0 = 1.0/800, gamma = 1.0/10000, omega_p = 0.001) \
							+ 1.0 * lorrentz_occilator_of_omega(omega, omega_0 = 1.0/250, gamma = 1.0/700, omega_p = 0.004) \
							+ 0.0 * drude_occilator_of_omega(omega, omega_p = 0.001, gamma = 1.0/2000)

	return epsilon



def epsilon_of_omega_test_5(omega): # mettallic conductivivty across the whole range, also hard
	epsilon =  epsilon_inf  + 0.0 * lorrentz_occilator_of_omega(omega, omega_0 = 1.0/800, gamma = 1.0/10000, omega_p = 0.001) \
							+ 0.0 * lorrentz_occilator_of_omega(omega, omega_0 = 1.0/400, gamma = 1.0/800, omega_p = 0.004) \
							+ 1.0 * drude_occilator_of_omega(omega, omega_p = 0.03, gamma = 1.0/250)

	return epsilon


epsilon_of_omega = epsilon_of_omega_test_4

omega = 1.0/lamda_list
epsilon = epsilon_of_omega(omega)
e1 = epsilon.real
e2 = epsilon.imag
emod = absolute(epsilon)
n = sqrt((emod+e1)/2)
k = sqrt((emod-e1)/2)




from TRANK import upper_bound_extrapolation_order_0, upper_bound_extrapolation_order_1




figure()
ub0 = upper_bound_extrapolation_order_0(lamda_list, k)
ub1 = upper_bound_extrapolation_order_1(lamda_list, k)
plot(lamda_list, ub0, label = 'ub0')
plot(lamda_list, ub1, label = 'ub01')

figure()

plot(lamda_list, n - n , label = 'n')

from TRANK import parallel_DKKT_n_from_lamda_k
from multiprocessing import Pool, cpu_count
my_compute_pool = Pool(cpu_count())
n_DKKT = parallel_DKKT_n_from_lamda_k(lamda_list, k, my_compute_pool)
shift_n = mean(n) - mean(n_DKKT)
#shift_n = n_0 - n_DKKT[0]
plot(lamda_list, n_DKKT + shift_n - n , label = 'n DKKT')



shift_n = mean(n) - mean(n_DKKT + ub0)
plot(lamda_list,  shift_n + n_DKKT + ub0 - n, label = 'ub0')


shift_n = mean(n) - mean(n_DKKT + ub1)
plot(lamda_list,  shift_n + n_DKKT + ub1 - n, label = 'ub1')

c0 = 1.0
c1 = 0.5

shift_n = mean(n) - mean(n_DKKT + c0*ub0 + c1*ub1)
plot(lamda_list , shift_n + n_DKKT + c0*ub0 + c1*ub1 - n, label = '%.2f*ub0+ %.2f*ub1'%(c0,c1) )

if force_mix:

	h = 1
	def diff_array(x):


		fc0      = x[0]
		fc1      = x[1]
		fshift_n = x[2]

		the_diff = n - (fshift_n + n_DKKT + fc0*ub0 + fc1*ub1)

		return the_diff[:-h]

	from scipy.optimize  import least_squares

	solution = least_squares(diff_array, x0 = (0.00, 0.00, 0.01)).x
	plot(lamda_list[:-h],  diff_array(solution) , label = 'force fit')

legend()
grid(True)


figure()

plot(lamda_list , n, label = 'n')

shift_n = mean(n) - mean(n_DKKT)
plot(lamda_list , n_DKKT +shift_n, label = 'n DKKT')


shift_n = mean(n) - mean(n_DKKT + ub0)
plot(lamda_list , shift_n + n_DKKT + ub0 , label = 'ub0')

shift_n = mean(n) - mean(n_DKKT + ub1)
plot(lamda_list , shift_n + n_DKKT + ub1 , label = 'ub1')

shift_n = mean(n) - mean(n_DKKT + c0*ub0 + c1*ub1)
plot(lamda_list , shift_n + n_DKKT + c0*ub0 + c1*ub1, label = '%.2f*ub0+ %.2f*ub1'%(c0,c1))
print (shift_n, c0, c1)


if force_mix:
	fc0 = solution[0]
	fc1 = solution[1]
	shift_n = solution[2]
	print (shift_n, fc0, fc1)
	plot(lamda_list , shift_n + n_DKKT + fc0*ub0 + fc1*ub1, label = 'force_fit')


plot(lamda_list , k, label = 'k')


grid(True)
legend()


my_compute_pool.terminate()
my_compute_pool.close()
show()
