

#def try_mkdir(direct):
#	from os import mkdir
#	try:
#		mkdir(direct)
#	except:
#		pass

def try_mkdir(direct):
	from os import mkdir
	from os.path import isdir
	if isdir(direct) == False:
		mkdir(direct)



def extrap(lamda, n, kind = 'linear'):
	'''Requires that lamda be in increasing order'''
	upper_value = n[-1]
	lower_value = n[0]
	from scipy.interpolate import interp1d
	interp1d.upper_bound = 0.0
	interp1d.lower_bound = 0.0
	def is_in_bounds(self,lamda):
		if (self.lower_bound <= lamda) and (lamda <= self.upper_bound):
			return True
		else:
			return False
	interp1d.is_in_bounds = is_in_bounds
	# now we instantiate
	func = interp1d(lamda, n, kind=kind, bounds_error = False, fill_value = (lower_value, upper_value))
	func.upper_bound = max(lamda)
	func.lower_bound = min(lamda)
	return func

def extrap_c(lamda, nk, kind = 'linear'):
	'''This is deprecated, old versions interp1d didn't always play nice with complex values and extrap is now the backend'''
	return extrap(lamda = lamda, n = nk, kind = kind)


def functionize_nk_file(file_name, skiprows = 0, kind = 'linear'):
	from numpy import loadtxt, argsort
	data = loadtxt(file_name, skiprows = skiprows).T

	indexing = argsort(data[0]) # creates lookup map of indices
	lamda = data[0][indexing] # this remaps to this order
	n = data[1][indexing]
	k = data[2][indexing]

	return extrap_c(lamda, n+1.0j*k, kind = kind)



def functionize_frequency_and_permittivity_file(file_name, skiprows = 0, kind = 'linear'):
	from numpy import loadtxt, absolute, sqrt, argsort
	data = loadtxt(file_name,skiprows = skiprows).T

	c = 299792458
	lamda = c/ data[0] * 1e9

	metal_epsilon = data[1]+1.0j*data[2]

	n = sqrt(( absolute(metal_epsilon) + metal_epsilon.real)/2.0 )
	k = sqrt(( absolute(metal_epsilon) - metal_epsilon.real)/2.0 )

	indexing = argsort(lamda)
	metal_nf = extrap(lamda[indexing], n[indexing], kind = kind)
	metal_kf = extrap(lamda[indexing], k[indexing], kind = kind)

	def metal_nk_f(lamda):
		return metal_nf(lamda) + 1.0j* metal_kf(lamda)

	return metal_nk_f




def nk_plot( nkf, lamda_fine, lamda_list = [], title_string = '', file_name = 'nk.pdf', show_nodes = False, zoom_window = [], show_plots = False):

	if show_plots == False:
		from matplotlib import use
		use('Agg')
		from matplotlib.pylab import ioff
		ioff()


	from matplotlib.pylab import figure, gca, subplots_adjust, tight_layout, show, savefig
	n_color = 'teal'
	k_color = 'orange'

	lamda_nodes = lamda_list
	lamda_smooth = lamda_fine

	fig = figure(figsize =(3.2,2.5), dpi = 220*2.0/3.0 )
	nax = gca()
	kax = nax.twinx()

	if show_nodes: nk_points =  nkf(lamda_nodes)
	nk_smooth =  nkf(lamda_smooth)


	if show_nodes: nax.plot(lamda_nodes, nk_points.real, linewidth = 2.0, color = 'k',   linestyle = '', marker = 'o', markersize = .4, zorder = 1)
	nax.plot(lamda_smooth, nk_smooth.real,   linewidth = 2.0, color = n_color, label ='n', linestyle = '-', zorder = 0.9)


	if show_nodes: kax.plot(lamda_nodes, nk_points.imag, linewidth = 2.0, color = 'k',   linestyle = '',   marker = 'o', markersize = .4, zorder = 1)
	kax.plot(lamda_smooth, nk_smooth.imag,   linewidth = 2.0, color = k_color, label ='k', linestyle = '-', zorder = 0.9)
	#subplots_adjust(left = 0.19, bottom = 0.18, right = 0.81, top = 0.96)

	nax.minorticks_on()


	kax.minorticks_on()
	nax.tick_params(axis='y',which = 'both', colors=n_color)
	kax.tick_params(axis='y',which = 'both', colors=k_color)
	nax.tick_params(axis='x',which = 'both',labelsize=10)

	nax.tick_params(axis='y',which = 'both',labelsize=10)
	kax.tick_params(axis='y',which = 'both',labelsize=10)

	nax.set_ylabel('n')
	kax.set_ylabel('k')

	nax.yaxis.label.set_color(n_color)
	kax.yaxis.label.set_color(k_color)

	nax.set_xlabel('Wavelength (nm)')

	nax.spines['left'].set_color(n_color)
	kax.spines['right'].set_color(k_color)
	kax.spines['left'].set_visible(False)
	nax.spines['right'].set_visible(False)


	nax.set_title(title_string , fontsize = 10)

	subplots_adjust(top = 0.99,
			bottom = 0.18,
			left = 0.15,
			right = 0.84)



	#fname = '%s_%s'%(simple_name,structure_name)
	### fix limits
	ylim =nax.get_ylim()
	nax.set_ylim(0,ylim[1])
	ylim =kax.get_ylim()
	kax.set_ylim(0,ylim[1])

	if zoom_window == []:
		nax.set_xlim(min(lamda_smooth),max(lamda_smooth))
	else:
		nax.set_xlim(zoom_window[0],zoom_window[1])

	tight_layout(pad = 0.1)
	savefig(file_name, dpi = 600,transparent = True)

	if show_plots:
		show()

	return fig




def error_plot(lamda_list, rms_spectrum,
				adaptation_threshold, adaptation_threshold_min, adaptation_threshold_max,
				reducible_error_spectrum = [],
				lamda_fine = [], rms_spectrum_fine = [], reducible_error_spectrum_fine = [],
				title_string = '',
				file_name = 'error_map.pdf', zoom_window = [], show_plots = False ):


	if show_plots == False:
		from matplotlib import use
		use('Agg')
		from matplotlib.pylab import ioff
		ioff()


	from matplotlib.pylab import figure, plot, axhline, ylabel, xlabel, title, minorticks_on, gca, gcf, savefig, show
	from numpy import array, mean, sqrt

	fig = figure(figsize=(3.2,2.5),dpi = 220*2/3)
	if len(rms_spectrum_fine)!=0:
		plot(lamda_fine, array(rms_spectrum_fine)*100,  color ='grey' )
		plot(lamda_list, array(rms_spectrum)*100,marker = 'o',markersize = 2, markerfacecolor = 'None', color ='grey', linestyle = '' )
	else:
		plot(lamda_list, array(rms_spectrum)*100,marker = 'o',markersize = 2, markerfacecolor = 'None', color ='grey', linestyle = '-' )

	if len(reducible_error_spectrum) != 0 : # check that it isn't null
		if len(reducible_error_spectrum_fine) != 0 : # lamda_fine can be null too...
			plot(lamda_fine, array(reducible_error_spectrum_fine)*100, color = 'royalblue' )
			plot(lamda_list, array(reducible_error_spectrum)*100, marker = 'o',markersize = 2, color = 'royalblue', linestyle = '' )
		else:
			plot(lamda_list, array(reducible_error_spectrum)*100, marker = 'o',markersize = 2, color = 'royalblue', linestyle = '-' )

	net_rms = sqrt( mean( array(rms_spectrum)**2 ) )
	axhline(net_rms*100,color='grey',linestyle=':', zorder = 4)

	axhline(adaptation_threshold*100,color='k',linestyle='--', zorder = 5)

	ylabel('RMS Error (%)', fontsize = 10)
	xlabel('Wavelength (nm)',fontsize = 10)
	title(title_string, loc = 'right', fontsize = 10)
	minorticks_on()
	gca().tick_params(axis='both', which='major', labelsize=10)

	ymax = gca().get_ylim()[1]
	max_rms = max(rms_spectrum)

	if max_rms > adaptation_threshold_max:
		axhline(adaptation_threshold_max*100,color='pink',linestyle=':', zorder = 4)
	axhline(adaptation_threshold_min*100,color='lightgreen',linestyle=':', zorder = 4)
	gca().set_ylim(0, ymax)

	if zoom_window == []:
		gca().set_xlim(min(lamda_list),max(lamda_list))
	else:
		gca().set_xlim(zoom_window[0],zoom_window[1])

	gcf().tight_layout(pad=0.1)
	savefig(file_name ,dpi=600, transparent = True)

	if show_plots:
		show()

	return fig
