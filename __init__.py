

from .TRANK_lib import (fit_TRA_nk_sqr,
						fit_TRA_nk_sqr_KK_compliant,
						rms_TRA_error_spectrum,
						reducible_rms_TRA_error_spectrum,
						gradient_at_lamda,
						hessian_at_lamda,
						pointwise_TRA_error_sum_wrapper,
						TR,
						TRA_spectra,
						single_lamda_TRA_error_map,
						find_min_indices_2d_array )

from .TRANK_lib import (fit_spectra_nk_sqr,
						fit_spectra_nk_sqr_KK_compliant,
						rms_error_spectrum,
						reducible_rms_error_spectrum,
						sqr_rms_gradient_at_lamda,
						sqr_rms_hessian_at_lamda,
						pointwise_TRA_error_sum_wrapper,
						TMM_spectrum_wrapper,
						TMM_spectra,
						single_lamda_rms_error_map,
						find_min_indices_2d_array )


from .helpers import (extrap, extrap_c, functionize_nk_file, functionize_frequency_and_permittivity_file, try_mkdir, nk_plot, error_plot)

from .iterator import error_adaptive_iterative_fit
from .iterator import error_adaptive_iterative_fit_spectra
