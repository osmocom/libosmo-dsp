/*
 * cxvec_math.h
 *
 * Complex vectors math and signal processing
 *
 * Copyright (C) 2011  Sylvain Munaut <tnt@246tNt.com>
 *
 * All Rights Reserved
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#ifndef __OSMO_DSP_CXVEC_MATH_H__
#define __OSMO_DSP_CXVEC_MATH_H__

/*! \defgroup cxvec_math Complex vectors math and signal processing
 *  \ingroup cxvec
 *  @{
 */

/*! \file cxvec_math.h
 *  \brief Osmocom Complex vectors math header
 */

#include <complex.h>
#include <math.h>

#include <osmocom/dsp/cxvec.h>


	/* Generic math stuff */

#define M_PIf (3.14159265358979323846264338327f) /*!< \brief PI value float */

/*! \brief Unnormalized sinc function
 *  \param[in] x Value for which to compute the sinc function.
 *  \returns The sinc(x) value
 *
 *  The function is defined as \f$\frac{\sin(x)}{x}\f$
 */
static inline float
osmo_sinc(float x)
{
	if ((x >= 0.01f) || (x <= -0.01f)) return (sinf(x)/x);
	return 1.0f;
}

/*! \brief Squared norm of a given complex
 *  \param[in] c Complex number for which to compute the squared norm
 *  \returns \f$|c|^2\f$
 */
static inline float
osmo_normsqf(float complex c)
{
	return crealf(c) * crealf(c) + cimagf(c) * cimagf(c);
}


	/* Complex vector math */

struct osmo_cxvec *
osmo_cxvec_scale(const struct osmo_cxvec *in, float complex scale,
                 struct osmo_cxvec *out);

struct osmo_cxvec *
osmo_cxvec_rotate(const struct osmo_cxvec *in, float freq_shift,
                  struct osmo_cxvec *out);

struct osmo_cxvec *
osmo_cxvec_delay(const struct osmo_cxvec *v, float delay,
                 struct osmo_cxvec *out);

/*! \brief Various possible types of convolution span */
enum osmo_cxvec_conv_type {
	/*! \brief Full span (every possible overlap of f onto g) */
	CONV_FULL_SPAN,	
	/*! \brief Every possible full overlap of f onto g */
	CONV_OVERLAP_ONLY,
	/*! \brief Center f sequence on every g sample */
	CONV_NO_DELAY,
};

struct osmo_cxvec *
osmo_cxvec_convolve(const struct osmo_cxvec *f, const struct osmo_cxvec *g,
                    enum osmo_cxvec_conv_type type, struct osmo_cxvec *out);

struct osmo_cxvec *
osmo_cxvec_correlate(const struct osmo_cxvec *f, const struct osmo_cxvec *g,
                     int g_corr_step, struct osmo_cxvec *out);

float complex
osmo_cxvec_interpolate_point(const struct osmo_cxvec *cv, float pos);

int
osmo_cxvec_peaks_scan(const struct osmo_cxvec *cv, int *peaks_idx, int N);

/*! \brief Various possible peak finding algorithms */
enum osmo_cxvec_peak_alg {
	/*! \brief Weigthed position for the max pwr window */
	PEAK_WEIGH_WIN,	
	/*! \brief Weighted position of the peak centered window */
	PEAK_WEIGH_WIN_CENTER,
	/*! \brief Early-Late balancing around peak */
	PEAK_EARLY_LATE,
};

float
osmo_cxvec_peak_energy_find(const struct osmo_cxvec *cv, int win_size,
                            enum osmo_cxvec_peak_alg alg,
                            float complex *peak_val_p);

struct osmo_cxvec *
osmo_cxvec_sig_normalize(const struct osmo_cxvec *sig,
                         int decim, float freq_shift,
                         struct osmo_cxvec *out);

/*! @} */

#endif /* __OSMO_DSP_CXVEC_MATH_H__ */
