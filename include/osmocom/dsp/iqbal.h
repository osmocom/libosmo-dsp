/*
 * iqbal.h
 *
 * IQ balance correction / estimation utilities
 *
 * Copyright (C) 2013  Sylvain Munaut <tnt@246tNt.com>
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

#ifndef __OSMO_DSP_IQBAL_H__
#define __OSMO_DSP_IQBAL_H__

/*! \defgroup iqbal IQ balance utilities
 *  @{
 */

/*! \file iqbal.h
 *  \brief Osmocom IQ balance utils header
 */

#include <complex.h>

#include <osmocom/dsp/cxvec.h>


/* IQ balance correction and estimation */

void osmo_iqbal_fix(float complex *out, float complex *in, unsigned int len,
                    float mag, float phase);

struct osmo_cxvec *
osmo_iqbal_cxvec_fix(const struct osmo_cxvec *in, float mag, float phase,
                     struct osmo_cxvec *out);

float
osmo_iqbal_estimate(const float complex *data,
                    int fft_size, int fft_count);

float
osmo_iqbal_cxvec_estimate(const struct osmo_cxvec *sig,
                          int fft_size, int fft_count);


/* IQ balance optimization */ 

/*! \brief Processing options for the IQ balance optimization algorithm */
struct osmo_iqbal_opts {
	int fft_size;   	/*!< \brief FFT size to use */
	int fft_count;  	/*!< \brief Number of FFT to use */
	int max_iter;   	/*!< \brief Max # iterations per pass */
	int start_at_prev;	/*!< \brief Use prev values as starting point */
};

extern const struct osmo_iqbal_opts osmo_iqbal_default_opts;

int
osmo_iqbal_cxvec_optimize(const struct osmo_cxvec *sig, float *mag, float *phase,
                          const struct osmo_iqbal_opts *opts);

/*! @} */

#endif /* __OSMO_DSP_IQBAL_H__ */
