/*
 * iqbal.c
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

/*! \addtogroup iqbal
 *  @{
 */

/*! \file iqbal.c
 *  \brief IQ balance utils implementation
 */

#include <complex.h>
#include <stdlib.h>
#include <string.h>

#include <fftw3.h>

#include <osmocom/dsp/cxvec.h>
#include <osmocom/dsp/cxvec_math.h>
#include <osmocom/dsp/iqbal.h>


/* ------------------------------------------------------------------------ */
/* IQ balance correction and estimation                                     */
/* ------------------------------------------------------------------------ */

/*! \brief Apply IQ balance correction to a given complex buffer
 *  \param[out] out Complex output buffer
 *  \param[in] in Complex input buffer
 *  \param[in] len Number of complex samples to process
 *  \param[in] mag Magnitude correction (approximated)
 *  \param[in] phase Phase correction (approximated)
 *
 *  The input and output buffers can be the same for in-place modification.
 *
 *  The applied transform is out[i] = (a * (1 + mag)) + (b + phase * a) * i
 *  (with in[i] = a+bi).
 */
void
osmo_iqbal_fix(float complex *out, float complex *in, unsigned int len,
               float mag, float phase)
{
	int i;

	for (i=0; i<len; i++) {
		float complex v = in[i];
		out[i] = (crealf(v) * (1.0f + mag)) +
		         (cimagf(v) + phase * crealf(v)) * I;
	}
}

/*! \brief Apply IQ balance correction to a given complex vector
 *  \param[in] in Complex input vector
 *  \param[in] mag Magnitude correction (approximated)
 *  \param[in] phase Phase correction (approximated)
 *  \param[out] out Complex output vector (can be NULL or equal to 'in')
 *  \returns The output complex vector (or NULL if error)
 *
 *  If the 'out' parameter is NULL, a new vector will be allocated
 *  See \ref osmo_iqbal_fix for details of the correction applied.
 */
struct osmo_cxvec *
osmo_iqbal_cxvec_fix(const struct osmo_cxvec *in, float mag, float phase,
                     struct osmo_cxvec *out)
{
	if (!out)
		out = osmo_cxvec_alloc(in->len);

	if (!out || out->max_len < in->len)
		return NULL;

	osmo_iqbal_fix(out->data, in->data, in->len, mag, phase);

	out->len = in->len;
	out->flags = in->flags;

	return out;
}


/*! \brief Cache for \ref _osmo_iqbal_estimate when doing lots of calls */
struct _iqbal_estimate_state {
	float complex *fft;	/*!< \brief Temporary memory for FFT */
	fftwf_plan fft_plan;	/*!< \brief FFTW plan */
};

/*! \brief Release a cache object created by \ref _osmo_iqbal_estimate */
static void
_osmo_iqbal_estimate_release(struct _iqbal_estimate_state *state)
{
	if (!state)
		return;

	fftwf_destroy_plan(state->fft_plan);
	free(state->fft);

	free(state);
}

/*! \brief Objectively estimate IQ balance in a given complex buffer
 *  \param[in] data Input complex buffer (at least fft_size * fft_count samples)
 *  \param[in] fft_size Size of the FFT to use internally
 *  \param[in] fft_count The number of consecutive FFT to use internally
 *  \param[out] state_p Cache object for multiple calls (can be NULL)
 *  \returns A number >= 0.0f estimating the IQ balance (the lower, the better)
 *
 *  The Cache object should only be used for multiple calls with the same parameters
 *  and the same size of input vector. Once you don't plan on using it anymore,
 *  you should call \ref _osmo_iqbal_estimate_release . The initial pointer value
 *  should also be initialized to NULL.
 */
static float
_osmo_iqbal_estimate(const float complex *data, int fft_size, int fft_count,
                     struct _iqbal_estimate_state **state_p)
{
	float complex *fft;
	float est = 0.0f;
	fftwf_plan fft_plan;
	int i, j;

	if (state_p && *state_p) {
		fft = (*state_p)->fft;
		fft_plan = (*state_p)->fft_plan;
	} else {
		fft = malloc(sizeof(float complex) * fft_size);
		fft_plan = fftwf_plan_dft_1d(fft_size, fft, fft, FFTW_FORWARD, FFTW_ESTIMATE);
	}

	for (i=0; i<fft_count; i++)
	{
		float complex corr = 0.0f;

		memcpy(fft, &data[i*fft_size], sizeof(float complex) * fft_size);
		fftwf_execute(fft_plan);

		for (j=1; j<fft_size/2; j++)
			corr += fft[fft_size-j] * conjf(fft[j]);

		est += osmo_normsqf(corr); /* / (fft_size / 2); */
	}

	/* est /= fft_count; */

	if (state_p && !*state_p) {
		*state_p = malloc(sizeof(struct _iqbal_estimate_state));
		(*state_p)->fft = fft;
		(*state_p)->fft_plan = fft_plan;
	} else if (!state_p) {
		fftwf_destroy_plan(fft_plan);
		free(fft);
	}

	return est;
}

/*! \brief Objectively estimate IQ balance in a given complex buffer
 *  \param[in] data Input complex buffer (at least fft_size * fft_count samples)
 *  \param[in] fft_size Size of the FFT to use internally
 *  \param[in] fft_count The number of consecutive FFT to use internally
 *  \returns A number >= 0.0f estimating the IQ balance (the lower, the better)
 */
float
osmo_iqbal_estimate(const float complex *data, int fft_size, int fft_count)
{
	return _osmo_iqbal_estimate(data, fft_size, fft_count, NULL);
}

/*! \brief Objectively estimate IQ balance in a given complex vector
 *  \param[in] sig Input complex vector (at least fft_size * fft_count samples)
 *  \param[in] fft_size Size of the FFT to use internally
 *  \param[in] fft_count The number of consecutive FFT to use internally
 *  \returns A number >= 0.0f estimating the IQ balance (the lower, the better)
 */
float
osmo_iqbal_cxvec_estimate(const struct osmo_cxvec *sig,
                          int fft_size, int fft_count)
{
	if (sig->len < fft_size * fft_count)
		return -1.0f;

	return osmo_iqbal_estimate(sig->data, fft_size, fft_count);
}


/* ------------------------------------------------------------------------ */
/* IQ balance optimization                                                  */
/* ------------------------------------------------------------------------ */

/*! \brief Default values for the optimization algorithm */
const struct osmo_iqbal_opts osmo_iqbal_default_opts = {
	.fft_size	= 1024,
	.fft_count	= 8,
	.max_iter	= 25,
	.start_at_prev	= 1,
};

/*! \brief Internal state structure for the IQ balance optimization algorithm */
struct _iqbal_state
{
	const struct osmo_iqbal_opts *opts; /*!< \brief Options */
	const struct osmo_cxvec *org;	/*!< \brief Original vector */
	struct osmo_cxvec *tmp;		/*!< \brief Temporary vector */
	int feval;			/*!< \brief # of function evaluation */
	struct _iqbal_estimate_state *cache; /*!< \brief Cache for estimate func */
};

/*! \brief Optimization objective function - Value
 *  \param[in] state Current state object of optimization loop
 *  \param[in] x An array of 2 float for (mag,phase) point to evaluate at
 *  \returns The value of the objective function at point 'x'
 */
static inline float
_iqbal_objfn_value(struct _iqbal_state *state, float x[2])
{
	state->feval++;
	osmo_iqbal_cxvec_fix(state->org, x[0], x[1], state->tmp);
	return _osmo_iqbal_estimate(state->tmp->data,
		state->opts->fft_size, state->opts->fft_count,
		&state->cache);
}

/*! \brief Optimization objective function - Gradient estimation
 *  \param[in] state Current state object of optimization loop
 *  \param[in] x An array of 2 float for (mag,phase) point to evaluate at
 *  \param[in] v The value of the objective function at point 'x'
 *  \param[out] grad An array of 2 float for the estimated gradient at point 'x'
 */
static void
_iqbal_objfn_gradient(struct _iqbal_state *state, float x[2], float v, float grad[2])
{
	const float GRAD_STEP = 1e-6f;
	float xd[2], vd[2];

	xd[0] = x[0] + GRAD_STEP; xd[1] = x[1];
	vd[0] = _iqbal_objfn_value(state, xd);

	xd[0] = x[0]; xd[1] = x[1] + GRAD_STEP;
	vd[1] = _iqbal_objfn_value(state, xd);

	grad[0] = (vd[0] - v) / GRAD_STEP;
	grad[1] = (vd[1] - v) / GRAD_STEP;
}

/*! \brief Optimization objective function - Value & Gradient estimation
 *  \param[in] state Current state object of optimization loop
 *  \param[in] x An array of 2 float for (mag,phase) point to evaluate at
 *  \param[out] grad An array of 2 float for the estimated gradient at point 'x'
 *  \returns The value of the objective function at point 'x'
 */
static inline float
_iqbal_objfn_val_gradient(struct _iqbal_state *state, float x[2], float grad[2])
{
	float v = _iqbal_objfn_value(state, x);
	_iqbal_objfn_gradient(state, x, v, grad);
	return v;
}


/*! \brief Finds the best IQ balance correction parameters for a given signal
 *  \param[in] sig The input signal to optimize for
 *  \param[in,out] mag Magnitude correction (See \ref osmo_iqbal_fix)
 *  \param[in,out] phase Phase correction (See \ref osmo_iqbal_fix)
 *  \param[in] opts Options of the optimization process (See \ref osmo_iqbal_opts)
 *
 *  The mag and phase parameters are pointers to float. If in the options,
 *  the 'start_at_prev' is enabled, the initial values of those will be used
 *  and so they should be initialized appropriately.
 */
int
osmo_iqbal_cxvec_optimize(const struct osmo_cxvec *sig, float *mag, float *phase,
                          const struct osmo_iqbal_opts *opts)
{
	struct _iqbal_state _state, *state = &_state;
	float cv, nv, step;
	float cx[2], nx[2];
	float cgrad[2];
	float p;
	int i;

	if (!opts)
		opts = &osmo_iqbal_default_opts;

	if (sig->len < (opts->fft_size * opts->fft_count))
		return -1;

	state->org = sig;
	state->tmp = osmo_cxvec_alloc(sig->len);
	state->opts = opts;
	state->feval = 0;
	state->cache = NULL;

	if (opts->start_at_prev) {
		cx[0] = *mag;
		cx[1] = *phase;
	} else {
		cx[0] = 0.0f;
		cx[1] = 0.0f;
	}

	cv = _iqbal_objfn_val_gradient(state, cx, cgrad);
	step = cv / (fabs(cgrad[0]) + fabs(cgrad[1]));

	for (i=0; i<opts->max_iter; i++)
	{
		nx[0] = cx[0] - step * (cgrad[0] / (fabs(cgrad[0]) + fabs(cgrad[1])));
		nx[1] = cx[1] - step * (cgrad[1] / (fabs(cgrad[0]) + fabs(cgrad[1])));

		nv = _iqbal_objfn_value(state, nx);

		if (nv <= cv) {
			p = (cv - nv) / cv;

			cx[0] = nx[0];
			cx[1] = nx[1];
			cv = nv;
			_iqbal_objfn_gradient(state, cx, cv, cgrad);

			if (p < 0.01f)
				break;
		} else {
			step /= 2.0 * (nv / cv);
		}
	}

	osmo_cxvec_free(state->tmp);
	_osmo_iqbal_estimate_release(state->cache);

	*mag   = cx[0];
	*phase = cx[1];

	return 0;
}

/*! @} */
