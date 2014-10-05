/*
 * cxvec_math.c
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

/*! \addtogroup cxvec_math
 *  @{
 */

/*! \file cxvec_math.c
 *  \brief Osmocom Complex vectors math implementation
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <osmocom/dsp/cxvec.h>
#include <osmocom/dsp/cxvec_math.h>

/*! \brief Scale a complex vector (multiply by a constant)
 *  \param[in] in Input complex vector
 *  \param[in] scale Factor to apply to each sample
 *  \param[out] out Output complex vector
 *  \returns The output complex vector (or NULL if error)
 *
 * \f$out(k) = in(k) \cdot scale\f$
 *
 * The output vector parameter 'out' can be NULL to allocate a new
 * vector, or can be equal to the 'in' input vector to perform the
 * transform in-place. If it's different, it must be long enough
 * to contain the result (i.e. in->len)
 */
struct osmo_cxvec *
osmo_cxvec_scale(const struct osmo_cxvec *in, float complex scale,
                 struct osmo_cxvec *out)
{
	int i;
	int seq_real;

	seq_real = !!(in->flags & CXVEC_FLG_REAL_ONLY);

	if (!out)
		out = osmo_cxvec_alloc(in->len);
	else if (out->max_len < in->len)
		return NULL;

	if (cimagf(scale) == 0.0f)
	{
		float scalef = crealf(scale);

		if (seq_real) {
			for (i=0; i<in->len; i++)
				out->data[i] = crealf(in->data[i]) * scalef;
		} else {
			for (i=0; i<in->len; i++)
				out->data[i] = in->data[i] * scalef;
		}
	}
	else
	{
		if (seq_real) {
			for (i=0; i<in->len; i++)
				out->data[i] = crealf(in->data[i]) * scale;
		} else {
			for (i=0; i<in->len; i++)
				out->data[i] = in->data[i] * scale;
		}

		out->flags &= ~CXVEC_FLG_REAL_ONLY;
	}

	out->len = in->len;

	return out;
}

/*! \brief Rotate a complex vector (frequency shift)
 *  \param[in] in Input complex vector
 *  \param[in] rps Rotation to apply in radian per sample
 *  \param[out] out Output complex vector
 *  \returns The output complex vector (or NULL if error)
 *
 * \f$out(k) = in(k) \cdot e^{j \cdot rps \cdot k}\f$
 *
 * The output vector parameter 'out' can be NULL to allocate a new
 * vector, or can be equal to the 'in' input vector to perform the
 * transform in-place. If it's different, it must be long enough
 * to contain the result (i.e. in->len)
 */
struct osmo_cxvec *
osmo_cxvec_rotate(const struct osmo_cxvec *in, float rps,
                  struct osmo_cxvec *out)
{
	int i;

	if (!out)
		out = osmo_cxvec_alloc(in->len);
	else if (out->max_len < in->len)
		return NULL;

	for (i=0; i<in->len; i++)
		out->data[i] = in->data[i] * cexpf(I*(rps*i));

	out->len = in->len;
	out->flags &= ~CXVEC_FLG_REAL_ONLY;

	return out;
}

/*! \brief Fractionally delay a vector while maintaining its length
 *  \param[in] in Input complex vector
 *  \param[in] delay The fractional delay to apply
 *  \param[out] out Output complex vector
 *  \returns The output complex vector (or NULL if error)
 *
 * The output always has the same length. Samples pushed out by
 * the delays are lost and new ones filled with zeroes are pushed in.
 *
 * The output vector parameter 'out' can be NULL to allocate a new
 * vector, or can be equal to the 'in' input vector to perform the
 * transform in-place. If it's different, it must be long enough
 * to contain the result (i.e. in->len)
 */
struct osmo_cxvec *
osmo_cxvec_delay(const struct osmo_cxvec *in, float delay,
                 struct osmo_cxvec *out)
{
	int ofs_int = (int) roundf(delay);
	float ofs_frac = delay - ofs_int;
	const struct osmo_cxvec *shifted_vect = NULL;
	int i, j;

	/* Get output vector */
	if (!out)
		out = osmo_cxvec_alloc(in->len);
	else if (out->max_len < in->len)
		return NULL;

	/* Set output length / flags */
	out->len = in->len;
	out->flags = in->flags;

	/* Fractional offset (if reasonable) */
	if (fabs(ofs_frac) > 0.05f) {
		float complex _d[21];
		struct osmo_cxvec _sinc_vect, *sinc_vect = &_sinc_vect;

		/* Create sinc vector */
		osmo_cxvec_init_from_data(sinc_vect, _d, 21);

		for (i=0; i<21; i++)
			sinc_vect->data[i] = osmo_sinc(M_PIf * (i - 10.0f - ofs_frac));

		sinc_vect->flags |= CXVEC_FLG_REAL_ONLY;

		/* Convolve */
		shifted_vect = osmo_cxvec_convolve(sinc_vect, in, CONV_NO_DELAY, NULL);
	}

	if (!shifted_vect)		/* Also covers failure of convolve ... */
		shifted_vect = in;

	/* Integer offset */
	if (ofs_int < 0) {
		ofs_int = - ofs_int;
		for (i=0; i<(shifted_vect->len-ofs_int); i++)
			out->data[i] = shifted_vect->data[i+ofs_int];
		for (; i<in->len; i++)
			out->data[i] = 0.0f;
	} else {
		for (i=in->len-1,j=shifted_vect->len-1; i>=ofs_int; i--, j--)
			out->data[i] = shifted_vect->data[j-ofs_int];
		for (; i>=0; i--)
			out->data[i] = 0.0f;
	}

	/* Release */
	if (in != shifted_vect)
		osmo_cxvec_free((struct osmo_cxvec *)shifted_vect);

	return out;
}

/*! \brief Convolve two complex vectors
 *  \param[in] f First input complex vector
 *  \param[in] g Second input complex vector
 *  \param[in] type The convolution span type
 *  \param[out] out Output complex vector
 *  \returns The output complex vector (or NULL if error)
 *
 * The convolution of discrete sequences is defined as :
 *
 * \f$(f * g)[n] = \sum_{m=-\infty}^{\infty} f[m] \; g[n-m]\f$
 *
 * Altough the mathematical operation is commutative, the constraint
 * of implementation limit this method. Depending on the type of span
 * chosen, it might not be and it's always recommended that 'g' be the longer
 * sequence. It should not be much of a limitation when this methos is used
 * for filtering or pulseshaping : use 'f' as the filter and 'g' as the
 * signal.
 *
 * The output vector parameter 'out' can be NULL to allocate a new
 * vector. If it's not NULL, it must be long enough to contain the result
 * (length depends on the exact convolution type)
 */
struct osmo_cxvec *
osmo_cxvec_convolve(const struct osmo_cxvec *f, const struct osmo_cxvec *g,
                    enum osmo_cxvec_conv_type type, struct osmo_cxvec *out)
{
	int Lf, Lg, Lo, si, ei, i, jf, jg;
	int f_real, g_real;

	if (!f || !g)
		return NULL;

	f_real = !!(f->flags & CXVEC_FLG_REAL_ONLY);
	g_real = !!(g->flags & CXVEC_FLG_REAL_ONLY);

	/* index / size */
	Lf = f->len;
	Lg = g->len;

	switch (type) {
	case CONV_FULL_SPAN:
		si = 0;
		Lo = Lf + Lg - 1;
		break;
	case CONV_OVERLAP_ONLY:
		si = Lf;
		Lo = abs(Lf-Lg) + 1;
		break;
	case CONV_NO_DELAY:
		si = (Lf >> 1) - ((Lf & 1) ^ 1);
		Lo = Lg;
		break;
	default:
		return NULL;
	}

	ei = si + Lo;

	/* Output vector */
	if (!out)
		out = osmo_cxvec_alloc(Lo);
	else if (out->max_len < Lo)
		return NULL;

	out->flags = 0;

	/* Do the math */
	if (f_real && g_real) {
		for (i=si; i<ei; i++) {
			float sum = 0.0f;
			for (jf=0,jg=i; jf<=Lf && jg>=0; jf++,jg--)
				if (jg < Lg)
					sum += crealf(f->data[jf]) * crealf(g->data[jg]);
			out->data[i-si] = sum;
		}
		out->flags |= CXVEC_FLG_REAL_ONLY;
	} else if (f_real) {
		for (i=si; i<ei; i++) {
			float complex sum = 0.0f;
			for (jf=0,jg=i; jf<Lf && jg>=0; jf++,jg--)
				if (jg < Lg)
					sum += crealf(f->data[jf]) * g->data[jg];
			out->data[i-si] = sum;
		}
	} else if (g_real) {
		for (i=si; i<ei; i++) {
			float complex sum = 0.0f;
			for (jf=0,jg=i; jf<Lf && jg>=0; jf++,jg--)
				if (jg < Lg)
					sum += f->data[jf] * crealf(g->data[jg]);
			out->data[i-si] = sum;
		}
	} else {
		for (i=si; i<ei; i++) {
			float complex sum = 0.0f;
			for (jf=0,jg=i; jf<Lf && jg>=0; jf++,jg--)
				if (jg < Lg)
					sum += f->data[jf] * g->data[jg];
			out->data[i-si] = sum;
		}
	}

	out->len = Lo;

	return out;
}

/*! \brief Cross-correlate two complex vectors
 *  \param[in] f First input complex vector
 *  \param[in] g Second input complex vector
 *  \param[in] g_corr_step Allow for oversampling of 'g' compared to 'f'
 *  \param[out] out Output complex vector
 *  \returns The output complex vector (or NULL if error)
 *
 * The cross-correlation of discrete sequences is defined as :
 *
 * \f$(f \star g)[n] = \sum_{m=-\infty}^{\infty} f^*[m] \; g[n+m]\f$
 *
 * In this implementation, the output vector will be for every n value
 * between 0 and (g->len - f->len + 1). This assumes that g is the longer
 * sequence and we 'fit' f at every positition inside it.
 *
 * With the parameter g_corr_step, it's also possible to have a g sequence
 * that is oversampled with regard to f. (if g_corr_step > 1)
 *
 * The output vector parameter 'out' can be NULL to allocate a new
 * vector. If it's not NULL, it must be long enough to contain the result
 * (i.e. g->len - f->len + 1)
 */
struct osmo_cxvec *
osmo_cxvec_correlate(const struct osmo_cxvec *f, const struct osmo_cxvec *g,
                     int g_corr_step, struct osmo_cxvec *out)
{
	int l, m, n, mn;
	int f_real, g_real;

	f_real = !!(f->flags & CXVEC_FLG_REAL_ONLY);
	g_real = !!(g->flags & CXVEC_FLG_REAL_ONLY);

	l = g->len - (f->len * g_corr_step) + 1;

	if (l < 0)
		return NULL;

	if (!out)
		out = osmo_cxvec_alloc(l);
	else if (out->max_len < l)
		return NULL;

	out->flags = 0;

	if (f_real && g_real) {
		for (m=0; m<l; m++) {
			float v = 0.0f;
			for (n=0,mn=m; n<f->len; n++,mn+=g_corr_step)
				v += crealf(f->data[n]) * crealf(g->data[mn]);
			out->data[m] = v;
		}
		out->flags |= CXVEC_FLG_REAL_ONLY;
	} else if (f_real) {
		for (m=0; m<l; m++) {
			complex float v = 0.0f;
			for (n=0,mn=m; n<f->len; n++,mn+=g_corr_step)
				v += crealf(f->data[n]) * g->data[mn];
			out->data[m] = v;
		}
	} else if (g_real) {
		for (m=0; m<l; m++) {
			complex float v = 0.0f;
			for (n=0,mn=m; n<f->len; n++,mn+=g_corr_step)
				v += f->data[n] * crealf(g->data[mn]);
			out->data[m] = conjf(v);
		}
	} else {
		for (m=0; m<l; m++) {
			complex float v = 0.0f;
			for (n=0,mn=m; n<f->len; n++,mn+=g_corr_step)
				v += conjf(f->data[n]) * g->data[mn];
			out->data[m] = v;
		}
	}

	out->len = l;

	return out;
}

/*! \brief Interpolate any fractional position in a vector using sinc filtering
 *  \param[in] cv Input complex vector
 *  \param[in] pos Position to interpolate
 *
 * pos must be >= 0 and < cv->len
 */
float complex
osmo_cxvec_interpolate_point(const struct osmo_cxvec *cv, float pos)
{
	const int N = 10;
	int b, e, i;
	float complex val;

	/* Index */
	i = (int)(floorf(pos));
	b = i - N;
	e = i + N + 1;

	if (b < 0)
		b = 0;

	if (e >= cv->len)
		e = cv->len - 1;

	/* Interpolate */
	if (cv->flags & CXVEC_FLG_REAL_ONLY) {
		float valf = 0.0f;
		for (i=b; i<e; i++)
			valf += crealf(cv->data[i]) * osmo_sinc(M_PIf * (i - pos));
		val = valf;
	} else {
		val = 0.0f;
		for (i=b; i<e; i++)
			val += cv->data[i] * osmo_sinc(M_PIf * (i - pos));
	}

	return val;
}

/*! \brief Find the maximum energy (\f$|x|^2\f$) peak in a sequence
 *  \param[in] cv Input complex vector
 *  \param[in] win_size Size of the window (for algorithms using windows)
 *  \param[in] alg Peak detection algorithm to use
 *  \param[out] peak_val_p Returns interpolated peak value if non-NULL
 *  \returns Peak position with sub-sample accuracy
 */
float
osmo_cxvec_peak_energy_find(const struct osmo_cxvec *cv, int win_size,
                            enum osmo_cxvec_peak_alg alg,
                            float complex *peak_val_p)
{
	float val, max_val;
	int idx, max_idx, hi;
	float he[win_size];
	float peak_pos = 0.0f;

	/* Safety */
	if (cv->len < win_size)
		win_size = cv->len;

	/* Scan for the window */
		/* State init */
	val = 0.0f;
	max_val = 0.0f;
	max_idx = 0;

		/* Prefill energy history array */
	for (hi=0; hi<win_size; hi++)
		he[hi] =  osmo_normsqf(cv->data[hi]);

		/* Main scan */
	for (idx=0; idx<cv->len-win_size; idx++)
	{
		hi = idx % win_size;

		val -= he[hi];
		he[hi] = osmo_normsqf(cv->data[idx+win_size]);
		val += he[hi];

		if (val > max_val) {
			max_val = val;
			max_idx = idx + 1;
		}
	}

	/* Find maximum peak within the window */
	/* (for PEAK_WEIGH_WIN_CENTER & PEAK_EARLY_LATE */
	if (alg == PEAK_WEIGH_WIN_CENTER || alg == PEAK_EARLY_LATE)
	{
		int mwi = 0;
		float mwv = 0.0f;

		for (idx=max_idx; idx<(max_idx+win_size); idx++) {
			val = osmo_normsqf(cv->data[idx]);
			if (val > mwv) {
				mwv = val;
				mwi = idx;
			}
		}

		if (alg == PEAK_WEIGH_WIN_CENTER) {
			max_idx = mwi - (win_size >> 1);

			if (max_idx < 0)
				max_idx = 0;
			if (max_idx > (cv->len - win_size - 1))
				max_idx = cv->len - win_size - 1;
		} else {
			max_idx = mwi;
		}
	}

	/* Find the fractional position */
	if (alg == PEAK_WEIGH_WIN || alg == PEAK_WEIGH_WIN_CENTER)
	{
		float wes = 0.0f;
		float es = 0.0f;

		for (idx=max_idx; idx<(max_idx+win_size); idx++) {
			val = osmo_normsqf(cv->data[idx]);
			wes += val * idx;
			es += val;
		}

		peak_pos = wes / es;
	}
	else if (alg == PEAK_EARLY_LATE)
	{
		float early_idx = max_idx - 1.0f;
		float late_idx  = max_idx + 1.0f;
		float complex early_pt;
		float complex late_pt;
		float incr = 0.5f;

		while (incr > (1.0f/1024.0f))
		{
			early_pt = osmo_cxvec_interpolate_point(cv, early_idx);
			late_pt  = osmo_cxvec_interpolate_point(cv, late_idx);

			if (osmo_normsqf(early_pt) < osmo_normsqf(late_pt))
				early_idx += incr;
			else if (osmo_normsqf(early_pt) > osmo_normsqf(late_pt))
				early_idx -= incr;
			else
				break;

			incr /= 2.0f;
			late_idx = early_idx + 2.0f;
		}

		peak_pos = early_idx + 1.0f;
	}

	/* Interpolate peak (if asked to) */
	if (peak_val_p)
		*peak_val_p = osmo_cxvec_interpolate_point(cv, peak_pos);

	return peak_pos;
}

/*! \brief 'Normalize' an IQ signal and apply decimation/frequency shift
 *  \param[in] sig Input complex signal
 *  \param[in] decim Decimation factor
 *  \param[in] freq_shift Frequency shift in radian per output sample
 *  \param[out] out Output complex vector
 *  \returns The output complex vector (or NULL if error)
 *
 * The operation performed are DC removal, amplitude normalization (divided
 * by the standard deviation), decimation, frequency shift.
 *
 * The output vector parameter 'out' can be NULL to allocate a new
 * vector, or can be equal to the 'in' input vector to perform the
 * transform in-place. If it's different, it must be long enough to contain
 * the result (i.e. (sig->len + decim - 1) / decim)
 */
struct osmo_cxvec *
osmo_cxvec_sig_normalize(const struct osmo_cxvec *sig,
                         int decim, float freq_shift,
                         struct osmo_cxvec *out)
{
	float complex avg = 0.0f;
	float sigma = 0.0f, stddev;
	int l, i, j;

	l = sig->len / decim;

	if (!out)
		out = osmo_cxvec_alloc(l);
	else if (out->max_len < l)
		return NULL;

	for (i=0; i<sig->len; i++)
		avg += sig->data[i];
	avg /= sig->len;

	for (i=0; i<sig->len; i++)
		sigma += osmo_normsqf(sig->data[i] - avg);
	sigma /= sig->len;

	stddev = sqrtf(sigma);
	if (stddev == 0.0f)
		stddev = 1.0f;	/* Safety against constant signals */

	for (i=0, j=0; i<l; i++,j+=decim)
		out->data[i] = (sig->data[j] - avg) / stddev;

	out->len = l;
	out->flags = sig->flags;

	if (freq_shift != 0.0f)
		for (i=0; i<out->len; i++)
			out->data[i] *= cexpf( I * (freq_shift * i) );

	return out;
}

/*! @} */
