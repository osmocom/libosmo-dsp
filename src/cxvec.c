/*
 * cxvec.c
 *
 * Complex vectors handling
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

/*! \addtogroup cxvec
 *  @{
 */

/*! \file cxvec.c
 *  \brief Osmocom Complex vectors implementation
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <osmocom/dsp/cxvec.h>

/*! \brief Initialize a vector structure with a given data array
 *  \param[out] cv The vector to be initialized
 *  \param[in] data Pointer to the complex data array
 *  \param[in] len Number of complex samples
 *
 *  The data is not copied, it is just referenced.
 */
void
osmo_cxvec_init_from_data(struct osmo_cxvec *cv,
                          float complex *data, int len)
{
	cv->len   = cv->max_len = len;
	cv->flags = 0;
	cv->data  = data;
}

/*! \brief Allocate a complex vector referencing a given data array
 *  \param[in] data Pointer to the complex data array
 *  \param[in] len Number of complex samples
 *
 *  The data is not copied, it is just referenced.
 */
struct osmo_cxvec *
osmo_cxvec_alloc_from_data(float complex *data, int len)
{
	struct osmo_cxvec *cv;

	cv = malloc(sizeof(struct osmo_cxvec));
	if (!cv)
		return NULL;

	osmo_cxvec_init_from_data(cv, data, len);

	return cv;
}

/*! \brief Allocate a complex vector of a given maximum length
 *  \param[in] max_len Maximum length of data
 *
 * Data array is allocated along with the structure, but is uninitialized.
 * Length is set to 0.
 */
struct osmo_cxvec *
osmo_cxvec_alloc(int max_len)
{
	struct osmo_cxvec *cv;

	cv = malloc(sizeof(struct osmo_cxvec) + max_len * sizeof(float complex));
	if (!cv)
		return NULL;

	cv->len = 0;
	cv->max_len = max_len;
	cv->flags = 0;
	cv->data = &cv->_data[0];

	return cv;
}

/*! \brief Free a complex vector (and possibly associated data)
 *  \param[in] cv Complex vector to free
 *
 * Notes: - Can be safely called with NULL
 *        - If the data was allocated with the vector using
 *          \ref osmo_cxvec_alloc , it will be free as well. If the
 *          data was pre-existing ( \ref osmo_cxvec_init_from_data or
 *          \ref osmo_cxvec_alloc_from_data ) it will not be free'd.
 */
void
osmo_cxvec_free(struct osmo_cxvec *cv)
{
	free(cv);
}

/*! \brief Save the data contained of a vector into a .cfile for debug
 *  \param[in] cv Complex vector to save
 *  \param[in] fname Filename to save the data to
 */
void
osmo_cxvec_dbg_dump(struct osmo_cxvec *cv, const char *fname)
{
	FILE *f = fopen(fname, "wb");
	int rv;
	if (!f)
		return;
	rv = fwrite(cv->data, sizeof(float complex), cv->len, f);
	if (rv != cv->len)
		fprintf(stderr, "[!] osmo_cxvec_dbg_dump: fwrite failed !\n");
	fclose(f);
}

/*! @} */
