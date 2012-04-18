/*
 * cxvec.h
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

#ifndef __OSMO_SDR_CXVEC_H__
#define __OSMO_SDR_CXVEC_H__

/*! \defgroup cxvec Complex vectors
 *  @{
 */

/*! \file cxvec.h
 *  \brief Osmocom Complex vectors header
 */

#include <complex.h>

#define CXVEC_FLG_REAL_ONLY	(1<<0)	/*!< \brief Real values only */

/*! \brief Complex vector */
struct osmo_cxvec {
	int len;		/*!< \brief Valid length */
	int max_len;		/*!< \brief Maximum length in data field */
	int flags;		/*!< \brief Flags, see CXVEC_FLG_xxx */
	float complex *data;	/*!< \brief Data field */
	float complex _data[0];	/*!< \brief Optional inline data array */
};

void
osmo_cxvec_init_from_data(struct osmo_cxvec *cv,
                          float complex *data, int len);

struct osmo_cxvec *
osmo_cxvec_alloc_from_data(float complex *data, int len);

struct osmo_cxvec *
osmo_cxvec_alloc(int max_len);

void
osmo_cxvec_free(struct osmo_cxvec *cv);

void
osmo_cxvec_dbg_dump(struct osmo_cxvec *cv, const char *fname);

/*! @} */

#endif /* __OSMO_SDR_CXVEC_H__ */
