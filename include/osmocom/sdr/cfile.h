/*
 * cfile.h
 *
 * Helpers to read .cfile (complex samples from gnuradio)
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

#ifndef __OSMO_SDR_CFILE_H__
#define __OSMO_SDR_CFILE_H__

/*! \defgroup cfile .cfile helpers
 *  @{
 */

/*! \file cfile.h
 *  \brief Osmocom .cfile helpers header
 */

#include <complex.h>

/*! \brief Structure representing a currently mapped .cfile */
struct cfile {
	float complex *data;	/*!< \brief Data array (read only !) */
	unsigned int len;	/*!< \brief Length (in samples) of the data */
	size_t _blen;		/*!< \brief Length (in bytes) of the data */
};

struct cfile *cfile_load(const char *filename);
void cfile_release(struct cfile *cf);

/*! }@ */

#endif /* __OSMO_SDR_CFILE_H__ */
