/*
 * cfile.c
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

/*! \addtogroup cfile
 *  @{
 */

/*! \file cfile.c
 *  \brief Osmocom .cfile helpers implementation
 */

#include <complex.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <osmocom/dsp/cfile.h>


/*! \brief .cfile loader: mmap() the data into memory (read-only)
 *  \param[in] filename Filename of the .cfile to map
 *  \returns A structure desribing the mapped data
 */
struct cfile *
cfile_load(const char *filename)
{
	struct cfile *cf;
	struct stat st;
	int rv, fd = -1;

	cf = calloc(1, sizeof(struct cfile));
	if (!cf) {
		perror("calloc");
		goto err;
	}

	fd = open(filename, O_RDONLY);
	if (fd < 0) {
		perror("open");
		goto err;
	}

	rv = fstat(fd, &st);
	if (rv) {
		perror("stat");
		goto err;
	}

	cf->_blen = st.st_size;
	cf->len = cf->_blen / sizeof(float complex);

	cf->data = mmap(NULL, cf->_blen, PROT_READ, MAP_SHARED, fd, 0);
	if (!cf->data) {
		perror("mmap");
		goto err;
	}

	close(fd);

	return cf;

err:
	if (cf) {
		if (fd >= 0)
			close(fd);

		free(cf);
	}

	return NULL;
}

/*! \brief Release all resources associated with a mapped .cfile
 *  \param[in] cf Structure describing the cfile to unmap
 */
void
cfile_release(struct cfile *cf)
{
	int rv;

	rv = munmap(cf->data, cf->_blen);
	if (rv)
		perror("munmap");

	free(cf);
}

/*! @} */
