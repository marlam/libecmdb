/*
 * Copyright (C) 2012
 * Computer Graphics Group, University of Siegen, Germany.
 * Written by Martin Lambers <martin.lambers@uni-siegen.de>.
 * See http://www.cg.informatik.uni-siegen.de/ for contact information.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef JPEGBLOB_H
#define JPEGBLOB_H

void jpegblob_enc(
        unsigned int width,
        unsigned int height,
        unsigned int components /* 1 or 3 */,
        const unsigned char* data,
        int quality /* 1 to 100 */,
        unsigned char* blob /* buffer at least as large as source buffer */,
        size_t* blob_size);

void jpegblob_dec(
        unsigned int width,
        unsigned int height,
        unsigned int components /* 1 or 3 */,
        const unsigned char* blob,
        size_t blob_size,
        unsigned char* data /* buffer large enough for image data */);

#endif
