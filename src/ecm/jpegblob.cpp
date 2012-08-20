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

#include "config.h"

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <stdexcept>

#include <setjmp.h>

#include <jpeglib.h>

#include "jpegblob.h"


/* libjpeg error callbacks */

struct my_error_mgr
{
    struct jpeg_error_mgr pub;
    jmp_buf setjmp_buffer;
};

static void my_error_exit(j_common_ptr cinfo)
{
    struct my_error_mgr* my_err = reinterpret_cast<struct my_error_mgr*>(cinfo->err);
    longjmp(my_err->setjmp_buffer, 1);
}

/* libjpeg compression memory management */

struct my_destination_mgr
{
    struct jpeg_destination_mgr pub;
    JOCTET* orig_buf;
    size_t orig_buf_size;
};

static void my_init_destination(j_compress_ptr cinfo)
{
    struct my_destination_mgr* my_dest = reinterpret_cast<struct my_destination_mgr*>(cinfo->dest);
    my_dest->pub.next_output_byte = my_dest->orig_buf;
    my_dest->pub.free_in_buffer = my_dest->orig_buf_size;
}

static boolean my_empty_output_buffer(j_compress_ptr /* cinfo */)
{
    // called when the buffer has no space left -- this must never happen.
    abort();
    return FALSE;
}

static void my_term_destination(j_compress_ptr /* cinfo */)
{
}

/* encoding */

void jpegblob_enc(
        unsigned int width,
        unsigned int height,
        unsigned int components /* 1 or 3 */,
        const unsigned char* data,
        int quality /* 1 to 100 */,
        unsigned char* blob /* buffer at least as large as source buffer */,
        size_t* blob_size)
{
    struct my_error_mgr jerr;
    struct jpeg_compress_struct cinfo;
    struct my_destination_mgr jdest;
    JSAMPROW jrow[1];

    // Setup error handling
    cinfo.err = jpeg_std_error(&jerr.pub);
    jerr.pub.error_exit = my_error_exit;
    if (setjmp(jerr.setjmp_buffer)) {
        /* libjpeg has signaled an error */
        char message[JMSG_LENGTH_MAX];
        (cinfo.err->format_message)(reinterpret_cast<jpeg_common_struct*>(&cinfo), message);
        jpeg_destroy_compress(&cinfo);
        throw std::runtime_error(message);
    }

    // Setup compression structure
    jpeg_create_compress(&cinfo);

    // Setup compression buffer management
    jdest.orig_buf = blob;
    jdest.orig_buf_size = width * height * components * sizeof(unsigned char);
    cinfo.dest = &jdest.pub;
    cinfo.dest->init_destination = my_init_destination;
    cinfo.dest->empty_output_buffer = my_empty_output_buffer;
    cinfo.dest->term_destination = my_term_destination;

    // Setup compression
    cinfo.image_width = width;
    cinfo.image_height = height;
    cinfo.input_components = components;
    cinfo.in_color_space = (components == 1 ? JCS_GRAYSCALE : JCS_RGB);
    jpeg_set_defaults(&cinfo);
    jpeg_set_quality(&cinfo, quality, TRUE);

    // Compress
    jpeg_start_compress(&cinfo, TRUE);
    jrow[0] = reinterpret_cast<JSAMPROW>(const_cast<unsigned char*>(data));
    size_t scanline_size = width * components * sizeof(unsigned char);
    while (cinfo.next_scanline < cinfo.image_height) {
        jpeg_write_scanlines(&cinfo, jrow, 1);
        jrow[0] += scanline_size;
    }
    jpeg_finish_compress(&cinfo);
    *blob_size = cinfo.dest->next_output_byte - blob;

    // Finish
    jpeg_destroy_compress(&cinfo);
}

/* libjpeg decompression memory management */

struct my_source_mgr
{
    struct jpeg_source_mgr pub;
    JOCTET* orig_buf;
    size_t orig_buf_size;
};

static void my_init_source(j_decompress_ptr cinfo)
{
    struct my_source_mgr* my_src = reinterpret_cast<struct my_source_mgr*>(cinfo->src);
    my_src->pub.next_input_byte = my_src->orig_buf;
    my_src->pub.bytes_in_buffer = my_src->orig_buf_size;
}

static boolean my_fill_input_buffer(j_decompress_ptr /* cinfo */)
{
    // called when buffer runs empty -- this must never happen
    abort();
    return FALSE;
}

static void my_skip_input_data(j_decompress_ptr cinfo, INT32 num_bytes)
{
    if (num_bytes > 0) {
        if (static_cast<size_t>(num_bytes) > cinfo->src->bytes_in_buffer) {
            cinfo->src->next_input_byte = NULL;
            cinfo->src->bytes_in_buffer = 0;
        } else {
            cinfo->src->next_input_byte += num_bytes;
            cinfo->src->bytes_in_buffer -= num_bytes;
        }
    }
}

static void my_term_source(j_decompress_ptr /* cinfo */)
{
}

/* decoding */

void jpegblob_dec(
        unsigned int width,
        unsigned int height,
        unsigned int components /* 1 or 3 */,
        const unsigned char* blob,
        size_t blob_size,
        unsigned char* data /* buffer large enough for image data */)
{
    struct my_error_mgr jerr;
    struct jpeg_decompress_struct cinfo;
    struct my_source_mgr jsrc;
    JSAMPROW jrow[1];

    // Setup error handling
    cinfo.err = jpeg_std_error(&jerr.pub);
    jerr.pub.error_exit = my_error_exit;
    if (setjmp(jerr.setjmp_buffer)) {
        /* libjpeg has signaled an error */
        char message[JMSG_LENGTH_MAX];
        (cinfo.err->format_message)(reinterpret_cast<jpeg_common_struct*>(&cinfo), message);
        jpeg_destroy_decompress(&cinfo);
        throw std::runtime_error(message);
    }

    // Setup decompression structure
    jpeg_create_decompress(&cinfo);

    // Setup decompression buffer management
    jsrc.orig_buf = const_cast<unsigned char*>(blob);
    jsrc.orig_buf_size = blob_size;
    cinfo.src = &jsrc.pub;
    cinfo.src->init_source = my_init_source;
    cinfo.src->fill_input_buffer = my_fill_input_buffer;
    cinfo.src->skip_input_data = my_skip_input_data;
    cinfo.src->resync_to_restart = jpeg_resync_to_restart;
    cinfo.src->term_source = my_term_source;

    // Decompress
    jpeg_read_header(&cinfo, TRUE);
    if (cinfo.image_width != width || cinfo.image_height != height
            || cinfo.num_components != static_cast<int>(components)) {
        throw std::runtime_error("unexpected image dimensions");
    }
    jpeg_start_decompress(&cinfo);
    unsigned char rowbuf[cinfo.output_width * cinfo.num_components];
    jrow[0] = rowbuf;
    unsigned char* dst = data;
    size_t scanline_size = cinfo.image_width * cinfo.num_components * sizeof(unsigned char);
    while (cinfo.output_scanline < cinfo.image_height) {
        jpeg_read_scanlines(&cinfo, jrow, 1);
        std::memcpy(dst, rowbuf, scanline_size);
        dst += scanline_size;
    }
    jpeg_finish_decompress(&cinfo);

    // Finish
    jpeg_destroy_decompress(&cinfo);
}
