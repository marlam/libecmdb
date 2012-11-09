/*
 * Copyright (C) 2011, 2012
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

#include <vector>
#include <cstring>
#include <cmath>

#include <gta/gta.hpp>

#include "dbg.h"
#include "exc.h"
#include "msg.h"
#include "str.h"
#include "fio.h"
#include "s11n.h"
#include "blob.h"

#include "ecmdb/ecmdb.h"
#include "jpegblob.h"


const int ecmdb::max_levels;
const int ecmdb::max_quad_size;
const int ecmdb::max_channels;
const int ecmdb::max_overlap;


void ecmdb::save(std::ostream& os) const
{
    s11n::save(os, semi_major_axis());
    s11n::save(os, semi_minor_axis());
    s11n::save(os, _quad_size);
    s11n::save(os, _levels);
    s11n::save(os, _extent_xmin, sizeof(_extent_xmin));
    s11n::save(os, _extent_xmax, sizeof(_extent_xmax));
    s11n::save(os, _extent_ymin, sizeof(_extent_ymin));
    s11n::save(os, _extent_ymax, sizeof(_extent_ymax));
    s11n::save(os, static_cast<int>(_category));
    s11n::save(os, static_cast<int>(_type));
    s11n::save(os, _channels);
    s11n::save(os, _overlap);
    s11n::save(os, _data_offset);
    s11n::save(os, _data_factor);
    s11n::save(os, _short_description);
    s11n::save(os, _description);
}

void ecmdb::load(std::istream& is)
{
    int x;
    double smaj, smin;
    s11n::load(is, smaj);
    s11n::load(is, smin);
    _ecm.reset(smaj, smin);
    s11n::load(is, _quad_size);
    s11n::load(is, _levels);
    s11n::load(is, _extent_xmin, sizeof(_extent_xmin));
    s11n::load(is, _extent_xmax, sizeof(_extent_xmax));
    s11n::load(is, _extent_ymin, sizeof(_extent_ymin));
    s11n::load(is, _extent_ymax, sizeof(_extent_ymax));
    s11n::load(is, x); _category = static_cast<category_t>(x);
    s11n::load(is, x); _type = static_cast<type_t>(x);
    s11n::load(is, _channels);
    s11n::load(is, _overlap);
    s11n::load(is, _data_offset);
    s11n::load(is, _data_factor);
    s11n::load(is, _short_description);
    s11n::load(is, _description);
}

std::string ecmdb::quad_filename(int side, int level, int x, int y)
{
    static const char hex[] = "0123456789abcdef";
    static const int bits = 12;         // Must be at least 4; should be a multiple of 4

    class _local
    {
    public:
        static size_t tile_path_component(const char* hex, const int bits, unsigned int val, std::string& str, size_t i)
        {
            for (int j = bits / 4 - 1; j >= 0; j--) {
                str[i++] = hex[(val >> (4 * j)) & 0xf];
            }
            return i;
        }
    };

    int64_t index = (static_cast<int64_t>(y) << level) | x;

    size_t length = 1 + 1 + 2 + 1 + ((2 * level - 1) / bits + 1) * (bits / 4 + 1) - 1 + 4;
    std::string path(length, '\0');

    size_t i = 0;
    // side and level
    path[i++] = '0' + side;
    path[i++] = '-';
    path[i++] = '0' + level / 10;
    path[i++] = '0' + level % 10;
    path[i++] = '/';
    // directory components
    for (int l = 2 * level; l > bits; l -= bits) {
        unsigned int component = (index >> ((l - 1) / bits * bits)) & ((1 << bits) - 1);
        i = _local::tile_path_component(hex, bits, component, path, i);
        path[i++] = '/';
    }
    // filename component
    unsigned int component = index & ((1 << bits) - 1);
    i = _local::tile_path_component(hex, bits, component, path, i);
    // extension component
    path[i++] = '.';
    path[i++] = 'g';
    path[i++] = 't';
    path[i++] = 'a';
    assert(i == length);
    return path;
}

void ecmdb::add_quad(int side, int x, int y)
{
    assert(side >= 0 && side <= 6);
    assert(x >= 0 && x < (1 << (_levels - 1)));
    assert(y >= 0 && y < (1 << (_levels - 1)));

    if (_extent_xmin[side] < 0) {
        _extent_xmin[side] = x;
        _extent_xmax[side] = x;
        _extent_ymin[side] = y;
        _extent_ymax[side] = y;
    } else {
        if (x < _extent_xmin[side]) {
            _extent_xmin[side] = x;
        } else if (x > _extent_xmax[side]) {
            _extent_xmax[side] = x;
        }
        if (y < _extent_ymin[side]) {
            _extent_ymin[side] = y;
        } else if (y > _extent_ymax[side]) {
            _extent_ymax[side] = y;
        }
    }
}

void ecmdb::check_validity(const std::string& source) const
{
    if (       !std::isfinite(semi_major_axis()) || semi_major_axis() <= 0.0
            || !std::isfinite(semi_minor_axis()) || semi_minor_axis() <= 0.0
            || semi_major_axis() < semi_minor_axis()
            || _quad_size < 0 || _quad_size > max_quad_size
            || _levels < 0 || _levels > max_levels
            || _overlap < 1 || _overlap > max_overlap || _overlap >= _quad_size / 2
            || !std::isfinite(_data_offset) || !std::isfinite(_data_factor)) {
        throw exc("Invalid database properties" + (source.empty() ? "" : " in " + source));
    }
    for (int i = 0; i < 6; i++) {
        if (_extent_xmin[i] == -1 && _extent_xmax[i] == -1 && _extent_ymin[i] == -1 && _extent_ymax[i] == -1) {
            continue;
        }
        if (_extent_xmin[i] < 0 || _extent_xmin[i] >= (1 << _levels)
                || _extent_xmax[i] < 0 || _extent_xmax[i] >= (1 << _levels)
                || _extent_xmax[i] < _extent_xmin[i]
                || _extent_ymin[i] < 0 || _extent_ymin[i] >= (1 << _levels)
                || _extent_ymax[i] < 0 || _extent_ymax[i] >= (1 << _levels)
                || _extent_ymax[i] < _extent_ymin[i]) {
            throw exc("Invalid database properties" + (source.empty() ? "" : " in " + source));
        }
    }
    if ((_category != category_elevation
                && _category != category_texture
                && _category != category_data
                && _category != category_sar_amplitude)
            || (_type != type_uint8
                && _type != type_int16
                && _type != type_float32)
            || (_category == category_elevation && (_channels != 1 || (_type != type_float32 && _type != type_int16)))
            || (_category == category_texture && ((_channels != 1 && _channels != 3) || _type != type_uint8))
            || (_category == category_data && (_channels != 1 || (_type != type_float32 && _type != type_int16)))
            || (_category == category_sar_amplitude && (_channels != 1 || _type != type_float32))) {
        throw exc("Invalid database properties" + (source.empty() ? "" : " in " + source));
    }
}

void ecmdb::create(
        double semi_major_axis, double semi_minor_axis, int quad_size, int levels,
        category_t category, type_t type, int channels, int overlap, float data_offset, float data_factor,
        const std::string& short_description, const std::vector<std::string>& description)
{
    assert(!is_valid());

    try {
        close();        // initializes all fields
        _ecm.reset(semi_major_axis, semi_minor_axis);
        _quad_size = quad_size;
        _levels = levels;
        _category = category;
        _type = type;
        _channels = channels;
        _overlap = overlap;
        _data_offset = data_offset;
        _data_factor = data_factor;
        _short_description = short_description;
        _description = description;
        _short_description = str::sanitize(_short_description);
        for (size_t i = 0; i < _description.size(); i++) {
            _description[i] = str::sanitize(_description[i]);
        }
        check_validity();
    }
    catch (exc& e) {
        close();
        throw e;
    }
    catch (std::exception& e) {
        close();
        throw exc(e);
    }
}

void ecmdb::open(const std::string& dirname, const std::string& original_url)
{
    assert(!is_valid());
    close();
    double smaj = -1.0, smin = -1.0;

    std::string orig_url = original_url.empty() ? dirname : original_url;

    try {
        FILE* info_file = fio::open(dirname + "/ecmdb.txt", "r");

        bool in_description_text = false;
        std::string line;
        while (!(line = fio::readline(info_file, dirname + "/ecmdb.txt")).empty()) {
            if (line[0] == '#') {
                continue;
            }
            if (in_description_text) {
                if (line[0] != ' ') {
                    in_description_text = false;
                } else {
                    _description.push_back(line.substr(1));
                    continue;
                }
            }
            if (       std::sscanf(line.c_str(), "Semi-Major-Axis: %lf", &smaj) == 1
                    || std::sscanf(line.c_str(), "Semi-Minor-Axis: %lf", &smin) == 1
                    || std::sscanf(line.c_str(), "Quad-Size: %d", &_quad_size) == 1
                    || std::sscanf(line.c_str(), "Levels: %d", &_levels) == 1
                    || std::sscanf(line.c_str(), "Extent-X-Min: %d %d %d %d %d %d",
                        _extent_xmin + 0, _extent_xmin + 1, _extent_xmin + 2, _extent_xmin + 3, _extent_xmin + 4, _extent_xmin + 5) == 6
                    || std::sscanf(line.c_str(), "Extent-X-Max: %d %d %d %d %d %d",
                        _extent_xmax + 0, _extent_xmax + 1, _extent_xmax + 2, _extent_xmax + 3, _extent_xmax + 4, _extent_xmax + 5) == 6
                    || std::sscanf(line.c_str(), "Extent-Y-Min: %d %d %d %d %d %d",
                        _extent_ymin + 0, _extent_ymin + 1, _extent_ymin + 2, _extent_ymin + 3, _extent_ymin + 4, _extent_ymin + 5) == 6
                    || std::sscanf(line.c_str(), "Extent-Y-Max: %d %d %d %d %d %d",
                        _extent_ymax + 0, _extent_ymax + 1, _extent_ymax + 2, _extent_ymax + 3, _extent_ymax + 4, _extent_ymax + 5) == 6
                    || std::sscanf(line.c_str(), "Channels: %d", &_channels) == 1
                    || std::sscanf(line.c_str(), "Overlap: %d", &_overlap) == 1
                    || std::sscanf(line.c_str(), "Data-Offset: %f", &_data_offset) == 1
                    || std::sscanf(line.c_str(), "Data-Factor: %f", &_data_factor) == 1)
            {
                continue;
            }
            if (line.compare(0, 10, "Category: ") == 0) {
                if (line.find("Elevation", 10) != std::string::npos) {
                    _category = category_elevation;
                } else if (line.find("Texture", 10) != std::string::npos) {
                    _category = category_texture;
                } else if (line.find("Data", 10) != std::string::npos) {
                    _category = category_data;
                } else if (line.find("SAR-Amplitude", 10) != std::string::npos) {
                    _category = category_sar_amplitude;
                } else {
                    throw exc("Invalid Category field in " + orig_url + "/ecmdb.txt");
                }
                continue;
            }
            if (line.compare(0, 6, "Type: ") == 0) {
                if (line.find("uint8", 6) != std::string::npos) {
                    _type = type_uint8;
                } else if (line.find("int16", 6) != std::string::npos) {
                    _type = type_int16;
                } else if (line.find("float32", 6) != std::string::npos) {
                    _type = type_float32;
                } else {
                    throw exc("Invalid Type field in " + orig_url + "/ecmdb.txt");
                }
                continue;
            }
            if (line.compare(0, 13, "Description: ") == 0) {
                _short_description = line.substr(13);
                in_description_text = true;
                continue;
            }
            msg::wrn(orig_url + "/ecmdb.txt" + ": ignoring unrecognized line '" + line + "'");
        }
        fio::close(info_file, dirname + "/ecmdb.txt");
        if (_short_description.empty())
            _short_description = orig_url;
        _short_description = str::sanitize(_short_description);
        for (size_t i = 0; i < _description.size(); i++)
            _description[i] = str::sanitize(_description[i]);
        if (smaj > 0.0 && smin > 0.0 && smaj >= smin)
            _ecm.reset(smaj, smin);

        /* Check validity */

        check_validity(orig_url + "/ecmdb.txt");
    }
    catch (exc& e) {
        close();
        throw e;
    }
    catch (std::exception& e) {
        close();
        throw exc(e);
    }
}

void ecmdb::write(const std::string& dirname) const
{
    assert(is_valid());

    fio::mkdir_p(dirname);
    FILE* info_file = fio::open(dirname + "/ecmdb.txt", "w");
    fprintf(info_file, "Semi-Major-Axis: %s\n", str::from(semi_major_axis()).c_str());
    fprintf(info_file, "Semi-Minor-Axis: %s\n", str::from(semi_minor_axis()).c_str());
    fprintf(info_file, "Quad-Size: %d\n", _quad_size);
    fprintf(info_file, "Levels: %d\n", _levels);
    fprintf(info_file, "Extent-X-Min: %d %d %d %d %d %d\n",
            _extent_xmin[0], _extent_xmin[1], _extent_xmin[2], _extent_xmin[3], _extent_xmin[4], _extent_xmin[5]);
    fprintf(info_file, "Extent-X-Max: %d %d %d %d %d %d\n",
            _extent_xmax[0], _extent_xmax[1], _extent_xmax[2], _extent_xmax[3], _extent_xmax[4], _extent_xmax[5]);
    fprintf(info_file, "Extent-Y-Min: %d %d %d %d %d %d\n",
            _extent_ymin[0], _extent_ymin[1], _extent_ymin[2], _extent_ymin[3], _extent_ymin[4], _extent_ymin[5]);
    fprintf(info_file, "Extent-Y-Max: %d %d %d %d %d %d\n",
            _extent_ymax[0], _extent_ymax[1], _extent_ymax[2], _extent_ymax[3], _extent_ymax[4], _extent_ymax[5]);
    fprintf(info_file, "Category: %s\n",
            _category == category_elevation ? "Elevation"
            : _category == category_texture ? "Texture"
            : _category == category_sar_amplitude ? "SAR-Amplitude"
            : "Data");
    fprintf(info_file, "Type: %s\n",
            _type == type_uint8 ? "uint8" : _type == type_int16 ? "int16" : "float32");
    fprintf(info_file, "Channels: %d\n", _channels);
    fprintf(info_file, "Overlap: %d\n", _overlap);
    fprintf(info_file, "Data-Offset: %s\n", str::from(_data_offset).c_str());
    fprintf(info_file, "Data-Factor: %s\n", str::from(_data_factor).c_str());
    fprintf(info_file, "Description: %s\n", _short_description.c_str());
    for (size_t i = 0; i < _description.size(); i++) {
        fprintf(info_file, " %s\n", _description[i].c_str());
    }
    fio::close(info_file, dirname + "/ecmdb.txt");
}

void ecmdb::close()
{
    _ecm = ecm();
    _quad_size = -1;
    _levels = -1;
    for (int i = 0; i < 6; i++) {
        _extent_xmin[i] = -1;
        _extent_xmax[i] = -1;
        _extent_ymin[i] = -1;
        _extent_ymax[i] = -1;
    }
    _category = category_invalid;
    _type = type_uint8;
    _channels = -1;
    _overlap = -1;
    _data_offset = 0.0f;
    _data_factor = 1.0f;
    _short_description.clear();
    _description.clear();
}

void ecmdb::metadata::create(ecmdb::category_t category)
{
    this->category = category;
    if (category == ecmdb::category_elevation) {
        elevation.min = 0.0f;
        elevation.max = 0.0f;
    } else if (category == ecmdb::category_sar_amplitude) {
        sar_amplitude.min = std::numeric_limits<float>::max();
        sar_amplitude.max = 0.0f;
        sar_amplitude.sum = 0.0;
        sar_amplitude.valid = 0;
    }
}

void ecmdb::metadata::check_validity(const std::string& source) const
{
    if ((category != ecmdb::category_elevation
                && category != ecmdb::category_texture
                && category != ecmdb::category_data
                && category != ecmdb::category_sar_amplitude)
            || (category == ecmdb::category_elevation
                && (!std::isfinite(elevation.min) || !std::isfinite(elevation.max)
                    || elevation.min > elevation.max))
            || (category == ecmdb::category_sar_amplitude
                && (!std::isfinite(sar_amplitude.min) || sar_amplitude.min < 0.0f
                    || !std::isfinite(sar_amplitude.max)
                    || !std::isfinite(sar_amplitude.sum)))
       ) {
        throw exc("Invalid metadata" + (source.empty() ? "" : " in " + source));
    }
}

void ecmdb::metadata::open(ecmdb::category_t category, const std::string& dirname, const std::string& original_url)
{
    create(category);
    std::string orig_url = original_url.empty() ? dirname : original_url;
    try {
        FILE* meta_file = fio::open(dirname + "/meta.txt", "r");
        std::string line;
        while (!(line = fio::readline(meta_file, dirname + "/meta.txt")).empty()) {
            if (line[0] == '#') {
                continue;
            }
            if (category == ecmdb::category_elevation) {
                if (       std::sscanf(line.c_str(), "Elevation-Min: %f", &elevation.min) == 1
                        || std::sscanf(line.c_str(), "Elevation-Max: %f", &elevation.max) == 1) {
                    continue;
                }
            } else if (category == ecmdb::category_sar_amplitude) {
                if (       std::sscanf(line.c_str(), "SAR-Amplitude-Min: %f", &sar_amplitude.min) == 1
                        || std::sscanf(line.c_str(), "SAR-Amplitude-Max: %f", &sar_amplitude.max) == 1
                        || std::sscanf(line.c_str(), "SAR-Amplitude-Sum: %lf", &sar_amplitude.sum) == 1
                        || std::sscanf(line.c_str(), "SAR-Amplitude-Valid: %Lu", &sar_amplitude.valid) == 1) {
                    continue;
                }
            }
            msg::wrn(orig_url + "/meta.txt" + ": ignoring unrecognized line '" + line + "'");
        }
        fio::close(meta_file, dirname + "/meta.txt");
        check_validity(orig_url + "/meta.txt");
    }
    catch (exc& e) {
        this->category = category_invalid;
        throw e;
    }
    catch (std::exception& e) {
        this->category = category_invalid;
        throw exc(e);
    }
}

void ecmdb::metadata::write(const std::string& dirname) const
{
    fio::mkdir_p(dirname);
    FILE* meta_file = fio::open(dirname + "/meta.txt", "w");
    if (category == ecmdb::category_elevation) {
        fprintf(meta_file, "Elevation-Min: %s\n", str::from(elevation.min).c_str());
        fprintf(meta_file, "Elevation-Max: %s\n", str::from(elevation.max).c_str());
    } else if (category == ecmdb::category_sar_amplitude) {
        fprintf(meta_file, "SAR-Amplitude-Min: %s\n", str::from(sar_amplitude.min).c_str());
        fprintf(meta_file, "SAR-Amplitude-Max: %s\n", str::from(sar_amplitude.max).c_str());
        fprintf(meta_file, "SAR-Amplitude-Sum: %s\n", str::from(sar_amplitude.sum).c_str());
        fprintf(meta_file, "SAR-Amplitude-Valid: %s\n", str::from(sar_amplitude.valid).c_str());
    }
    fio::close(meta_file, dirname + "/meta.txt");
}

void ecmdb::metadata::save(std::ostream& os) const
{
    s11n::save(os, static_cast<int>(category));
    if (category == ecmdb::category_elevation) {
        s11n::save(os, elevation.min);
        s11n::save(os, elevation.max);
    } else if (category == ecmdb::category_texture) {
    } else if (category == ecmdb::category_data) {
    } else if (category == ecmdb::category_sar_amplitude) {
        s11n::save(os, sar_amplitude.min);
        s11n::save(os, sar_amplitude.max);
        s11n::save(os, sar_amplitude.sum);
        s11n::save(os, sar_amplitude.valid);
    }
}

void ecmdb::metadata::load(std::istream& is)
{
    int x;
    s11n::load(is, x); category = static_cast<ecmdb::category_t>(x);
    if (category == ecmdb::category_elevation) {
        s11n::load(is, elevation.min);
        s11n::load(is, elevation.max);
    } else if (category == ecmdb::category_texture) {
    } else if (category == ecmdb::category_data) {
    } else if (category == ecmdb::category_sar_amplitude) {
        s11n::load(is, sar_amplitude.min);
        s11n::load(is, sar_amplitude.max);
        s11n::load(is, sar_amplitude.sum);
        s11n::load(is, sar_amplitude.valid);
    } else {
        category = category_invalid;
    }
}

static void get_meta_from_gta(ecmdb::category_t category, const std::string& filename, const gta::header& data_hdr, ecmdb::metadata* meta)
{
    meta->create(category);
    if (category == ecmdb::category_elevation) {
        const char* val;
        val = data_hdr.component_taglist(0).get("X-ECM-ELEVATION-MIN");
        if (!val)
            val = data_hdr.component_taglist(0).get("X-ECM-MIN");
        if (!val)
            val = data_hdr.component_taglist(0).get("X-PLEW-MIN");      // backward compat; remove soon
        if (!val)
            val = data_hdr.component_taglist(0).get("MIN");             // backward compat; remove soon
        if (!val)
            throw exc(std::string("Missing meta information in ") + filename);
        meta->elevation.min = str::to<float>(val);
        val = data_hdr.component_taglist(0).get("X-ECM-ELEVATION-MAX");
        if (!val)
            val = data_hdr.component_taglist(0).get("X-ECM-MAX");
        if (!val)
            val = data_hdr.component_taglist(0).get("X-PLEW-MAX");      // backward compat; remove soon
        if (!val)
            val = data_hdr.component_taglist(0).get("MAX");             // backward compat; remove soon
        if (!val)
            throw exc(std::string("Missing meta information in ") + filename);
        meta->elevation.max = str::to<float>(val);
        if (!std::isfinite(meta->elevation.min) || !std::isfinite(meta->elevation.max)
                || !(meta->elevation.min <= meta->elevation.max)) {
            throw exc(std::string("Invalid meta information in ") + filename);
        }
    } else if (category == ecmdb::category_sar_amplitude) {
        const char* val;
        val = data_hdr.component_taglist(0).get("X-ECM-SAR-AMPLITUDE-MIN");
        if (!val)
            throw exc(std::string("Missing meta information in ") + filename);
        meta->sar_amplitude.min = str::to<float>(val);
        val = data_hdr.component_taglist(0).get("X-ECM-SAR-AMPLITUDE-MAX");
        if (!val)
            throw exc(std::string("Missing meta information in ") + filename);
        meta->sar_amplitude.max = str::to<float>(val);
        val = data_hdr.component_taglist(0).get("X-ECM-SAR-AMPLITUDE-SUM");
        if (!val)
            throw exc(std::string("Missing meta information in ") + filename);
        meta->sar_amplitude.sum = str::to<float>(val);
        val = data_hdr.component_taglist(0).get("X-ECM-SAR-AMPLITUDE-VALID");
        if (!val)
            throw exc(std::string("Missing meta information in ") + filename);
        meta->sar_amplitude.valid = str::to<unsigned long long>(val);
        if (!std::isfinite(meta->sar_amplitude.min) || !std::isfinite(meta->sar_amplitude.max)
                || !(meta->sar_amplitude.min <= meta->sar_amplitude.max)
                || !std::isfinite(meta->sar_amplitude.sum)) {
            throw exc(std::string("Invalid meta information in ") + filename);
        }
    }
}

void ecmdb::load_quad_meta(const std::string& filename, metadata* meta) const
{
    if (category() == category_elevation || category() == category_sar_amplitude) {
        FILE* file = fio::open(filename, "r", O_NOATIME);
        gta::header data_hdr;
        data_hdr.read_from(file);
        fio::close(file, filename);
        get_meta_from_gta(category(), filename, data_hdr, meta);
    }
}

void ecmdb::load_quad(const std::string& filename, void* data, uint8_t* mask, bool* all_valid, metadata* meta) const
{
    FILE* f = fio::open(filename, "r", O_NOATIME);
    gta::header data_hdr;
    data_hdr.read_from(f);
    if (category() == category_texture
            && data_hdr.components() == 1
            && data_hdr.component_type(0) == gta::blob
            && data_hdr.dimensions() == 1
            && data_hdr.dimension_size(0) == 1) {
        blob jpegblob(data_hdr.data_size());
        data_hdr.read_data(f, jpegblob.ptr());
        try {
            jpegblob_dec(total_quad_size(), total_quad_size(), channels(),
                    jpegblob.ptr<const unsigned char>(), data_hdr.data_size(),
                    static_cast<unsigned char*>(data));
        }
        catch (std::exception& e) {
            throw exc(filename + ": " + e.what());
        }
    } else {
        if (data_hdr.data_size() != data_size())
            throw exc(std::string("Invalid content in ") + filename);
        data_hdr.read_data(f, data);
    }
    const char* validity_tag = data_hdr.global_taglist().get("X-ECM-VALIDITY");
    if (!validity_tag)
        validity_tag = data_hdr.global_taglist().get("X-PLEW-VALIDITY");        // backward compat; remove soon
    *all_valid = (validity_tag && std::strcmp(validity_tag, "ALL") == 0);
    if (!(*all_valid)) {
        gta::header mask_hdr;
        mask_hdr.read_from(f);
        if (mask_hdr.data_size() != mask_size())
            throw exc(std::string("Invalid content in ") + filename);
        mask_hdr.read_data(f, mask);
    }
    fio::advise(f, POSIX_FADV_DONTNEED, filename);
    fio::close(f);
    get_meta_from_gta(category(), filename, data_hdr, meta);
}

void ecmdb::save_quad(const std::string& filename, const void* data, const uint8_t* mask, bool all_valid,
        const metadata* meta,
        int compression /* 0=none, 1=lossless, 2=lossy */, int lossy_compression_quality /* 1-100*/) const
{
    FILE* quadfile = fio::open(filename, "w");
    gta::header data_hdr;
    data_hdr.global_taglist().set("X-ECM-VALIDITY", all_valid ? "ALL" : "MASK");
    if (compression == 2 && category() == category_texture) {
        blob jpegblob(data_size());
        size_t jpegblob_size;
        try {
            jpegblob_enc(total_quad_size(), total_quad_size(), channels(),
                    static_cast<const unsigned char*>(data),
                    lossy_compression_quality,
                    jpegblob.ptr<unsigned char>(), &jpegblob_size);
        }
        catch (std::exception& e) {
            throw exc(filename + ": " + e.what());
        }
        data_hdr.set_dimensions(1);
        data_hdr.set_components(gta::blob, jpegblob_size);
        data_hdr.write_to(quadfile);
        data_hdr.write_data(quadfile, jpegblob.ptr());
    } else {
        data_hdr.set_dimensions(total_quad_size(), total_quad_size());
        gta::type data_type = (type() == type_uint8 ? gta::uint8
                : type() == type_int16 ? gta::int16 : gta::float32);
        std::vector<gta::type> component_types(channels(), data_type);
        data_hdr.set_components(component_types.size(), &component_types[0]);
        if (compression)
            data_hdr.set_compression(type() == type_uint8 ? gta::bzip2 : gta::xz);
        if (category() == category_elevation) {
            data_hdr.component_taglist(0).set("X-ECM-ELEVATION-MIN", str::from(meta->elevation.min).c_str());
            data_hdr.component_taglist(0).set("X-ECM-ELEVATION-MAX", str::from(meta->elevation.max).c_str());
        } else if (category() == category_sar_amplitude) {
            data_hdr.component_taglist(0).set("X-ECM-SAR-AMPLITUDE-MIN", str::from(meta->sar_amplitude.min).c_str());
            data_hdr.component_taglist(0).set("X-ECM-SAR-AMPLITUDE-MAX", str::from(meta->sar_amplitude.max).c_str());
            data_hdr.component_taglist(0).set("X-ECM-SAR-AMPLITUDE-SUM", str::from(meta->sar_amplitude.sum).c_str());
            data_hdr.component_taglist(0).set("X-ECM-SAR-AMPLITUDE-VALID", str::from(meta->sar_amplitude.valid).c_str());
        }
        data_hdr.write_to(quadfile);
        data_hdr.write_data(quadfile, data);
    }
    if (!all_valid) {
        gta::header mask_hdr;
        mask_hdr.set_dimensions(total_quad_size(), total_quad_size());
        mask_hdr.set_components(gta::uint8);
        if (compression)
            mask_hdr.set_compression(gta::bzip2);
        mask_hdr.write_to(quadfile);
        mask_hdr.write_data(quadfile, mask);
    }
    fio::advise(quadfile, POSIX_FADV_DONTNEED, filename);
    fio::close(quadfile, filename);
}
