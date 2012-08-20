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

#ifndef ECM_DB_H
#define ECM_DB_H

#include <vector>
#include <string>
#include <cstdint>

#include <ecm/ecm.h>

/**
 * \file ecmdb.h
 * \brief The ECM DB class.
 *
 * The ECM DB class manages Ellipsoid Cube Map data bases.
 */


/**
 * Ellipsoid Cube Map (ECM) data bases
 *
 * An ECM data base consists of a directory that is filled with data files, organized in subdirectories.
 * This class helps to manage such a directory structure.
 *
 * The data stored in an ECM data base has one of three types: elevation data, texture data (in the SRGB
 * color space), or undefined data that is later converted to elevation or texture data by a rendering
 * system. The data can have several channels (between 1 and ecmdb::max_channels), e.g. 3 for SRGB data.
 * The data is stored using one of three data types: uint8, int16, or float32. Elevation data has one
 * channel and data type int16 or float32, SRGB texture data has one or three channels and data type uint8.
 *
 * The data is managed in six quadtree hierarchies (one for each side of the circumscribing cube).
 * One quad provides a fixed number of samples (quad_size x quad_size). Additionally, quads may overlap.
 * This is useful to avoid artefacts when the rendering system needs to access elevation or texture samples
 * at the quad borders.
 *
 * Each quad stores one data array, containing the data samples for all the channels, and one mask array that
 * determines which of the data samples are valid. This allows to manage holes in data sets, which is common
 * for remote sensing data.
 *
 * For elevation data, a rendering system typically needs meta data for level of detail purposes. This meta data
 * consists of the global minimum and maximum elevation, and of the minimum and maximum elevation per quad.
 * The global meta data is read when opening a data base, the per-quad metadata is stored with each quad data file.
 *
 * Each data set is only valid for a specific ellipsoid defined by its semi-major and semi-minor axes.
 *
 * A rendering system can typically combine data bases that use the same ellipsoid and the same quad size.
 * Such data bases are called compatible data bases.
 */
class ecmdb
{
public:
    /** \name Limits **/
    /*@{*/

    static const int max_levels = 31;           /**< \brief Maximum number of quadtree levels; required to prevent overflow of quad coordinates */
    static const int max_channels = 16;         /**< \brief Arbitrary limit, but keep it low to avoid overhead */
    static const int max_quad_size = 2048;      /**< \brief Arbitrary limit, but keep it low to avoid overhead */
    static const int max_overlap = 16;          /**< \brief Arbitrary limit, but keep it low to avoid overhead */

    /*@}*/
    /** \name Constants **/
    /*@{*/

    /** Data category */
    enum category_t {
        category_invalid = 0,           /**< \brief Invalid data base */
        category_elevation = 1,         /**< \brief Elevation data (one channel, type int16 or float32) */
        category_texture = 2,           /**< \brief Texture data SRGB (three channels, type uint8) */
        category_data = 3               /**< \brief Unspecified data */
    };

    /** Data type */
    enum type_t {
        type_uint8 = 0,                 /**< \brief Data type uint8 */
        type_int16 = 1,                 /**< \brief Data type int16 */
        type_float32 = 2                /**< \brief Data type float32 */
    };

    /*@}*/
    /** \name Meta data **/
    /*@{*/

    /** Global meta data for the whole data base */
    class global_metadata
    {
    public:
        ecmdb::category_t category;   /**< \brief ecmdb::category */
        union {
            /** Meta data for ecmdb::category_elevation */
            struct {
                float min;      /**< Minimum elevation in this data base */
                float max;      /**< Maximum elevation in this data base */
            } elevation;
            /* Meta data for ecmdb::category_texture */
            // empty
            /* Meta data for ecmdb::category_data */
            // empty
        };
        /** Constructor for invalid meta data */
        global_metadata() : category(category_invalid) { }
        /** Constructor for valid meta data of the given category */
        global_metadata(ecmdb::category_t category) { create(category); }
        /** Create valid meta data for the given category */
        void create(ecmdb::category_t category);
        /** Return whether this meta data is valid */
        bool is_valid() const { return category >= 0; }
        /** Serialization: save to stream */
        void save(std::ostream& os) const;
        /** Serialization: load from stream */
        void load(std::istream& is);
        /** Open a metadata file */
        void open(ecmdb::category_t category, const std::string& dirname, const std::string& original_url = "");
        /** Save a metadata file */
        void write(const std::string& dirname) const;
    private:
        void check_validity(const std::string& source = "") const;
    };

    /** Per-quad meta data */
    class quad_metadata
    {
    public:
        ecmdb::category_t category;   /**< \brief ecmdb::category */
        union {
            /** Meta data for ecmdb::category_elevation */
            struct {
                float min;      /**< Minimum elevation in this data base */
                float max;      /**< Maximum elevation in this data base */
            } elevation;
            /* Meta data for ecmdb::category_texture */
            // empty
            /* Meta data for ecmdb::category_data */
            // empty
        };
        /** Constructor for invalid meta data */
        quad_metadata() : category(category_invalid) { }
        /** Constructor for valid meta data of the given category */
        quad_metadata(ecmdb::category_t category) { create(category); }
        /** Create valid meta data for the given category */
        void create(ecmdb::category_t category);
        /** Return whether this meta data is valid */
        bool is_valid() const { return category >= 0; }
    };

private:
    // Properties of the database. Only databases with matching properties can be combined.
    class ecm _ecm;             // The ECM, defined by semi-major axis and semi-minor axis of the ellipsoid.
    int _quad_size;             // Base quad size (without overlap), e.g. 512 for 512x512 quads

    // Number of levels in this database.
    int _levels;

    // Extent of this database at the maximum level, for the six sides of the cube,
    // as a bounding box.
    int _extent_xmin[6];
    int _extent_xmax[6];
    int _extent_ymin[6];
    int _extent_ymax[6];

    // Properties of the contents of this database.
    category_t _category;       // Category
    type_t _type;               // Data type
    int _channels;              // Number of channels (e.g. 3 for RGB)
    int _overlap;               // Quad overlap
    float _data_offset;         // Offset and factor for the data:
    float _data_factor;         // real_value = factor * value + offset
    std::string _short_description;             // Short description (one line)
    std::vector<std::string> _description;      // Long description (multiple lines)

    void check_validity(const std::string& source = "") const;

public:
    /** Constructor for an invalid (empty) database */
    ecmdb() : _levels(0) {};

    /** Serialization: save this object to a stream */
    void save(std::ostream& os) const;
    /** Serialization: load this object from a stream */
    void load(std::istream& is);

    /** Check if this database is valid (non-empty). */
    bool is_valid() const
    {
        return _levels > 0;
    }

    /**
     * Check if this database is compatible with another database,
     * i.e. data from both databases can be combined into a single view.
     */
    bool is_compatible(const ecmdb& database)
    {
        return (_ecm == database._ecm && _quad_size == database._quad_size);
    }

    /** Return the ECM. */
    const class ecm& ecm() const
    {
        return _ecm;
    }

    /** Return the semi-major axis of the ellipsoid. */
    double semi_major_axis() const
    {
        return _ecm.semi_major_axis();
    }

    /** Return the semi-minor axis of the ellipsoid. */
    double semi_minor_axis() const
    {
        return _ecm.semi_minor_axis();
    }

    /** Return the category of data stored in this database. */
    category_t category() const
    {
        return _category;
    }

    /** Return the data type used by this database. */
    type_t type() const
    {
        return _type;
    }

    /** Return the size of the data type used by this database. */
    size_t type_size() const
    {
        return (type() == type_uint8 ? sizeof(uint8_t)
                : type() == type_int16 ? sizeof(int16_t)
                : sizeof(float));
    }

    /** Return the number of channels of the data, e.g. 3 for RGB texture data. */
    int channels() const
    {
        return _channels;
    }

    /** Return the size of one data element. */
    size_t element_size() const
    {
        return channels() * type_size();
    }

    /** Return the base quad size, without overlap, e.g. 512 for quads of size 512x512. */
    int quad_size() const
    {
        return _quad_size;
    }

    /**
     * Return the overlap size. For example, a 512x512 quad with an overlap of 2
     * leads to 516x516 data samples.
     */
    int overlap() const
    {
        return _overlap;
    }

    /** Return the total quad size, including overlap areas. */
    int total_quad_size() const
    {
        return quad_size() + 2 * overlap();
    }

    /** Return the size of one quad worth of data, including overlap. */
    size_t data_size() const
    {
        return total_quad_size() * total_quad_size() * element_size();
    }

    /** Return the size of the validity mask for one quad, including overlap. */
    size_t mask_size() const
    {
        return total_quad_size() * total_quad_size() * sizeof(uint8_t);
    }

    /**
     * Return the data offset that needs to be applied when interpreting values
     * from this database: real_value = factor * value + offset.
     */
    float data_offset() const
    {
        return _data_offset;
    }

    /**
     * Return the data factor that needs to be applied when interpreting values
     * from this database: real_value = factor * value + offset.
     */
    float data_factor() const
    {
        return _data_factor;
    }

    /** Return the short description (one line) for this database. */
    const std::string& short_description() const
    {
        return _short_description;
    }

    /** Return the long description (a vector of lines) for this database. */
    const std::vector<std::string>& description() const
    {
        return _description;
    }

    /** Return the number of quadtree levels available in this database. */
    int levels() const
    {
        return _levels;
    }

    /**
     * \param side      ECM cube side.
     * \param level     Cube side quadtree level.
     * \param x         X coordinate of quad.
     * \param y         Y coordinate of quad.
     *
     * Check if this database potentially contains a given quad.
     *
     * If this returns false, then you know the database does not contain the quad.
     * If this returns true, then the database contains the quad, but note that this
     * does not necessarily mean that the corresponding quad file exists. Only quads
     * with valid data are stored on disk. If the quad file does not exist, then the
     * quad does not contain any valid data (for example if the data set has a hole).
     */
    bool has_quad(int side, int level, int x, int y) const
    {
        if (level >= _levels || _extent_xmin[side] < 0) {
            return false;
        } else {
            int div = (1 << (_levels - 1 - level));
            int sxmin = _extent_xmin[side] / div;
            int sxmax = _extent_xmax[side] / div;
            int symin = _extent_ymin[side] / div;
            int symax = _extent_ymax[side] / div;
            return (x >= sxmin && x <= sxmax && y >= symin && y <= symax);
        }
    }

    /**
     * \param dirname           Database directory.
     * \param original_url      Original URL.
     *
     * Open a database. The original URL is useful for error messages if this database
     * actually lives somewhere on the network and the directory only refers to a locally
     * cached copy.
     */
    void open(const std::string& dirname, const std::string& original_url = "");

    /**
     * Return the name of the quad file that stores data and mask for the given quad.
     * Note that this quad file only exists if the quad contains valid data.
     * The quad file name is relative to the database directory:
     * abs_quadfilename = dbdir + '/' + quad_filename();
     */
    static std::string quad_filename(int side, int level, int x, int y);

    /**
     * Load quad data and validity mask from a quad file.
     * If the all_valid flag is set afterwards, then all data samples are valid, and the mask buffer was not touched.
     * Otherwise, the mask buffer stores the value 255 for all valid data samples and 0 for all invalid data samples.
     */
    void load_quad(const std::string& filename, void* data, uint8_t* mask, bool* all_valid,
            quad_metadata* meta) const;

    /* The rest of the interface is only interesting for creating new databases. This is what the ecmdb tool does.
     * Refer to its documentation and source code. */
    void create(
            double semi_major_axis, double semi_minor_axis, int quad_size, int levels,
            category_t category, type_t type, int channels, int overlap, float data_offset, float data_factor,
            const std::string& short_description, const std::vector<std::string>& description);
    void write(const std::string& dirname) const;
    void load_quad_meta(const std::string& filename, quad_metadata* meta) const;
    void save_quad(const std::string& filename, const void* data, const uint8_t* mask, bool all_valid,
            const quad_metadata* meta,
            int compression /* 0=none, 1=lossless, 2=lossy */, int lossy_compression_quality = 80 /* 1-100*/) const;
    void add_quad(int side, int x, int y);
    void close();
};

#endif
