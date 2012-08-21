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

#ifndef ECM_H
#define ECM_H

/**
 * \file ecm.h
 * \brief The ECM class.
 *
 * The ECM class provides the basic Ellipsoid Cube Map model.
 */

/**
 * \mainpage Ellipsoid Cube Map (ECM)
 */

#include <cstddef>
#include <ecmdb/version.h>

/**
 * Ellipsoid Cube Map (ECM):
 *
 * An Ellipsoid Cube Map (ECM) maps an ellipsoid to a cube.
 *
 * The ellipsoid is defined by its semi-major axis (a) and semi-minor axis (b). The cube side
 * length is 2a.
 *
 * By convention, the cube is axis-aligned with the origin at its center. Its six sides are
 * numbered as follows: side 0 is at x=+a (front), 1 is at y=+a (right), 2 is at x=-a (back),
 * 3 is at y=-a (left), 4 is at z=+a (top), and 5 is at z=-a (bottom).
 *
 * Also by convention, the four directions are numbered as follows:
 * up = 0, right = 1, bottom = 2, left = 3.
 *
 * Ellipsoid Cube Map coordinates (ecx,ecy) use the following coordinate system:
 *
 * \verbatim
 ecy
  
 +3  +-------+
  |  |       |
 +2  |   4   |
  |  |       |
 +1  +-------+-------+-------+-------+
  |  |       |       |       |       |
  0  |   0   |   1   |   2   |   3   |
  |  |       |       |       |       |
 -1  +-------+-------+-------+-------+
  |  |       |
 -2  |   5   |
  |  |       |
 -3  +-------+

    -1---0---1---2---3---4---5---6---7  ecx
\endverbatim
 * 
 * For sides 0-3, the following applies: the lower bound for ecx,ecy for each side is
 * inclusive, the higher bound is exclusive. Thus, the cartesian point (-a,-a,0) belongs
 * to area 0, and (+a,-a,0) belongs to area 1.
 * For side 4, both lower bounds and upper bounds of ecx,ecy are exclusive (i.e. side 4
 * includes none of its edges).
 * For side 5, both lower bounds and upper bounds of ecx,ecy are inclusive (i.e. side 5
 * includes all of its edges).
 *
 * Inside each cube side, the Quadrilateralized Spherical Cube projection is used to map
 * the unit sphere to the unit cube. This is an equal-area projection and therefore
 * minimizes sampling problems (both in terms of sampling errors and in terms of data
 * storage efficiency). This projection is described in this paper:
 * E.M. O'Neill and R.E. Laubscher, "Extended Studies of a Quadrilateralized Spherical Cube Earth Data Base",
 * Naval Environmental Prediction Research Facility Tech. Report NEPRF 3-76 (CSC), May 1976.
 *
 * Ellipsoid Cube Map Quadtree:
 *
 * A quadtree on an Ellipsoid Cube Map is a group of 6 quadtrees (i.e. it has 6 root
 * nodes), one for each side of the cube, with proper neighboring relations as
 * given by the cube edges. For example, the level 0 quad that represents the
 * top side of the cuboid is a neighbor of the level 0 quads on the front, back,
 * left, and right sides of the cuboid; similar for higher quadtree levels.
 *
 * A quadtree node always either has zero children or four children. By convention,
 * the four children of a quad are numbered as follows:
 * \verbatim
+---+---+
| 0 | 1 |
+---+---+
| 2 | 3 |
+---+---+
\endverbatim
 */
class ecm
{
public:
    /** \name Constants **/
    /*@{*/

    /** Cube sides */
    enum side {
        side_front = 0,         /**< \brief Cube side: front */
        side_right = 1,         /**< \brief Cube side: right */
        side_back = 2,          /**< \brief Cube side: back */
        side_left = 3,          /**< \brief Cube side: left */
        side_top = 4,           /**< \brief Cube side: top */
        side_bottom = 5         /**< \brief Cube side: bottom */
    };

    /** Directions */
    enum direction {
        dir_up = 0,             /**< \brief Direction: up */
        dir_right = 1,          /**< \brief Direction: right */
        dir_down = 2,           /**< \brief Direction: down */
        dir_left = 3            /**< \brief Direction: left */
    };

    /** Quad corners */
    enum corner {
        corner_tl = 0,          /**< \brief Corner: top left */
        corner_tr = 1,          /**< \brief Corner: top right */
        corner_bl = 2,          /**< \brief Corner: bottom left */
        corner_br = 3           /**< \brief Corner: bottom right */
    };

    static const double semi_major_axis_earth_wgs84;    /**< \brief Semi-major axis of Earth (WGS84) */
    static const double semi_minor_axis_earth_wgs84;    /**< \brief Semi-minor axis of Earth (WGS84) */
    static const double radius_moon_nasa;               /**< \brief Radius of Moon (NASA) */
    static const double semi_major_axis_mars_nasa;      /**< \brief Semi-major axis of Mars (NASA) */
    static const double semi_minor_axis_mars_nasa;      /**< \brief Semi-minor axis of Mars (NASA) */
    /*@}*/

private:
    /** \name Properties **/
    /*@{*/
    double _semi_major_axis;    /**< Ellipsoid definition: semi-major axis */
    double _semi_minor_axis;    /**< Ellipsoid definition: semi-minor axis */
    /*@}*/

public:
    /** \name Constructor/Destructor, basic operations **/
    /*@{*/

    /**
     * \param semi_major_axis   Semi-major axis of the ellipsoid.
     * \param semi_minor_axis   Semi-minor axis of the ellipsoid.
     *
     * Constructor.
     */
    ecm(double semi_major_axis, double semi_minor_axis);

    /**
     * Constructor. This will create an *invalid* ECM!
     * Only the constructor with two arguments creates valid ECMs.
     */
    ecm();

    /**
     * \param semi_major_axis   Semi-major axis of the ellipsoid.
     * \param semi_minor_axis   Semi-minor axis of the ellipsoid.
     *
     * Reset.
     */
    void reset(double semi_major_axis, double semi_minor_axis);

    /**
     * Return true if two ECMs are identical.
     */
    bool operator==(const ecm& e) const {
        return (_semi_major_axis <= e._semi_major_axis && _semi_major_axis >= e._semi_major_axis
                && _semi_minor_axis <= e._semi_minor_axis && _semi_minor_axis >= e._semi_minor_axis);
    }

    /**
     * Return true if two ECMs are not identical.
     */
    bool operator!=(const ecm& e) const {
        return !(*this == e);
    }

    /**
     * Return true if this ECM is valid.
     */
    bool is_valid() const
    {
        return (_semi_major_axis > 0.0);
    }

    /**
     * Return the semi-major axis.
     */
    double semi_major_axis() const
    {
        return _semi_major_axis;
    }

    /**
     * Return the semi-major axis.
     */
    double semi_minor_axis() const
    {
        return _semi_minor_axis;
    }

    /*@}*/
    /** \name Coordinate conversions **/
    /*@{*/

    /**
     * \param ecm_x     ECM x coordinate.
     * \param ecm_y     ECM y coordinate.
     * \return          The cube side.
     *
     * Extract the cube side from ecm coordinates.
     */
    static int ecm_to_sidenumber(double ecm_x, double ecm_y);

    /**
     * \param ecm_x     ECM x coordinate.
     * \param ecm_y     ECM y coordinate.
     * \param side_x    Pointer to the resulting side x coordinate.
     * \param side_y    Pointer to the resulting side y coordinate.
     * \param sidenumber Pointer to the resulting side number (optional; may be zero).
     *
     * Return the side-relative coordinates.
     */
    static void ecm_to_side(double ecm_x, double ecm_y, double* side_x, double* side_y, int* sidenumber = NULL);

    /**
     * \param sidenumber Side number.
     * \param side_x    Side x coordinate.
     * \param side_y    Side y coordinate.
     * \param ecm_x     Pointer to the resulting ECM x coordinate.
     * \param ecm_y     Pointer to the resulting ECM y coordinate.
     *
     * Return ECM coordinates for side-relative coordinates.
     */
    static void side_to_ecm(int sidenumber, double side_x, double side_y, double* ecm_x, double* ecm_y);

    /**
     * \param ecm_x     ECM x coordinate.
     * \param ecm_y     ECM y coordinate.
     * \param geoc_lat  Pointer to the resulting geocentric latitude.
     * \param lon       Pointer to the resulting longitude.
     *
     * Convert ECM coordinates to geocentric latitude and longitude.
     */
    void ecm_to_geocentric(double ecm_x, double ecm_y, double* geoc_lat, double* lon) const;

    /**
     * \param ecm_x     ECM x coordinate.
     * \param ecm_y     ECM y coordinate.
     * \param geod_lat  Pointer to the resulting geodetic latitude.
     * \param lon       Pointer to the resulting longitude.
     *
     * Convert ECM coordinates to geodetic latitude and longitude.
     */
    void ecm_to_geodetic(double ecm_x, double ecm_y, double* geod_lat, double* lon) const;

    /**
     * \param ecm_x     ECM x coordinate.
     * \param ecm_y     ECM y coordinate.
     * \param cart      Pointer to the resulting cartesian coordinates (3 components).
     *
     * Convert ECM coordinates to planetocentric cartesian coordinates.
     */
    void ecm_to_cartesian(double ecm_x, double ecm_y, double* cart) const;

    /**
     * \param cart      Cartesian coordinates (3 components).
     * \param geoc_lat  Pointer to the resulting geocentric latitude.
     * \param lon       Pointer to the resulting longitude.
     *
     * Convert planetocentric cartesian coordinates to geocentric latitude and longitude.
     */
    void cartesian_to_geocentric(const double* cart, double* geoc_lat, double* lon) const;

    /**
     * \param geoc_lat  Geocentric latitude.
     * \param lon       Longitude.
     * \param geoc_alt  Geocentric altitude.
     * \param cart      Pointer to the resulting cartesian coordinates (3 components).
     *
     * Convert geocentric coordinates to planetocentric cartesian coordinates.
     */
    void geocentric_to_cartesian(double geoc_lat, double lon, double geoc_alt, double* cart) const;

    /**
     * \param geoc_lat  Geocentric latitude.
     * \return          Geodetic latitude.
     *
     * Convert geocentric latitude to geodetic latitude.
     * (Note that there is no difference between geocentric and geodetic longitude).
     */
    double geocentric_to_geodetic(double geoc_lat) const;

    /**
     * \param geod_lat  Geodetic latitude.
     * \return          Geocentric latitude.
     *
     * Convert geodetic latitude to geocentric latitude.
     * (Note that there is no difference between geocentric and geodetic longitude).
     */
    double geodetic_to_geocentric(double geod_lat) const;

    /**
     * \param geod_lat  Geodetic latitude.
     * \param lon       Longitude.
     * \param geod_alt  Geodetic altitude.
     * \param cart      Pointer to the resulting cartesian coordinates (3 components).
     *
     * Convert geodetic coordinates to planetocentric cartesian coordinates.
     */
    void geodetic_to_cartesian(double geod_lat, double lon, double geod_alt, double* cart) const;

    /**
     * \param cart      Cartesian coordinates (3 components).
     * \param geod_lat  Pointer to the resulting geodetic latitude.
     * \param lon       Pointer to the resulting longitude.
     * \param geod_alt  Pointer to the resulting geodetic altitude.
     *
     * Convert planetocentric cartesian coordinates to geodetic coordinates.
     */
    void cartesian_to_geodetic(const double* cart, double* geod_lat, double* lon, double* geod_alt) const;

    /**
     * \param geod_lat  Geodetic latitude.
     * \param lon       Longitude.
     * \param normal    Pointer to the resulting normal (3 components).
     *
     * Return the ellipsoid surface normal at the geodetic latitude and longitude.
     * The returned normal has length 1.
     */
    static void geodetic_normal(double geod_lat, double lon, double* normal);

    /*@}*/
    /** \name ECM quadtree quads **/
    /*@{*/
    
    /**
     * \param quad_side Quadtree side.
     * \param quad_level Quadtree level.
     * \param quad_x    Quad x coordinate.
     * \param quad_y    Quad y coordinate.
     * \param qx        Quad-relative x position.
     * \param qy        Quad-relative y position.
     * \param ecm_x     Pointer to the resulting ECM x coordinate.
     * \param ecm_y     Pointer to the resulting ECM y coordinate.
     *
     * Compute the ECM coordinates of a point on a quad.
     */
    static void quad_to_ecm(int quad_side, int quad_level, int quad_x, int quad_y, double qx, double qy,
            double* ecm_x, double* ecm_y);
    
    /**
     * \param quad_side Quadtree side.
     * \param quad_level Quadtree level.
     * \param quad_x    Quad x coordinate.
     * \param quad_y    Quad y coordinate.
     * \param corner    The corner of the quad.
     * \param ecm_x     Pointer to the resulting ECM x coordinate.
     * \param ecm_y     Pointer to the resulting ECM y coordinate.
     *
     * Compute the ECM coordinates of a corner of a quad.
     */
    static void quad_to_ecm(int quad_side, int quad_level, int quad_x, int quad_y, int corner,
            double* ecm_x, double* ecm_y);

    /**
     * \param quad_side         Quadtree side.
     * \param quad_level        Quadtree level.
     * \param quad_x            Quad x coordinate.
     * \param quad_y            Quad y coordinate.
     * \param quad_tl_cart      Quad corner top left: cartesian coordinates (3 components).
     * \param quad_tr_cart      Quad corner top right: cartesian coordinates (3 components).
     * \param quad_bl_cart      Quad corner bottom left: cartesian coordinates (3 components).
     * \param quad_br_cart      Quad corner bottom right: cartesian coordinates (3 components).
     * \param plane_normal      Pointer to the normalized quad plane normal (3 components).
     * \param plane_distance    Pointer to the quad plane distance.
     *
     * Returns the quad plane E for the given quad in Hesse normal form:
     * E: plane_normal * x = plane_distance
     */
    void quad_plane(int quad_side, int quad_level, int quad_x, int quad_y,
            const double* quad_tl_cart, const double* quad_tr_cart, const double* quad_bl_cart, const double* quad_br_cart,
            double* plane_normal, double* plane_distance) const;

    /**
     * \param quad_side         Quadtree side.
     * \param quad_level        Quadtree level.
     * \param quad_x            Quad x coordinate.
     * \param quad_y            Quad y coordinate.
     * \param quad_plane_normal Normalized quad plane normal (3 components).
     * \param quad_plane_distance Quad plane distance.
     *
     * Estimate the maximum distance of the quad plane to the ellipsoid surface.
     * This estimation might be bad, but can be computed quickly.
     * The correct value is computed by quad_base_data().
     */
    double max_quad_plane_distance_estimation(int quad_side, int quad_level, int quad_x, int quad_y,
            const double* quad_plane_normal, double quad_plane_distance) const;

    /**
     * \param quad_side         Quadtree side.
     * \param quad_level        Quadtree level.
     * \param quad_x            Quad x coordinate.
     * \param quad_y            Quad y coordinate.
     * \param sym_quad_side     Pointer to the symmetry quadtree side.
     * \param sym_quad_level    Pointer to the symmetry quadtree level.
     * \param sym_quad_x        Pointer to the symmetry quad x coordinate.
     * \param sym_quad_y        Pointer to the symmetry quad y coordinate.
     * \param mirror_x          Pointer to the mirror-horizontally flag, or NULL.
     * \param mirror_y          Pointer to the mirror-vertically flag, or NULL.
     * \param matrix            The transformation matrix, or NULL.
     *
     * Quads in ECM quadtrees are symmetrical: the top and bottom sides of the cube are mirrors of each other,
     * and the four remaining sides can be rotated into each other. Furthermore, the four quarters
     * of each side can be mirrored into each other. This can be used to reduce the number of quads
     * for which base data needs to be computed with quad_base_data().\n
     * This function returns a symmetry quad for an original quad, and tells you how to transform
     * the base data computed by quad_base_data() for this symmetry quad so that it applies to the
     * original quad: mirror the texture coordinates according to the mirror flags, and apply the
     * column-major transformation matrix.
     */
    static void symmetry_quad(int quad_side, int quad_level, int quad_x, int quad_y,
            int* sym_quad_side, int* sym_quad_level, int* sym_quad_x, int* sym_quad_y,
            bool* mirror_x, bool* mirror_y, float matrix[9]);

    /**
     * \param quad_side         Quadtree side.
     * \param quad_level        Quadtree level.
     * \param quad_x            Quad x coordinate.
     * \param quad_y            Quad y coordinate.
     * \param quad_tl_cart      Quad corner top left: cartesian coordinates (3 components).
     * \param quad_tr_cart      Quad corner top right: cartesian coordinates (3 components).
     * \param quad_bl_cart      Quad corner bottom left: cartesian coordinates (3 components).
     * \param quad_br_cart      Quad corner bottom right: cartesian coordinates (3 components).
     * \param quad_plane_normal Normalized quad plane normal (3 components).
     * \param quad_plane_distance Quad plane distance.
     * \param quad_size         The quad size in samples.
     * \param quad_overlap      The quad overlap size in samples.
     * \param offsets           Storage space for (quad_size + 2*quad_overlap)^2 offset vectors (3 floats).
     * \param normals           Storage space for (quad_size + 2*quad_overlap)^2 normal vectors (2 floats, only x and y).
     * \param max_quad_plane_distance Pointer to the maximum distance of the quad plane to the ellipsoid surface.
     *
     * Compute base data for a quad.\n
     * It is recommended that this expensive function is only called for symmetry quads (see summetry_quad()).
     * Offsets between positions interpolated from the quad corners and the real ellipsoid surface are stored in the
     * offsets array. Ellipsoid surface normals are stored in the normals array (but only the x and y components; z can
     * be computed from them since the normals have length 1).
     */
    void quad_base_data(int quad_side, int quad_level, int quad_x, int quad_y,
            const double* quad_tl_cart, const double* quad_tr_cart, const double* quad_bl_cart, const double* quad_br_cart,
            const double* quad_plane_normal, double quad_plane_distance,
            int quad_size, int quad_overlap,
            float* offsets, float* normals, double* max_quad_plane_distance) const;

    /*@}*/
};

/** \name Library version **/
/*@{*/

/**
 * \brief       Get the libecmdb version.
 * \param major Buffer for the major version number, or NULL.
 * \param minor Buffer for the minor version number, or NULL.
 * \param patch Buffer for the patch version number, or NULL.
 * \return      The libecmdb version string.
 *
 * This function returns the version string "MAJOR.MINOR.PATCH".
 * If the pointers \a major, \a minor, \a patch are not NULL,
 * the requested version number will be stored there.
 */
const char* ecmdb_version(int* major, int* minor, int* patch);

/*@}*/

#endif
