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

#include <cstring>

#include "dbg.h"
#include "compiler.h"

#include "glvm.h"

#include "ecm/ecm.h"

using namespace glvm;

#ifndef HAVE_SINCOS
static void sincos(double x, double* s, double* c)
{
    *s = sin(x);
    *c = cos(x);
}
#endif

const double ecm::semi_major_axis_earth_wgs84 = 6378137.0;
const double ecm::semi_minor_axis_earth_wgs84 = 6356752.314245;
const double ecm::radius_moon_nasa = 1737400.0;
const double ecm::semi_major_axis_mars_nasa = 3396000.0;
const double ecm::semi_minor_axis_mars_nasa = 3376800.0;

ecm::ecm(double semi_major_axis, double semi_minor_axis) :
    _semi_major_axis(semi_major_axis),
    _semi_minor_axis(semi_minor_axis)
{
    assert(isfinite(semi_major_axis));
    assert(isfinite(semi_minor_axis));
    assert(semi_minor_axis > 0.0);
    assert(semi_major_axis >= semi_minor_axis);
}

ecm::ecm() :
    _semi_major_axis(-1.0),
    _semi_minor_axis(-1.0)
{
}

void ecm::reset(double semi_major_axis, double semi_minor_axis)
{
    _semi_major_axis = semi_major_axis;
    _semi_minor_axis = semi_minor_axis;
    assert(isfinite(semi_major_axis));
    assert(isfinite(semi_minor_axis));
    assert(semi_minor_axis > 0.0);
    assert(semi_major_axis >= semi_minor_axis);
}

#ifdef NDEBUG
# define CHECK_ECM(ecm_x, ecm_y)
# define CHECK_SIDE(sn, sx, sy)
#else
# define CHECK_ECM(ecm_x, ecm_y) \
    assert(ecm_x >= -1.0 && ecm_x <= +7.0); \
    assert(ecm_y >= -3.0 && ecm_y <= +3.0); \
    assert(ecm_x <= +1.0 || (ecm_y >= -1.0 && ecm_y <= +1.0));
# define CHECK_SIDE(sn, sx, sy) \
    assert(sn >= 0 && sn <= 5); \
    assert(sx >= -1.0 && sx <= +1.0 && sy >= -1.0 && sy <= +1.0);
#endif

int ecm::ecm_to_sidenumber(double ecm_x, double ecm_y)
{
    CHECK_ECM(ecm_x, ecm_y);
    int sidenumber;
    if (ecm_y > +1.0)
        sidenumber = side_top;
    else if (ecm_y <= -1.0)
        sidenumber = side_bottom;
    else if (ecm_x < +1.0)
        sidenumber = side_front;
    else if (ecm_x < +2.0)
        sidenumber = side_right;
    else if (ecm_x < +3.0)
        sidenumber = side_back;
    else
        sidenumber = side_left;
    return sidenumber;
}

void ecm::ecm_to_side(double ecm_x, double ecm_y, double* side_x, double* side_y, int* sidenumber)
{
    CHECK_ECM(ecm_x, ecm_y);
    int sn;
    double sx, sy;
    if (ecm_y > +1.0) {
        sn = side_top;
        sx = ecm_x;
        sy = ecm_y - 2.0;
    } else if (ecm_y <= -1.0) {
        sn = side_bottom;
        if (ecm_y >= -1.0) {
            // special case: borders of sides 0..4
            if (ecm_x < 1.0) {
                sx = ecm_x;
                sy = +1.0;
            } else if (ecm_x < 3.0) {
                sx = +1.0;
                sy = -(ecm_x - 2.0);
            } else if (ecm_x < 5.0) {
                sx = -(ecm_x - 4.0);
                sy = -1.0;
            } else /* ecm_x <= 7.0 */ {
                sx = -1.0;
                sy = ecm_x - 6.0;
            }
        } else {
            sx = ecm_x;
            sy = ecm_y + 2.0;
        }
    } else {
        if (ecm_x < +1.0) {
            sn = side_front;
        } else if (ecm_x < +3.0) {
            sn = side_right;
        } else if (ecm_x < +5.0) {
            sn = side_back;
        } else {
            sn = side_left;
        }
        sx = ecm_x - (2 * sn);
        sy = ecm_y;
    }
    *side_x = sx;
    *side_y = sy;
    if (sidenumber)
        *sidenumber = sn;
    CHECK_SIDE(sn, sx, sy);
}

void ecm::side_to_ecm(int sidenumber, double side_x, double side_y, double* ecm_x, double* ecm_y)
{
    CHECK_SIDE(sidenumber, side_x, side_y);
    double ex, ey;
    if (sidenumber == side_top) {
        ex = side_x;
        ey = side_y + 2.0f;
    } else if (sidenumber == side_bottom) {
        ex = side_x;
        ey = side_y - 2.0f;
    } else {
        ex = side_x + (2 * sidenumber);
        ey = side_y;
    }
    *ecm_x = ex;
    *ecm_y = ey;
    CHECK_ECM(ex, ey);
}

void ecm::ecm_to_geocentric(double ecm_x, double ecm_y, double* geoc_lat, double* lon) const
{
    CHECK_ECM(ecm_x, ecm_y);

    int sidenumber;
    double side_x, side_y;
    ecm_to_side(ecm_x, ecm_y, &side_x, &side_y, &sidenumber);

    /* Convert the side coordinates to the mu and nu angles as used by QSC.
     * This depends on the area of the cube face. */
    double nu = atan(length(dvec2(side_x, side_y)));
    double mu = atan2(side_y, side_x);
    int area;
    if (side_x >= 0.0 && side_x >= abs(side_y)) {
        area = 0;
    } else if (side_y >= 0.0 && side_y >= abs(side_x)) {
        area = 1;
        mu -= const_pi_2<double>();
    } else if (side_x < 0.0 && -side_x >= abs(side_y)) {
        area = 2;
        mu = (mu < 0.0 ? mu + const_pi<double>() : mu - const_pi<double>());
    } else {
        area = 3;
        mu += const_pi_2<double>();
    }
    assert(abs(mu) <= const_pi_4<double>());

    /* Compute phi and theta for the area of definition.
     * The inverse projection is not described in the original paper, but some
     * good hints can be found here (as of 2011-12-14):
     * http://fits.gsfc.nasa.gov/fitsbits/saf.93/saf.9302
     * (search for "Message-Id: <9302181759.AA25477 at fits.cv.nrao.edu>") */
    double t = (const_pi<double>() / 12.0) * tan(mu);
    double tantheta = sin(t) / (cos(t) - (1.0 / const_sqrt2<double>()));
    double theta = atan(tantheta);
    double cosmu = cos(mu);
    double tannu = tan(nu);
    double cosphi = 1.0 - cosmu * cosmu * tannu * tannu * (1.0 - cos(atan(1.0 / cos(theta))));
    cosphi = clamp(cosphi, -1.0, +1.0);

    /* Apply the result to the real area on the cube face.
     * For the top and bottom face, we can compute phi and lam directly.
     * For the other faces, we must use unit sphere cartesian coordinates
     * as an intermediate step. */
    dvec2 geoc;
    if (sidenumber == side_top) {
        geoc[0] = const_pi_2<double>() - acos(cosphi);
        if (area == 0) {
            geoc[1] = theta + const_pi_2<double>();
        } else if (area == 1) {
            geoc[1] = (theta < 0.0 ? theta + const_pi<double>() : theta - const_pi<double>());
        } else if (area == 2) {
            geoc[1] = theta - const_pi_2<double>();
        } else /* area == 3 */ {
            geoc[1] = theta;
        }
    } else if (sidenumber == side_bottom) {
        geoc[0] = acos(cosphi) - const_pi_2<double>();
        if (area == 0) {
            geoc[1] = -theta + const_pi_2<double>();
        } else if (area == 1) {
            geoc[1] = -theta;
        } else if (area == 2) {
            geoc[1] = -theta - const_pi_2<double>();
        } else /* area == 3 */ {
            geoc[1] = (theta < 0.0 ? -theta - const_pi<double>() : -theta + const_pi<double>());
        }
    } else {
        /* Compute phi and lambda via cartesian unit sphere coordinates q,r,s. */
        double q, r, s, t;
        q = cosphi;
        t = q * q;
        if (t >= 1.0) {
            s = 0.0;
        } else {
            s = sqrt(1.0 - t) * sin(theta);
        }
        t += s * s;
        if (t >= 1.0) {
            r = 0.0;
        } else {
            r = sqrt(1.0 - t);
        }
        /* Rotate q,r,s into the correct area. */
        if (area == 1) {
            t = r;
            r = -s;
            s = t;
        } else if (area == 2) {
            r = -r;
            s = -s;
        } else if (area == 3) {
            t = r;
            r = s;
            s = -t;
        }
        /* Rotate q,r,s into the correct cube face. */
        if (sidenumber == side_right) {
            t = q;
            q = -r;
            r = t;
        } else if (sidenumber == side_back) {
            q = -q;
            r = -r;
        } else if (sidenumber == side_left) {
            t = q;
            q = r;
            r = -t;
        }
        /* Now compute phi and lam from the unit sphere coordinates q,r,s. */
        geoc[0] = acos(-s) - const_pi_2<double>();
        geoc[1] = atan2(r, q);
    }
    assert(isfinite(geoc[0]));
    assert(geoc[0] >= -const_pi_2<double>() && geoc[0] <= +const_pi_2<double>());
    assert(isfinite(geoc[1]));
    assert(geoc[1] >= -const_pi<double>() && geoc[1] <= +const_pi<double>());
    *geoc_lat = geoc[0];
    *lon = geoc[1];
}

void ecm::ecm_to_geodetic(double ecm_x, double ecm_y, double* geod_lat, double* lon) const
{
    CHECK_ECM(ecm_x, ecm_y);
    double geoc_lat;
    ecm_to_geocentric(ecm_x, ecm_y, &geoc_lat, lon);
    *geod_lat = geocentric_to_geodetic(geoc_lat);
}

void ecm::ecm_to_cartesian(double ecm_x, double ecm_y, double* cart) const
{
    CHECK_ECM(ecm_x, ecm_y);
    double geoc_lat, lon;
    ecm_to_geocentric(ecm_x, ecm_y, &geoc_lat, &lon);
    geocentric_to_cartesian(geoc_lat, lon, 0.0, cart);
}

void ecm::cartesian_to_geocentric(const double* cart, double* geoc_lat, double* lon) const
{
    assert(isfinite(cart[0]) && isfinite(cart[1]) && isfinite(cart[2]));
    *geoc_lat = acos(-cart[2] / length(dvec3(cart))) - const_pi_2<double>();
    *lon = atan2(cart[1], cart[0]);
}

void ecm::geocentric_to_cartesian(double geoc_lat, double lon, double geoc_alt, double* cart) const
{
    double sinlat, coslat;
    double sinlon, coslon;
    ::sincos(geoc_lat, &sinlat, &coslat);
    ::sincos(lon, &sinlon, &coslon);

    // From geocentric to unit-sphere-cartesian:
    dvec3 c;
    c.x = coslat * coslon;
    c.y = coslat * sinlon;
    c.z = sinlat;

    // From sphere-cartesian to ellipsoid-cartesian:
    double a2 = semi_major_axis() * semi_major_axis();
    double b2 = semi_minor_axis() * semi_minor_axis();
    c /= sqrt(c.x * c.x / a2 + c.y * c.y / a2 + c.z * c.z / b2);

    // Apply altitude:
    if (geoc_alt < 0.0 || geoc_alt > 0.0)
        c += geoc_alt * normalize(c);

    cart[0] = c.x;
    cart[1] = c.y;
    cart[2] = c.z;
}

double ecm::geocentric_to_geodetic(double geoc_lat) const
{
    double t = semi_minor_axis() / semi_major_axis();
    double tanlat = tan(geoc_lat);
    double xa = semi_minor_axis() / sqrt(tanlat * tanlat + t * t);
    double geod_lat = atan(sqrt(semi_major_axis() * semi_major_axis() - xa * xa) / (t * xa));
    if (geoc_lat < 0.0)
        geod_lat = -geod_lat;
    return geod_lat;
}

double ecm::geodetic_to_geocentric(double geod_lat) const
{
    double t = semi_minor_axis() / semi_major_axis();
    return atan((t * t) * tan(geod_lat));
}

void ecm::geodetic_to_cartesian(double geod_lat, double lon, double geod_alt, double* cart) const
{
    double sinlat, coslat;
    double sinlon, coslon;
    ::sincos(geod_lat, &sinlat, &coslat);
    ::sincos(lon, &sinlon, &coslon);

    double t = semi_minor_axis() / semi_major_axis();
    double e2 = (1.0 - t) * (1.0 + t);
    double N = semi_major_axis() / sqrt(1.0 - e2 * sinlat * sinlat);

    cart[0] = (N + geod_alt) * coslat * coslon;
    cart[1] = (N + geod_alt) * coslat * sinlon;
    cart[2]= (N * (1.0 - e2) + geod_alt) * sinlat;
}

void ecm::cartesian_to_geodetic(const double* cart, double* geod_lat, double* lon, double* geod_alt) const
{
    // See http://en.wikipedia.org/wiki/Geodetic_system
    // in its version from 2009-05-13 under 'Convert ECEF to WGS-84'.
    // XXX: It is unclear if this single-step algorithm is useful for anything other than WGS84!
    const double a2 = semi_major_axis() * semi_major_axis();
    const double b2 = semi_minor_axis() * semi_minor_axis();
    const double z2 = cart[2] * cart[2];

    const double t = semi_minor_axis() / semi_major_axis();
    const double e2 = (1.0 - t) * (1.0 + t);
    const double ep2 = e2 / (t * t);
    const double r2 = cart[0] * cart[0] + cart[1] * cart[1];
    const double r = sqrt(r2);
    const double E2 = a2 - b2;
    const double F = 54.0 * b2 * z2;
    const double G = r2 + (1.0 - e2) * z2 - e2 * E2;
    const double c = (e2 * e2 * F * r2) / (G * G * G);
    const double s = cbrt(1.0 + c + sqrt(c * c + 2.0 * c));
    const double P = F / (3.0 * (s + 1.0 / s + 1.0) * (s + 1.0 / s + 1.0) * G * G);
    const double Q = sqrt(1.0 + 2.0 * e2 * e2 * P);
    const double ro_radicand = (a2 / 2.0) * (1.0 + 1.0 / Q) - ((1.0 - e2) * P * z2) / (Q * (1.0 + Q)) - P * r2 / 2.0;
    const double ro = -(e2 * P * r) / (1.0 + Q) + (ro_radicand >= 0.0 ? sqrt(ro_radicand) : 0.0);
    const double tmp = (r - e2 * ro) * (r - e2 * ro);
    const double U = sqrt(tmp + z2);
    const double V = sqrt(tmp + (1.0 - e2) * z2);
    const double zo = (b2 * cart[2]) / (semi_major_axis() * V);

    *geod_alt = U * (1.0 - b2 / (semi_major_axis() * V));
    *geod_lat = (r > 0.0 ? atan((cart[2] + ep2 * zo) / r) : const_pi_2<double>());
    *lon = atan2(cart[1], cart[0]);
}

void ecm::geodetic_normal(double geod_lat, double lon, double* normal)
{
    double sinlat, coslat;
    double sinlon, coslon;
    ::sincos(geod_lat, &sinlat, &coslat);
    ::sincos(lon, &sinlon, &coslon);
    dvec3 n = normalize(dvec3(coslat * coslon, coslat * sinlon, sinlat));
    normal[0] = n.x;
    normal[1] = n.y;
    normal[2] = n.z;
}

void ecm::quad_to_ecm(int quad_side, int quad_level, int quad_x, int quad_y, double qx, double qy,
        double* ecm_x, double* ecm_y)
{
    assert(qx >= 0.0 && qx <= 1.0);
    assert(qy >= 0.0 && qy <= 1.0);
    int quads_in_level = (1 << quad_level);
    assert(quad_x >= 0 && quad_x < quads_in_level);
    assert(quad_y >= 0 && quad_y < quads_in_level);
    assert(quad_side >= 0 && quad_side < 6);
    side_to_ecm(quad_side,
            2.0 * ((quad_x                        + qx) / quads_in_level) - 1.0,
            2.0 * (((quads_in_level - 1 - quad_y) + qy) / quads_in_level) - 1.0,
            ecm_x, ecm_y);
}

void ecm::quad_to_ecm(int quad_side, int quad_level, int quad_x, int quad_y, int corner,
        double* ecm_x, double* ecm_y)
{
    dvec2 qxy;
    if (corner == corner_tl)
        qxy = dvec2(0.0, 1.0);
    else if (corner == corner_tr)
        qxy = dvec2(1.0, 1.0);
    else if (corner == corner_br)
        qxy = dvec2(1.0, 0.0);
    else
        qxy = dvec2(0.0, 0.0);
    quad_to_ecm(quad_side, quad_level, quad_x, quad_y, qxy.x, qxy.y, ecm_x, ecm_y);
}

void ecm::quad_plane(int quad_side, int quad_level, int quad_x, int quad_y,
        const double* quad_tl_cart, const double* quad_tr_cart, const double* quad_bl_cart, const double* quad_br_cart,
        double* plane_normal, double* plane_distance) const
{
    // The quad plane normal is the ellipsoid normal at the quad center.
    double quad_center_ecm_x, quad_center_ecm_y;
    quad_to_ecm(quad_side, quad_level, quad_x, quad_y, 0.5, 0.5, &quad_center_ecm_x, &quad_center_ecm_y);
    double quad_center_geod_lat, quad_center_lon;
    ecm_to_geodetic(quad_center_ecm_x, quad_center_ecm_y, &quad_center_geod_lat, &quad_center_lon);
    dvec3 quad_plane_normal;
    geodetic_normal(quad_center_geod_lat, quad_center_lon, quad_plane_normal.vl);
    // The quad plane distance is the smallest distance of the four quad corners to the origin.
    double quad_tl_distance = dot(quad_plane_normal, dvec3(quad_tl_cart));
    double quad_tr_distance = dot(quad_plane_normal, dvec3(quad_tr_cart));
    double quad_bl_distance = dot(quad_plane_normal, dvec3(quad_bl_cart));
    double quad_br_distance = dot(quad_plane_normal, dvec3(quad_br_cart));
    double quad_plane_distance = min(min(min(quad_tl_distance, quad_tr_distance), quad_bl_distance), quad_br_distance);
    assert(quad_plane_distance > 0.0);
    // Return the result
    plane_normal[0] = quad_plane_normal.x;
    plane_normal[1] = quad_plane_normal.y;
    plane_normal[2] = quad_plane_normal.z;
    *plane_distance = quad_plane_distance;
}

double ecm::max_quad_plane_distance_estimation(int quad_side, int quad_level, int quad_x, int quad_y,
        const double* quad_plane_normal, double quad_plane_distance) const
{
    // Use the distance at the quad center as the max distance estimation.
    double ecm_x, ecm_y;
    quad_to_ecm(quad_side, quad_level, quad_x, quad_y, 0.5, 0.5, &ecm_x, &ecm_y);
    dvec3 cart;
    ecm_to_cartesian(ecm_x, ecm_y, cart.vl);
    double estimation = dot(dvec3(quad_plane_normal), cart) - quad_plane_distance;
    assert(estimation >= 0.0);
    return estimation;
}

void ecm::symmetry_quad(int quad_side, int quad_level, int quad_x, int quad_y,
        int* sym_quad_side, int* sym_quad_level, int* sym_quad_x, int* sym_quad_y,
        bool* mirror_x, bool* mirror_y, float matrix[9])
{
    int sq_side = quad_side;
    int sq_level = quad_level;
    int sq_x = quad_x;
    int sq_y = quad_y;
    bool mx = false;
    bool my = false;
    mat3 M = mat3(1.0f);

    if (sq_side == ecm::side_bottom) {
        sq_side = ecm::side_top;
    } else if (sq_side != ecm::side_top && sq_side != ecm::side_front) {
        sq_side = ecm::side_front;
    }
    int quads_in_level = (1 << sq_level);
    if (sq_x >= quads_in_level / 2) {
        sq_x = quads_in_level - 1 - sq_x;
        mx = true;
        M[1][1] = -1.0f;
    }
    if (sq_y >= quads_in_level / 2) {
        sq_y = quads_in_level - 1 - sq_y;
        my = true;
        if (quad_side <= 3)
            M[2][2] = -1.0f;
        else
            M[0][0] = -1.0f;
    }
    // correction for the destination side
    if (quad_side == ecm::side_right) {
        mat3 m = mat3(
                 0.0f, -1.0f,  0.0f,
                +1.0f,  0.0f,  0.0f,
                 0.0f,  0.0f, +1.0f);
        M *= m;
    } else if (quad_side == ecm::side_back) {
        mat3 m = mat3(
                -1.0f,  0.0f,  0.0f,
                 0.0f, -1.0f,  0.0f,
                 0.0f,  0.0f, +1.0f);
        M *= m;
    } else if (quad_side == ecm::side_left) {
        mat3 m = mat3(
                 0.0f, +1.0f,  0.0f,
                -1.0f,  0.0f,  0.0f,
                 0.0f,  0.0f, +1.0f);
        M *= m;
    } else if (quad_side == ecm::side_bottom) {
        mat3 m = mat3(
                -1.0f,  0.0f,  0.0f,
                 0.0f, +1.0f,  0.0f,
                 0.0f,  0.0f, -1.0f);
        M *= m;
    }

    *sym_quad_side = sq_side;
    *sym_quad_level = sq_level;
    *sym_quad_x = sq_x;
    *sym_quad_y = sq_y;
    if (mirror_x)
        *mirror_x = mx;
    if (mirror_y)
        *mirror_y = my;
    if (matrix)
        std::memcpy(matrix, M.vl, 9 * sizeof(float));
}

static void go_left(int from_side, dvec2 from_sqxy, double d, int& side, dvec2& sqxy)
{
    if (from_side == ecm::side_top) {
        side = ecm::side_left;
        sqxy.x = 1.0 - from_sqxy.y;
        sqxy.y = 1.0 - d;
    } else if (from_side == ecm::side_bottom) {
        side = ecm::side_left;
        sqxy.x = from_sqxy.x;
        sqxy.y = d;
    } else {
        side = (from_side == ecm::side_front ? ecm::side_left : from_side - 1);
        sqxy.x = 1.0 - d;
        sqxy.y = from_sqxy.y;
    }
}

static void go_right(int from_side, dvec2 from_sqxy, double d, int& side, dvec2& sqxy)
{
    if (from_side == ecm::side_top) {
        side = ecm::side_right;
        sqxy.y = 1.0 - d;
        sqxy.x = from_sqxy.x;
    } else if (from_side == ecm::side_bottom) {
        side = ecm::side_right;
        sqxy.y = d;
        sqxy.x = 1.0 - from_sqxy.y;
    } else {
        side = (from_side == ecm::side_left ? ecm::side_front : from_side + 1);
        sqxy.x = d;
        sqxy.y = from_sqxy.y;
    }
}

static void go_bottom(int from_side, dvec2 from_sqxy, double d, int& side, dvec2& sqxy)
{
    if (from_side == ecm::side_top) {
        side = ecm::side_front;
        sqxy.x = from_sqxy.x;
        sqxy.y = 1.0 - d;
    } else if (from_side == ecm::side_bottom) {
        side = ecm::side_back;
        sqxy.x = 1.0 - from_sqxy.x;
        sqxy.y = d;
    } else if (from_side == ecm::side_front) {
        side = ecm::side_bottom;
        sqxy.x = from_sqxy.x;
        sqxy.y = 1.0 - d;
    } else if (from_side == ecm::side_right) {
        side = ecm::side_bottom;
        sqxy.x = 1.0 - d;
        sqxy.y = 1.0 - from_sqxy.x;
    } else if (from_side == ecm::side_back) {
        side = ecm::side_bottom;
        sqxy.x = 1.0 - from_sqxy.x;
        sqxy.y = d;
    } else {
        assert(from_side == ecm::side_left);
        side = ecm::side_bottom;
        sqxy.x = d;
        sqxy.y = from_sqxy.x;
    }
}

static void go_top(int from_side, dvec2 from_sqxy, double d, int& side, dvec2& sqxy)
{
    if (from_side == ecm::side_top) {
        side = ecm::side_back;
        sqxy.x = 1.0 - from_sqxy.x;
        sqxy.y = 1.0 - d;
    } else if (from_side == ecm::side_bottom) {
        side = ecm::side_front;
        sqxy.x = from_sqxy.x;
        sqxy.y = d;
    } else if (from_side == ecm::side_front) {
        side = ecm::side_top;
        sqxy.x = from_sqxy.x;
        sqxy.y = d;
    } else if (from_side == ecm::side_right) {
        side = ecm::side_top;
        sqxy.x = 1.0 - d;
        sqxy.y = from_sqxy.x;
    } else if (from_side == ecm::side_back) {
        side = ecm::side_top;
        sqxy.x = 1.0 - from_sqxy.x;
        sqxy.y = 1.0 - d;
    } else {
        assert(from_side == ecm::side_left);
        side = ecm::side_top;
        sqxy.x = d;
        sqxy.y = 1.0 - from_sqxy.x;
    }
}

void ecm::quad_base_data(int quad_side, int quad_level, int quad_x, int quad_y,
        const double* quad_tl_cart, const double* quad_tr_cart, const double* quad_bl_cart, const double* quad_br_cart,
        const double* quad_plane_normal, double quad_plane_distance,
        int quad_size, int quad_overlap,
        float* offsets, float* normals, double* max_quad_plane_distance) const
{
    const int quads_in_level = (1 << quad_level);

    assert(quad_level > 0);
#if 1
    // Force use of symmetry quads
    assert(quad_side == ecm::side_top || quad_side == ecm::side_front);
    assert(quad_x < quads_in_level / 2);
    assert(quad_y < quads_in_level / 2);
#endif

    *max_quad_plane_distance = -1.0;
    for (int y = 0; y < quad_size + 2 * quad_overlap; y++) {
        for (int x = 0; x < quad_size + 2 * quad_overlap; x++) {
            // Quad-relative coordinates.
            dvec2 qxy = dvec2(x - quad_overlap + 0.5, y - quad_overlap + 0.5) / static_cast<double>(quad_size);
            // Interpolated position (computed in the same way as in the vertex shader,
            // but in planetocentric coordinates with full precision).
            dvec3 p_interp =
                  dvec3(quad_bl_cart) * (1.0 - qxy.x) * (1.0 - qxy.y)
                + dvec3(quad_br_cart) * (      qxy.x) * (1.0 - qxy.y)
                + dvec3(quad_tr_cart) * (      qxy.x) * (      qxy.y)
                + dvec3(quad_tl_cart) * (1.0 - qxy.x) * (      qxy.y);
            // Side coordinates. These need to be adjusted for samples in the border area
            // since these belong to neighboring sides. Note that for the diagonal
            // directions (e.g. x==0,y==0 for bottom-left), we chose only one way to go
            // into the neighbor since the result of diagonal movement is ambiguous.
            int side = quad_side;
            dvec2 sqxy = dvec2(
                    (quad_x + qxy.x) / quads_in_level,
                    ((quads_in_level - 1 - quad_y) + qxy.y) / quads_in_level);
            if (unlikely(sqxy.x < 0.0)) {
                go_left(side, clamp(sqxy, 0.0, 1.0), -sqxy.x, side, sqxy);
            } else if (unlikely(sqxy.x > 1.0)) {
                go_right(side, clamp(sqxy, 0.0, 1.0), sqxy.x - 1.0, side, sqxy);
            } else if (unlikely(sqxy.y < 0.0)) {
                go_bottom(side, clamp(sqxy, 0.0, 1.0), -sqxy.y, side, sqxy);
            } else if (unlikely(sqxy.y > 1.0)) {
                go_top(side, clamp(sqxy, 0.0, 1.0), sqxy.y - 1.0, side, sqxy);
            }
            assert(sqxy.x >= 0.0 && sqxy.x <= 1.0);
            assert(sqxy.y >= 0.0 && sqxy.y <= 1.0);
            /*
            if (qxy.x >= 0 && qxy.x <= 1.0 && qxy.y >= 0.0 && qxy.y <= 1.0) {
                assert(sqxy.x >= 0.0 && sqxy.x <= 0.5);
                assert(sqxy.y >= 0.5 && sqxy.y <= 1.0);
                assert(sqxy.x >= quad[2] / static_cast<double>((1 << quad[1])));
                assert(sqxy.x <= (quad[2] + 1) / static_cast<double>((1 << quad[1])));
                assert(sqxy.y <= 1.0 - quad[2] / static_cast<double>((1 << quad[1])));
                assert(sqxy.y >= (quad[2] + 1) / static_cast<double>((1 << quad[1])));
            }
            */
            dvec2 sxy = 2.0 * sqxy - dvec2(1.0, 1.0);
            // Cartesian coordinates in full precision.
            dvec2 p_ecm;
            side_to_ecm(side, sxy.x, sxy.y, &(p_ecm.x), &(p_ecm.y));
            dvec3 p_cart;
            ecm_to_cartesian(p_ecm.x, p_ecm.y, p_cart.vl);
            if (x >= quad_overlap && y >= quad_overlap && x < quad_size + quad_overlap && y < quad_size + quad_overlap) {
                // assert(length(p_cart) >= length(p_interp));
            }
            // The offset between interpolated and full precision coordinates.
            dvec3 offset = p_cart - p_interp;
            // The ellipsoid normal.
            double geod_lat, lon;
            ecm_to_geodetic(p_ecm.x, p_ecm.y, &geod_lat, &lon);
            dvec3 normal;
            geodetic_normal(geod_lat, lon, normal.vl);
            // The quad plane distance.
            double qpd = dot(dvec3(quad_plane_normal), p_cart) - quad_plane_distance;
            if (x >= quad_overlap && y >= quad_overlap && x < quad_size + quad_overlap && y < quad_size + quad_overlap) {
                assert(qpd >= 0.0);
            }
            // Store the results.
            offsets[3 * (y * (quad_size + 2 * quad_overlap) + x) + 0] = offset.x;
            offsets[3 * (y * (quad_size + 2 * quad_overlap) + x) + 1] = offset.y;
            offsets[3 * (y * (quad_size + 2 * quad_overlap) + x) + 2] = offset.z;
            normals[2 * (y * (quad_size + 2 * quad_overlap) + x) + 0] = normal.x;
            normals[2 * (y * (quad_size + 2 * quad_overlap) + x) + 1] = normal.y;
            if (qpd > *max_quad_plane_distance)
                *max_quad_plane_distance = qpd;
        }
    }
}

const char* ecm_version(int* major, int* minor, int* patch)
{
    if (major)
        *major = ECM_VERSION_MAJOR;
    if (minor)
        *minor = ECM_VERSION_MINOR;
    if (patch)
        *patch = ECM_VERSION_PATCH;
    return ECM_VERSION;
}
