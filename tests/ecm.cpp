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

#include <cstdio>

#include "msg.h"

#include "ecm/ecm.h"
#include "glvm-str.h"

using namespace glvm;

#define CHECK(condition) \
    if (!(condition)) \
    { \
        msg::err("%s:%d: %s: Check '%s' failed.", \
                __FILE__, __LINE__, __PRETTY_FUNCTION__, #condition); \
        exit(1); \
    }

int main(void)
{
    const double eps = 1e-5;
    class ecm ecm(3.0, 2.0);

    // cube face centers
    {
        dvec3 a, b, c;
        ecm.geodetic_to_cartesian(0.0, 0.0, 0.0, a.vl);
        ecm.ecm_to_cartesian(0.0, 0.0, b.vl);
        c = dvec3(3.0, 0.0, 0.0);
        CHECK(length(a - b) <= eps);
        CHECK(length(a - c) <= eps);
        ecm.geodetic_to_cartesian(0.0, radians(90.0), 0.0, a.vl);
        ecm.ecm_to_cartesian(2.0, 0.0, b.vl);
        c = dvec3(0.0, 3.0, 0.0);
        CHECK(length(a - b) <= eps);
        CHECK(length(a - c) <= eps);
        ecm.geodetic_to_cartesian(0.0, radians(180.0), 0.0, a.vl);
        ecm.ecm_to_cartesian(4.0, 0.0, b.vl);
        c = dvec3(-3.0, 0.0, 0.0);
        CHECK(length(a - b) <= eps);
        CHECK(length(a - c) <= eps);
        ecm.geodetic_to_cartesian(0.0, radians(-90.0), 0.0, a.vl);
        ecm.ecm_to_cartesian(6.0, 0.0, b.vl);
        c = dvec3(0.0, -3.0, 0.0);
        CHECK(length(a - b) <= eps);
        CHECK(length(a - c) <= eps);
        ecm.geodetic_to_cartesian(radians(90.0), 0.0, 0.0, a.vl);
        ecm.ecm_to_cartesian(0.0, +2.0, b.vl);
        c = dvec3(0.0, 0.0, 2.0);
        CHECK(length(a - b) <= eps);
        CHECK(length(a - c) <= eps);
        ecm.geodetic_to_cartesian(radians(-90.0), 0.0, 0.0, a.vl);
        ecm.ecm_to_cartesian(0.0, -2.0, b.vl);
        c = dvec3(0.0, 0.0, -2.0);
        CHECK(length(a - b) <= eps);
        CHECK(length(a - c) <= eps);
    }

    // cube corners: check if coordinates are identical from all three cube sides
    {
        dvec3 a, b, c;
        // front corners: tl, tr, br, bl
        ecm.ecm_to_cartesian(-1.0 + eps * eps, +1.0 - eps * eps, a.vl);    // front
        ecm.ecm_to_cartesian(-1.0 + eps * eps, +1.0 + eps * eps, b.vl);    // top
        ecm.ecm_to_cartesian(+7.0 - eps * eps, +1.0 - eps * eps, c.vl);    // left
        CHECK(length(a - b) <= eps);
        CHECK(length(b - c) <= eps);
        ecm.ecm_to_cartesian(+1.0 - eps * eps, +1.0 - eps * eps, a.vl);    // front
        ecm.ecm_to_cartesian(+1.0 - eps * eps, +1.0 + eps * eps, b.vl);    // top
        ecm.ecm_to_cartesian(+1.0 + eps * eps, +1.0 - eps * eps, c.vl);    // right
        CHECK(length(a - b) <= eps);
        CHECK(length(b - c) <= eps);
        ecm.ecm_to_cartesian(+1.0 - eps * eps, -1.0 + eps * eps, a.vl);    // front
        ecm.ecm_to_cartesian(+1.0 - eps * eps, -1.0 - eps * eps, b.vl);    // bottom
        ecm.ecm_to_cartesian(+1.0 + eps * eps, -1.0 + eps * eps, c.vl);    // right
        CHECK(length(a - b) <= eps);
        CHECK(length(b - c) <= eps);
        ecm.ecm_to_cartesian(-1.0 + eps * eps, +1.0 - eps * eps, a.vl);    // front
        ecm.ecm_to_cartesian(-1.0 + eps * eps, +1.0 + eps * eps, b.vl);    // bottom
        ecm.ecm_to_cartesian(+7.0 - eps * eps, +1.0 - eps * eps, c.vl);    // left
        CHECK(length(a - b) <= eps);
        CHECK(length(b - c) <= eps);
        // back corners: tl, tr, br, bl
        ecm.ecm_to_cartesian(+3.0 + eps * eps, +1.0 - eps * eps, a.vl);    // back
        ecm.ecm_to_cartesian(+1.0 - eps * eps, +3.0 - eps * eps, b.vl);    // top
        ecm.ecm_to_cartesian(+3.0 - eps * eps, +1.0 - eps * eps, c.vl);    // right
        CHECK(length(a - b) <= eps);
        CHECK(length(b - c) <= eps);
        ecm.ecm_to_cartesian(+5.0 - eps * eps, +1.0 - eps * eps, a.vl);    // back
        ecm.ecm_to_cartesian(-1.0 + eps * eps, +3.0 - eps * eps, b.vl);    // top
        ecm.ecm_to_cartesian(+5.0 + eps * eps, +1.0 - eps * eps, c.vl);    // left
        CHECK(length(a - b) <= eps);
        CHECK(length(b - c) <= eps);
        ecm.ecm_to_cartesian(+5.0 - eps * eps, -1.0 + eps * eps, a.vl);    // back
        ecm.ecm_to_cartesian(-1.0 + eps * eps, -3.0 + eps * eps, b.vl);    // bottom
        ecm.ecm_to_cartesian(+5.0 + eps * eps, -1.0 + eps * eps, c.vl);    // left
        CHECK(length(a - b) <= eps);
        CHECK(length(b - c) <= eps);
        ecm.ecm_to_cartesian(+3.0 + eps * eps, -1.0 + eps * eps, a.vl);    // back
        ecm.ecm_to_cartesian(+1.0 - eps * eps, -3.0 + eps * eps, b.vl);    // bottom
        ecm.ecm_to_cartesian(+3.0 - eps * eps, -1.0 + eps * eps, c.vl);    // right
        CHECK(length(a - b) <= eps);
        CHECK(length(b - c) <= eps);
    }

    // quad corners level 1
    {

    }

    return 0;
}
