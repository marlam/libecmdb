/*
 * C++ vector and matrix classes that resemble GLSL style.
 *
 * Copyright (C) 2011, 2012
 * Martin Lambers <marlam@marlam.de>
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

#ifndef GLVM_STR_H
#define GLVM_STR_H

#include "str.h"

#include "glvm.h"

namespace str
{
    template<typename T, int rows, int cols>
    std::string from(const glvm::array<T, rows, cols>& x)
    {
        std::string s = "[ ";
        for (int i = 0; i < rows * cols; i++) {
            s += str::from(x.vl[i]);
            if (i < rows * cols - 1)
                s += ' ';
        }
        s += ']';
        return s;
    }
}

#endif
