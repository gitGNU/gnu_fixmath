#  Automake file for the Fixmath src dir.

#  Copyright (C) 2005-2011, Axis Communications AB, LUND, SWEDEN
#
#  This file is part of Fixmath.
#
#  Fixmath is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published
#  by the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  You can use the comments under either the terms of the GNU Lesser General
#  Public License version 3 as published by the Free Software Foundation,
#  either version 3 of the License or (at your option) any later version, or
#  the GNU Free Documentation License version 1.3 or any later version
#  published by the Free Software Foundation; with no Invariant Sections, no
#  Front-Cover Texts, and no Back-Cover Texts.
#  A copy of the license is included in the documentation section entitled
#  "GNU Free Documentation License".
#
#  Fixmath is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License and a copy of the GNU Free Documentation License along
#  with Fixmath. If not, see <http://www.gnu.org/licenses/>.

# Include the exported API headers
AM_CPPFLAGS = -I$(top_srcdir)/include -I../include

# The library to build
lib_LTLIBRARIES = libfixmath.la

# Sources
libfixmath_la_SOURCES = fixmath.c      \
                        fixmath_ver.c  \
                        fixmath_impl.c \
                        fixmath_str.c

# Set the library version
libfixmath_la_LDFLAGS = -no-undefined \
                        -version-info @FX_LIB_CUR@:@FX_LIB_REV@:@FX_LIB_AGE@
