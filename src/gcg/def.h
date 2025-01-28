/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2025 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/* This program is free software; you can redistribute it and/or             */
/* modify it under the terms of the GNU Lesser General Public License        */
/* as published by the Free Software Foundation; either version 3            */
/* of the License, or (at your option) any later version.                    */
/*                                                                           */
/* This program is distributed in the hope that it will be useful,           */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with this program; if not, write to the Free Software               */
/* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.*/
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file    def.h
 * @brief   common defines and data types used in all packages of GCG
 * @author  Stefanie Koß
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_DEF_H_
#define GCG_DEF_H_

/*
 * Add some macros for differing functions on Windows
 */
#if defined(_WIN32) || defined(_WIN64)
#define realpath(x,y) _fullpath(y,x,SCIP_MAXSTRLEN)
#endif

#if defined(_WIN32) || defined(_WIN64)
#define fileno(x) _fileno(x)
#endif

/*
 * Define the macro GCG_EXPORT if it is not included from the generated header
 */
#ifndef GCG_EXPORT
#if defined(_WIN32) || defined(_WIN64)
#define GCG_EXPORT __declspec(dllexport)
#elif defined(__GNUC__) && __GNUC__ >= 4
#define GCG_EXPORT __attribute__((__visibility__("default")))
#else
#define GCG_EXPORT
#endif
#endif

#define GCG_VERSION         371 /**< GCG version number (multiplied by 100 to get integer number) */
#define GCG_SUBVERSION        0 /**< GCG sub version number */

#endif
