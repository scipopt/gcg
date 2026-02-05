/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2026 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with GCG; see the file LICENSE. If not visit gcg.or.rwth-aachen.de.*/
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

#define GCG_VERSION         410 /**< GCG version number (multiplied by 100 to get integer number) */
#define GCG_SUBVERSION        0 /**< GCG sub version number */

#endif
