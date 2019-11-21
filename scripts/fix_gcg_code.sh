#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *#
#*                                                                           *#
#*                  This file is part of the program                         *#
#*          GCG --- Generic Column Generation                                *#
#*                  a Dantzig-Wolfe decomposition based extension            *#
#*                  of the branch-cut-and-price framework                    *#
#*         SCIP --- Solving Constraint Integer Programs                      *#
#*                                                                           *#
#* Copyright (C) 2010-2019 Operations Research, RWTH Aachen University       *#
#*                         Zuse Institute Berlin (ZIB)                       *#
#*                                                                           *#
#* This program is free software; you can redistribute it and#or             *#
#* modify it under the terms of the GNU Lesser General Public License        *#
#* as published by the Free Software Foundation; either version 3            *#
#* of the License, or (at your option) any later version.                    *#
#*                                                                           *#
#* This program is distributed in the hope that it will be useful,           *#
#* but WITHOUT ANY WARRANTY; without even the implied warranty of            *#
#* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *#
#* GNU Lesser General Public License for more details.                       *#
#*                                                                           *#
#* You should have received a copy of the GNU Lesser General Public License  *#
#* along with this program; if not, write to the Free Software               *#
#* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.*#
#*                                                                           *#
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *#

# lowercase first character in doxygen comment
sed -i 's,\/\*\*\(\ *\)\([A-Z]\)\([a-z]\),\/\*\*\1\L\2\3,' src/*.{c,h}

# fix parameter comments
sed -i 's,\(\ *\)\([A-Za-z].*\)\(\ *\)\([A-Za-z].*\)\(\ *\)\/\*\*\ \([A-Za-z]\),\1\2\3\4\5\/\*\*<\ \6,' src/*.{c,h}

# fix for,while loops and if, switch clauses
sed -i 's/for(\ *\(.*[^\ ]\)\ *)$/for(\ \1\ )/' src/*.c
sed -i 's/while(\ *\(.*[^\ ]\)\ *)$/while(\ \1\ )/' src/*.c
sed -i 's/switch(\ *\(.*[^\ ]\)\ *)$/switch(\ \1\ )/' src/*.c
sed -i 's/if(\ *\(.*[^\ ]\)\ *)$/if(\ \1\ )/' src/*.c

# wrap all allocs in SCIP_CALL
sed -i 's/^\(\ *\)SCIPalloc\(.*\)(\(.*\));$/\1SCIP_CALL( SCIPalloc\2(\3) );/' src/*.c

# fix SCIP calls
sed -i 's/SCIP_CALL(\ *\(.*\)(\(.*\))\ *);$/SCIP_CALL( \1(\2) );/' src/*.c

# delete all trailing whitespaces
sed -i 's,\([^ ]\)\ *$,\1,' Makefile src/*.{c,cpp,h} check/*.{awk,sh}
