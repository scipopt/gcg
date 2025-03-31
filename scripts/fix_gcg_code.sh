#!/bin/bash

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
