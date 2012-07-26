# this file is licensed under the GNU Lesser General Public License and part of the program GCG
# Copyright 2012 Martin Bergner

# lowercase first character in doxygen comment
sed -i 's,\/\*\*\(\ *\)\([A-Z]\)\([a-z]\),\/\*\*\1\L\2\3,' src/*.{c,h}

# fix parameter comments
sed -i 's,\(\ *\)\([A-Za-z].*\)\(\ *\)\([A-Za-z].*\)\(\ *\)\/\*\*\ \([A-Za-z]\),\1\2\3\4\5\/\*\*<\ \6,' src/*.{c,h}
