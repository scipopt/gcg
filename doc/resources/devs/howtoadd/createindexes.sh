# Adds new .md and .dox pages to the index on the "How to ..." page
# and therefore makes them subpages of the aforementioned

cd $(dirname $0)
# Remove ending and subpages of .dox file
sed -i'' '/\*\//d' howtoadd.dox
sed -i'' '/\\subpage /d' howtoadd.dox

# Get index list and append to .dox
ls | egrep '\.md$|\.dox$' | sed 's/.dox//' | sed 's/.md//' | sed 's/^/- \\subpage /' | sed '/howto/d' >> howtoadd.dox
echo "*/" >> howtoadd.dox

echo "Subpage indexing for devs/howtoadd/ built sucessfully."
