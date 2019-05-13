# Adds new .md and .dox pages to the index on the "How to ..." page
# and therefore makes them subpages of the aforementioned

cd $(dirname $0)
# Remove ending and subpages of .dox file
sed -i '/\*\//d' ./howto.dox
sed -i '/\\subpage /d' ./howto.dox

# Get index list and append to .dox
ls | egrep '\.md$|\.dox$' | sed 's/.dox//' | sed 's/.md//' | sed 's/^/- \\subpage /' | sed '/howto/d' >> howto.dox
echo "*/" >> howto.dox

echo "Subpage indexing for devs/howtos/ built sucessfully."
