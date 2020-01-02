# Adds new .md and .dox pages to the index on the "How to ..." page
# and therefore makes them subpages of the aforementioned

cd $(dirname $0)
# Remove ending and subpages of .dox file
sed -i '/\*\//d' howtouse.dox
sed -i '/\\subpage /d' howtouse.dox

# Get index list and append to .dox
ls | egrep '\.md$|\.dox$' | sed 's/.dox//' | sed 's/.md//' | sed 's/^/- \\subpage /' | sed '/howto/d' >> howtouse.dox
echo "*/" >> howtouse.dox

echo "Subpage indexing for devs/howtouse/ built sucessfully."
