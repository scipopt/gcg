# Adds new .md and .dox pages to the index on the "How to ..." page
# and therefore makes them subpages of the aforementioned

cd $(dirname $0)
# Remove ending and subpages of .dox file
sed -i '/\*\//d' ./detectors.dox
sed -i '/\\subpage /d' ./detectors.dox

# Get index list and append to .dox
ls | egrep '\.md$|\.dox$' | sed 's/.dox//' | sed 's/.md//' | sed 's/^/- \\subpage /' | sed '/detectors/d' >> detectors.dox
echo "*/" >> detectors.dox

echo "Subpage indexing for devs/detectors/ built sucessfully."
