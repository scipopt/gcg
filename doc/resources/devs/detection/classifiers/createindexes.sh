# Adds new .md and .dox pages to the index on the "How to ..." page
# and therefore makes them subpages of the aforementioned

cd $(dirname $0)

cd consclass
echo "# Constraint Classifier {#clscons}" > clscons.md
ls | egrep '\.md$|\.dox$' | sed 's/clscons\.md//' | sed 's/.md//' | sed 's/^/- @subpage /' | sed '/howto/d' >> clscons.md
cd ..

cd varclass
echo "# Variable Classifier {#clsvar}" > clsvar.md
ls | egrep '\.md$|\.dox$' | sed 's/clsvar\.md//' | sed 's/.md//' | sed 's/^/- @subpage /' | sed '/howto/d' >> clsvar.md
cd ..

echo "Subpage indexing for devs/detection/ built sucessfully."
