# Adds new .md and .dox pages to the index on the "How to ..." page
# and therefore makes them subpages of the aforementioned

cd $(dirname $0)

cd consclass
echo "# Constraint Classifier {#clscons}" > clscons.md
ls | egrep '\.md$|\.dox$' | sed '/clscons\.md/d' | sed 's/.md//' | sed 's/^/- @subpage /' >> clscons.md
cd ..

cd varclass
echo "# Variable Classifier {#clsvar}" > clsvar.md
ls | egrep '\.md$|\.dox$' | sed '/clsvar\.md/d' | sed 's/.md//' | sed 's/^/- @subpage /' >> clsvar.md
cd ..

echo "Subpage indexing for devs/classifiers/ built sucessfully."
