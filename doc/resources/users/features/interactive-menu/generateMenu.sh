# Script that creates menu documentation
cd resources/users/features/interactive-menu
rm menu.html

echo "Getting GCG interactive menu structure."
python3 getMenu.py

echo "Compiling interactive menu documentation."
cat menu_start.html.in  > menu.html
cat menu.txt            >> menu.html
cat menu_end.html.in    >> menu.html

rm menu.txt
