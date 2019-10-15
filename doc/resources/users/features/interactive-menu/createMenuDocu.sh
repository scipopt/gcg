# Script that creates menu documentation
# Params: menu.txt
cd resources/users/features/interactive-menu
cat menu.html.in > menu.html

sed '/^$/d' menu.txt | sed -e "s/</\&lt;/g" | sed -e "s/>/\&gt;/g" | awk '{$2="</code>:" OFS $2} 1' | sed -e "s/\ <\/code>/<\/code>/" | sed -e "s/^/<li><a href=\"#\"><code>/" | sed -e "s/$/<\/a><\/li>/" >> menu.html
#sed '/^$/d' menu.txt | sed -e "s/</\&lt;/g" | sed -e "s/>/\&gt;/g" | sed -e "s/[a-z+<+>]*/<\/code>/" | sed -e "s/\ <\/code>/<\/code>/" | sed -e "s/^/<li><a href=\"#\"><code>/" | sed -e "s/$/<\/a><\/li>/" >> menu.html
#sed '/^$/d' menu.txt | sed -e "s/</\&lt;/g" | sed -e "s/>/\&gt;/g" | sed -e "s/^/<li><a href=\"#\">/" | sed -e "s/$/<\/a><\/li>/" >> menu.html

echo "</ul>" >> menu.html
echo "</div>" >> menu.html
