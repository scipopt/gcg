# The Interactive Menu {#interactive-menu}

Once you started GCG, you can navigate through the menu. All expandable points of the menu are written inside brackets (`<display>`, `<set>`, ...) and everything else is just a keyword (`read`, `quit`, ...).
Switching back to a higher level inside the menu can be done by simply typing `..`, going to the top level by just pressing enter. Finally, there are also real submenus (@ref master-menu "the master menu", @ref explore-menu "the explore menu"), which you have to leave with `quit`.

Most often, you will need to change settings. This can be done using `set` or `fix` (parameters for `fix` are the same as for `set`, thus they are omitted in the list below).
If you use `fix`, the parameters will not be changeable throughout the execution of GCG. If you use `set`, GCG will
sometimes modify some parameters (e.g. for the @ref presolving). This also results in the `set diffsave` command to not only write the parameters you changed (they appear at the bottom), but also
some that GCG will always change by default.

> You can click on an item to show its description. When searching, the command as well as its description will be searched for.

\htmlinclude menu.html
