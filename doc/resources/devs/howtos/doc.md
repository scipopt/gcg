# document in this documentation {#doc}
<!-- The very first line of the .md document should be the page title and {# name of the site}
     #sandbox is the name of this site. Needed when making a link to this site. -->

# Documentation Readme #

Note that this page will be moved to somewhere in the Developer's Guide before
publication of the new documentation.

**Updates:**
- landingpage cleanup
  - sections that are unnecessary for most users removed or made less obvious
  - updated some stuff (e.g. removed "student assistant") and made some things
    consistent with the rest of the pages on the SCIP site
  - renaming to "User's Guide"/"Developer's Guide" with link to documentation
- documentation cleanup/preparation
  - creation of the basic structures so that new subpages can be added easily
  - script to index howto- and detectors- pages automatically (also added to builddoc.sh)

**Why's there a mix of .md and .dox files?**
- markdown does not seem to support hierarchical functions of doxygen,
therefore .dox files are used for all non-leaf pages

**How can I add new pages?**
- If you don't expect your page to get any subpages:
Create an .md file with the first line being
```
# your title # {#identifier}
```
in the right folder.
IMPORTANT: The filename must match the page identifier (xy.md <-> "{#xy}")
- Otherwise, use a .dox file.

**How can I structure it folderwise?**
- You can move and create folders as you like inside 'resources',
but remember to add your folder path to INPUT in gcg.dxy (line 787).

**How can I add pictures to my page?**
- Copy your picture to doc/pictures/. If used in a md/dox file, it will automatically
be moved to the html/ folder.
- Add a tag `![Pic Title](filename.png)`
- If you need to control the image size, use HTML: `<img src="../../img/rout-cC-5.png" alt="alt text" width="50%" />`

# Markdown Basics #
This is a piece of code. <!-- C indicates that C language syntax is recognized. -->
```C
int main()
```
This is a piece of inline code.
`cd ..`

**This is the title of a section**
- \ref doc "Title of the link to this page"
- [This is the title of the link to Google](https://www.google.com)
- <a href="https://www.google.com">Another Link to Google</a>

__For more information on Basic Markdown Syntax, visit [this page](https://help.github.com/en/articles/basic-writing-and-formatting-syntax).__
