# How to build this documentation {#doc}
<!-- The very first line of the .md document should be the page title and {# name of the site}
     #doc is the name of this site. Needed when making a link to this site. -->

# Building this Documentation
## Requirements
- Rights to write inside the `doc` directory
- common bash tools (`egrep`, `sed`, ...)
- Doxygen and graphviz (for graph visualizations)
- PHP (for FAQ generation)
- Python 2.7 or newer with `re`, `sys` and `argparse` importable (for FAQ generation)
- a working LaTeX installation (for formulas inside the use cases)
- a fully linked and functional default installation of the most recent GCG
(used to read the current interactive menu and settings file)
- a javascript-enabled browser

### Additional Notes
**For macOS Users (tested with Catalina):**<br>
To compile the SCIP documentation, please install gnu-sed and add it to your path. 
(For the GCG documentation, this is not needed.)
```
brew install gnu-sed
PATH="/usr/local/opt/gnu-sed/libexec/gnubin:$PATH"
```

## Instructions
The documentation for your installed GCG version including its code can be
generated using a simple

    make doc

inside the GCG root directory. Alternatively, one can also navigate to the
folder `doc` and execute the corresponding script

    ./builddoc.sh

Please always make sure that the Requirements (as mentioned above) are met for your
installation.

## Important Information
### Subpage Indexing
Subpage indexing for Howtos (use and add), as well as all classifiers
and detectors is being built automatically. If one is added, it will
be made a subpage of the respective parent page and automatically linked
on it upon the next `make doc`.
### FAQ
The FAQ are generated using PHP. The questions can be added in the file
`/doc/resources/misc/faq/faqtext.txt`.


## How to add new pages
### Folder structure
You can move and create folders and pages (see below) as you like inside the folder `resources`,
but remember to add new folder paths to INPUT in gcg.dxy (vicinity of line 787).
Note that the pages will not automatically be subpages of the containing folder's main page.
To make this the case, you have to add a tag `\@subpage` to the containing page.

### .md files
The first line has to be of the form
```
# your title {#identifier}
```
in the right folder.\n
**Important**: For all pages that are referenced by the subpage indexing scripts
(see above), the filename has to match the page identifier (xy.md <-> `\{\#xy\}`).

Some of the main features of .md files are explained in the following.
#### HTML compatibility
The Markdown format is based on HTML, with more restricted, but better readable functionalities.
Thus, most HTML tags (such as `iframe`s, anchors, embeds, etc.) can be used in the .md file.
If you intend on only adding HTML code, use the tag `\htmlinclude HTMLFILE`.

#### Referencing
You can reference other pages using `@ref`.

#### Adding Pictures
Copy your picture to doc/pictures/. If used in a md file, it will automatically
be moved to the html/ folder. Add a tag `![Pic Title](filename.png)`.
If you need to control the image size, use HTML: `\image html picture.png "caption width=50% `

# Markdown Basics #
In the following, you will be presented code snippets and their result after compiling
the markdown code.

Code:

    ```C
    int main()
    ```

Result:

```C
int main()
```
\n
Code:

    `cd ..`

Result: `cd ..`

\n
Code:

    ## Title 1
    ### Subtitle 1
    #### Subsubtitle 1

Result:

## Title 1
### Subtitle 1
#### Subsubtitle 1

\n
Code:

    **This is the title of a section**

Result: **This is the title of a section**

\n
Code:

    @ref doc "Title of the link to this page"

Result: @ref doc "Title of the link to this page"

\n
Code:

    [This is the title of the link to Google](https://www.google.com)

Result: [This is the title of the link to Google](https://www.google.com)

\n
Code:

    <a href="https://www.google.com">Another Link to Google</a>

Result: <a href="https://www.google.com">Another Link to Google</a>


__For more information on Basic Markdown Syntax, visit [this page](https://help.github.com/en/articles/basic-writing-and-formatting-syntax).__
