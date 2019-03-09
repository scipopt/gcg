# Section 1 {#section}
-\subpage subpage_1

# This is the title of the Sandbox page {#subpage_1}
<!-- The very first line of the .md document should be the page title and {# name of the site}
     #sandbox is the name of this site. Needed when making a link to this site. -->

<!-- This is a piece of code. C indicates that C language syntax is recognized. -->
```C
int main()
```
<!-- This is a piece of inline code. -->
`cd ..`

**This is the title of a section**
- \ref sandbox      "Title of the link to sandbox"
- \ref assdgda      "This is the title of a non-existing link"
- [**Title of the link to sandbox**](sandbox.html)
- [Link to install](INSTALL.html)
- [This is the title of the link to Google](https://www.google.com)
- <a href="https://www.google.com">Google</a>

images can be easily included with the `!` command:

![instance MIPLIB2003/rout with 5 blocks](../../img/rout-cC-5.png)

we need to discuss where images are located; meanwhile I use the `../../img/` folder

if you need to control the image size, the only way appears to use `HTML` code:

```
<img src="../../img/rout-cC-5.png" alt="alt text" width="50%"></img>
```
\subpage subpage_1_1
