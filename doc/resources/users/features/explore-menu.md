# The Explore Menu {#explore-menu}

## Enter the explore menu
First, you have to `read` and `optimize`, or at least `detect` your problem. Then, like with all commands, just type `explore` in the interactive shell. GCG will show
an overview of all found decompositions:
```
GCG> explore

==================================================================================
Summary              presolved       original
                     ---------       --------
detected                     6              0
==================================================================================
   nr   id nbloc nmacon nlivar nmavar nstlva spfwh history  pre nopcon nopvar  sel
 ---- ---- ----- ------ ------ ------ ------ ----- ------- ---- ------ ------ ----
    0   13    50     50      0      0      0 0.745      cC  yes      0      0   no
    1   12    50     50      0     50      0 0.245      cC  yes      0      0   no
    2   15   100      0   2500      0      0 0.009      vC  yes      0      0   no
    3   16     1      0     50      0      0 0.006      vC  yes      0      0   no
    4   11     1      0      0      0      0 0.000       C  yes      0      0   no
    5   14     1      1      0      0      0 0.000      nC  yes      0      0   no
Please enter command or decomposition id to select (or "h" for help) :
```

When typing `h`, you will get a list of possible commands:
```
GCG/explore> help
==================================================================================

List of selection commands

                       command     description
                       -------     -----------
                        select     selects/unselects decomposition with given id
                      previous     displays the preceding decompositions (if there are any)
                          next     displays the subsequent decompositions (if there are any)
                           top     displays the first decompositions
                           end     displays the last decompositions
                        legend     displays the legend for table header and history abbreviations
                          help     displays this help
                number_entries     modifies the number of displayed decompositions
                     visualize     visualizes the specified decomposition (requires gnuplot)
                       inspect     displays detailed information for the specified decomposition
                     set_score     sets the score by which the "goodness" of decompositions is evaluated
                      sort_asc     sets whether to sort (by score) in ascending or descending order
                          quit     return to main menu

==================================================================================
```

### Select structures for solving with `select`
You can choose one of the decompositions to
be used for solving with the `select <id>` command.
```
GCG/explore> select 2
==================================================================================
Please specify the nr of the decomposition to be selected:

==================================================================================
Summary              presolved       original
                     ---------       --------
detected                     6              0
==================================================================================
   nr   id nbloc nmacon nlivar nmavar nstlva spfwh history  pre nopcon nopvar  sel
 ---- ---- ----- ------ ------ ------ ------ ----- ------- ---- ------ ------ ----
    0   13    50     50      0      0      0 0.745      cC  yes      0      0   no
    1   12    50     50      0     50      0 0.245      cC  yes      0      0   no
    2   15   100      0   2500      0      0 0.009      vC  yes      0      0  yes
    3   16     1      0     50      0      0 0.006      vC  yes      0      0   no
    4   11     1      0      0      0      0 0.000       C  yes      0      0   no
    5   14     1      1      0      0      0 0.000      nC  yes      0      0   no
Please enter command or decomposition id to select (or "h" for help) :

```
After that, you can visualize (see below) the decomposition or also leave the
explore menu (`..`) and optimize the problem. If you don't select a decomposition,
GCG will always use the 'best' one, i.e. nr. 0 in the example above.

### Visualize your decomposition with `visualize`
With the `visualize` command, you can export a visualization of the selected
decomposition. Prior to that, you will have to install "gnuplot" by typing
`sudo apt-get install gnuplot` into your console.

### Use a different score with `set_score`
By default, GCG uses the so-called max white score. This means that the decomposition with the most zero entries is used. Since this might not always be the best measure, it is possible to choose a different score:
```
GCG/explore> set_score
==================================================================================
Please specify the new score:
0: max white,
1: border area,
2: classic,
3: max foreseeing white,
4: ppc-max-white,
5: max foreseeing white with aggregation info,
6: ppc-max-white with aggregation info,
7: experimental benders score
8: strong decomposition score
Note: Sets the detection/scoretype parameter to the given score.
```
