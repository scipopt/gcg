# Detection Scores {#detection-scores}
> **This page is still in development and may be incomplete. Please excuse any inconveniences.**

# Scoring of Decompositions in GCG
As explained in @ref detection, GCG is decomposing problems, reformulating each to finally solve the whole 
problem more efficiently. The algorithmics are generating a **high number of different decompositions** 
(see @ref detection-process).
The challenge is to choose **which one of them is the most promising**, i.e. that can be reformulated and solved
best. In many cases, this desired decomposition coincides with the model structure that users probably already
had in mind when they created their problem. Without that knowledge (e.g. coming form one of our interfaces to
modeling languages), we deploy **scores to rate decompositions** for their suitability.

## List of Scores
Take the scores implemented in GCG with a **big grain of salt**: They all have existing theory in mind, however,
it is just theory. None of our scores is currently well-tested and thus all scores can be declared as rather
experimental. Our default score is the set partitioning foreseeing max white score (`spfwh`).\n
If you want to implement your own score, have a look at the guide @ref own-score.

### Border Area (`border`)
Minimizes the area that the borders take:
\f$ \text{s}_\text{borderarea} = 1 - \frac{\text{borderarea}}{\text{totalarea}} \f$.

### Classic (`classi`)
```
1 - (alphaborderarea * ( borderscore ) + alphalinking * ( linkingscore ) + alphadensity * ( densityscore ) )
```

### Max White Score
In the following, we present all variants of the **white area score** that, if we maximize it, will 
**maximize the non-block and non-border area** (see @ref structure-types for more information on those).
It can be seen as if the white space (zeroes in the coefficient matrix) is maximized, which is illustrated very 
nicely when visualizing decompositions using @ref explore-menu (for an example, see @ref structure-types).

#### Foreseeing Max White Score (`forswh`/`fawh` [with aggregation info])
This score is the unchanged variant of the max white score as described above. Stair-linking variables count 
as usual linking variables.

#### Set Partitioning Foreseeing Max White Score (`spfwh`/`spfawh` [with aggregation info])
This score encourages set partitioning master constraints through combining the usual foreseeing max white score with a 
boolean score rewarding a master containing only set partitioning and cardinality constraints.

### Experimental Benders Score (`bender`)
This score should always and exlusively be used when using the @ref benders "Benders mode". Its calculation is as follows:

\f{align}{ 
    \text{s}_\text{benders} =& \max(0, 1- (\text{s}_\text{blockarea} + (\text{s}_\text{borderarea} - \text{s}_\text{bendersarea})))\\
    \text{with} \\
    \text{s}_\text{blockarea} =& 1- \frac{\text{blockarea}}{\text{totalarea}}\\
    \text{s}_\text{borderarea} =& 1 - \frac{\text{borderarea}}{\text{totalarea}}\\
    \text{s}_\text{bendersarea} =& \frac{\text{bendersarea}}{\text{totalarea}}\\
    \text{bendersarea} =& \text{nmasterconshittingonlyblockvars} \cdot \text{nblockvarshittingNOmasterconss} + \\
                        & \text{nlinkingvarshittingonlyblockconss} \cdot \text{nblockconsshittingonlyblockvars} - p\\
    p=&\sum_{b=1}^{\text{nblocks}} \sum_{\text{blockvars }b_v\text{ of block }b\text{ hitting a master constraint}} \sum_{\text{all blocks }b_2 != b} \text{nblockcons}(b_2)
\f}

### Strong Decomposition Score (`strode`)
The strong decomposition score is not implemented for decompositions belonging to the original (non-presolved) problem.