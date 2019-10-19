# Colors Detector {#det-colors}
> **This page is still in development and may be incomplete. Please excuse any inconveniences.**

The detector assigns color classes to constraints and tries combinations of colors in the master.

### Details

To detect decompositions the detector uses the following steps:

1. Each constraints gets assigned to a color class. Constraints having the same lhs, rhs and conshdlr get assigned to the same color class.
2. For each set X of color classes with cardinality k with 2 <= k <= nbits, where nbits is set to 2, a decomposition with constraints that belong to a color class in X in the master is created.

### Parameters
```
no parameters; we should introduce a parameter to determine the value of nbits
```
### Future work
Test the detector with nbits > 2. Try different types of color classes (maybe also depending on number of non-zeros or something else)

### Links
 * Documentation: http://www.or.rwth-aachen.de/gcg/doc/dec__colors_8cpp.html
