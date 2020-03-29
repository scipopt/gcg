# Structure Types {#structure-types}

The following structure types are being detected in the current release (c.f. DEC_DECTYPE):

- DEC_DECTYPE_UNKNOWN: unknown structure (used for initialization)
- DEC_DECTYPE_ARROWHEAD: arrowhead structure (linking variables and constraints)
- DEC_DECTYPE_STAIRCASE: staircase structure (linking variables between consecutive blocks)
- DEC_DECTYPE_DIAGONAL: block diagonal structure (no linking variables and constraints)
- DEC_DECTYPE_BORDERED: bordered block diagonal structure (linking constraints only)


# Arrowhead {#arrowhead}
The arrowhead structure is a double bordered decomposition. It consists of
master constraints and linking variables.

@todo add decomp picture here

# Staircase {#staircase}
So-called staircase structures define decompositions where each block is connected
with the next one via linking variables. This often happens for programs that
possess a time dependency, where one time step influences the following.

@todo add decomp picture here

# Block-Diagonal {#block-diagonal}

# Single-Bordered {#single-bordered}
This type is a subtype of the arrowhead structure, but either without linking variables, or without linking constraints.

@todo add decomp picture here
