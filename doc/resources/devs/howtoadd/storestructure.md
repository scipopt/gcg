# Storing structure information (deprecated) {#storestructure}
> **This page is currently being refactored. Some things might still be outdated.**

struct_decomp.h is responsible for storing structure information. The memory has to be allocated by caller and is freed
later. You can either create the decomposition by calling DECcreateDecompFromMasterconss() or DECfilloutDecompFromConstoblock()
or use the getter and setter functions ins pub_decomp.h to fill this structure.
 *
Very quick until more elaborate, these are the relevant fields of the structure and its basic usage:
 - <code>subscipcons   </code>
  - an array of array of constraints in each block
  - Usage: <code>subscipcons[b][c]</code> is constraint <em>c</em> in block <em>b</em>
 - <code>nsubscipconss </code>
  - an array of the number of constraints in each block
  - Usage: <code>nsubscipcons[b]</code> gives you the number of constraints of block <em>b</em>
 - <code>subscipvars   </code>
  - an array of arrays of variables in each block
  - Usage: <code>subscipvars[b][v]</code> is variable <em>v</em> in block <em>b</em>
 - <code>nsubscipvars  </code>
  - an array of the number of variables in each block
  - Usage: <code>nsubscipvars[b]</code> gives you the number of variables of block <em>b</em>
 - <code>nblocks       </code>
  - number of blocks/pricing problems in the matrix resp. reformulation
 - <code>type          </code>
  - Type of the decomposition (DEC_DECTYPE_ARROWHEAD is the most general)
  - Current supported types: DEC_DECTYPE_DIAGONAL, DEC_DECTYPE_BORDERED, DEC_DECTYPE_ARROWHEAD, DEC_DECTYPE_STAIRCASE,
    DEC_DECTYPE_UNKNOWN
 - <code>constoblock   </code>
  - SCIP_HASHMAP linking constraints to blocks
  - Usage: <code>SCIPhashmapGetImage(constoblock, cons)</code> returns <em>b+1</em>, where <em>b</em> is the block of
    constraint cons. This map is somehow inverse of the subscipconss array.
 - <code>vartoblock    </code>
  - SCIP_HASHMAP linking variables to blocks
  - Usage: <code>SCIPhashmapGetImage(vartoblock, var)</code> returns <em>b+1</em>, where <em>b</em> is the block of
    variable var. This map is somehow inverse of the subscipcvars array.
 - <code>linkingvars   </code>
  - array of linking variables (to be in the master)
  - Usage: <code>linkingvars[v]</code> is linking variable number <em>v</em>
 - <code>nlinkingvars  </code>
  - number of linking variables
 - <code>linkingconss  </code>
  - array of linking constraints (to be in the master)
  - Usage: <code>linkingcons[c]</code> is linking constraint number <em>c</em>
 - <code>nlinkingconss </code>
  - number of linking constraints
