# Visualization arguments {#visu-args}
On this page, you find all arguments (e.g. the overall time or the pricing calls made)
that are parsed by the parsers of the visualization suite.
# General Parser {#general-args}

The following arguments are contained inside a dataframe. \n

### Times
| **time** | **type** | 
|---|---|
|<code>OVERALL TIME</code> | duration |
|<code>READING TIME</code> | duration |
|<code>COPYING TIME</code> | duration |
|<code>PRESOLVING TIME</code> | duration |
|<code>DETECTION TIME</code> | duration |
|<code>ORIGINAL LP TIME</code> | duration |
|<code>RMP LP TIME</code> | duration |
|<code>PRICING TIME</code> | duration |
|<code>PRICING SOLVER TIME</code> | duration |
|<code>FARKAS TIME</code> | duration |
|<code>HEURISTICS TIME</code> | duration |
|<code>CUTS TIME</code> | duration |
|<code>ROOT NODE TIME</code> | duration |
|<code>SOLUTIONS FIRST</code> | point in time | 
|<code>SOLUTIONS BEST</code> | point in time | 

### General Data
| **data** | **type** | **possible values** | **fallback** |
|---|---|---|---|
|<code>MASTER TIME</code> | float | \f$[0, \infty[\f$ | NaN|
|<code>STATUS</code> | enum | \f$\{0, \ldots, 5\}\f$ | \f$-1\f$|
|<code>HEURISTICS CALLS</code> | int | \f$[0, \infty[\f$ | \f$-1\f$|
|<code>HEURISTICS FOUND</code> | int | \f$[0, \infty[\f$ | \f$-1\f$|
|<code>CUTS CALLS</code> | int | \f$[0, \infty[\f$ | \f$-1\f$|
|<code>CUTS FOUND</code> | int | \f$[0, \infty[\f$ | \f$-1\f$|
|<code>CUTS APPLIED</code> | int | \f$[0, \infty[\f$ | \f$-1\f$|
|<code>MIP OR KNAPSACK</code> | int | \f$[1, 3]\f$ | \f$-1\f$|
|<code>CONS LINEAR</code> | int | \f$[0, \infty[\f$ | \f$-1\f$|
|<code>CONS KNAPSACK</code> | int | \f$[0, \infty[\f$ | \f$-1\f$|
|<code>CONS LOGICOR</code> | int | \f$[0, \infty[\f$ | \f$-1\f$|
|<code>CONS SETPPC</code> | int | \f$[0, \infty[\f$ | \f$-1\f$|
|<code>CONS VARBOUND</code> | int | \f$[0, \infty[\f$ | \f$-1\f$|
|<code>LINKING VARS</code> | int | \f$[0, \infty[\f$ | \f$-1\f$|
|<code>NBLOCKS</code> | int | \f$[0, \infty[\f$ | \f$-1\f$|
|<code>AGGREGATIONS</code> | int | \f$[0, \infty[\f$ | \f$-1\f$|
|<code>SOLUTIONS FOUND</code> | int | \f$[0, \infty[\f$ | \f$-1\f$|
|<code>PD INTEGRAL</code> | float | \f$[0, \infty[\f$ | NaN
|<code>MASTER NCONSS</code> | int | \f$[0, \infty[\f$ | \f$-1\f$|
|<code>MASTER NVARS</code> | int | \f$[0, \infty[\f$ | \f$-1\f$|
|<code>BNB TREE NODES</code> | int | \f$[0, \infty[\f$ | \f$-1\f$|
|<code>BNB TREE DEPTH</code> | int | \f$[0, \infty[\f$ | \f$-1\f$|
|<code>BNB TREE LEFT</code> | int | \f$[0, \infty[\f$ | \f$-1\f$|
|<code>MASTER LP ITERATIONS</code> | int | \f$[0, \infty[\f$ | \f$-1\f$|
|<code>MASTER LP CALLS</code> | int | \f$[0, \infty[\f$ | \f$-1\f$|
|<code>ORIGINAL LP CALLS</code> | int | \f$[0, \infty[\f$ | \f$-1\f$|
|<code>ORIGINAL LP ITERATIONS</code> | int | \f$[0, \infty[\f$ | \f$-1\f$|
|<code>LP FILE</code> | string | | \f$-1\f$|
|<code>DEC FILE</code> | string | | \f$-1\f$|

# Bounds Parser {#bounds-args}

The following arguments are contained inside a dataframe. \n
| **data** | **type** | **possible values** | **fallback** |
|---|---|---|---|
|<code>iter</code> | integer | \f$[0, \infty[\f$ | -|
|<code>pb</code> | float | \f$]-\infty, \infty[\f$ | NaN|
|<code>db</code> | float | \f$]-\infty, \infty[\f$ | NaN|
|<code>time</code> | float | \f$[0, \infty[\f$ | NaN|
|<code>dualdiff</code> | float | \f$]-\infty, \infty[\f$ | NaN|
|<code>dualoptdiff</code> | float | \f$]-\infty, \infty[\f$ | NaN|
|<code>gap</code> | float | \f$]-\infty, \infty[\f$ | NaN|
|<code>nlpvars</code> | float | \f$[0, \infty[\f$ | NaN|
|<code>nlpvars\_cum</code> | float | \f$[0, \infty[\f$ | NaN|
|<code>lpvars</code> | float | \f$[0, 1]\f$ | NaN|
|<code>nipvars</code> | float | \f$[0, \infty[\f$ | NaN|
|<code>nipvars\_cum</code> | float | \f$[0, \infty[\f$ | NaN|
|<code>ipvars</code> | float | \f$[0, 1]\f$ | NaN|
|<code>time\_count</code> | integer | \f$[1, \infty[\f$ | NaN|
|<code>time\_first</code> | integer | \f$[0, \infty[\f$ | NaN|
|<code>time\_diff</code> | float | \f$[0, \infty[\f$ | NaN|
|<code>db\_ma</code> | float | \f$]-\infty, \infty[\f$ | NaN|


# Pricing Parser {#pricing-args}

The following arguments are contained inside a dataframe. \n
| **data** | **type** | **possible values** | **fallback** |
|---|---|---|---|
|<code>incumbent_times</code> | float[] | \f$]-\infty, \infty[\f$ | -|
|<code>rootlpsol_times</code> | float[] | \f$]-\infty, \infty[\f$ | -|
|<code>info</code> |dictionary| | |
|<code>pricing_data</code> |dictionary| | |
|<code>root_bounds</code> |dictionary| | |

Contents of the dictionary <code>info</code>:\n
| **data** | **type** | **possible values** | **fallback** |
|---|---|---|---|
|<code>instance</code> | string | | "unknown"|
|<code>settings</code> | string | | "unknown"|
|<code>status</code> | string | | "unknown"|

Contents of the dictionary <code>pricing_data</code>:\n
| **data** | **type** | **possible values** | **fallback** |
|---|---|---|---|
|<code>node</code> | integer | \f$[0, \infty[\f$ | (empty)|
|<code>pricing_round</code> | integer | \f$[0, \infty[\f$ | NaN|
|<code>stab_round</code> | integer | \f$[0, \infty[\f$ | NaN|
|<code>round</code> | integer | \f$[0, \infty[\f$ | NaN|
|<code>pricing_prob</code> | integer | \f$[0, \infty[\f$ | NaN|
|<code>time</code> | float | \f$[0, \infty[\f$ | NaN|
|<code>nVars</code> | integer | \f$[0, \infty[\f$ | NaN|
|<code>farkas</code> | bool | \f$[0, 1]\f$ | False|

Contents of the dictionary <code>root_bounds</code>:\n
| **data** | **type** | **possible values** | **fallback** |
|---|---|---|---|
|<code>iter</code> | integer | \f$[0, \infty[\f$ | -|
|<code>pb</code> | float | \f$]-\infty, \infty[\f$ | NaN|
|<code>db</code> | float | \f$]-\infty, \infty[\f$ | NaN|
|<code>time</code> | float | \f$[0, \infty[\f$ | NaN|
|<code>dualdiff</code> | float | \f$]-\infty, \infty[\f$ | NaN|
|<code>dualoptdiff</code> | float | \f$]-\infty, \infty[\f$ | NaN|
