# Visualization arguments {#visu-args}
On this page, you find all arguments (e.g. the overall time or the pricing calls made)
that are parsed by the parsers of the visualization suite.
# General Parser {#general-args}
The general parser is the most commonly used one and is e.g. used for the [strIPlib](https://striplib.or.rwth-aachen.de) pipeline.

In the following, we first present all time arguments and their inclusion. Afterwards, a list of all runtime features contained inside a dataframe. They are also exported when executing the general parser as `feature_descriptions.json`.\n

### Time Inclusion
\image html VisualizationArguments.png "Illustration of all plottable times and their overlaps." width=90%

### Data
\htmlinclude visu_args_table.html


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
