# Settings for Reports {#report-settings}
## Giving a Report Settings File
In every case, the settings file, e.g. `settings.vset`, should simply look like this:
```
# GCG script settings file
DRAFT=true
TREE=false
...
```
Comments are possible.

## Test Set Report Settings
### General settings
**For variables with boolean type, we require them to be set to either `false` or `true`. No capital letters or numbers are supported.**
|Variable|Type|Default Value|Description|
|---|---|---|---|
|`DRAFT`|Boolean|`false`|Draft mode (reduced report length)|
|`REPORTDIR`|Path|`check/reports/report_<TEST>_<SET>_<TIME>`|Output folder for report|
|`TREE`|Boolean|`true`|Plot toggle|
|`TIME`|Boolean|`true`|Plot toggle|
|`DETECTION`|Boolean|`true`|Plot toggle|
|`BOUNDS`|Boolean|`true`|Plot toggle|
|`PRICING`|Boolean|`true`|Plot toggle|
|`TREEARGS`|String|-|Script arguments (see @ref tree-plotter "here")|
|`TIMEARGS`|String|-|Script arguments (see @ref general-plotter "here")|
|`DETECTIONARGS`|String|-|Script arguments (see @ref detection-plotter "here")|
|`BOUNDSARGS`|String|-|Script arguments (see @ref bounds-plotter "here")|
|`PRICINGARGS`|String|-|Script arguments (see @ref pricing-plotter "here")|
|`CUSTOMCLASSIFIER`|String|`nonzeros`|Classifier to plot detection visualizations for|
|`CUSTOMDETECTOR`|String|`SetPartMaster`|Detector to plot detection visualizations for (currently not implemented)|

\n
### Additional settings for manual execution
**Parameters that only have to be defined when generating the report manually.** \n
X: mandatory for all build systems\n
C: should be given (otherwise "unknown") when using CMake

|Required|Variable|Type|Description|
|:-:|---|---|---|
|C|`BINARY_ID`|File|Name of your GCG Binary (e.g. `gcg-3.1.0.linux.x86_64.gnu.opt.spx2`).|
|C|`VERSION`|String|GCG version to be printed on front page.|
|C|`LPS`|String|LP solver used during your test run to be printed on front page.|
|C|`THREADS`|Number|Threads used during your test run to be printed on front page.|
|C|`MODE`|String|Test mode (e.g. `readdec`) to be printed on front page.|
|C|`FEASTOL`|String|Feasibility tolerance setting (default: `default`) to be printed on front page.|
|C*|`LAST_STATISTICS`|Boolean|Flag to indicate whether you compiled GCG with `STATISTICS=true`.|
|X|`OUTFILE`|File|Path (absolute or relative to GCG root directory) to your `out` file.|
|X|`RESFILE`|File|Path (absolute or relative to GCG root directory) to your `res` file.|
|X|`VBCFILES`|Folder|Path (absolute or relative to GCG root directory) to your `vbc` files folder for this run.|
| |`TSTNAME`|String|Name of your test set to be printed on front page.|
| |`SETNAME`|String|Name of your settings to be printed on front page.|
| |`TIMELIMIT`|Number|Time limit of your test run to be printed on front page.|
| |`MEMLIMIT`|Number|Memory limit of your test run to be printed on front page.|
| |`NODELIMIT`|Number|Node limit of your test run to be printed on front page.|

*Set to true by default for CMake.

## Comparison Report Settings {#comparison-report-settings}
**For variables with boolean type, we require them to be set to either `false` or `true`. No capital letters or numbers are supported.**

### General Settings
|Variable|Type|Default Value|Description|
|---|---|---|---|
|`DATADIR`|Directory|-|Must be the directory including your `.out` and `.res` files (will be asked for if not given).|
|`DRAFT`|Boolean|`false`|Draft mode (reduced report length)|
|`REPORTDIR`|Path|`check/reports/report_<TEST>_<SET>_<TIME>`|Output folder for report|
|`TABLE`|Boolean|`true`|Comparison Table toggle|
|`PERFPROF`|Boolean|`true`|Plot toggle|
|`GENERAL`|Boolean|`true`|Plot toggle|
|`TIME`|Boolean|`true`|Plot toggle|
|`DETECTION`|Boolean|`true`|Plot toggle|
|`BOUNDS`|Boolean|`true`|Plot toggle|
|`PERFPROFARGS`|Boolean|-|Script arguments (see @ref performance-profile-plotter "here")|
|`TIMEARGS`|String|-|Script arguments (see @ref general-plotter "here")|
|`DETECTIONARGS`|String|-|Script arguments (see @ref detection-plotter "here")|
|`BOUNDSARGS`|String|-|Script arguments (see @ref bounds-plotter "here")|

### Advanced Settings
|Variable|Type|Default Value|Description|
|---|---|---|---|
|`MAXNINSTANCES_FRONTPAGE`|Integer|30|Max number of instances s.t. table printed on front page|
|`MAXNINSTANCES_SECONDPAGE`|Integer|50|Max number of instances s.t. table printed on second page|
|`CUSTOMCLASSIFIER`|String|`nonzeros`|Classifier to plot detection visualizations for|
|`CUSTOMDETECTOR`|String|`SetPartMaster`|Detector to plot detection visualizations for (currently not implemented)|
|`LAST_STATISTICS`|Boolean|as compiled|Flag indicating whether you compiled with `STATISTICS` flag (if false, no bounds visus will be generated)|