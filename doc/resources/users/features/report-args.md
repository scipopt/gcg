# Arguments for Reports {#report-args}
## Test Set Report
In every case, the settings file, e.g. `settings.scset`, should simply look like this:
```
# GCG script settings file
DRAFT=true
TREE=false
...
```
Comments are possible.

> For variables with binary type, we require them to be set to either `false` or `true`. No capital letters or numbers are supported. 

|Variable|Type|Default Value|Description|
|---|---|---|---|
|`DRAFT`|Binary|`false`|Draft mode (reduced report length)|
|`TREE`|Binary|`true`|Plot toggle|
|`TIME`|Binary|`true`|Plot toggle|
|`DETECTION`|Binary|`true`|Plot toggle|
|`BOUNDS`|Binary|`true`|Plot toggle|
|`PRICING`|Binary|`true`|Plot toggle|
|`TREE`|String|-|Script arguments (see tree-plotter)|
|`TIME`|String|-|Script arguments (see @ref general-plotter "here")|
|`DETECTION`|String|-|Script arguments (see @ref detection-plotter "here")|
|`BOUNDS`|String|-|Script arguments (see @ref bounds-plotter "here")|
|`PRICING`|String|-|Script arguments (see @ref pricing-plotter "here")|
|`CUSTOMCLASSIFIER`|String|`nonzeros`|Classifier to plot detection visualizations for|
|`CUSTOMDETECTOR`|String|`SetPartMaster`|Detector to plot detection visualizations for (currently not implemented)|

> Parameters that only have to be defined when generating the report manually. All of them are not initialized otherwise!

|Required|Variable|Type|Description|
|:-:|---|---|---|
| |`TSTNAME`|String|Name of your test set to be printed on front page.|
| |`SETNAME`|String|Name of your settings to be printed on front page.|
|X|`OUTFILE`|File|Path (absolute or relative to GCG root directory) to your `out` file.|
|X|`OUTFILE`|File|Path (absolute or relative to GCG root directory) to your `out` file.|
|X|`RESFILE`|File|Path (absolute or relative to GCG root directory) to your `res` file.|
|X|`VBCFILES`|Folder|Path (absolute or relative to GCG root directory) to your `vbc` files folder for this run.|
|X|`LAST_STATISTICS`|Binary|Flag to indicate whether you compiled GCG with `STATISTICS=true`.|
|-|`BINARY_ID`|File|Name of your GCG binary (e.g. `gcg-3.1.0.linux.x86_64.gnu.opt.spx2`).|
|-|`TIMELIMIT`|Number|Time limit of your test run to be printed on front page.|
|-|`MEMLIMIT`|Number|Memory limit of your test run to be printed on front page.|
|-|`THREADS`|Number|Threads used during your test run to be printed on front page.|
|-|`FEASTOL`|String|Feasibility tolerance setting (default: `default`) to be printed on front page.|
|-|`VERSION`|String|GCG version to be printed on front page.|
|-|`MODE`|String|Test mode (e.g. `readdec`) to be printed on front page.|

## Comparison Report