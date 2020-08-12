# Detection Parameters {#detection-params}
# Modifying Parameters of the GCG Detection 
> GCG performs its detection using @ref classifiers and @ref detectors that
> classify different variable and constraint types, respectively detect structures each.

There are multiple reasons to deactivate and activate classifiers and detectors.
Most often, you want to...
  * **Activate** classifiers/detectors that are deactivated by default (and deactivate the respective other ones)
  Activating and deactivating classifiers and detectors can be particularly helpful not only to find specific structures 
  but also to test the behavior of your own detector with all others deactivated.
  * **Change priorities** classifiers/detectors to find a specific structure in your program
  (see @ref own-detector "How to add a detector" and @ref own-classifier "How to add a classifier" for more information)

### Activate and Deactivate Classifiers and Detectors
Sometimes you might be tempted to deactivate all detectors except for one. Note that you can always choose a decomposition
using @ref explore-menu, but then there might be decompositions where different detectors (e.g. also ones you don't want to use) were used. Thus, before detecting and solving your problem, you will have to type `set detection detectors` to completely
deactivate detectors. Once you are in that submenu, you can select detectors, e.g. `consclass` and disable them with `enabled FALSE`.
You will see the following feedback:
```
GCG/set/detection/detectors/consclass> enabled FALSE
detection/detectors/consclass/enabled = FALSE
```
You can also selectively deactivate one of the two more functionalities except for propagating 
(finishing/postprocessing, see @ref detectors "here"), then you will have to use
`{finishing,postprocessing}enabled FALSE`. \n

If you want to directly deactivate _all_ detectors (e.g. to activate one single one again), 
you can also do that using the `off` function of the emphasis setting (see "Emphasis Settings").

### Change Priorities of Classifiers and Detectors
We use the same commands as above to enter the settings for a detector, again e.g. `consclass`. Then, instead of disabling or enabling
them, we can also tell GCG to change the priority of that detector, e.g. with `priority 42`. You will be able to observe the following:
```
GCG/set/detection/detectors/consclass> priority
current value: 0, new value [-2147483648,2147483647]: 42
detection/detectors/consclass/priority = 42
```
This detector will now be called before others (see @ref detection).

### Emphasis Settings
Another method to change parameters of the detection is by using our predefined emphasis settings. 
There are three different settings that we deliver and an `off` setting that will deactivate all functionalities
of all detectors.
```
  aggressive            sets detection <aggressive>
  default               sets detection <default>
  fast                  sets detection <fast>
  off                   turns <off> all detectors
```
The different emphasis settings will modify parameters from the activation of specified detectors to the number of rounds they are called (see @ref detection).

### Change Parameters of the Detection
You can also modify other parameters of the detection. 
You can set them by entering `set` and then `detection`. 
In the following, we give a list of all parameters that can then be changed.
\n
A list of all detection parameters can be seen (and searched) when inserting `set detection` into the search box on the page of @ref interactive-menu.