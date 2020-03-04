# Your own Classifier {#own-classifier}
> **This page is still in development and may be incomplete. Please excuse any inconveniences.**

# Properties of a Classifier {#CLS_PROPERTIES}
\par DEC_CLASSIFIERNAME: name of classifier

\par DEC_DESC: short description of classification

\par DEC_PRIORITY: priority of classifier

\par DEC_ENABLEDORIG: classify on original problem

\par DEC_ENABLEDPRESOLVED: classify on presolved problem

# Interface Methods {#CLS_INTERFACE}
In order to make use of the classifier, you have to include it in SCIP.
```C++
SCIP_RETCODE SCIPincludeConsClassifierXYZ(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_CALL( DECincludeConsClassifier(scip, DEC_CLASSIFIERNAME, DEC_DESC, DEC_PRIORITY, DEC_ENABLEDORIG, DEC_ENABLEDPRESOLVED, classifierInit, classifierFree, classifierClassify) );

   return SCIP_OKAY;
}
```

# Fundamental Callback Methods of a Detector {#DEC_FUNDAMENTALCALLBACKS}
These methods have to be static.

## DEC_DECL_CONSCLASSIFY


# Additional Callback Methods of a Classifier {#CLS_ADDITIONALCALLBACKS}

## classifierFree
## classifierInit
