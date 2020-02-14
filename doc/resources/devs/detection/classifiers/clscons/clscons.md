# Constraint Classifiers {#clscons}
Constraint classifiers can store in `ConsClassDecompInfo` how the constraints of a class should be handled by a detector:
 * BOTH -> not specified
 * ONLY_MASTER -> assigned to master problem
 * ONLY_PRICING -> assigned to pricing problem

 The following constraint classifiers are available:
 - @subpage clscons_scipconstypes
 - @subpage clscons_miplibconstypes
 - @subpage clscons_nnonzeros
 - @subpage clscons_consnamenonumbers
 - @subpage clscons_consnamelevenshtein
