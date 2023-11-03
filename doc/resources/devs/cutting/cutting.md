# The GCG Cutting {#cutting}
> The **generation and application of cuts** is, among with the pricing, the most important component in solving when
> using a Branch-and-Price-and-Cut approach, since it can decrease the size of the feasible region significantly and
> lead to faster finding of the optimal solution.

During the **execution of Branch-and-Price-and-Cut**, GCG, as the name suggests, also **applies cuts**.
In the following, we will present what separators (which generate cutting planes for your problem) are implemented in 
GCG, along with some theory and implementational details.

@subpage sepa_basis "Basic Separator" \n
@subpage sepa_master "Master Separator" \n
@subpage sepa_subsetrow "Subset Row Cuts Separator" \n

<hr>

In order to **create your own separator**, please consider SCIP's "How to add" for that. There should be no bigger
differences when doing it for GCG specifically.

[**How to add Separators (SCIP guide)**](https://scipopt.org/doc/html/SEPA.php) \n
