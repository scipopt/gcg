<div id="faq">

    <h4><a name="Section0">Sections</a></h4>
    <ul>
     <li><a href="#Section1">General Questions about GCG</a></li>
     <li><a href="#Section2">Using SCIP as a standalone MIP-Solver</a></li>
     <li><a href="#Section3">Using SCIP included in another source code</a></li>
     <li><a href="#Section4">Using SCIP as a Branch-Cut-And-Price-Framework</a></li>
     <li><a href="#Section5">Specific questions about Column Generation and Branch-And-Price with SCIP</a></li>
     <li><a href="#Section6">Known Bugs</a></li>
    </ul>

    <h4><a name="Section1">General Questions about GCG</a></h4>
    <p class="q">
      <a name="Q1">What is GCG?</a>
    </p>

    <p class="a">
      GCG is a generic solver for Dantzig-Wolfe reformulation
      using <b>branch-cut-and-price</b>. It contains all code to serve as
      a <b>standalone solver</b> for structured mixed integer programs (MIPs)
      and the accompanying structure information.

      This FAQ contains common questions as well as starting information and known bugs.
   </p>

    <p class="q">
      <a name="Q2">What MIPs can GCG solve?</a>
    </p>
    <p class="a">
      GCG can solve all problems for which a structure for Dantzig-Wolfe
      reformulation is known. The structure is either automatically detected by
      GCG (see  <a href="#Q2">next question</a>) or supplied as a decomposition
      file (see <a href="group__FILEREADERS.html">filereaders</a>).
   </p>

    <p class="q">
      <a name="Q2">What are the detected structures?</a>
    </p>
    <p class="a">
      GCG has built in detection code for block diagonal problems and for
      classical set partitioning or setcovering master problems. In order to
      benefit from the automatic detection, your problem must contain only set
      partitioning or set covering constraints in the master and non set
      covering or set partitioning constraints in each of the pricing
      subproblems. <br>
      After removing all set partitioning or covering constraints from the
      problem, the pricing problems need to be disjoint, e.g. they share no
      common variables. In these cases, GCG can detect the structure
      automatically. <br>
      This works for bin packing, coloring, capacitated p-median, generalized
      assignment problems as well as for cutting stock, in which the set
      partitioning is relaxed to set partitioning with integer right hand side.
   </p>

    <h4><a name="Section2">Decomposition file related questions</a></h4>
    <p class="q">
      <a name="Q1">Which filereaders can read decompositions</a>
    </p>
    <p class="a">
      Currently GCG reads three different decomposition structure information:
      <ul>
	<li>dec (<a href="reader__dec_8h.html">documentation</a>)</li>
	<li>blk (<a href="reader__blk_8h.html">documentation</a>)</li>
	<li>ref (<a href="reader__ref_8h.html">documentation</a>)</li>
      </ul>
    </p>

    <p class="q">
      <a name="Q1">Why is presolving important for the decomposition</a>
    </p>
    <p class="a">
      As GCG uses several presolving methods from SCIP, the transformed problem
      (see also SCIP FAQ) may change significantly from the original
      problem. Variables or constraints may be deleted or added which renders
      the current decomposition invalid or unsuitable in some cases. GCG does
      some basic sanity checks, however, it doesn't handle all problems and may
      crash if the decomposition is read at the wrong time (e.g. a
      decomposition found after presolving is read before the problem is
      presolved).<br>
      Use only non presolved decompositions and disable presolving if you are
      in doubt!
    </p>

    <p class="q">
      <a name="Q1">After reading the decomposition file, GCG tells me that the
	behaviour is undefined</a>
    </p>
    <p class="a">
      If GCG can not detect whether your decomposition is for the presolved or
      the original problem, it cannot guard you against errors. See the
      <a href="#Q??"previous question</a> for more information.
    </p>

    <h4><a name="Section6">Known Bugs</a></h4>
    <p class="q">
      <a name="Q1">Known Bugs</a>
    </p>
    <p class="a">
      GCG contains bugs which should be reported
      per <a href="mailto:gcg-bugs@or.rwth-aachen.de?subject=New GCG bug&body=Please describe your bug here.">mail</a>.<br>
      The following bugs are known to the developers and will be fixed in
      future bugfix versions:
      <ul>
	<li>Reading in the wrong decomposition leads to a crash</li>
	<li>Presolved decompositions are unstable with respect to SCIP
	  improvements</li>
	<li>Reading presolved decompositions should cause GCG to presolve the problem</li>
	<li>Unsafe presolving reductions should be disabled if a decomposition
	  is read before presolving</li>
	<li>GCG can not recover from some numerical errors and may not fail
	  gracefully</li>
	<li>GCG does not fail gracefully on most bogus inputs</li>
	<li>GCG may consume a lot of memory</li>
      </ul>
    </p>


</div>

