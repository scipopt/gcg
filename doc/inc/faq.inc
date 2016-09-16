<div id="faq">

    <h3><a name="Section0">Sections</a></h3>
    <ul>
     <li><a href="#Section1">General Questions about GCG</a></li>
     <li><a href="#Section2">Decomposition file related questions</a></li>
     <li><a href="#Section3">Miscellaneous</a></li>
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
      GCG (see  <a href="#Q3">next question</a>) or supplied as a decomposition
      file (see <a href="group__FILEREADERS.html">filereaders</a>).
   </p>

    <p class="q">
      <a name="Q3">What are the detected structures?</a>
    </p>
    <p class="a">
      GCG contains several detectors to find a block diagonal structures in
      the constraint matrix, pure ones as well as such with a
      border. Additionally, some detectors are able to detect staircase
      structures where no linking constraints are present and linking variables
      only occur between two consecutive blocks. <br />
      With applications such as bin packing, capacitated p-median and
      generalized assignment problems in mind, we furthermore included a
      detector that only assigns set partitioning and set covering constraints
      to the master problem.
   </p>

    <p class="q">
      <a name="Q3a">Why are some detectors switch off by default?</a>
    </p>
    <p class="a">
      While several detectors are only experimental, others need external
      software.
      Some detectors rely
      on <a href="http://glaros.dtc.umn.edu/gkhome/metis/hmetis/overview">hMETIS</a>
      and only work under Unix based systems. To use them, you must either
      install hMETIS from the package sources (if available there) or download it
      and put the executable in a directory contained in the <tt>PATH</tt>
      environment variable.
      The Isomorphism detector needs Bliss (see <a href="#Q8">here</a>).
   </p>

    <h4><a name="Section2">Decomposition file related questions</a></h4>
    <p class="q">
      <a name="Q4">Which filereaders can read decompositions?</a>
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
      <a name="Q5">Why is presolving important for the decomposition?</a>
    </p>
    <p class="a">
      As GCG uses several presolving methods from SCIP, the transformed problem
      (see also SCIP FAQ) may change significantly from the original
      problem. Variables or constraints may be deleted or added which renders
      the current decomposition invalid or unsuitable in some cases. GCG does
      some basic sanity checks, however, it doesn't handle all problems and may
      crash if the decomposition is read at the wrong time (e.g. a
      decomposition found after presolving is read before the problem is
      presolved).<br />
      Use only non presolved decompositions and disable presolving if you are
      in doubt!
    </p>

    <p class="q">
      <a name="Q6">After reading the decomposition file, GCG tells me that the
	behaviour is undefined.</a>
    </p>
    <p class="a">
      If GCG can not detect whether your decomposition is for the presolved or
      the original problem, it cannot guard you against errors. See the
      <a href="#Q5">previous question</a> for more information.
    </p>

    <h4><a name="Section3">Miscellaneous</a></h4>
    <p class="q">
      <a name="Q7">Why is CTRL-C unsafe to use?</a>
    </p>
    <p class="a">
       In its current version, SCIP will finish the current node when pressing
    CTRL-C. This is a problem if the master problem is not completely finished
    solving. The current node will be marked as finished and two branching
    children are created. Depending on when you press CTRL-C, this may and
    will lead to different solving processes. If an optimal solution is found,
    they will be the same in all runs, however, it might take a substantially
    different amount of time to do so.
    </p>

    <p class="q">
      <a name="Q8">How does aggregation work in GCG?</a>
    </p>
    <p class="a">
       Detecting whether (some) pricing problems can be aggregated is done using
      <a href="http://www.tcs.hut.fi/Software/bliss/">Bliss</a>, an external
      tool for computing automorphism groups of graphs. After a decomposition has been
      selected, we check for each pair of pricing problems if they are identical
      (and their variables have the same coefficients in the same master constraints)
      by creating an auxiliary graph and trying to find a suitable automorphism
      using Bliss. Identical pricing problems are then aggregated.
      See the <a href="INSTALL.html">installation information </a>
      for an instruction on how to include bliss.<br />

       If Bliss is not included, aggregation in GCG is done very simple. If the
      same variables appear in the same order with the same coefficients in the same
      constraints for all pricing problems (and with the same coefficients in the
      master constraints), GCG will detect that and aggregate the pricing problems.
    </p>
</div>

