<!-- this file is used to generate the local doxygen documentation -->
<!-- using make doc from within scip or soplex -->

<style>
.reveal:hover {
    text-shadow: 1px 1px 1px #777;
}
.answer {
    background-color: #fff;
    padding-left:   1em;
}
</style>

<h3><a class="reveal_faq" href="#faq_generalquestionsaboutgcg"><span class="fa fa-caret-right"></span> General Questions about GCG</a></h3><ol>  <li>
    <a class="reveal_faq" href="#whatisgcg">
      What is GCG?    </a>
  </li>
    <li>
    <a class="reveal_faq" href="#whenusegcg">
      When should I use GCG?    </a>
  </li>
    <li>
    <a class="reveal_faq" href="#licensefaq">
      I heard something about licenses. Do I have to pay for using GCG?    </a>
  </li>
    <li>
    <a class="reveal_faq" href="#howtogetstarted">
      How do I get started?    </a>
  </li>
    <li>
    <a class="reveal_faq" href="#createsettingsfile">
      How do I create my own settings file?    </a>
  </li>
  </ol>
<hr />
<h3 id="faq_generalquestionsaboutgcg" class="anchor"><span class="fa fa-caret-right"></span> General Questions about GCG<a href="#" class="pull-right"><span title="go to top" class="fa fa-caret-up"></span></a></h3><ol>  <li id="whatisgcg" class="anchor">
    <h4>
      <a class="reveal_faq" href="#whatisgcg">What is GCG?</a>

      <a href="#" class="pull-right"><span class="fa fa-caret-up" title="go to top"></span></a>
    </h4>
    <div id="whatisgcg_ans" class="answer">
      <p>
GCG is...
</p>
<ul>
<li>You can install GCG using the <a href="./install.html">installation guide</a> and solve your own problem. To read your problem, learn more about the <a href="./input-formats.html">file formats supported by GCG</a>,
and get <a href="./WHATPROBLEMS.shtml">an overview</a> of the supported problem classes and additional recommendations for solving them.
</li>
<li>You can use GCG for solving MINLPs and more general constraint integer programs from your own source code.</li>
<li>You can use GCG as a framework in which you implement your own plugins.</li>
<li>You can use GCG in any combination of the three purposes above.</li>
</ul>
<p>
This FAQ contains separate sections covering each of these usages of SCIP. It further considers specific questions for some features.
</p>
    </div>
  </li>
    <li id="whenusegcg" class="anchor">
    <h4>
      <a class="reveal_faq" href="#whenusegcg">When should I use GCG?</a>

      <a href="#" class="pull-right"><span class="fa fa-caret-up" title="go to top"></span></a>
    </h4>
    <div id="whenusegcg_ans" class="answer">
      <p>
If you are looking for a speedup on many instances while using \SCIP, a fast non-commercial MIP/MINLP-solver.
GCG allows you to implement your own methods and take full control of the solving process.
</p>
    </div>
  </li>
    <li id="licensefaq" class="anchor">
    <h4>
      <a class="reveal_faq" href="#licensefaq">I heard something about licenses. Do I have to pay for using GCG?</a>

      <a href="#" class="pull-right"><span class="fa fa-caret-up" title="go to top"></span></a>
    </h4>
    <div id="licensefaq_ans" class="answer">
      <p>
Just like for SCIP, as long as you use it for academic, non-commercial purposes: No.
This will not change. For the other cases, check the explanation of the
<a href="http://scip.zib.de/#license">ZIB academic license</a> and always feel free to ask us.
</p>
    </div>
  </li>
    <li id="howtogetstarted" class="anchor">
    <h4>
      <a class="reveal_faq" href="#howtogetstarted">How do I get started?</a>

      <a href="#" class="pull-right"><span class="fa fa-caret-up" title="go to top"></span></a>
    </h4>
    <div id="howtogetstarted_ans" class="answer">
      <p>
An easy way is to use the SCIP-binaries and call SCIP from a shell, see <a href="./SHELL.shtml">here</a> for a tutorial.
For that, you just have to download one of the precompiled binaries from the
<a href="http://scip.zib.de/#download">download section</a>, or the zipped source code and compile it
with your favorite settings. This is described in detail in the <code><a href="./INSTALL.shtml">INSTALL</a></code> file in the SCIP main directory.
</p>
<p>
Another way is to use SCIP as a solver integrated into your own program source code.
See the directories &quot;examples/MIPsolver/&quot; and &quot;examples/Queens/&quot;
for simple examples and <a href="#howtocreateproblem">this point</a>.
</p>
<p>
A third way is to implement your own plugins into SCIP.
This is explained in the HowTos for all plugin types, which you can find in the
<a href="./index.shtml">doxygen documentation</a>.
See also <a href="./START.shtml">How to start a new project</a>.
</p>
    </div>
  </li>
    <li id="createsettingsfile" class="anchor">
    <h4>
      <a class="reveal_faq" href="#createsettingsfile">How do I create my own settings file?</a>

      <a href="#" class="pull-right"><span class="fa fa-caret-up" title="go to top"></span></a>
    </h4>
    <div id="createsettingsfile_ans" class="answer">
      <p>
Hm.
</p>    </div>
  </li>
  </ol>
<!--WE don't use this stuff because it is uncontrollable
<script>
$(".reveal").click(function() {
    var questionId = "#"+this.id+"_ans";
    $(questionId).toggle("fast");
});
</script>-->
