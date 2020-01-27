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
    <a class="reveal_faq" href="#createprogramfile">
      How do I create a file out of my pen-and-paper program?    </a>
  </li>
    <li>
    <a class="reveal_faq" href="#createsettingsfile">
      How do I create my own settings file?    </a>
  </li>
    <li>
    <a class="reveal_faq" href="#createdecfile">
      How do I create my own decomposition file?    </a>
  </li>
    <li>
    <a class="reveal_faq" href="#decvisu">
      How do I generate decomposition visualizations?    </a>
  </li>
    <li>
    <a class="reveal_faq" href="#runvisu">
      How do I generate runtime visualizations?    </a>
  </li>
    <li>
    <a class="reveal_faq" href="#gams">
      How do I use the GAMS interface?    </a>
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
GCG is a <b>branch-price-and-cut solver</b> that also enables easy inclusion of own methods.
With its different techniques such as Dantzig-Wolfe and Benders decomposition,
GCG is able to solve some <b>instances much quicker</b> than other solvers, so
it is always worth to let GCG solve your instance.
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
You should use GCG if you want to...
<ul>
<li>solve a problem that might have a <b>structure in it</b> (whether known or unknown) much faster (see <a href="gcg-vs-scip.html">GCG vs SCIP</a>)</li>
<li>find out about how your problem is structured (see <a href="explore-menu.html">Explore Menu</a>)</li>
<li>implement your own rules, heuristics and methods <b>without having to implement branch-and-price (see <a href="example-projects.html">Example Projects</a>)</b></li>
</ul>
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
For that question, we made a dedicated page: <a href="MANUAL.html">The Manual</a>.
</p>
    </div>
  </li>
    <li id="createprogramfile" class="anchor">
    <h4>
      <a class="reveal_faq" href="#createprogramfile">How do I create a file out of my pen-and-paper program?</a>

      <a href="#" class="pull-right"><span class="fa fa-caret-up" title="go to top"></span></a>
    </h4>
    <div id="createprogramfile_ans" class="answer">
      <p>
You will first have to create a program file that can be read by GCG.
We offer readers for multiple input file formats, along which are the most frequently
used ones (.lp, .mps). Additionally, a ZIMPL interface is already present and a GAMS
interface is in development. With all these possible file formats, you can choose
which one to go with. For new users, we recommend using <a href="https://zimpl.zib.de/download/zimpl.pdf">ZIMPL</a>
to convert your paper-and-pen-program to a computer file.
After you obtained the required file, get started with our <a href="getting-started.html">guide</a> and read it.
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
A settings file should be located in the settings folder inside the root directory.
The most common way to generate a settings file is through the <code>set diffsave</code>
command inside GCG that exports a <code>.set</code>-file with all parameters that
you changed.\n
Otherwise, you can also define a settings file by yourself. In each line, a parameter
can be set. This parameter is of the same form that it is in GCG,
for example "detection/detectors/connected/enabled". Often, these parameters are bools,
so turning on the connected detector could be done by "detection/detectors/connected/enabled = TRUE".
All possible parameters of GCG's current version can be found <a href="PARAMETERS.html">here</a>.
Note that the #-symbol will start comments.
</p>
    </div>
  </li>
    <li id="createdecfile" class="anchor">
    <h4>
      <a class="reveal_faq" href="#createdecfile">How do I create my own decomposition file?</a>

      <a href="#" class="pull-right"><span class="fa fa-caret-up" title="go to top"></span></a>
    </h4>
    <div id="createdecfile_ans" class="answer">
      <p>
First off; we do not recommend to write any decomposition files. Instead, you
should make GCG detect your desired decomposition by using the settings. For example,
you could look up which detector decomposes your instance in the way you want it to be
and then only activate this detector. To find out which detector you need, look up
which structure each detector finds <a href="detectors.html">here</a>. Then, create your
<a href"#createsettingsfile">own settings file</a>, where only this detector is enabled.
</p>
<p>
If you still want to write your .dec-file by hand, you can find the required syntax <a href="reader__dec_8h.html">here</a>.
</p>
    </div>
  </li>
    <li id="decvisu" class="anchor">
    <h4>
      <a class="reveal_faq" href="#decvisu">How do I generate decomposition visualizations?</a>

      <a href="#" class="pull-right"><span class="fa fa-caret-up" title="go to top"></span></a>
    </h4>
    <div id="decvisu_ans" class="answer">
      <p>
Just like the picture on the GCG landing page, you can export images of how GCG decomposed
your program. This can be done in the <a href="explore-menu.html">Explore Menu</a>.
</p>
<p>
To find out which decomposition stems from what other, GCG can also generate a decomposition "family tree".
This is done with "write familytree".
</p>
    </div>
  </li>
    <li id="runvisu" class="anchor">
    <h4>
      <a class="reveal_faq" href="#runvisu">How do I generate runtime visualizations?</a>

      <a href="#" class="pull-right"><span class="fa fa-caret-up" title="go to top"></span></a>
    </h4>
    <div id="runvisu_ans" class="answer">
      <p>
Apart from the detection visualizations, GCG also comes with some python scripts that
allow to make graphics showing the pricing process, time distribution, bounds development and more.
A guide on how to use those scripts can be found <a href="generatevisu.html">here</a>
</p>
    </div>
  </li>
    <li id="gams" class="anchor">
    <h4>
      <a class="reveal_faq" href="#gams">How do I use the GAMS interface?</a>

      <a href="#" class="pull-right"><span class="fa fa-caret-up" title="go to top"></span></a>
    </h4>
    <div id="gams_ans" class="answer">
      <p>
For the GAMS interface to run with GCG, you have to make sure that the dependencies were compiled with SHARED=true and READLINE=false.
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
