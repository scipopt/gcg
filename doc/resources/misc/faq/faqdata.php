<?php 

   $faq = array(
        
          array(
             'title'=>'General Questions about GCG',
             'label'=>'faq_generalquestionsaboutgcg',
             'content'=>array(
          array(
              'question'=>'What is GCG?',
              'answer'=>'<p>
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
',
              'label'=>'whatisgcg'
             ),

          array(
              'question'=>'When should I use GCG?',
              'answer'=>'<p>
If you are looking for a speedup on many instances while using \SCIP, a fast non-commercial MIP/MINLP-solver.
GCG allows you to implement your own methods and take full control of the solving process.
</p>
',
              'label'=>'whenusegcg'
             ),

          array(
              'question'=>'I heard something about licenses. Do I have to pay for using GCG?',
              'answer'=>'<p>
Just like for SCIP, as long as you use it for academic, non-commercial purposes: No.
This will not change. For the other cases, check the explanation of the
<a href="http://scip.zib.de/#license">ZIB academic license</a> and always feel free to ask us.
If you want to use SCIP commercially, please write an e-mail to koch@zib.de.
</p>
',
              'label'=>'licensefaq'
             ),

          array(
              'question'=>'How do I get started?',
              'answer'=>'<p>
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
</p>',
              'label'=>'howtogetstarted'
             )),
             )
        );
 ?>