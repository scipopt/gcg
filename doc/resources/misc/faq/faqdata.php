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
',
              'label'=>'whatisgcg'
             ),

          array(
              'question'=>'When should I use GCG?',
              'answer'=>'<p>
If you are looking for a speedup on many instances while using SCIP, a fast non-commercial MIP/MINLP-solver.
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
</p>
',
              'label'=>'licensefaq'
             ),

          array(
              'question'=>'How do I get started?',
              'answer'=>'<p>
For that question, we made a dedicated page: <a href="MANUAL.html">The Manual</a>.
</p>
',
              'label'=>'howtogetstarted'
             ),

          array(
              'question'=>'How do I create a file out of my pen-and-paper program?',
              'answer'=>'<p>
You will first have to create a program file that can be read by GCG.
We offer readers for multiple input file formats, along which are the most frequently
used ones (.lp, .mps). Additionally, a ZIMPL interface is already present and a GAMS
interface is in development. With all these possible file formats, you can choose
which one to go with. For new users, we recommend using <a href="https://zimpl.zib.de/download/zimpl.pdf">ZIMPL</a>
to convert your paper-and-pen-program to a computer file.
After you obtained the required file, get started with our <a href="getting-started.html">guide</a> and read it.
</p>
',
              'label'=>'createprogramfile'
             ),

          array(
              'question'=>'How do I create my own settings file?',
              'answer'=>'<p>
A settings file should be located in the settings folder inside the root directory.
In each line, a parameter can be set. This parameter is of the same form that it is in GCG,
for example "detection/detectors/connected/enabled". Often, these parameters are bools,
so turning on the connected detector could be done by "detection/detectors/connected/enabled = TRUE".
All possible parameters of GCG\'s current version can be found <a href="PARAMETERS.html">here</a>.
Note that the #-symbol will start comments.
</p>
',
              'label'=>'createsettingsfile'
             ),

          array(
              'question'=>'How do I create my own decomposition file?',
              'answer'=>'<p>
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
',
              'label'=>'createdecfile'
             ),

          array(
              'question'=>'How do I generate decomposition visualizations?',
              'answer'=>'<p>
Just like the picture on the GCG landing page, you can export images of how GCG decomposed
your program. This can be done in the <a href="explore-menu.html">Explore Menu</a>.
</p>
<p>
To find out which decomposition stems from what other, GCG can also generate a decomposition "family tree".
This is done with "write familytree".
</p>
',
              'label'=>'decvisu'
             ),

          array(
              'question'=>'How do I generate runtime visualizations?',
              'answer'=>'<p>
Apart from the detection visualizations, GCG also comes with some python scripts that
allow to make graphics showing the pricing process, time distribution, bounds development and more.
A guide on how to use those scripts can be found <a href="generatevisu.html">here</a>
</p>',
              'label'=>'runvisu'
             )),
             )
        );
 ?>