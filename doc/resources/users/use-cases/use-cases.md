# Use Cases of GCG {#use-cases}
> <span style="color:gray;">Some guides are not complete yet, they are colored in gray.</span>

In the following, you find short **descriptions of all the use cases that
we offer guides for**, which are ordered according to their complexity.\n
If you're new here, we recommend starting with the first one.
\n

> <div style="width:80%">
> <img src="user.png" style="vertical-align:middle; height:3%; position:absolute; right:40px; margin-top: -1em;">
> <h3> Solving existing structured mixed-integer programs</h3>
> If you already have **existing instances** and have probably also tried them in other
> solvers already, we will show you how to get started with GCG for the first time 
> **without much of technical detail** and explain the output superficially.\n
> ⇨ @subpage u1
> \n\n
> </div>

> <div style="width:80%">
> <img src="user.png" style="vertical-align:middle; height:3%; position:absolute; right:40px; margin-top: -1em;">
> <h3> Solving freshly created instances and start to understand GCG</h3>
> If you first want to **create your instance**, for example using ZIMPL,
> this is the right guide for you. We will use GCG to solve a freshly created problem
> optimally and afterwards **explain the theory** behind the most important
> steps GCG does to do that.\n 
> ⇨ @subpage u2
> \n\n
> </div>

> <div style="width:80%">
> <img src="user.png" style="vertical-align:middle; height:3%; position:absolute; right:40px; margin-top: -1em;">
> <h3> Making GCG use more information from your model</h3>
> Using GCG as an external server for GAMS, GCG is able to read more information from your model than
> with a usual LP-File. We assume that you **already created your instance** using GAMS. 
> We will elaborate on how to set up GCG with GAMS to solve your problem.
> Contrary to the previous use case, we will **not go into details concerning the theory**
> behind what GCG does.\n 
> ⇨ @subpage u3
> \n\n
> </div>

> <div style="width:80%">
> <img src="expert.png" style="vertical-align:middle; height:3.4%; position:absolute; right:40px; margin-top: -1em;">
> <h3> Gain knowledge about the structure of your model</h3>
> To get to know about **how GCG sees your instance** and why, we check out the 
> @ref explore-menu "Explore Menu". With it, we can make GCG show the structure that its
> @ref detectors "detectors" found and also influence and tweak that a bit without 
> too much technical knowledge. \n
> ⇨ @subpage u4
> \n\n
> </div>

> <div style="width:80%">
> <img src="expert.png" style="vertical-align:middle; height:3.4%; position:absolute; right:40px; margin-top: -1em; opacity: 0.5;">
> <h3 style="color:gray;">Tell GCG how to decompose your problem using your own knowledge</h3>
> <span style="color:gray;">
> It can happen that you are **able to fully or at least partially grasp how your problem should be decomposed**,
> e.g. with the help of the explore menu (see above). In this guide, we will teach you how you can use that knowledge.\n
> ⇨ @subpage u5 
> </span>
> \n\n
> </div>

> <div style="width:80%">
> <img src="expert.png" style="vertical-align:middle; height:3.4%; position:absolute; right:40px; margin-top: -1em; opacity: 0.5;">
> <h3 style="color:gray;"> Tell GCG how to solve your problem using GCG settings</h3>
> <span style="color:gray;">
> Apart from the detection, we have more settings, in particular the pricing settings, which can be 
> **changed to achieve speedups** in solving your instance. In this guide, we will also give you 
> **some hints on what could be good parameters**.\n
> ⇨ @subpage u6
> </span>
> \n\n
> </div>

> <div style="width:80%">
> <img src="scientist.png" style="vertical-align:middle; height:3.2%; position:absolute; right:40px; margin-top: -1em; opacity: 0.5;">
> <h3 style="color:gray;"> Gain knowledge about how GCG solves your instance</h3>
> <span style="color:gray;">
> To see how GCG ran with your instance, for example how long it detected or how performant
> the pricing was, you should check out the @ref visu-suite "Visualization Suite". It shows not only
> timings, but also **illustrates algorithmic behaviour**.\n
> ⇨ @subpage u7
> </span>
> \n\n
> </div>

> <div style="width:80%">
> <img src="scientist.png" style="vertical-align:middle; height:3.2%; position:absolute; right:40px; margin-top: -1em; opacity: 0.5;">
> <h3 style="color:gray;"> Modify GCG's algorithmic behaviour</h3>
> <span style="color:gray;">
> If you have to modify more than just parameters, it might be necessary to **directly add or change parts of code** 
> within GCG. You will need the guides on @ref howtoadd "adding your own plugins" and it will also
> be helpful to read through the @ref example-projects. If you are doing your own project with, 
> we would be very thankful if you contacted us such that we can add it to our documentation!\n
> ⇨ @subpage u8
> </span>
> \n\n
> </div>

> <div style="width:80%">
> <img src="developer.png" style="vertical-align:middle; height:3.2%; position:absolute; right:40px; margin-top: -1em;">
> <h3> Develop completely new features for GCG</h3>
> If it is not enough to add a pricing solver here or a branching rule there, then you might need to
> **add completely new features**. Please contact us before and, if possible, while doing that.
> Note that this guide is also presenting some of the **general internal workflow** for which Git access might
> sometimes be required.\n
> ⇨ @ref new-gcg-dev
> \n\n
> </div>
