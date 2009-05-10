Modiphy
=======

Library for phylogenetic analysis using Scala.

Why?
====

Because I want to be able to rapidly develop evolutionary models and test them.
Scala is fast (enough) to run, fast to develop, and flexible enough that models 
can be declared by mixing in parameters.

For instance, to declare models like:

//standard model
val model = new Model[DNA] with GTR
//gamma model
val model2 = new Model[DNA] with GTR with GammaRates
//model where nodes 2,3 have different gamma rate
val model3 = new Model[DNA] with GTR with GammaRates with AltModel{val altModel=new GammaRates;val altNodes=List(2,3)}

It won't quite be that simple - but you get the idea.

Status
======

Currently this software is in the testing phase. Documentation is sparse.

Future
======

I want to include HMM-based models (i.e. tree-HMM or phylo-HMM type features).
I would like to develop this into a full domain specific language.
It hasn't had any optimisation yet.

Problems
========

I periodically get a load of build errors, where scala can't tell that my model classes
inherit from the base trait Model. Removing target/ and rebuilding seems to fix things.

Building:
=========

Dependencies:

Put colt-1.2.0.jar in lib/ (can't get from MVN repo as they are mistagged!)
You also need need a [commons-math-2.0-snapshot](http://commons.apache.org/math/).jar
(I use the optimisation routines), as well as
[pal-1.5.1](http://www.cebl.auckland.ac.nz/pal-project/download.html).jar (I use
the Gamma Distribution routines) in lib/.

You need [Simple Built Tool](http://code.google.com/p/simple-build-tool/)
You also need Apache Ivy in the classpath for SBT:

    sbt update # to pull in latest jars
    sbt compile


