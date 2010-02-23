Modiphy
=======

Library for phylogenetic analysis using Scala.

Why?
====

Because I want to be able to rapidly develop evolutionary models and test them.
Scala is fast (enough) to run, fast to develop, and has threadless parallelism with 
Actors.

Status
======

Currently this software is in the testing phase. Documentation is sparse.

Future
======

I want to include HMM-based models (i.e. tree-HMM or phylo-HMM type features).
I would like to develop this into a full domain specific language.
It hasn't had any optimisation yet.

License
=======

Modiphy is licensed under the MIT/X license.

Building:
=========

Dependencies:

You need [Simple Build Tool](http://code.google.com/p/simple-build-tool/)
Dependencies should be automatically handled by SBT.
You also need Apache Ivy in the classpath for SBT:

    sbt update # to pull in latest jars
    sbt compile


