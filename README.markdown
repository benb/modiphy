Building:
=========

Put colt-1.2.0.jar in lib/ (can't get from MVN repo as they are mistagged!)

You need (Simple Built Tool)[http://code.google.com/p/simple-build-tool/]
You also need Apache Ivy in the classpath for SBT

run:

  sbt update # to pull in latest jars
  sbt compile

