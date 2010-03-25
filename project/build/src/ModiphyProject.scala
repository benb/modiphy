import sbt._

class ModiphyProject(info: ProjectInfo) extends DefaultProject(info)
{
  override def compileOptions = super.compileOptions ++ Seq(Unchecked)

  val newReleaseToolsRepository = "Scala Tools Repository" at "http://nexus.scala-tools.org/content/repositories/snapshots/" 
  val scalatest = "org.scalatest" % "scalatest" % "1.0.1-for-scala-2.8.0.Beta1-with-test-interfaces-0.3-SNAPSHOT" % "test" 



  val logspace = "Logspace maven repo" at "http://www.logspace.co.uk/maven/"
  val tlf = "uk.co.logspace.tlf" % "tlf_2.8.0.Beta1"  % "1.0.1"
  val colt = "coltjar" % "colt" % "1.2.0" from "http://repo1.maven.org/maven2/colt/colt/1.2.0/colt-1.2.0.jar"
  val commonsMath = "org.apache.commons" % "commons-math" % "2.0"
  val commonsCollections = "commons-collections" % "commons-collections" % "3.2.1"
  val beastBeauti = "dr.math" % "beauti" % "1.5.3" from "http://www.logspace.co.uk/jar/beauti-1.5.3.jar"
 // val sbt = "sbt" % "sbt_2.7.7" % "0.5.7" from "http://simple-build-tool.googlecode.com/svn-history/r1125/artifacts/0.5.7-p1/jars/sbt_2.7.7.jar"
  
}
