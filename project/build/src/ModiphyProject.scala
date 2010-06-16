import sbt._
import java.io.File


class ModiphyProject(info: ProjectInfo) extends ProguardProject(info)
{
  override def compileOptions = super.compileOptions ++ Seq(Unchecked)

  val scalaToolsSnapshots = "Scala-Tools Maven2 Snapshots Repository" at "http://scala-tools.org/repo-snapshots"

  val scalatest = "org.scalatest" % "scalatest" % "1.2-for-scala-2.8.0.RC5-SNAPSHOT"
  val logspace = "Logspace maven repo" at "http://www.logspace.co.uk/maven/"
  val tlf = "uk.co.logspace.tlf" %% "tlf" % "1.2.0"
  val colt = "coltjar" % "colt" % "1.2.0" from "http://repo1.maven.org/maven2/colt/colt/1.2.0/colt-1.2.0.jar"
  val commonsMath = "org.apache.commons" % "commons-math" % "2.0"
  val commonsCollections = "commons-collections" % "commons-collections" % "3.2.1"
  val beastBeauti = "dr.math" % "beauti" % "1.5.3" from "http://www.logspace.co.uk/jar/beauti-1.5.3.jar"
 // val sbt = "sbt" % "sbt_2.7.7" % "0.5.7" from "http://simple-build-tool.googlecode.com/svn-history/r1125/artifacts/0.5.7-p1/jars/sbt_2.7.7.jar"
  
}
