import sbt._

class ModiphyProject(info: ProjectInfo) extends DefaultProject(info)
{
  override def compileOptions = super.compileOptions ++ Seq(Unchecked)

  val logspace = "Logspace maven repo" at "http://www.logspace.co.uk/maven/"
  val tlf = "uk.co.logspace.tlf" % "tlf" % "1.0.1"
  val colt = "coltjar" % "colt" % "1.2.0" from "http://repo1.maven.org/maven2/colt/colt/1.2.0/colt-1.2.0.jar"
  val commonsMath = "org.apache.commons" % "commons-math" % "2.0"
  
}
