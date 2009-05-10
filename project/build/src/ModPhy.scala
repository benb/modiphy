import sbt._

class ModPhyProject(info: ProjectInfo) extends DefaultProject(info)
{
    //override def mainClass = Some("Test")

    val scalatest = "org.scala-tools.testing" % "scalatest" % "0.9.5"
    //val commonsMath = "commons-math" % "commons-math" % "1.2"
    //use 2.0-SNAPSHOT from svn
    override def ivyXML =
    <dependencies>
      <dependency org="sanger" name="tlf" rev="1.0" conf="default">
        <artifact name="tlf" url="http://www.sanger.ac.uk/~bb4/jar/tlf-1.0.jar" />
      </dependency>
    </dependencies>
}

