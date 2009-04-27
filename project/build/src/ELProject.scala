import sbt._

class ELProject(info: ProjectInfo) extends DefaultProject(info)
{
    //override def mainClass = Some("Test")

    val scalatest = "org.scala-tools.testing" % "scalatest" % "0.9.5"
}

