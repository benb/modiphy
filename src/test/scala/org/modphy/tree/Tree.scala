package org.modphy.tree
import org.scalatest.Suite
import org.modphy.sequence.BioSeq._

class TreeSuite extends Suite{
  def testNode{
    import org.modphy.sequence.DNA
    val l1=new Leaf("l1","AGCT",DNA,0.1)
    val l2=new Leaf("l2","ACCT",DNA,0.2)
    val root = new INode(List(l1,l2),DNA,0.2D)
    assert((root child 0).get == l1)
    assert((root child 1).get == l2)
    assert((root length 0) == 0.1D )
    assert((root length 1) == 0.2D )
    assert((root length l2) == 0.2D )
  }
}

