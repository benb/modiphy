package org.modphy.tree
import org.scalatest.Suite
import org.modphy.sequence.BioSeq._
import org.modphy.sequence._
import org.modphy.tree._

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
  def testTreeParse{
    val treeTxt="""(((((((((((((((((hg18:0.003731,panTro2:0.005501):0.013010,gorGor1:0.01):0.01,
      ponAbe2:0.020000):0.020000,rheMac2:0.031571):0.01,calJac1:0.01):0.020000,
    tarSyr1:0.100000):0.050000,(micMur1:0.084110,
    otoGar1:0.145437):0.033956):0.020000,tupBel1:0.203975):0.020000,
  (((((mm9:0.104920,rn4:0.109421):0.020000,dipOrd1:0.200000):0.050000,
    cavPor3:0.150000):0.050000,speTri1:0.150000):0.020000,(oryCun1:0.100000,
  ochPri2:0.200000):0.050000):0.020000):0.020000,(((vicPac1:0.100000,
  (turTru1:0.120000,bosTau4:0.162368):0.050000):0.050000,((equCab2:0.150000,
    (felCat3:0.098674,canFam2:0.114682):0.002783):0.050000,(myoLuc1:0.142600,
    pteVam1:0.142460):0.003381):0.007170):0.030000,(eriEur1:0.279121,
  sorAra1:0.309867):0.023929):0.040000):0.01,(((loxAfr2:0.110021,
  proCap1:0.150000):0.003000,echTel1:0.265218):0.066339,(dasNov2:0.178799,
choHof1:0.200000):0.090000):0.053170):0.213469,monDom4:0.320721):0.088647,
ornAna1:0.488110):0.118797,((galGal3:0.23,taeGut1:0.16):0.16,
anoCar1:0.513962):0.093688):0.151358,xenTro2:0.778272):0.174596,
(((tetNig1:0.203933,fr2:0.239587):0.203949,(gasAcu1:0.314162,
  oryLat2:0.501915):0.055354):0.346008,danRer5:0.730028):0.174596):0.1,
petMar1:0.1);"""

  
  val tree = DataParse(treeTxt,DNA)

  assert(tree.descendents contains "cavPor3")
  assert(tree.descendents contains "oryLat2")
  println(tree.descendents)
  
  true

  val treeTxt2="((one:0.1,two:0.1):0.2,three:0.1);"
  val tree2 = DataParse(treeTxt2,DNA)
  val tree3 = tree2.setBranchLengths({1 to 4}.map{i=>0.0}.toList)


  println(tree3)
  assert (tree3.child(0).get.lengthTo==0.0D)
  



  }

}

