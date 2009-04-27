package org.modphy.math
import  org.modphy.math.EnhancedMatrix._
import cern.colt.matrix.DoubleFactory2D._
import org.scalatest.Suite
import org.modphy.tree._
import org.modphy.sequence.DNA

class MathSuite extends Suite{
    def testMatrix{
      val size=4
      val x = dense.make(size,size)
      (0 to size-1) foreach{i=> x(i,i)=i+1}
      println(x)
      assert(x(3,3)==4)
    }

    def testLikelihood{
      import DNA._ 
      val sMat=dense.make(4,4)
      sMat.assign(0.1)
      val pi = Vector(4)
      pi(A.id)=0.3
      pi(C.id)=0.2
      pi(G.id)=0.3
      pi(T.id)=0.2
      val model = (sMat sToQ pi,pi)
      println(model)

      val fasta = List(">one","AGGT-A",">two","AGGTCC",">three","ATGT-A").elements
      val tree="((one:0.1,two:0.2):0.12625,three:0.01281);"     

      val (rawtre,aln)=DataParse(tree,fasta,DNA)
      val tre:INode[DNA.type,Node[DNA.type]] with LikelihoodNode[DNA.type] = rawtre.mkLkl(model)

      val leaf=tre.child(0).get.child(0).get
      assert(leaf.name=="one")
      assert(leaf.lengthTo<0.1001 && leaf.lengthTo > 0.0999)
      assert(leaf.isInstanceOf[LikelihoodNode[DNA.type]])

      println(leaf.likelihoods)
      assert(leaf.likelihoods.length == "AGGT-A".length)

      println(leaf)

      val lkl = tre.mkLkl(model)
      println("HELLO")
      println(lkl.likelihoods)
      true
      
    }
}


