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
      assert(x(3,3)==4)
    }

    def testLikelihood{
      import DNA._ 
      val sMat=dense.make(4,4)
      sMat.assign(0.1)
      val pi = Vector(4)
      pi(A)=0.25
      pi(C)=0.25
      pi(G)=0.25
      pi(T)=0.25

      val fasta = List(">one","AGGT-A-",">two","AGGTCC-",">three","ATGT-A-").elements
      //val tree="((one:0.1,two:0.2):0.12625,three:0.01281);"     
      val tree="(one:0.1,two:0.2);"     

      val (rawtre,aln)=DataParse(tree,fasta,DNA)
      val model = new EnhancedModel(pi,sMat,rawtre)
      assert(aln("one")=="AGGT-A-")
      val tre = rawtre.mkLkl(model)

      val leaf=tre.child(0).get.asInstanceOf[Leaf[DNA.type] with LikelihoodNode[DNA.type]]
      assert(leaf.name=="one")
      assert(leaf.lengthTo<0.1001 && leaf.lengthTo > 0.0999)
      assert(leaf.isInstanceOf[LikelihoodNode[DNA.type]])

      println("LEAF LKL " + leaf.likelihoods)
      assert(leaf.likelihoods.length == "AGGT-A-".length)

      val lkl = model.getTree
      println("HELLO")
      println("IS ROOT NODE? " + lkl.isRoot)
      assert(lkl.isInstanceOf[RootNode[DNA.type]])
      assert(lkl.realLikelihoods.last < 1.00001 && lkl.realLikelihoods.last > 0.9999)
      println("LOG LKL " + lkl.logLikelihood)

    
      ModelOptimiser.nelderMead(model) 
      println("Opt mat " + model.sMat)
      println("Opt lkl " + model.logLikelihood)

      //in alignment above, there are A->T and A->C but no A->G transitions
      assert(model.sMat(A,C)>model.sMat(A,G))
      assert(model.sMat(A,T)>model.sMat(A,G))


      true
      
    }
}


