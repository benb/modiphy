package org.modiphy.tree
import org.modiphy.sequence._
import org.modiphy.math.EnhancedMatrix._
import cern.colt.matrix.DoubleFactory1D
import scala.collection.Set
import org.modiphy.math._
import tlf.Logging


import scala.util.parsing.combinator._
/**
  Used internally for parsing trees
*/
class TreeParser[A <: BioEnum](aln:Alignment[A]) extends JavaTokenParsers{
  val treegen = new TreeGen[A]
  import treegen._
  def root: Parser[INode[A]] = "("~node~rep(","~>node)~");" ^^ {
    case "("~node1~nodeList~");" => getINode(node1::nodeList,aln,0.0)
  }
  def node: Parser[Node[A]] = leaf | "("~node~rep(","~>node)~"):"~floatingPointNumber ^^ {
    case "("~node1~nodeList~"):"~length => getINode(node1::nodeList,aln,length.toDouble)
  } 
  def seqName: Parser[String] = regex(new scala.util.matching.Regex("[a-zA-Z0-9_.+-]+"))
  def leaf: Parser[Leaf[A]] = seqName~":"~floatingPointNumber ^^ {case name~":"~length => getLeaf(name,aln,length.toDouble)}
}

class TreeGen[A <: BioEnum]{
  var nodeCount:Int= -1
  def getINode(children:List[Node[A]],aln:Alignment[A], lengthTo:Double)={nodeCount += 1;new INode[A](children,aln,lengthTo,nodeCount)}
  def getLeaf(name:String,aln:Alignment[A], lengthTo:Double)={nodeCount += 1;new Leaf[A](name,aln,lengthTo,nodeCount)}
}

/**
 Parses Newick format trees
*/
object DataParse{
  type Tree[A <: BioEnum]=Node[A] with RootNode[A]
  import sequence.Alignment

  private def cleanTree(t:String)=t.split("\\s+").mkString("").split("\n").mkString("")

  /**
   Produce a parsed Tree and alignment from a newick file and a fasta file
  */
  def apply[A <: BioEnum](tree:String,alignment:Iterator[String],alphabet:A):(Tree[A],Alignment[A])={
    val aln = new Fasta(alignment)
    val seqMap = aln.foldLeft(Map[String,String]()){_+_}
    apply(tree,new Alignment(seqMap,alphabet))
  }

  
  def apply[A <: BioEnum](tree:String,aln:Alignment[A]):(Tree[A],Alignment[A])={
    val root = new TreeParser[A](aln){def parseAll=parse(root,cleanTree(tree))}.parseAll.get.iNode.get.setRoot
    (root,aln)
  }

  def dropNodes[A <: BioEnum](tree:String,aln:Alignment[A]):(Tree[A],Alignment[A])={
    val t1:Tree[A] = apply(tree,aln)._1
    val t2 = t1.restrictTo(aln.map.keySet).iNode.get.setRoot
    //println(t2)
    (t2.setAlign(aln).iNode.get.setRoot,aln)
  }
}



import DataParse._


trait Node[A <: BioEnum] extends Logging{
  val id:Int
  val isRoot=false
  val aln:Alignment[A]
  val alphabet:A=aln.alphabet
  
  def lengthTo:Double
  def child(i:Int):Option[Node[A]]
  def children:List[Node[A]]
  def name:String
  def descendents:List[String]
  def numChildren:Int
  def removeUseless:Node[A]
  def resize(bl:Double):Node[A]
  def setAlign[B <: BioEnum](a2:Alignment[B]):Node[B]
  def restrictTo(allowed:Set[String]):Node[A]
  def descendentNodes:List[Node[A]]
  def setBranchLengths(l:List[Double]):Node[A]
  def getBranchLengths:List[Double]
  def branchTo:String

  def nodes=this::descendentNodes
  def likelihoods(m:Model[A]):List[Vector]=m.likelihoods(this)

  def iNode:Option[INode[A]]=None
}

trait RootNode[A <: BioEnum] extends INode[A]{
  val children:List[Node[A]]
  override val isRoot=true
  override def toString="("+children.mkString(",")+");"

  override def factory[B <: BioEnum](c:List[Node[B]],aln:Alignment[B],len:Double):INode[B] with RootNode[B]=new INode[B](c,aln,len,id) with RootNode[B]

  override def getBranchLengths={
    children.map{c=>c.getBranchLengths}.flatten[Double]
  }
  override def setBranchLengths(l:List[Double]):INode[A] with RootNode[A]={
    var listPtr = l
    val newchildren = children.map{c=>val c2 = c.setBranchLengths(listPtr); listPtr = listPtr.drop(c.descendentNodes.size + 1);c2}
    factory(newchildren,aln,l.head)
  }
  override def removeUseless:INode[A] with RootNode[A]=super.removeUseless.iNode.get.setRoot
  override def restrictTo(allowed:Set[String]):INode[A] with RootNode[A]=super.restrictTo(allowed).iNode.get.setRoot
}

class INode[A <: BioEnum](val children:List[Node[A]],val aln:Alignment[A],val lengthTo:Double,val id:Int) extends Node[A]{ 
  override def iNode=Some(this)
  def realLikelihoods(m:Model[A]) = {
    likelihoods(m).map{vec=>
      val ans = m.piVals.toArray.elements.zipWithIndex.map{t=>
        val(p,i)=t
        vec(i)*p
      }.foldLeft(0.0D){_+_}
      ans
    }
  }
  def logLikelihood(m:Model[A])={
    val lkl = realLikelihoods(m)
    lkl.zip(aln.pCount).foldLeft(0.0D){(i,j)=>i+Math.log(j._1)*j._2}
  }
  val name=""

  def factory[B<:BioEnum](c:List[Node[B]],aln:Alignment[B],len:Double)=if (isRoot){new INode[B](c,aln,len,id) with RootNode[B]}else {new INode[B](c,aln,len,id)}

  def numChildren = children.size
  def childElements:Iterator[Node[A]] = children.elements
  def child(i:Int)=if (i < children.length){Some(children(i))}else{None}
  def length(i:Int):Double=children(i).lengthTo
  def length(n:Node[A]):Double=children.find{node=>node==n}.get.lengthTo
  def restrictTo(allowed:Set[String])={
    val newChildren=children.filter{child=>child.descendents.exists{name=>allowed contains name}}map{_.restrictTo(allowed)}
    //println(children.toString + " => " + newChildren.toString)
    factory(newChildren,aln,lengthTo).removeUseless
  }
  def descendents:List[String]=childElements.map{i=>i.descendents}.toList.flatten[String]
  def descendentNodes = {children ++ children.map{c=>c.descendentNodes}.flatten[Node[A]]}
  def removeUseless:INode[A]={
    val newChildren = children.map{child=>
      if (child.numChildren >1){
        child
      }else if (child.numChildren==1){
        var endchild = child
        var bl = child.lengthTo
        while (endchild.numChildren==1){
          endchild=endchild.child(0).get
          bl= bl + endchild.lengthTo
        }
        endchild.resize(bl)
      
      }else {//if (child.numChildren==0){
        child
        //maybe need to remove if useless internal
      }
    }.map{_.removeUseless}
    if (newChildren.length==1){// This check can't be made because of erasure: && newChildren(0).isInstanceOf[INode[A]]){
        newChildren(0).asInstanceOf[INode[A]]
    }else{
      factory(newChildren,aln,lengthTo)
    }
  }

  def setAlign[B<:BioEnum](a2:Alignment[B]):Node[B]=factory(children.map{i=>i.setAlign(a2)},a2,lengthTo)
  def resize(bl:Double)=factory(children,aln,bl)
          
  override def toString={
    "("+children.mkString(",")+"):"+lengthTo
  }

  def setRoot=new INode[A](children,aln,0.0D,id) with RootNode[A]

  def setBranchLengths(l:List[Double])={
    var listPtr=l.tail
    val newchildren = children.map{c=>val c2 = c.setBranchLengths(listPtr); listPtr = listPtr.drop(c.descendentNodes.size + 1);c2}
    factory(newchildren,aln,l.head)
  }

  def getBranchLengths={
    lengthTo :: children.map{c=>c.getBranchLengths}.flatten[Double]
  }

  def branchTo=descendents.mkString(",")
}

class Leaf[A <: BioEnum](val name:String,val aln:Alignment[A],val lengthTo:Double, val id:Int) extends Node[A]{

  def children:List[Node[A]]=Nil

  val sequence:List[alphabet.Value]=aln.getPatterns(name).asInstanceOf[List[alphabet.Value]]
  lazy val likelihoods:List[Vector]={
    sequence.map{a:alphabet.Value=>
    
      val vec = Vector(alphabet.matLength)
        alphabet.getNums(a).foreach{i=>
        vec(i)=1.0D}
        vec
    }.toList
  }
  override def likelihoods(m:Model[A]):List[Vector]=likelihoods

  def descendentNodes=List()
  def setBranchLengths(l:List[Double])={
  //   println("Node " + this + " setting branch length " + l.head)
     new Leaf[A](name,aln,l.head,id)
 }
  def child(i:Int)=None
  def descendents=List(name)
  def resize(bl:Double) = new Leaf[A](name,aln,bl,id)
  def numChildren=0
  def removeUseless=this
  def setAlign[B<:BioEnum](a2:Alignment[B]):Node[B]=new Leaf[B](name,a2,lengthTo,id)
  def restrictTo(allowed:Set[String])=this
  override def toString=name + ":" + lengthTo
  def getBranchLengths=List(lengthTo)
  def factory[B <: BioEnum](n:String,a:Alignment[B],len:Double)=new Leaf[B](n,a,len,id)
  def branchTo=name
}


