package org.modphy.tree
import org.modphy.sequence._
import org.modphy.math.EnhancedMatrix._
import cern.colt.matrix.DoubleFactory1D
import scala.collection.Set
import org.modphy.math._
import tlf.Logging


import scala.util.parsing.combinator._
/**
  Used internally for parsing trees
*/
class TreeParser[A <: BioEnum](m:String=>String,alphabet:A) extends JavaTokenParsers{
  val treegen = new TreeGen[A]
  import treegen._
  def root: Parser[INode[A]] = "("~node~rep(","~>node)~");" ^^ {
    case "("~node1~nodeList~");" => getINode(node1::nodeList,alphabet,0.0)
  }
  def node: Parser[Node[A]] = leaf | "("~node~rep(","~>node)~"):"~floatingPointNumber ^^ {
    case "("~node1~nodeList~"):"~length => getINode(node1::nodeList,alphabet,length.toDouble)
  } 
  def seqName: Parser[String] = regex(new scala.util.matching.Regex("[a-zA-Z0-9_.+-]+"))
  def leaf: Parser[Leaf[A]] = seqName~":"~floatingPointNumber ^^ {case name~":"~length => getLeaf(name,m(name),alphabet,length.toDouble)}
}

class TreeGen[A <: BioEnum]{
  var nodeCount:Int= -1
  def getINode(children:List[Node[A]],alphabet:A, lengthTo:Double)={nodeCount += 1;new INode[A](children,alphabet,lengthTo,nodeCount)}
  def getLeaf(name:String,seq:String, alphabet:A, lengthTo:Double)={nodeCount += 1;new Leaf[A](name,seq,alphabet,lengthTo,nodeCount)}
}

/**
 Parses Newick format trees
*/
object DataParse{
  type Tree[A <: BioEnum]=Node[A] with RootNode[A]
  type Alignment=Map[String,String]

  private def cleanTree(t:String)=t.split("\\s+").mkString("").split("\n").mkString("")

  /**
   Produce a parsed Tree and alignment from a newick file and a fasta file
  */
  def apply[A <: BioEnum](tree:String,alignment:Iterator[String],alphabet:A):(Tree[A],Alignment)={
    val aln = new Fasta(alignment)
    val seqMap = aln.foldLeft(Map[String,String]()){_+_}
    apply(tree,seqMap,alphabet)
  }

  
  def apply[A <: BioEnum](tree:String,aln:Alignment,alphabet:A):(Tree[A],Alignment)={
    val root = new TreeParser[A](aln,alphabet){def parseAll=parse(root,cleanTree(tree))}.parseAll.get.setRoot
    (root,aln)
  }

  def dropNodes[A <: BioEnum](tree:String,aln:Alignment,alphabet:A):(Tree[A],Alignment)={
    val t1:Tree[A] = apply(tree,alphabet)
    val t2 = t1.restrictTo(aln.keySet).setRoot
    //println(t2)
    (t2.setAlign(aln).setRoot,aln)
  }

  def apply[A <: BioEnum](tree:String,alphabet:A):Tree[A]={
    //println("TREE= " + cleanTree(tree))
    new TreeParser[A]({s:String=>""},alphabet){def parseAll=parse(root,cleanTree(tree))}.parseAll.get.setRoot
  }
}



import DataParse._


trait Node[A <: BioEnum] extends Logging{
  val id:Int
  val isRoot=false
  def mkLkl(mod:Model[A]):LikelihoodNode[A]
  def lengthTo:Double
  def child(i:Int):Option[Node[A]]
  def name:String
  def descendents:List[String]
  def numChildren:Int
  def removeUseless:Node[A]
  def resize(bl:Double):Node[A]
  def setAlign(m:Map[String,String]):Node[A]
  def restrictTo(allowed:Set[String]):Node[A]
  def descendentNodes:List[Node[A]]
  def setBranchLengths(l:List[Double]):Node[A]
  def getBranchLengths:List[Double]
  def branchTo:String
  def setNewDataType[B <: BioEnum](alphabet:B):Node[B]

  def nodes=this::descendentNodes
}

trait LikelihoodNode[A <: BioEnum] extends Node[A]{
  val alphabet:A
  def likelihoods:List[Vector]
}

class CalcLikelihoodNode[A <: BioEnum](children:List[LikelihoodNode[A]],alphabet:A,lengthTo:Double,val model:Model[A],id:Int) extends INode[A](children,alphabet,lengthTo,id) with LikelihoodNode[A]{
  //return a list of the lists of probabilties for each site



  override def childElements:Iterator[LikelihoodNode[A]] = children.elements
  
  def realLikelihoods = {

    likelihoods.map{vec=>
      val ans = model.piVals.toArray.elements.zipWithIndex.map{t=>
        val(p,i)=t
        vec(i)*p
      }.foldLeft(0.0D){_+_}
      ans
    }

  }
  def logLikelihood={
    val lkl = realLikelihoods
    lkl.foldLeft(0.0D){(i,j)=>i+Math.log(j)}
  }
  override lazy val likelihoods:List[Vector]={
    model.likelihoods(this)
  }
}

class LeafLikelihoodNode[A <: BioEnum](name:String,seq:String,alpha:A,lengthTo:Double,id:Int) extends Leaf[A](name,seq,alpha,lengthTo,id) with LikelihoodNode[A]{
  val sequence=alphabet.parseString(seq.toUpperCase)
  override lazy val likelihoods:List[Vector]={
    val elements = alphabet.matLength
    sequence.map{a=>
      val vec = DoubleFactory1D.dense.make(alphabet.matLength)
        alphabet.getNums(a).foreach{i=>
        vec(i)=1.0D}
        vec
    }.toList
  }
}

trait RootNode[A <: BioEnum] extends INode[A]{
  val children:List[Node[A]]
  override val isRoot=true
  override def toString="("+children.mkString(",")+");"

  override def setNewDataType[B <: BioEnum](alpha:B)=new INode[B](children.map{i=> i.setNewDataType(alpha)},alpha,lengthTo,id) with RootNode[B]
  override def factory(c:List[Node[A]],al:A,len:Double):INode[A] with RootNode[A]=new INode[A](c,al,len,id) with RootNode[A]

  override def getBranchLengths={
    children.map{c=>c.getBranchLengths}.flatten[Double]
  }
  override def setBranchLengths(l:List[Double]):INode[A] with RootNode[A]={
    var listPtr = l
    val newchildren = children.map{c=>val c2 = c.setBranchLengths(listPtr); listPtr = listPtr.drop(c.descendentNodes.size + 1);c2}
    factory(newchildren,alphabet,l.head)
  }
  override def removeUseless:INode[A] with RootNode[A]=super.removeUseless.setRoot
  override def restrictTo(allowed:Set[String]):INode[A] with RootNode[A]=super.restrictTo(allowed).setRoot
}

class INode[A <: BioEnum](val children:List[Node[A]],val alphabet:A,val lengthTo:Double,val id:Int) extends Node[A]{ 
  val name=""

  def setNewDataType[B <: BioEnum](alpha:B)=new INode[B](children.map{i=> i.setNewDataType(alpha)},alpha,lengthTo,id)

  def factory(c:List[Node[A]],al:A,len:Double)=if (isRoot){new INode[A](c,al,len,id) with RootNode[A]}else {new INode[A](c,al,len,id)}

  def numChildren = children.size
  def childElements:Iterator[Node[A]] = children.elements
  def child(i:Int)=if (i < children.length){Some(children(i))}else{None}
  def length(i:Int):Double=children(i).lengthTo
  def length(n:Node[A]):Double=children.find{node=>node==n}.get.lengthTo
  def mkLkl(mod:Model[A]):CalcLikelihoodNode[A]=if (isRoot){
    new CalcLikelihoodNode[A](children.map{t=>t.mkLkl(mod)}.toList,alphabet,lengthTo,mod,id) with RootNode[A]
  }else { 
    new CalcLikelihoodNode[A](children.map{t=>t.mkLkl(mod)}.toList,alphabet,lengthTo,mod,id) 
  }
  def restrictTo(allowed:Set[String])={
    val newChildren=children.filter{child=>child.descendents.exists{name=>allowed contains name}}map{_.restrictTo(allowed)}
    //println(children.toString + " => " + newChildren.toString)
    factory(newChildren,alphabet,lengthTo).removeUseless
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
    if (newChildren.length==1 && newChildren(0).isInstanceOf[INode[A]]){
        newChildren(0).asInstanceOf[INode[A]]
    }else{
      factory(newChildren,alphabet,lengthTo)
    }
  }

  def setAlign(m:Map[String,String])=factory(children.map{i=>i.setAlign(m)},alphabet,lengthTo)
  def resize(bl:Double)=factory(children,alphabet,bl)
          
  override def toString={
    "("+children.mkString(",")+"):"+lengthTo
  }

  def setRoot=new INode[A](children,alphabet,0.0D,id) with RootNode[A]

  def setBranchLengths(l:List[Double])={
    var listPtr=l.tail
    val newchildren = children.map{c=>val c2 = c.setBranchLengths(listPtr); listPtr = listPtr.drop(c.descendentNodes.size + 1);c2}
    factory(newchildren,alphabet,l.head)
  }

  def getBranchLengths={
    lengthTo :: children.map{c=>c.getBranchLengths}.flatten[Double]
  }

  def branchTo=descendents.mkString(",")
}

class Leaf[A <: BioEnum](val name:String,val seq:String,val alphabet:A,val lengthTo:Double, val id:Int) extends Node[A]{

  def setNewDataType[B <: BioEnum](alpha:B)=new Leaf[B](name,seq,alpha,lengthTo,id)

  def descendentNodes=List()
  def setBranchLengths(l:List[Double])={
  //   println("Node " + this + " setting branch length " + l.head)
     new Leaf[A](name,seq,alphabet,l.head,id)
 }
  def mkLkl(mod:Model[A])=new LeafLikelihoodNode[A](name,seq,alphabet,lengthTo,id)
  def child(i:Int)=None
  def descendents=List(name)
  def resize(bl:Double) = new Leaf[A](name,seq,alphabet,bl,id)
  def numChildren=0
  def removeUseless=this
  def setAlign(m:Map[String,String])=new Leaf[A](name,m(name),alphabet,lengthTo,id)
  def restrictTo(allowed:Set[String])=this
  override def toString=name + ":" + lengthTo
  def getBranchLengths=List(lengthTo)
  def factory(n:String,s:String,a:A,len:Double)=new Leaf[A](n,s,a,len,id)
  def branchTo=name
}


