package org.modiphy.tree
import org.modiphy.sequence._
import org.modiphy.math.EnhancedMatrix._
import cern.colt.matrix.DoubleFactory1D
import scala.collection.Set
import org.modiphy.math._
import tlf.Logging
import scala.actors._


import scala.util.parsing.combinator._
/**
  Used internally for parsing trees
*/

class TreeParser[A <: BioEnum](aln:Alignment[A]) extends JavaTokenParsers{
  val treegen = new TreeGen[A]
  import treegen._
  def root: Parser[INode[A]] = "("~node~rep(","~>node)~");" ^^ {
    case "("~node1~nodeList~");" => getINode(node1::nodeList,aln,0.0,None)
  }
  def node: Parser[Node[A]] = leaf | "("~node~rep(","~>node)~"):"~branchLength ^^ {
    case "("~node1~nodeList~"):"~length => getINode(node1::nodeList,aln,length._1,length._2)
  }
  def seqName: Parser[String] = regex(new scala.util.matching.Regex("[a-zA-Z0-9_.+-]+"))
  def leaf: Parser[Leaf[A]] = seqName~":"~branchLength ^^ {case name~":"~branchLength => getLeaf(name,aln,branchLength._1,branchLength._2)}
  def branchLength:Parser[(Double,Option[Int])] = {floatingPointNumber~"#"~wholeNumber ^^{
     case length~"#"~id => (length.toDouble,Some(id.toInt))
  } | floatingPointNumber  ^^ {
     case length => (length.toDouble,None)
  }
  }
}

class TreeGen[A <: BioEnum]{
  var nodeCount:Int= -1
  def getINode(children:List[Node[A]],aln:Alignment[A], lengthTo:Double, label:Option[Int]):INode[A]={nodeCount += 1;new INode[A](children,aln,lengthTo,nodeCount,label)}
  def getLeaf(name:String,aln:Alignment[A], lengthTo:Double,label:Option[Int]):Leaf[A]={nodeCount += 1;new Leaf[A](name,aln,lengthTo,nodeCount,label)}
}

/**
 Parses Newick format trees
*/
object DataParse{
  type Tree[A <: BioEnum]=Node[A] with RootNode[A]
  import org.modiphy.sequence.Alignment

  private def cleanTree(t:String)=t.split("\\s+").mkString("").split("\n").mkString("")

  /**
   Produce a parsed Tree and alignment from a newick file and a fasta file
  */
  def apply[A <: BioEnum](tree:String,alignment:Iterator[String],alphabet:A):(Tree[A],SimpleAlignment[A])={
    val aln = GenAlnParser(alignment)
    val seqMap = aln.foldLeft(Map[String,String]()){_+_}
    apply(tree,SimpleAlignment(seqMap,alphabet))
  }
  def fromFiles[A<:BioEnum](treeFile:String,alnFile:String,alphabet:A):(Tree[A],SimpleAlignment[A])={
    apply(
      scala.io.Source.fromFile(treeFile).getLines().map{_.trim}.mkString(""),
      scala.io.Source.fromFile(alnFile).getLines(),
      alphabet
    )
  }

  
  def apply[A <: BioEnum](tree:String,aln:SimpleAlignment[A]):(Tree[A],SimpleAlignment[A])={
    val root = new TreeParser[A](aln){def parseAll=parse(root,cleanTree(tree))}.parseAll.get.iNode.get.setRoot
    (root,aln)
  }

  def dropNodes[A <: BioEnum](tree:String,aln:SimpleAlignment[A]):(Tree[A],SimpleAlignment[A])={
    val t1:Tree[A] = apply(tree,aln)._1
    val t2 = t1.restrictTo(aln.sequenceNames.toSet).iNode.get.setRoot
    //println(t2)
    (t2.setAlign(aln).iNode.get.setRoot,aln)
  }
}



import DataParse._

import scala.actors.Actor
import scala.actors.Actor._

trait Node[A <: BioEnum] extends Actor with Logging{
  val label:Option[Int]
  val me = this
  def isLeaf=false
  val id:Int
  val isRoot=false
  val aln:Alignment[A]
  val alphabet:A=aln.alphabet
  
  def toLabelledString:String
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
  def setBranchLengths(m:Map[Int,Double]):Node[A]
  def getBranchLengths:List[Double]
  def branchTo:String
  def addLabels:Node[A]

  def nodes=this::descendentNodes

  def iNode:Option[INode[A]]=None

  def cromulent:Boolean= lengthTo > -scala.Double.Epsilon && children.foldLeft(true){(a,b)=>a && b.cromulent}

  def splitAln(i:Int):List[Node[A]]

  override def start={
   children.foreach{_.start}
   super.start
  }
  
}

object ReadTree{
  def fromFiles[A <: BioEnum](tree:String,aln:String,alphabet:A):Tree[A]={
    DataParse.fromFiles( tree,aln,alphabet)._1
  }
}



trait RootNode[A <: BioEnum] extends INode[A]{
  override def addLabels = factory(children.map{_.addLabels},aln,lengthTo,Some(id))
  val children:List[Node[A]]
  import scala.collection.immutable.IntMap
  val nodeMap=(this::descendentNodes).foldLeft[Map[Int,Node[A]]](IntMap[Node[A]]()){(m,n)=> m + ((n.id,n))}
  def apply(i:Int)=nodeMap(i)

  override val isRoot=true
  override def toString="("+children.mkString(",")+");"
  override def toLabelledString="("+children.map{_.toLabelledString}.mkString(",")+");"


  override def factory[B <: BioEnum](c:List[Node[B]],aln:Alignment[B],len:Double):Tree[B]=new INode[B](c,aln,len,id,label) with RootNode[B]
  override def factory[B <: BioEnum](c:List[Node[B]],aln:Alignment[B],len:Double,lab:Option[Int]):Tree[B]=new INode[B](c,aln,len,id,lab) with RootNode[B]

  override def getBranchLengths={
    descendentNodes.map{_.lengthTo}
  }
  override def setBranchLengths(l:List[Double]):Tree[A]={
    val nMap = descendentNodes.zip(l).foldLeft(Map[Int,Double]()){(m,t)=>
      m + ((t._1.id,t._2))
    }
    setBranchLengths(nMap)
    /*
    var listPtr = l
    val newchildren = children.map{c=>
      val c2 = c.setBranchLengths(listPtr)
      listPtr = listPtr.drop(c.descendentNodes.size + 1)
      c2
    }
    factory(newchildren,aln,l.head)
    */
  }
  override def setBranchLengths(m:Map[Int,Double]):Tree[A]={
    val newChildren = children.map{c=>c.setBranchLengths(m)}
    factory(newChildren,aln,0.0)
  }

  override def removeUseless:INode[A] with RootNode[A]=super.removeUseless.iNode.get.setRoot
  override def restrictTo(allowed:Set[String]):INode[A] with RootNode[A]=super.restrictTo(allowed).iNode.get.setRoot
  override val cromulent = children.foldLeft(true){(a,b) => a && b.cromulent}
  override def splitAln(i:Int):List[RootNode[A]]=super.splitAln(i).map{_.asInstanceOf[INode[A]].setRoot}

  def fPi = aln.getFPi

  override def act{
    react {
      case LogLikelihoodCalc(model,actor) =>
        children.foreach{c => c ! LikelihoodCalc(model)}
        model ! NewMatReq(self.asInstanceOf[Node[A]]) 
        act2(Nil,0,actor,None)
      case a:Any =>
        println(this + " received unexpected message: " + a)
    }
  }
  def act2(pl:List[List[Vector]],done:Int,requester:Actor,pi:Option[Vector]){
    //println("ROOT ACT2")
    if (done < numChildren+1){
      //println("0")
    react{
      case PartialLikelihoods(pl2,eMat)=>
        //println("1")
        plCalc ! PartialLikelihoods(pl2,eMat)
        act2(pl,done,requester,pi)
      case CalculatedPartialLikelihoods(pl2)=>
        //println("2")
        act2(pl2::pl,done+1,requester,pi)
      case MatReq(m,a,p)=>
        //println("3")
        if (m==me){
          act2(pl,done+1,requester,p)
        }else {
          println(this + " received unexpected message: " +  MatReq(m,a,p))
        }
      case a:Any => 
        println(this + " received unexpected message: " + a)
    }
  }
  else {
    act3(pl,requester,pi.get)
  }
}


  def act3(pl:List[List[Vector]],requester:Actor,pi:Vector){
        //println("ROOT ACT3")
        assert (pl.length==numChildren)
    val myPl = BasicLikelihoodCalc.combinePartialLikelihoods(pl)
    val likelihoods = myPl.map{vec=>
      val ans = pi.elements.zipWithIndex.map{t=>
        val(p,i)=t
        vec(i)*p
      }.foldLeft(0.0D){_+_}
      ans
    }
   // println("parallel likelihoods " + likelihoods)
    val lnL = likelihoods.zip(aln.pCount).foldLeft(0.0D){(i,j)=>i+math.log(j._1)*j._2}
    //println("LnL = " + lnL)
    //println("ROOT SENDING LOG LIKELIHOOD " + requester)
    requester ! LogLikelihood(lnL)
    act
  }




}


case class PartialLikelihoods(pl:List[Vector],eMat:Matrix)
case class CalculatedPartialLikelihoods(pl:List[Vector])
case class LogLikelihood(lnl:Double)
case class LikelihoodCalc(m:ActorModelComponent)
case class LogLikelihoodCalc(m:ActorModelComponent,receiver:Actor)

class PartialLikelihoodCalc extends Actor{
  def act{
     loop{
       react{
         case PartialLikelihoods(pl2,eMat) => {
           sender ! CalculatedPartialLikelihoods(BasicLikelihoodCalc.partialLikelihoodCalc(pl2,eMat))
         }
       }
     }
  }
}

abstract class MatExpActor extends Actor
class INode[A <: BioEnum](val children:List[Node[A]],val aln:Alignment[A],val lengthTo:Double,val id:Int,val label:Option[Int]) extends Node[A]{ 
  val plCalc = new PartialLikelihoodCalc
  plCalc.start

  def addLabels = factory(children.map{_.addLabels},aln,lengthTo,Some(id))
  def toLabelledString="("+children.map{_.toLabelledString}.mkString(",")+"):" + lengthTo + {if (label.isDefined){"#" + label.get}else {""}}


  def act{
    react {
      case LikelihoodCalc(model) => 
        //println("INode got")
        children.foreach{c => c ! LikelihoodCalc(model)}
        model ! NewMatReq(self.asInstanceOf[Node[A]])
        qMat(Nil,0,sender,model,None)
      case a:Any => 
        println(this + " WTF " + a)
    }
  }
  def qMat(pl:List[List[Vector]],done:Int,requester:OutputChannel[Any],matExp:ActorModelComponent,myEMat:Option[Matrix]){
    //println("DONE " + done)
    if (done < numChildren+1){
    react{
      case MatReq(m,mat,pi) =>
        //println("CASE 1")
        if (m==me){
          qMat(pl,done+1,requester,matExp,mat)
        }else {
          println(this + " received unexpected message: " + MatReq(m,mat,pi))
        }
      case PartialLikelihoods(pl2,eMat)=>
        //println("CASE 2")
        plCalc ! PartialLikelihoods(pl2,eMat)
        qMat(pl,done,requester,matExp,myEMat)
      case CalculatedPartialLikelihoods(pl2)=>
        //println("CASE 3")
        qMat(pl2::pl,done+1,requester,matExp,myEMat)
      case a:Any => 
        println(this + " received unexpected message: " + a)
      
    }
  }
  else {
    answer(pl,requester,myEMat.get)
  }
}


  def answer(pl:List[List[Vector]],requester:OutputChannel[Any],eMat:Matrix){
    assert(pl.length==numChildren)
    val myPl = BasicLikelihoodCalc.combinePartialLikelihoods(pl)
    requester ! PartialLikelihoods(myPl,eMat) 
    act
  }

  override def iNode=Some(this)
        

  val name=""

  def factory[B<:BioEnum](c:List[Node[B]],aln:Alignment[B],len:Double)=if (isRoot){new INode[B](c,aln,len,id,label) with RootNode[B]}else {new INode[B](c,aln,len,id,label)}
  def factory[B<:BioEnum](c:List[Node[B]],aln:Alignment[B],len:Double,lab:Option[Int])=if (isRoot){new INode[B](c,aln,len,id,lab) with RootNode[B]}else {new INode[B](c,aln,len,id,lab)}


  def numChildren = children.size
  def childElements:Iterator[Node[A]] = children.iterator
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
  def removeUseless:Node[A]={
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
        newChildren(0)
    }else{
      factory(newChildren,aln,lengthTo)
    }
  }

  def setAlign[B<:BioEnum](a2:Alignment[B]):Node[B]=factory(children.map{i=>i.setAlign(a2)},a2,lengthTo)
  def resize(bl:Double)=factory(children,aln,bl)
          
  override def toString={
    "("+children.mkString(",")+"):"+lengthTo
  }

  def setRoot=new INode[A](children,aln,0.0D,id,label) with RootNode[A]

  override def setBranchLengths(l:List[Double])={
    val nMap = descendentNodes.zip(l).foldLeft(Map[Int,Double]()){(m,t)=>
      m + ((t._1.id,t._2))
    }
    setBranchLengths(nMap)
  }
 
  def setBranchLengths(m:Map[Int,Double])={
    val newChildren = children.map{c=>c.setBranchLengths(m)}
    factory(newChildren,aln,m(id))
  }

  def getBranchLengths={
    //lengthTo :: children.map{c=>c.getBranchLengths}.flatten[Double]
    lengthTo :: descendentNodes.map{_.lengthTo}
  }

  def branchTo=descendents.mkString(",")
  def splitAln(i:Int)={
    val splitChildren = children.map{_.splitAln(i)}
    def makeNode(chi:List[Node[A]],nextChi:List[List[Node[A]]]):List[Node[A]]= nextChi.head match{
      case Nil => List(factory(chi,chi.head.aln,lengthTo))
      case _ => factory(chi,chi.head.aln,lengthTo)::makeNode(nextChi.map{_.head},nextChi.map{_.tail})
    }
    makeNode(splitChildren.map{_.head},splitChildren.map{_.tail})
  }
}

class Leaf[A <: BioEnum](val name:String,val aln:Alignment[A],val lengthTo:Double, val id:Int, val label:Option[Int]) extends Node[A]{

  override def isLeaf=true
  def children:List[Node[A]]=Nil
  def addLabels = factory(name,aln,lengthTo,Some(id))

  def toLabelledString = toString + {if (label.isDefined){"#" + label.get}else {""}}
  val sequence:Seq[alphabet.Value]=aln.getPatterns(name).asInstanceOf[Seq[alphabet.Value]]
  lazy val likelihoods:List[Vector]={
    sequence.map{a:alphabet.Value=>
    
      val vec = Vector(alphabet.matLength)
        alphabet.getNums(a).foreach{i=>
        vec(i)=1.0D}
        vec
    }.toList
  }

  def descendentNodes=List()
  def setBranchLengths(l:List[Double])={
  //   println("Node " + this + " setting branch length " + l.head)
     new Leaf[A](name,aln,l.head,id,label)
 }
  def setBranchLengths(m:Map[Int,Double])={
    new Leaf[A](name,aln,m(id),id,label)
  }
  def child(i:Int)=None
  def descendents=List(name)
  def resize(bl:Double) = new Leaf[A](name,aln,bl,id,label)
  def numChildren=0
  def removeUseless=this
  def setAlign[B<:BioEnum](a2:Alignment[B]):Node[B]=new Leaf[B](name,a2,lengthTo,id,label)
  def restrictTo(allowed:Set[String])=this
  override def toString=name + ":" + lengthTo
  def getBranchLengths=List(lengthTo)
  def factory[B <: BioEnum](n:String,a:Alignment[B],len:Double)=new Leaf[B](n,a,len,id,label)
  def factory[B <: BioEnum](n:String,a:Alignment[B],len:Double,lab:Option[Int])=new Leaf[B](n,a,len,id,lab)
  def branchTo=name

  def splitAln(i:Int)={
    val alnSplit = aln.split(i)
    val ans = alnSplit.map{a => new Leaf(name,a,lengthTo,id,label)}
    //println(ans.length)
    assert(ans.length==i)
    ans
  }
  def act{
    react {
      case LikelihoodCalc(model) => 
        //println("LeafNode Got")
        model ! NewMatReq(self.asInstanceOf[Node[A]])
        qMat(sender,model)
      case a:Any => 
        println(this + " received unexpected message: " + a)
      
    }
  }
  def qMat(requester:OutputChannel[Any],matExp:ActorModelComponent){
    //println("LeafNode Got2")
    react{
      case MatReq(m,eMat,pi) =>
        if (m==me){
          requester ! PartialLikelihoods(likelihoods,eMat.get)
          act
        }else{
          println(this + " received unexpected message: " + MatReq(m,eMat,pi))
        }
      case a:Any => 
        println(this + " received unexpected message: " + a)
    }
  }
}

