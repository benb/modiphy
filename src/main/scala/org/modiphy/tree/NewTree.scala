package org.modiphy.tree
import org.modiphy.math._
import scala.actors.Actor._ 
import scala.actors._ 
import org.modiphy.math.EnhancedMatrix._
import org.modiphy.sequence._
import scala.util.parsing.combinator._

case object BranchLength
case class UpdateDist(d:Double)
case class Unclean[A <: BioEnum](direction:DirBranch[A])
case class UncleanNode[A <: BioEnum](direction:Node[A])
case object UpdateMat
case class LikelihoodCalcDir[A <: BioEnum](model:ActorModelComponent,b:DirBranch[A])
case class LikelihoodCalc(model:ActorModelComponent)
case class LogLikelihood(d:Double)
case class LogLikelihoodCalc(model:ActorModelComponent,pi:Vector)
case class CalculatedPartialLikelihoods(pl:List[Vector])
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
  var branchCount:Int= -1
  def getINode(children:List[Node[A]],aln:Alignment[A], lengthTo:Double)={
    nodeCount += 1;
    val iNode =  new INode[A](nodeCount,lengthTo,aln)
    children.foreach{c =>
      branchCount += 1;
      new Branch(c,iNode,c.initialLengthTo,branchCount)
    }
    iNode
  }
  def getLeaf(name:String,aln:Alignment[A], lengthTo:Double)={nodeCount += 1;new Leaf[A](nodeCount,aln,name,lengthTo)}
}

/**
 Parses Newick format trees
*/
object DataParse{
  type Tree[A <: BioEnum]=INode[A]
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
    val root = new TreeParser[A](aln){def parseAll=parse(root,cleanTree(tree))}.parseAll.get
    (root,aln)
  }
}


class DirBranch[A <: BioEnum](val down:Node[A],var dist:Double,val up:Node[A],val myBranch:Branch[A]) extends Actor{
  def id = myBranch.id.toString + "->" + up.id
  override def toString={
    println(super.toString + " " + dist)
    down.dirToString(this)+":"+dist
  }


  def act{
    main(None)
  }
  //partial likelihoods go up the tree - so down -> to
  def main(pl:Option[List[Vector]]){
    react{
      case BranchLength=>
        reply(dist)
        main(pl)
      case UpdateDist(d)=>
        dist=d
        reply('ok)
        main(None)
      case UpdateMat=>
        up !? Unclean(this)
        reply('ok)
        main(None)
      case LikelihoodCalc(model)=>
        println("LikelihoodCalc " + id)
        if (!pl.isDefined){
          println("LikelihoodCalc Update " + id)
          //get PL up to 'descendent' node
          println("Sending LikelihoodCalcDir to " + down.id)
          down ! LikelihoodCalcDir(model,this)
          //get Matrix
          model ! NewMatReq(myBranch)
          getAns1(sender)
        }else {
          sender ! CalculatedPartialLikelihoods(pl.get)
          main(pl)
        }
       case UncleanNode(node)=>
         up !? Unclean(this)
         reply('ok)
      }
    }
  def getAns1(replyTo:OutputChannel[Any]){
    println("branch " + id + " waiting on MatReq")
    react{
      case MatReq(_,Some(mat),Some(pi))=>
        getAns2(replyTo,mat,pi)
    }
  }
  def getAns2(replyTo:OutputChannel[Any],mat:Matrix,pi:Vector){
    println("branch " + id + " waiting on PL")
    react{
      case CalculatedPartialLikelihoods(pl)=>
        val ans = BasicLikelihoodCalc.partialLikelihoodCalc(pl,mat)
        replyTo ! CalculatedPartialLikelihoods(ans)
        main(Some(ans))
    }
  }
}

class Branch[A <: BioEnum](val a:Node[A],val b:Node[A],var dist:Double,val id:Int) extends Actor{
  override def start={
    println("Branch " + id + " starting")
    endA.start
    endB.start
    super.start
  }

  val endA = new DirBranch(a,dist,b,this)
  val endB = new DirBranch(b,dist,a,this)
  b addBranch endA
  a addBranch endB 
  def act{
    loop{  
      react{
        case BranchLength=>
          println("Requested branch length " + dist)
          reply(dist)
        case UpdateDist(d)=>
          endA !? UpdateDist(d)
          endB !? UpdateDist(d)
          reply('ok)
      }
    }
  }
}


abstract class Node[A <: BioEnum] extends Actor{
  def dirToString(dir:DirBranch[A]):String
  def initialLengthTo:Double
  def addBranch(b:DirBranch[A])
  def bList:List[DirBranch[A]]
  def id:Int
  private var started=false
 def descendentBranches(dir:DirBranch[A]):List[Branch[A]]
}

class INode[A <: BioEnum](val id:Int,val initialLengthTo:Double,val aln:Alignment[A]) extends Node[A]{
 
 def startTree={
   println("Branches " + descendentBranches.map{_.id})
   descendentBranches.foreach{_.start}
   println("Nodes " + nodeList.map{_.id})
   nodeList.filter{_ != this}.foreach{n=>
     println("STARTING " + n.id)
     n.start
   }
   super.start
 }
 override def start={
   println("node " + id + " starting")
   super.start
 }
 
 def alphabet = aln.alphabet
 override def toString={
   "(" + branchEnds.reverse.mkString(",") + ");" //reverse to match original tree ordering
 }
 def dirToString(branch:DirBranch[A])={
   println("node(" + id + ") " + branch.id + " " + branchX(branch).map{_.id})
   "(" + branchX(branch).reverse.mkString(",") + ")"
 }
 def copy=DataParse(this.toString,aln)._1
 def descendentBranches(dir:DirBranch[A]) = {
   println("node " + this + "Coming from " + dir)
   branchX(dir).map{_.myBranch} ++ branchX(dir).map{b=> b.down.descendentBranches(b)}.flatten[Branch[A]]
 }
 lazy val nodeList={
   (descendentBranches.map{_.a} ++ descendentBranches.map{_.b}).removeDuplicates.sort{_.id < _.id}
 }
 lazy val descendentBranches={
   (branchEnds.map{_.myBranch} ++ branchEnds.map{b=> b.down.descendentBranches(b)}.flatten[Branch[A]]).sort{(a,b)=> a.id < b.id}
 }
 def getBranchLengths={
   import scala.collection.immutable.IntMap
   descendentBranches.map{b=>(b.id,b.dist)}.foldLeft[Map[Int,Double]](IntMap[Double]()){_+_}
 }
 def setBranchLengths(bl:List[Double])={
   descendentBranches.zip(bl).foreach{t=>
     t._1.dist = t._2
   }
   this
 }
 def setBranchLengths(f:Int=>Double)={
   descendentBranches.foreach{b=>
     b.dist = f(b.id)
   }
 }
 def branchX(dir:DirBranch[_])=branchEnds.filter{_.myBranch != dir.myBranch}
 
 private var branchEnds:List[DirBranch[A]]=Nil  
 def bList = branchEnds
  def addBranch(branch:DirBranch[A]){
    branchEnds = branch::branchEnds
  }

  def act{
    main
  }
  type Cache=Map[Option[DirBranch[A]],List[Vector]] 
  def main{
    react{
      case LogLikelihoodCalc(model,pi) =>
        println("Calc Log Likelihood " + id)
        var i = 0
        branchEnds.foreach{c => c ! LikelihoodCalc(model);i=i+1}
        logLikelihood(Nil,i,sender,pi)
      case LikelihoodCalcDir(model,branch)=>
        println("LikelihoodCalcDir " + id)
        //request partial likelihoods along each branch
        var i=0
        branchX(branch).foreach{c => c ! LikelihoodCalc(model);i=i+1}
        partialLikelihoods(branch,Nil,i,sender)
      case a:Any => 
        println("Node WTF " + a)
      }
    }
    
    def logLikelihood(plList:List[List[Vector]],toDo:Int,replyTo:OutputChannel[Any],pi:Vector){
      if (toDo==0){
        println("Combining " + plList.length + " partialLikelihoods")
        val ans = BasicLikelihoodCalc.combinePartialLikelihoods(plList)
        replyTo ! LogLikelihood(BasicLikelihoodCalc.logLikelihood(ans,pi))
        main
      }else {
        react{
          case CalculatedPartialLikelihoods(pl)=>
            logLikelihood(pl::plList,toDo-1,replyTo,pi)
        }
      }
    }

    def partialLikelihoods(branch:DirBranch[_],plList:List[List[Vector]],toDo:Int,replyTo:OutputChannel[Any]){
      if (toDo==0){
        println("Combining " + plList.length + " partialLikelihoods")
        val ans = BasicLikelihoodCalc.combinePartialLikelihoods(plList)
        replyTo ! CalculatedPartialLikelihoods(ans)
        main
      }else {
        react{
          case CalculatedPartialLikelihoods(pl) =>
            partialLikelihoods(branch,pl::plList,toDo-1,replyTo)
        }
      }
    }

}
class Leaf[A <: BioEnum](val id:Int,aln:Alignment[A],val name:String,val initialLengthTo:Double) extends Node[A]{
  override def start={
    println("Leaf " + id + " starting")
    super.start
  }
  def descendentBranches(n:DirBranch[A])=Nil
  
  override def toString=name
  def dirToString(d:DirBranch[A])=name
  var branchEnd:Option[DirBranch[A]]=None
  def bList=branchEnd.toList
  def addBranch(d:DirBranch[A]){
    branchEnd=Some(d)
  }
  val sequence:List[alphabet.Value]=aln.getPatterns(name).asInstanceOf[List[alphabet.Value]]
  val alphabet = aln.alphabet

  lazy val likelihoods:List[Vector]={
    sequence.map{a:alphabet.Value=>
      val vec = Vector(alphabet.matLength)
        alphabet.getNums(a).foreach{i=>
        vec(i)=1.0D}
        vec
    }.toList
  }

  def act{
    loop{
      react{
        case LikelihoodCalcDir(model,dir)=>
          println("Leaf " + id + " got LikelihoodCalcDir")
          sender ! CalculatedPartialLikelihoods(likelihoods)          
          println("Leaf " + id + " sent CalculatedPartialLikelihoods")
        case a:Any =>
          println("Leaf Node WTF " + a)
      }
    }
  } 
}
