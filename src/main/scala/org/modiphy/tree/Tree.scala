package org.modiphy.tree
import org.modiphy.math._
import scala.actors.Actor._ 
import scala.actors._ 
import org.modiphy.math.EnhancedMatrix._
import org.modiphy.sequence._
import scala.util.parsing.combinator._
import tlf._

case object BranchLength
case class UpdateDist(d:Double)
case class Unclean[A <: BioEnum](direction:DirBranch[A])
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
  def id = down.id + "--" + myBranch.id.toString + "->" + up.id
  override def toString={
    down.dirToString(this)+":"+dist
  }

  lazy val reverse=myBranch.other(this)


  def act{
    main(None)
  }
  //partial likelihoods go up the tree - so down -> to
  def main(pl:Option[List[Vector]]){
//    if (myBranch.id==0){println("node " + id + " " + pl.isDefined)}
    react{
      case BranchLength=>
      if (myBranch.id==0){println("node " + id + " BranchLength")}
        reply(dist)
        main(pl)
      case UpdateDist(d)=>
      if (myBranch.id==0){println("node " + id + " " + UpdateDist(d))}
        println(id + " unclean")
        dist=d
        if (pl.isDefined){
          up !? Unclean(this)
        }else {
          println(id + " already unclean so stopping")
        }
        reply('ok)
        main(None)
      case UpdateMat=>
      println("node " + id + " " + UpdateMat)
        if (pl.isDefined){
          println(id + " Waiting")
          up !? Unclean(this)
          println(id + " Got")
        }
        reply('ok)
        main(None)
      case LikelihoodCalc(model)=>
      if (myBranch.id==0){println("node " + id + " LikelihoodCalc")}
      if (myBranch.id==2){println("node " + id + " LikelihoodCalc")}
        
        if (pl.isEmpty){
          //if(myBranch.id==0){println("Branch " + this.id + " recalculating")}
          println("Branch " + this.id + " recalculating")
          //get PL up to 'descendent' node
          down ! LikelihoodCalcDir(model,this)
          //get Matrix
          model ! NewMatReq(myBranch)
          getAns1(sender)
        }else {
          println("Branch " + this.id + " cached")
          sender ! CalculatedPartialLikelihoods(pl.get)
          main(pl)
        }
       case a:Any =>
         println("dirbranch " + id + " unexpected msg " + a)
      }
    }
  def getAns1(replyTo:OutputChannel[Any]){
    println(id + " getAns1")
    react{
      case MatReq(_,Some(mat),Some(pi))=>
        getAns2(replyTo,mat,pi)
    }
  }
  def getAns2(replyTo:OutputChannel[Any],mat:Matrix,pi:Vector){
    println(id + " getAns2")
    react{
      case CalculatedPartialLikelihoods(pl)=>
        val ans = BasicLikelihoodCalc.partialLikelihoodCalc(pl,mat)
        replyTo ! CalculatedPartialLikelihoods(ans)
        //main(Some(ans))
    println(id + " getAns3")
        main(None)
    }
  }
}

class Branch[A <: BioEnum](val a:Node[A],val b:Node[A],var dist:Double,val id:Int) extends Actor{
  override def start={
    endA.start
    endB.start
    super.start
  }


  var chainedBranches:List[Branch[_]]=Nil
  def chain(b:Branch[_]){
    chainedBranches=b::chainedBranches
  }

  val endA = new DirBranch(a,dist,b,this)
  val endB = new DirBranch(b,dist,a,this)

  def other(b:DirBranch[A])={
    if (b==endA){endB}else{endA}
  }

  b addBranch endA
  a addBranch endB 
  def act{
    loop{  
      react{
        case BranchLength=>
          if (id==0){ println("Branch " + id + " length " + dist)}
          reply(dist)
        case UpdateDist(d)=>
          endA !? UpdateDist(d)
          endB !? UpdateDist(d)
          chainedBranches.foreach{b=>
            println("Sending to chained " + id)
            b !? UpdateDist(d)
          }
          dist=d
          reply('ok)
        case UpdateMat =>
          endA !? UpdateMat
          endB !? UpdateMat
          chainedBranches.foreach{b=>
            b!? UpdateMat
          }
          reply('ok)
       case a:Any =>
         println("branch " + id + " unexpected msg " + a)
      }
    }
  }
}


abstract class Node[A <: BioEnum] extends Actor with Logging{
  def dirToString(dir:DirBranch[A]):String
  def initialLengthTo:Double
  def addBranch(b:DirBranch[A])
  def bList:List[DirBranch[A]]
  def id:Int
  private var started=false
 def descendentBranches(dir:DirBranch[A]):List[Branch[A]]
}

class INode[A <: BioEnum](val id:Int,val initialLengthTo:Double,val aln:Alignment[A]) extends Node[A]{

 var started=false 
 override def start={
   if (started==true){
     println("STARTED TWICE")
   }
   started=true
   super.start
 }
 
 def splitAln(i:Int)={
   aln.split(i).map{this.copy(_).chainFrom(this)}
 }
 def chainFrom[B<:BioEnum] (other:INode[B])={
   descendentBranches.zip(other.descendentBranches).foreach{t=>
     t._2 chain t._1
   }
   this
 }
 def startTree={
   descendentBranches.foreach{_.start}
   nodeList.filter{_ != this}.foreach{n=>
     n.start
   }
   this.start
 }
 
 def alphabet = aln.alphabet
 override def toString={
   "(" + branchEnds.reverse.mkString(",") + ");" //reverse to match original tree ordering
 }
 def dirToString(branch:DirBranch[A])={
   "(" + branchX(branch).reverse.mkString(",") + ")"
 }
 def copy=DataParse(this.toString,aln)._1
 def copy[B <: BioEnum](newAln:Alignment[B])=DataParse(this.toString,newAln)._1
 def descendentBranches(dir:DirBranch[A]) = {
   branchX(dir).map{_.myBranch} ++ branchX(dir).map{b=> b.down.descendentBranches(b)}.flatten[Branch[A]]
 }
 lazy val nodeList={
   (descendentBranches.map{_.a} ++ descendentBranches.map{_.b}).removeDuplicates.sort{_.id < _.id}
 }
 lazy val iNodeList={
   nodeList.filter{_.isInstanceOf[INode[A]]}.map{_.asInstanceOf[INode[A]]}
 }
 lazy val descendentBranches={
   (branchEnds.map{_.myBranch} ++ branchEnds.map{b=> b.down.descendentBranches(b)}.flatten[Branch[A]]).sort{(a,b)=> a.id < b.id}
 }
 def getBranchLengths={
   descendentBranches.map{_.dist}.toArray//.foldLeft[Map[Int,Double]](IntMap[Double]()){_+_}
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
   this
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
        debug{"Calc Log Likelihood " + id}
        var i = 0
        branchEnds.foreach{c => c ! LikelihoodCalc(model);i=i+1}
        logLikelihood(Nil,i,sender,pi)
      case LikelihoodCalcDir(model,branch)=>
        debug{"LikelihoodCalcDir " + id}
        //request partial likelihoods along each branch
        var i=0
        branchX(branch).foreach{c => c ! LikelihoodCalc(model);i=i+1}
        partialLikelihoods(branch,Nil,i,sender)
      case Unclean(dir)=>
        branchX(dir).map{_.reverse}.foreach{ _ !? UpdateMat}
        reply('ok)
        main
      case a:Any => 
        println("Node received unexpected message " + a)
      }
    }
    
    def logLikelihood(plList:List[List[Vector]],toDo:Int,replyTo:OutputChannel[Any],pi:Vector){
      if (toDo==0){
        val ans = BasicLikelihoodCalc.combinePartialLikelihoods(plList)
        val ans2 = BasicLikelihoodCalc.logLikelihood(ans,pi,aln)
        replyTo ! LogLikelihood(ans2)
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
          sender ! CalculatedPartialLikelihoods(likelihoods)          
        case Unclean(dir)=>
          reply('ok)
        case a:Any =>
          println("Leaf Node received unexpected message " + a)
      }
    }
  } 
}
