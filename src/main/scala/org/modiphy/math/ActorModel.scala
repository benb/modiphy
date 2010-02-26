package org.modiphy.math
import EnhancedMatrix._
import org.modiphy.tree._
import org.modiphy.tree.DataParse.Tree
import org.modiphy.sequence._
import scala.collection.immutable.{Map,IntMap}
import scala.actors.Actor
import scala.actors.Actor._
import scala.actors.OutputChannel
import tlf.Logging

abstract class ActorModelComponent extends Actor with Logging{
}
case class Params(list:List[ActorParamComponent])
object GetParams extends Params(Nil)
class MatBuilder(m:Option[Map[Node[_],Matrix]])
case class MatReq(n:Node[_],m:Option[Matrix],pi:Option[Vector])
object NewMatReq{
  def apply(n:Node[_])=MatReq(n,None,None)
}
case class Unclean(sender:Actor)
case object RequestOpt
object RequestParam 
abstract class ParamName

case class Pi(i:Int) extends ParamName
case class S(i:Int) extends ParamName
case class BranchLengths(i:Int) extends ParamName
case class Alpha(i:Int) extends ParamName
case class InvarPrior(i:Int) extends ParamName
case class Sigma(i:Int) extends ParamName
case class SingleParam(p:ParamName) extends ParamName
//case class JoinedParam(p:List[ParamName]) extends ParamName

object LazyP{
  type PFact={def apply(i:Int):ParamName}
  implicit def lazymaker(p:PFact)=p(0)
}

case class ParamUpdate[A](d:A)
case class ParamChanged[A](p:ParamName,param:A){
  def toParamUpdate = ParamUpdate(param)
}

object Lower
object Upper


/**
 Optimiser may have different view of params (log of real params, for example)
*/
case class OptUpdate(d:Array[Double])

class BasicActorModel(piParam:ActorPiComponent,sParam:ActorSComponent,rec1:Actor) extends ActorModelComponent{
  val piName = piParam.name
  val sName = sParam.name
  override def start={
    piParam.start
    sParam.start
    piParam addActor this
    sParam addActor this
    rec1.start
    super.start
  }
  def act{
    finest{this + " STARTING"}
    piParam ! RequestParam
    sParam ! RequestParam
    initialise(None,None)
  }
  def initialise(s:Option[Matrix],pi:Option[Vector]){
    finest{"BasicActor initialise " + s + " " + pi}
    if (s.isDefined && pi.isDefined){
      finest{this + " STARTED"}
      main(s.get,pi.get,None)
    }
    react{
      case ParamChanged(piName,v:Vector)=>initialise(s,Some(v))
      case ParamChanged(sName,m:Matrix)=>initialise(Some(m),pi)
    }
  }
  def main(s:Matrix,pi:Vector,mat:Option[Matrix]){
    react{
      case Params(l)=>
        rec1 forward Params(piParam :: sParam :: l)
        main(s,pi,mat)
      case MatReq(n,oldM,oldPi)=>
        val myMat = if (mat.isDefined){
          mat.get
        }else {
          s.sToQ(pi)
        }
        rec1 forward MatReq(n,Some(myMat),Some(pi))
        main(s,pi,Some(myMat))
        case ParamChanged(piName,v:Vector)=>
          rec1 forward Unclean(piParam)
          main(s,v.copy,None)
        case ParamChanged(sName,m:Matrix)=>
          rec1 forward Unclean(sParam)
          main(m.copy,pi,None)
     }
  }
}

class THMMActorModel(sigmaParam:ActorFullSComponent,numClasses:Int,rec1:Actor) extends ActorModelComponent{
  val sigmaName = sigmaParam.name
  override def start={
    sigmaParam.start
    sigmaParam addActor this
    rec1.start
    super.start
  }
  /**
   pi must be the new length, m is pi.length - numAA
  */
  def applyMat(m:Matrix,pi:Vector,sigma:Matrix):Matrix={

   val qStart = m.copy
   val numAlpha = m.rows / numClasses
   for (i <- 0 until numClasses){
        for (j <- i+1 until numClasses){
          for (x <- 0 until numAlpha){
            qStart(i * numAlpha + x, j * numAlpha + x) = sigma(i,j) *  pi(j * numAlpha + x) // i->j transition
            qStart(j * numAlpha + x, i * numAlpha + x) = sigma(i,j) * pi(i * numAlpha + x) // j->i transition
          }
        }
      }
      qStart.fixDiag
      qStart
  }
  def act{
    finest{this + " STARTING"}
   // sigmaParam ! RequestParam
  //  initialise(None)
    main(sigmaParam.myParam,None)
  }
  def initialise(sigma:Option[Matrix]){
    finest{"Sigma Initialize " + sigma}
    react{
      case ParamChanged(sigmaName,a:Matrix)=>
        finest{this + " STARTED"}
        main(a,None)
    }
    initialise(None)

  }
  def main(sigma:Matrix,mat:Option[Matrix]){
    react{
      case Params(l)=> 
        rec1 forward Params(sigmaParam :: l)
        main(sigma,mat)
      case ParamChanged(sigmaName,array:Matrix)=>
        rec1 forward Unclean(sigmaParam)
        main(array,None)
      case MatReq(n,Some(m),Some(p))=>
        val myMat = if (mat.isDefined){
          mat.get
        }else {
          applyMat(m,p,sigma)
        }
        rec1 forward MatReq(n,Some(myMat),Some(p))
        main(sigma,Some(myMat))
      case Unclean(a) => 
        rec1 forward Unclean(a)
        main(sigma,None)
    }
  }
}




class InvarActorModel(priorParam:ActorProbComponent,piParam:ActorPiComponent,numClasses:Int,rec1:Actor) extends ActorModelComponent{
  val priorName=priorParam.name
  val piName = piParam.name
  override def start={
    priorParam.start
    priorParam addActor this
    piParam.start
    piParam addActor this
    rec1.start
    super.start
  }
  /**
   pi must be the new length, m is pi.length - numAA
  */
  def applyMat(m:Matrix,pi:Vector):Matrix={
    val newMatSize = pi.size
    val mat = Matrix(newMatSize,newMatSize)
    mat.viewPart(0,0,m.rows,m.rows).assign(m).normalize(pi)
    mat
  }
  def applyPi(pi:Vector,oldPi:Vector,prior:Double):Vector={
    val newPi = Vector(oldPi.size + pi.size)
    newPi.viewPart(0,oldPi.size).assign(oldPi * (1.0D-prior))
    newPi.viewPart(oldPi.size,pi.size).assign(pi * prior)
    newPi
  }
  def act{
    finest{this + " STARTING"}
    piParam ! RequestParam
    priorParam ! RequestParam
    initialise(None,None)
  }
  def initialise(prior:Option[Double],pi:Option[Vector]){
    finest{"Invar Initialise " + prior + " " + pi}
    if (prior.isDefined && pi.isDefined){
      finest{this + " STARTED"}
      main(prior.get,pi.get,None,None)
    }
    react{
      case ParamChanged(priorName,d:Double)=>initialise(Some(d),pi)
      case ParamChanged(piName,v:Vector)=>initialise(prior,Some(v))
    }
  }
  def main(prior:Double,rawpi:Vector,processedPi:Option[Vector],mat:Option[Matrix]){
    react{
      case Params(l)=> 
        rec1 forward Params(priorParam :: piParam :: l)
        main(prior,rawpi,processedPi,mat)
      case ParamChanged(priorName,d:Double)=>
        rec1 forward Unclean(priorParam)
        main(d,rawpi,None,None)
      case ParamChanged(piName,v:Vector)=>
        rec1 forward Unclean(piParam)
        main(prior,v,None,None)
      case MatReq(n,Some(m),Some(p))=>
        val myPi = if (processedPi.isDefined){
          processedPi.get
        }else{
          applyPi(rawpi,p,prior)
        }
        val myMat = if (mat.isDefined){
          mat.get
        }else {
          applyMat(m,myPi)
        }
        rec1 forward MatReq(n,Some(myMat),Some(myPi))
        main(prior,rawpi,Some(myPi),Some(myMat))
      case Unclean(a) => 
        rec1 forward Unclean(a)
        main(prior,rawpi,None,None)
    }
  }
}

class GammaActorModel(shape:ActorGammaComponent,numClasses:Int,rec1:Actor) extends ActorModelComponent{
  val shapeName = shape.name
  override def start={
    shape.start
    shape addActor this
    rec1.start
    super.start
  }
  val gamma = new Gamma(numClasses)
  
  def applyMat(m:Matrix,alpha:Double,pi:Vector):Matrix={
    val r = gamma(alpha)
    debug{"GAMMA " + alpha + " " + r.mkString(" ")}
    val myMat = Matrix(m.rows*numClasses,m.columns * numClasses)
    val numAlpha = m.rows 
    (0 until numClasses).foreach{i=>
      myMat.viewPart(i * numAlpha,i * numAlpha,numAlpha,numAlpha).assign(m)* r(i) //(m.normalize(pi,r(i))))
    }
    myMat
  }
  def applyPi(pi:Vector):Vector={
    val numAlpha = pi.size
    val myPi = Vector(numAlpha * numClasses)
    val postVals = pi / numClasses.toDouble
    (0 until numClasses).foreach{i=>
      myPi.viewPart(i*numAlpha,numAlpha).assign(postVals)
    }
    myPi
  }
  def act{
    finest{ this + " STARTING"}
    shape ! RequestParam
    react{
      case ParamChanged(shapeName,alpha:Double)=>
        finest{ this + " STARTED"}
        main(alpha,None,None)
    }
  }
  def main(alpha:Double,mat:Option[Matrix],pi:Option[Vector]){
    react{
      case Params(l)=>
        rec1 forward Params(shape :: l)
        main(alpha,mat,pi)
      case MatReq(n,m,p)=>
        val myMat = if (mat.isDefined){
          mat.get
        }else {
          applyMat(m.get,alpha,p.get)
        } 
        val myPi = if(pi.isDefined){
          pi.get
        }else{
          applyPi(p.get)
        }
        rec1 forward MatReq(n,Some(myMat),Some(myPi))
        main(alpha,Some(myMat),Some(myPi))
      case ParamChanged(shapeName,d:Double)=>
        rec1 forward Unclean(shape)
        debug{"NEWGAMMA " + d}
        main(d,None,pi)
      case Unclean(a) =>
        rec1 forward Unclean(a)
        main(alpha,None,None)
      case a:Any => 
        warning{"Gamma WTF " + a}
        main(alpha,mat,pi)
    }
  }
}

class ForkActor[A <: BioEnum](tree:Tree[A],rec1Map:Map[Int,Actor]) extends ActorModelComponent{
  val rec = rec1Map.values.toList.removeDuplicates
  val numRec = rec.length
  override def start={
    rec1Map.values.foreach{v=>
      v.start
    }
    super.start
  }

  def act{
    loop{
      react{
        case Params(l)=>
          rec.foreach{v=>
            v ! Params(l)             
          }
          getListReplies(sender,numRec,Nil)
        case Unclean(a) =>
          rec.foreach{v=>
            v ! Unclean(self)
          }
          getCleanReplies(a,numRec,self)
        case MatReq(n,m,p)=>
          if(rec1Map contains n.id){
            rec1Map(n.id) forward MatReq(n,m,p)
          }else {
            if (! n.isRoot){warning{"Map does not contain " + n.id + " " + n}}
            sender ! MatReq(n,m,p) 
          }
      }
    }
  }

  def getListReplies(output:OutputChannel[Any],toDo:Int,l:List[ActorParamComponent]){
    if (toDo==0){
      output ! Params(l) 
      act
    }
    react{
      case Params(l2) => 
        getListReplies(output,toDo-1,l2++l)
    }
  }

  def getCleanReplies(output:Actor,toDo:Int,me:Actor){
    if (toDo==0){
      output ! Unclean(output)
      act
    }
    react{
      case Unclean(`me`) => 
        getCleanReplies(output,toDo-1,me)
    }
  }

}

class BasicSingleExpActorModel[A <: BioEnum](tree:Tree[A],branchLengthParams:ActorTreeComponent[A],rec1:Option[Actor]) extends ActorModelComponent{
  val blName = branchLengthParams.name
  val nodes = tree.descendentNodes.toArray //root has no length

  override def start={
    branchLengthParams.start
    branchLengthParams addActor this
    if (rec1.isDefined){rec1.get.start}
    super.start
  }
  var eigen:MatrixExponential=null 

  def act{
    main(None,nodes.foldLeft(Map[Int,Double]()){(m,n)=>m+((n.id,n.lengthTo))})
  }
  def main(eigen:Option[MatrixExponential],lengths:Map[Int,Double]){
    case class ExpReq(n:Node[_],e:MatrixExponential,lengthTo:Double,pi:Vector)
      react{
        case Params(l)=>
          if (rec1.isDefined){rec1.get forward Params(branchLengthParams :: l)}
          else {sender ! Params(branchLengthParams :: l)}
          main(eigen,lengths)
        case MatReq(n,Some(m),Some(pi)) => 
          if (n.isRoot){//don't need exp(qt)
            if (rec1.isDefined){rec1.get forward MatReq(n,None,Some(pi))}
            else { sender ! MatReq(n,None,Some(pi))}
            main(eigen,lengths)
          }
          val myEigen = if (eigen.isEmpty){
            val ans = new MatExpNormal(m,pi)
            ans
          }else {
            eigen.get
          }
          new Actor{
            def act{
              react{
                case ExpReq(n,e,lengthTo,pi)=> 
                  val ans = e.exp(lengthTo)
                  if (rec1.isDefined){rec1.get forward MatReq(n,Some(ans),Some(pi))}
                  else { sender ! MatReq(n,Some(ans),Some(pi)) }
                  exit
              }
            }
          }.start forward ExpReq(n,myEigen,lengths(n.id),pi)
          main(Some(myEigen),lengths)
        case Unclean(a)=>
          if (rec1.isDefined){rec1.get forward Unclean(a)} else {
            a ! Unclean(a)
          }
          main(None,lengths)
        case ParamChanged(blName,a:Array[Double])=>
          if (rec1.isDefined){rec1.get forward Unclean(branchLengthParams)} else {branchLengthParams ! Unclean(branchLengthParams)}
          val nodeMap = nodes.zip(a).foldLeft(Map[Int,Double]()){(m,t)=>m + ((t._1.id,t._2))}
          main(eigen,nodeMap)
      }
  }
}

trait ActorParamComponent extends Actor{
  def name:ParamName
  var modelComp:List[Actor]=Nil
  def addActor(l:Actor){
    modelComp = l::modelComp
  }
}


case class SingleParamWrapper(p:ActorParamComponent) extends ActorParamComponent{
  start
  val name = SingleParam(p.name)
  def act{
    val mainParam = (p !? RequestOpt).asInstanceOf[Array[Double]]
    main(Array.make(mainParam.length,mainParam(0)))
  }
  def lower = (p !? Lower).asInstanceOf[Array[Double]](0)
  def upper = (p !? Upper).asInstanceOf[Array[Double]](0)
  def main(param:Array[Double]){
    react{
    case RequestParam=>
      reply(ParamChanged(name,param(0)))
      main(param)
    case ParamUpdate(x:Array[Double])=>
      val newParam = Array.make(param.length,x(0))
      p forward OptUpdate(newParam)
      main(newParam)
    case OptUpdate(x:Array[Double]) =>
      val newParam = Array.make(param.length,x(0))
      p forward OptUpdate(newParam)
      main(newParam)
    case RequestOpt =>
      reply(Array(param(0)))
      main(param)
    case Lower =>
      sender ! Array(lower)
      main(param)
    case Upper =>
      sender ! Array(upper)
      main(param)
    }
  }
}
/*
case class JoinedParamWrapper(p:List[ActorParamComponent]) extends ActorParamComponent{
 start
 val name = JoinedParam(p)
 def act{
   val current = p.map{p2 => (p2 !? RequestOpt).asInstanceOf[Array[Double]].toList}.flatten[Double].toArray
 }
 
}*/

abstract class AbstractActorParam[A] extends ActorParamComponent with Logging{
  type Param = {
    def getParams:Array[Double]
    def setParams(a:Array[Double]):Unit
  }
  def internal:Param
  def myParam:A
  override def toString={
    name.toString + " " + internal.getParams.toList
  }

  def lower:Double
  def upper:Double
  def lowerArray = Array.make(internal.getParams.length,lower)
  def upperArray = Array.make(internal.getParams.length,upper)
  //override this if raw params are different!
  def setRaw(p:Array[Double]){internal setParams p}
  def waitOn(i:Int,o:OutputChannel[Any]){
    finest{super.toString + " WAITING ON " + i + ":::" + this}
    react{
      case Unclean(_)=>
        finest{"GOT " + i}
        if (i==1){
          o ! 'ok
          act
        }else {
          waitOn(i-1,o)
        }
    }
  }
  private val myHandler:PartialFunction[Any,Unit] = {
        case RequestParam=>
         reply(ParamChanged(name,myParam))
        case ParamUpdate(x:Array[Double])=> 
          setRaw(x)
          modelComp.foreach{c=>
            finest{"SENDING to " + c}
            c ! ParamChanged(name,myParam)
          }
          waitOn(modelComp.length,sender)
        case ParamUpdate(v:Vector)=>
          setRaw(v.toArray)
          modelComp.foreach{c=>
            c ! ParamChanged(name,myParam)
          }
          waitOn(modelComp.length,sender)
        case OptUpdate(x)=>
          internal setParams x
          modelComp.foreach{c=>
            c ! ParamChanged(name,myParam)
          }
          waitOn(modelComp.length,sender)
        case RequestOpt=>
          reply(internal.getParams)
        case Lower=>
          sender ! lowerArray
        case Upper=>
          sender ! upperArray
  }
  def handler:PartialFunction[Any,Unit]= {
     myHandler 
  }

  def act{
    finest{"ACT"}
    loop{
      receive(handler)
    }
  }
}
class ActorPiComponent(pi:Vector,val name:ParamName) extends AbstractActorParam[Vector]{
  val internal = new PiParam(pi)
  def lower = -10 // this is log-odds scaled prob
  def upper = 10
  def myParam = pi
  override def setRaw(p:Array[Double]){
    internal setPi p
  }
  override def toString=name.toString + " " + internal.view
}
class ActorSComponent(s:Matrix,val name:ParamName) extends AbstractActorParam[Matrix] with SMatUtil{
  val internal = new SMatParam(s,name.toString)
  def myParam = s
  def lower = 0.0
  def upper = 200
  private val myHandler:PartialFunction[Any,Unit]={
      case ParamUpdate(m:Matrix)=>
        setRaw(linearSMatFull(m).toArray)
        modelComp.foreach{c=>
          c ! ParamChanged(name,myParam)
        }
        waitOn(modelComp.length,sender)
  }
  override def handler={
    myHandler orElse (super.handler)
  }

}
class ActorFullSComponent(s:Matrix,val name:ParamName) extends AbstractActorParam[Matrix] with SMatUtil{
  val internal = new FullSMatParam(s,name.toString)
  def lower = 0.0D
  def upper = 200
  def myParam = s
  private val myHandler:PartialFunction[Any,Unit]={
      case RequestParam=>
          reply(ParamChanged(name,myParam))
      case ParamUpdate(m:Matrix)=>
        setRaw(linearSMatFull(m).toArray)
        modelComp.foreach{c=>
          c ! ParamChanged(name,myParam)
        }
        waitOn(modelComp.length,sender)
  }
  override def handler={
    myHandler orElse (super.handler)
  }
}


class ActorTreeComponent[B <: BioEnum](tree:Tree[B],name:ParamName) extends ActorArrayComponent((tree.descendentNodes.map(_.lengthTo).toArray),name,0.0,100.0)
class ActorDoubleComponent(param:Double,val name:ParamName,val lower:Double,val upper:Double) extends AbstractActorParam[Double]{
  class DoubleParam{
    var myP:Double=param
    def getParams=Array(myP)
    def setParams(a:Array[Double]){myP=a(0)}
  }
  val internal = new DoubleParam
  def myParam = internal.getParams(0)

  private val myHandler:PartialFunction[Any,Unit]={
    case ParamUpdate(x:Double)=>
      internal.myP=x
      modelComp.foreach{c=>
        c ! ParamChanged(name,myParam)
      }
      waitOn(modelComp.length,sender)
      reply('ok)
  }
  override def handler={
    myHandler orElse (super.handler)
  }
}
class ActorArrayComponent(param:Array[Double],val name:ParamName,val lower:Double, val upper:Double) extends AbstractActorParam[Array[Double]]{
  class ArrayParam{
    var myP=param
    def getParams = myP.toArray
    def setParams(a:Array[Double]){Array.copy(a,0,myP,0,Math.min(a.length,myP.length))}
  }
  val internal = new ArrayParam
  def myParam = internal.getParams
}
class ActorProbComponent(prob:Double,name:ParamName,min:Double,max:Double) extends ActorDoubleComponent(prob,name,min,max){
  def this(prob:Double,name:ParamName)=this(prob,name,0,0.9)
}
class ActorGammaComponent(alpha:Double,name:ParamName) extends ActorDoubleComponent(alpha,name,0.01,1000)
object SimpleModel{
  def apply[A <: BioEnum](tree:Tree[A])={
    val pi = new ActorPiComponent(WAG.pi,Pi(0))
    val s = new ActorSComponent(WAG.S,S(0))
    val branchLength = new ActorTreeComponent(tree,BranchLengths(0))
    val components = new BasicActorModel(pi,s, new BasicSingleExpActorModel(tree,branchLength,None))
    val pList = List(pi,s,branchLength)
    val pMap = pList.map{p => (p.name,p)}.foldLeft(Map[ParamName,ActorParamComponent]()){_+_}
    new ActorModel(tree,components,pMap)
  }
}
object GammaModel{
  def apply[A <: BioEnum](tree:Tree[A])={
    val pi = new ActorPiComponent(WAG.pi,Pi(0))
    val s = new ActorSComponent(WAG.S,S(0))
    val branchLength = new ActorTreeComponent(tree,BranchLengths(0))
    val alpha = new ActorGammaComponent(0.5D,Alpha(0))
    val components = new BasicActorModel(pi,s, 
      new GammaActorModel(alpha,tree.alphabet.numClasses,
       new BasicSingleExpActorModel(tree,branchLength,None)))
    val pList = List(pi,s,branchLength,alpha)
    val pMap = pList.map{p => (p.name,p)}.foldLeft(Map[ParamName,ActorParamComponent]()){_+_}
    new ActorModel(tree,components,pMap)
  }
}
object InvarGammaModel{
  def apply[A <: BioEnum](tree:Tree[A])={
    val pi = new ActorPiComponent(WAG.pi,Pi(0))
    val s = new ActorSComponent(WAG.S,S(0))
    val branchLength = new ActorTreeComponent(tree,BranchLengths(0))
    val alpha = new ActorGammaComponent(0.5D,Alpha(0))
    val invarPrior = new ActorProbComponent(0.2D,InvarPrior(0))
    val components = new BasicActorModel(pi,s,
      new GammaActorModel(alpha,tree.alphabet.numClasses-1,
        new InvarActorModel(invarPrior,pi,tree.alphabet.numClasses,
          new BasicSingleExpActorModel(tree,branchLength,None))))
    val pList = List(pi,s,branchLength,alpha,invarPrior)
    val pMap = pList.map{p => (p.name,p)}.foldLeft(Map[ParamName,ActorParamComponent]()){_+_}
     new ActorModel(tree,components,pMap)
  }
}

object InvarThmmModel{
def apply[A <: BioEnum](tree:Tree[A])={
  val numClasses = tree.alphabet.numClasses
    val pi = new ActorPiComponent(WAG.pi,Pi(0))
    val s = new ActorSComponent(WAG.S,S(0))
    val branchLength = new ActorTreeComponent(tree,BranchLengths(0))
    val alpha = new ActorGammaComponent(0.5D,Alpha(0))
    val invarPrior = new ActorProbComponent(0.2D,InvarPrior(0))
    val sigma = new ActorFullSComponent(Matrix(numClasses,numClasses),Sigma(0))
    val components = new BasicActorModel(pi,s,
      new GammaActorModel(alpha,numClasses-1,
        new InvarActorModel(invarPrior,pi,numClasses,
          new THMMActorModel(sigma,numClasses,
            new BasicSingleExpActorModel(tree,branchLength,None)))))
    val pList = List(pi,s,branchLength,alpha,invarPrior,sigma)
    val pMap = pList.map{p => (p.name,p)}.foldLeft(Map[ParamName,ActorParamComponent]()){_+_}
     new ActorModel(tree,components,pMap)
  }
}

object BranchSpecificThmmModel{
  def apply[A <: BioEnum](tree:Tree[A]):ActorModel={
    var i = -1
    val nodeMap=tree.descendentNodes.foldLeft[Map[Int,Int]](IntMap[Int]()){(m,n)=>
        i=i+1
        m + ((n.id,i))
    }
    apply(tree,nodeMap)
  }
  def apply[A <: BioEnum](tree:Tree[A],map:Map[Int,Int]):ActorModel={
    val numClasses = tree.alphabet.numClasses
    val pi = new ActorPiComponent(WAG.pi,Pi(0))
    val s = new ActorSComponent(WAG.S,S(0))
    val branchLength = new ActorTreeComponent(tree,BranchLengths(0))
    val alpha = new ActorGammaComponent(0.5D,Alpha(0))
    val invarPrior = new ActorProbComponent(0.2D,InvarPrior(0))
    var pList:List[ActorParamComponent] = List(pi,s,branchLength,alpha,invarPrior)
    val modelNums = map.values.toList.removeDuplicates
    val modelMapTmp = modelNums.map{a=>
      val sigma = new ActorFullSComponent(Matrix(numClasses,numClasses),Sigma(a))
      pList = sigma :: pList
      val ans = new THMMActorModel(sigma,numClasses, new BasicSingleExpActorModel(tree,branchLength,None))
      (a,ans)
    }.foldLeft[Map[Int,Actor]](IntMap[Actor]()){_+_}

    val pMap = pList.map{p => (p.name,p)}.foldLeft(Map[ParamName,ActorParamComponent]()){_+_}

      
    val modelMap = (map.keys).map{id=> (id,modelMapTmp(map(id)))}.foldLeft[Map[Int,Actor]](IntMap[Actor]()){_+_}
    val components = new BasicActorModel(pi,s,
      new GammaActorModel(alpha,numClasses-1,
        new InvarActorModel(invarPrior,pi,numClasses,
          new ForkActor(tree,modelMap))))
          //modelMap(1))))


     new ActorModel(tree,components,pMap)

   
  }
  def apply[A <: BioEnum](tree:Tree[A],l:List[Int]):ActorModel={
    apply(tree,(tree::tree.descendentNodes).map{n => (n.id,if (l contains n.id){1}else{0})}.foldLeft(Map[Int,Int]()){_+_})
  }
}
class BadParameterException(s:String) extends Exception(s)

trait PSetter[A]{
  def internal:ActorParamComponent
  def update(i:Int,x:Double) { throw new BadParameterException("Can't accept 1D parameter location")}
  def update(i:Int,j:Int,x:Double) { throw new BadParameterException("Can't accept 2D parameter location")}
  def apply(i:Int):Double = throw new BadParameterException("Can't accept 1D parameter location")
  def apply(i:Int,j:Int):Double = throw new BadParameterException("Can't accept 2D parameter location")
  def apply():Double = throw new BadParameterException("Can't accept 0D parameter location")
  def update(x:Double) { throw new BadParameterException("Can't accept 0D parameter location")}
  def <<[B](other:PSetter[B]){apply[B](other.getP)}
  def apply[B](a:B)={internal !? ParamUpdate(a)}
  def getP:A
  def name:ParamName
  override def toString = name + " : " + paramString
  def paramString:String
}
trait OptPSetter {
  def lower:Array[Double]
  def upper:Array[Double]
  def numArguments:Int
  def apply(d:Array[Double])
  def currentArgs:Array[Double]
  def latestArgs:Array[Double]
}



class ActorModel(t:Tree[_],components:ActorModelComponent,val paramMap:Map[ParamName,ActorParamComponent]) extends Logging{
  components.start
  val tree = t.splitAln(4)
  tree.foreach{_.start}

   val params = paramMap.values.toList.removeDuplicates.map{_.name}
  debug{"Param Map " + paramMap}
  def getParam(p:ParamName)={
    p match {
      case SingleParam(p2)=>SingleParamWrapper(paramMap(p2))
      case p => paramMap(p)
    }
  }
  
  def setParamsFrom(other:ActorModel){
    other.paramMap.foreach{t=> val (k,v)=t
      if (paramMap.contains(k) && !k.isInstanceOf[BranchLengths]){
      debug{"Setting " + k}
          this(k)=(v !? RequestParam).asInstanceOf[ParamChanged[_]]
      }
    }
  }
  def <<(other:ActorModel){setParamsFrom(other)}

  val paramLengthMap = paramMap.keys.map{t=>(t,optGet(t).length)}.foldLeft(Map[ParamName,Int]()){(m,p)=>m+((p._1,p._2))}
  def update(p:ParamName,x:Array[Double]) {getParam(p) !? ParamUpdate[Array[Double]](x)}
  def update(p:ParamName,x:Double) {getParam(p) !? ParamUpdate(Array(x))}
  def update(p:ParamName,x:Vector) {getParam(p) !? ParamUpdate(x.toArray)}
  def update(p:ParamName,x:Matrix) {getParam(p) !? ParamUpdate(x)}
  def update[B](p:ParamName,pc:ParamChanged[B]){
    getParam(p) !? pc.toParamUpdate
  }
  def optUpdate(p:ParamName,x:Array[Double]){getParam(p) !? OptUpdate(x)}
  def optGet(p:ParamName):Array[Double]={getParam(p) !? RequestOpt}.asInstanceOf[Array[Double]]
  def optUpdate(pList:List[ParamName])(x:Array[Double]){
    val iter = x.elements
    pList.foreach{p=>
      optUpdate(p,iter.take(paramLengthMap(p)).toList.toArray)
    }
  }
  def optGet(p:List[ParamName]):Array[Double]={
    //this is a crappy implementation
    //but how often is this really going to be called?
    p.map{optGet(_).toList}.flatten[Double].toArray
  }
  def optSetter(p:ParamName):OptPSetter={
    new OptPSetter{
      val param = getParam(p) 
      val lower = (param !? Lower).asInstanceOf[Array[Double]]
      val upper = (param !? Upper).asInstanceOf[Array[Double]]
      val numArguments=lower.length
      val currentArgs=(param !? RequestOpt).asInstanceOf[Array[Double]]
      def latestArgs={
        val latestArgs=(param !? RequestOpt).asInstanceOf[Array[Double]]
        Array.copy(latestArgs,0,currentArgs,0,currentArgs.length)
        latestArgs
      }
      def apply(d:Array[Double])={
        if (d.zip(currentArgs).foldLeft(false){(bool,t)=> bool || t._1!=t._2}){
          Array.copy(d,0,currentArgs,0,d.length)
          param !? OptUpdate(d)
        }
      }
      def paramString=latestArgs.mkString(" ")
      override def toString = p + " : (" +paramString+")"
    }
  }
  def optSetter(p:List[ParamName]):OptPSetter={
    new OptPSetter{
      val numSub=p.length
      val sub = p.map{optSetter}
      val lower = sub.map{_.lower.toList}.flatten[Double].toArray
      val upper = sub.map{_.upper.toList}.flatten[Double].toArray
      val subLengths = sub.map{_.numArguments}
      val numArguments=subLengths.foldLeft(0){_+_}

      val currentArgs = sub.map{_.currentArgs.toList}.flatten[Double].toArray
      def latestArgs={
        val latestArgs = sub.map{_.latestArgs.toList}.flatten[Double].toArray
        Array.copy(latestArgs,0,currentArgs,0,currentArgs.length)
        latestArgs

      }
      def apply(d:Array[Double])={
        var pointer =0
        val subIter = sub.elements
        val lenIter=subLengths.elements
        for (i<- 0 until numSub){
          val len = lenIter.next
          val s =  subIter.next
          s(d.slice(pointer,pointer + len).force)
          pointer = pointer + len
        }
      }
      override def toString = sub.map{_.toString}.mkString(" ")
    }
  }


  type SensibleParam1D = {
    def update(i:Int,d:Double)
    def apply(i:Int):Double
  }
  type SensibleParam2D = {
    def update(i:Int,j:Int,d:Double)
    def apply(i:Int,j:Int):Double
  }

  def paramSetterVector(p:ActorParamComponent)(startParam:Vector)={
    class APSetter extends PSetter[Vector]{
      val internal=p
      def getP={(p !? RequestParam).asInstanceOf[ParamChanged[Vector]].param}
      override def update(i:Int,x:Double){
        val array=getP
        array(i)=x
        p !? ParamUpdate(array)
      }
      override def apply(i:Int):Double={
        val array = getP
        array(i)
      }
      def name = p.name
      def paramString = getP.toString
    }
    
    new APSetter
  }
  //Separate array type necessary otherwise casting back to [A] seems to give a
  //java array rather than scala array - maybe this is something that changes
  //in Scala 2.8?
  def paramSetterArray(p:ActorParamComponent)(startParam:Array[Double])={
    class APSetter extends PSetter[Array[Double]]{
      def getP={(p !? RequestParam).asInstanceOf[ParamChanged[Array[Double]]].param}
      val internal=p
      override def update(i:Int,x:Double){
        val array = getP
        array(i)=x
        p !? ParamUpdate(array)
      }
      override def apply(i:Int):Double={
        val array = getP
        array(i)
      }
      def name = p.name
      def paramString =  getP.mkString(" ") 
    }
    new APSetter
  }

 
  def paramSetterMatrix(p:ActorParamComponent)(startParam:Matrix)={
    class APSetter extends PSetter[Matrix]{
      val internal=p
      def getP=(p !? RequestParam).asInstanceOf[ParamChanged[Matrix]].param
      override def update(i:Int,j:Int,x:Double){
        val array = getP
        array(i,j)=x
        p !? ParamUpdate(array)
      }
      override def apply(i:Int,j:Int):Double={
        getP(i,j)
      }
      def name = p.name
      def paramString = getP.toString
    }
    new APSetter
  }
  def paramSetterDouble(p:ActorParamComponent)(startParam:Double)={
    class APSetter extends PSetter[Double]{
      val internal=p
      def getP=(p !? RequestParam).asInstanceOf[ParamChanged[Double]].param
      override def update(x:Double){
        p !? ParamUpdate(x)
      }
      override def apply():Double={
        getP
      }
      def name = p.name
      def paramString = getP.toString
    }
    new APSetter
  }


  //model(Pi(0))(S)=0.03
  def apply(p:ParamName)={
    val param = getParam(p)
    val start = param !? RequestParam
    start match {
      case ParamChanged(_,a:Array[Double])=>paramSetterArray(param)(a)
      case ParamChanged(_,a:Vector)=>paramSetterVector(param)(a)
      case ParamChanged(_,a:Matrix)=>paramSetterMatrix(param)(a)
      case ParamChanged(_,a:Double)=>paramSetterDouble(param)(a)
      case _ => null//TODO error handling
    }
  }


  def logLikelihood = {val ans = 
    {
    object Send
    Actor.actor{
      receive{
      case Send => 
        for (t<-tree){(t ! LogLikelihoodCalc(components,self))}
        var ans:Double=0.0D
        var returned = 0
        while (returned < tree.length){
          receive{
            case LogLikelihood(d) => 
              returned=returned+1
              ans=ans+d
          }
        }
        reply{ans}
      }
      exit
    } !? Send
  }.asInstanceOf[Double]
  if (ans.isNaN){-1E100}else{ans}
}

  def paramString=paramMap.map{t=> t._2.toString}.mkString("\n")

  override def toString = paramString + "\nlog-likelihood: " + logLikelihood

  def optimise(params:ParamName*):Double={
    optimise(params.toList)
  }
  def optimise(params:List[ParamName]):Double={
    import ModelOptimiser._
    ModelOptimiser.optimise(getConjugateDirection,params,this)
  }
  

  
}

