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

object ParamScheduler extends scala.actors.scheduler.ForkJoinScheduler(10,400,true,false)

trait ParamActor extends Actor{
  override def scheduler = ParamScheduler
}
abstract class ActorModelComponent extends ParamActor with Logging

trait SimpleActorModelComponent extends ActorModelComponent{
  def params:List[ActorParamComponent]
  def rec1:Option[Actor]
  lazy val numParams = params.length
  lazy val paramNames=params.map{_.name}
  lazy val pActorMap=paramNames.zip(params).foldLeft(Map[ParamName,ActorParamComponent]()){_+_}
  var pMap:Map[ParamName,ParamChanged[_]]=null
  override def start={
    debug{"Setup start" + this}
    params.foreach{p=>
      p addActor this
    }
    pMap = params.map{p=>
    var x = (p !? (4000,RequestParam))
    while (x.isEmpty){
       debug{"LOOP"}
       x = (p !? (4000,RequestParam))
    }
    x.get
    }.map{_.asInstanceOf[ParamChanged[_]]}.map{p=>(p.name,p)}.foldLeft(Map[ParamName,ParamChanged[_]]()){_+_}
    pMap.map{_._2}.foreach{updateParam}
    if (rec1.isDefined){rec1.get.start}
    debug{"Setup " + this}
    super.start
  }
  def act{
    loop{
      react{myHandler.orElse { 
          case a:Any => 
            if (rec1.isDefined){
              rec1.get forward a
            }else {
              println(self + " received unexpected message" + a)
            }}
          }
    }
  }

 def matReq[A <: BioEnum](m:MatReq[A]):MatReq[A]
 def updateParam(p:ParamChanged[_])

 def sendUnclean(p:ParamName)={
    if (rec1.isDefined){ rec1.get forward Unclean(pActorMap(p))}
    else {reply(Unclean(pActorMap(p)))}
 }
 /**
   called when a downstream component has a parameter changed
 */
 def unclean:Unit
 def myHandler = handlers
  lazy val handlers = {

    val mainHandler:PartialFunction[Any,Unit]={
      case m:QMatReq[BioEnum]=>
        val ans = matReq(m.toMatReq).toQMatReq
        if (rec1.isDefined){
          rec1.get forward ans
        }
        else {sender ! ans}
      case m:MatReq[BioEnum]=>
        val ans = matReq(m)
        if (rec1.isDefined){rec1.get forward ans}
        else {sender ! ans}
      case Unclean(a)=>
        unclean
        if (rec1.isDefined){rec1.get forward Unclean(a)}
        else {reply(Unclean(a))}
      }
      paramNames.map{p=>
      {
        case ParamChanged(`p`,param)=>
        sendUnclean(p)
        updateParam(ParamChanged(`p`,param))
      }:PartialFunction[Any,Unit]
    }.foldLeft(mainHandler){_ orElse _}
  }
  
  def problem(s:String)={
    warning{s}
  }
  

}

sealed trait NMessage{
  def n:Int
}

case class SwitchingRate(val n:Int) extends NMessage



case class Params(list:List[ActorParamComponent])
object GetParams extends Params(Nil)
class MatBuilder(m:Option[Map[Node[_],Matrix]])
case class Unclean(sender:Actor)
case object RequestOpt
object RequestParam 
sealed trait ParamName
trait ConcreteParamName extends ParamName{
  def i:Int
}
case class Pi(i:Int) extends ConcreteParamName
case class S(i:Int) extends ConcreteParamName
case class BranchLengths(i:Int) extends ConcreteParamName
case class Alpha(i:Int) extends ConcreteParamName
case class InvarPrior(i:Int) extends ConcreteParamName
case class Sigma(i:Int) extends ConcreteParamName
case class SingleParam(p:ParamName) extends ParamName
case class JoinedParam(p:List[ParamName]) extends ParamName


trait ParamMatcher extends ParamName{
  def matches(p2:ParamName):Boolean
}

case class All(p:ConcreteParamName) extends ParamMatcher{
  def matches(p2:ParamName)={
    p2!=null && (p2.getClass==p.getClass) 
  }
}
case class AllMatches(p:ConcreteParamName,f:Int=>Boolean) extends ParamMatcher{
  def matches(p2:ParamName)={
    p2!=null && (p2.getClass==p.getClass) && f(p2.asInstanceOf[ConcreteParamName].i)
  }
}
//This doesn't seem ideal and may be improved...
case class AllSingle(p:ConcreteParamName) extends ParamMatcher{
  def matches(p2:ParamName)={
    p2!=null && (p2.getClass==p.getClass) 
  }
}

object LazyP{
  type PFact={def apply(i:Int):ConcreteParamName}
  implicit def lazymaker(p:PFact)=p(0)
}

case class ParamUpdate[A](d:A)
case class ParamChanged[A](name:ParamName,param:A){
  def toParamUpdate = ParamUpdate(param)
}

object Lower
object Upper


/**
 Optimiser may have different view of params (log of real params, for example)
*/
case class OptUpdate(d:Array[Double])

class BasicActorModel(piParam:ActorPiComponent,sParam:ActorSComponent,rec:Actor) extends ActorModelComponent with SimpleActorModelComponent{
  val params = piParam :: sParam :: Nil
  val piName = piParam.name
  val sName = sParam.name
  val rec1 = Some(rec)
  var pi:Vector = null
  var sMat:Matrix = null
  var mat:Option[Matrix]=None
  def updateParam(p:ParamChanged[_]){
    p match {
      case ParamChanged(`piName`,v:Vector)=>
        pi = v
        unclean
      case ParamChanged(`sName`,m:Matrix) =>
        sMat = m
        unclean
    }
  }
  def unclean{mat=None}
  def matReq[A <: BioEnum](m:MatReq[A])={
    val myMat = if (mat.isDefined){
      mat.get
    }else {
      mat = Some(sMat.sToQ(pi))
        mat.get
      }
    MatReq(m.n,Some(myMat),Some(pi))
  }
}
class NormaliserActorModel(rec:Actor) extends SimpleActorModelComponent{
  val params = Nil
  val rec1 = Some(rec)
  var mat:Option[Matrix]=None
  def updateParam(p:ParamChanged[_]){
  }
  def unclean{mat=None}
  def matReq[A <: BioEnum](m:MatReq[A])={
    if (mat.isDefined){
      MatReq(m.n,mat,m.pi)
    }else {
      mat = Some(m.m.get.normalize(m.pi.get))
      MatReq(m.n,mat,m.pi)
    }
  }
}
   
class THMMActorModel(sigmaParam:ActorFullSComponent,numClasses:Int,rec:Actor) extends SimpleActorModelComponent{
  val sigmaName = sigmaParam.name
  val params = sigmaParam :: Nil
  val rec1 = Some(rec)
  var mat:Option[Matrix]=None
  var switchingRate:Option[Double]=None
  var sigma:Matrix = null
  def updateParam(p:ParamChanged[_]){
    p match {
      case ParamChanged(`sigmaName`,m:Matrix)=>
        sigma = m
        unclean
    }
  }
  def unclean{ mat=None }
  /**
   pi must be the new length, m is pi.length - numAA
  */
  def applyMat(m:Matrix,pi:Vector,sigma:Matrix):(Matrix,Double)={
   val qStart = m.copy
   val numAlpha = m.rows / numClasses
   var switching = 0.0
   for (i <- 0 until numClasses){
        for (j <- i+1 until numClasses){
          for (x <- 0 until numAlpha){
            qStart(i * numAlpha + x, j * numAlpha + x) = sigma(i,j) * pi(j * numAlpha + x) // i->j transition
            qStart(j * numAlpha + x, i * numAlpha + x) = sigma(i,j) * pi(i * numAlpha + x) // j->i transition
            switching = switching +  sigma(i,j) * pi(j * numAlpha + x) *  pi(i * numAlpha + x) * 2
          }
        }
      }
      qStart.fixDiag
      (qStart,switching)
  }

  def matReq[A <: BioEnum](m:MatReq[A])={
    m match {
      case MatReq(n,Some(m),Some(pi))=>
        if (mat.isDefined){mat}
        else {
          val (myMat,switch) = applyMat(m,pi,sigma)
          mat = Some(myMat)
          switchingRate=Some(switch)
        }
        MatReq(n,mat,Some(pi))
      case m:MatReq[_]=>
        problem(m + " is not defined")
        m
    }
  }
  override def myHandler = {
    val h:PartialFunction[Any,Unit]={
      case SwitchingRate(n) =>
        reply(switchingRate)
    }
    h orElse (super.myHandler)
  }

}




class InvarActorModel(priorParam:ActorProbComponent,piParam:ActorPiComponent,numClasses:Int,rec:Actor) extends SimpleActorModelComponent{
  val priorName=priorParam.name
  val piName = piParam.name
  val params = priorParam :: piParam :: Nil
  var rec1=Some(rec)
  var processedPi:Option[Vector]=None
  var mat:Option[Matrix]=None
  var prior:Double = 0.0D
  var pi:Vector = null

  def updateParam(p:ParamChanged[_]){
    p match {
      case ParamChanged(`priorName`,d:Double)=>
        prior = d
        unclean
      case ParamChanged(`piName`,v:Vector)=>
        pi=v
        unclean
    }
  }
  def unclean{
    processedPi = None
    mat = None
  }
  /**
   pi must be the new length, m is pi.length - numAA
  */
  def applyMat(m:Matrix,pi:Vector):Matrix={
    val newMatSize = pi.size
    val mat = Matrix(newMatSize,newMatSize)
    mat.viewPart(0,0,m.rows,m.rows).assign(m)
    mat
  }
  def applyPi(oldPi:Vector):Vector={
    val newPi = Vector(oldPi.size + pi.size)
    newPi.viewPart(0,oldPi.size).assign(oldPi * (1.0D-prior))
    newPi.viewPart(oldPi.size,pi.size).assign(pi * prior)
    newPi
  }

   def matReq[A <: BioEnum](m:MatReq[A])={
     m match{
      case MatReq(n,Some(m),Some(p))=>
        val myPi = if (processedPi.isDefined){
          processedPi.get
        }else{
          processedPi=Some(applyPi(p))
          processedPi.get
        }
        val myMat = if (mat.isDefined){
          mat.get
        }else {
          mat = Some(applyMat(m,myPi))
          mat.get
        }
        MatReq(n,mat,processedPi)
      case m:MatReq[_] =>
        problem(m + " not completely defined")
        m
      }
    }
  }

class GammaActorModel(shape:ActorGammaComponent,numClasses:Int,rec:Actor) extends SimpleActorModelComponent{
  val rec1 = Some(rec)
  val shapeName = shape.name
  val params = shape :: Nil
  val gamma = new Gamma(numClasses)
  var alpha = 1.0
  var mat:Option[Matrix] = None
  var processedPi:Option[Vector] = None
  
  def applyMat(m:Matrix,alpha:Double,pi:Vector):Matrix={
    val r = gamma(alpha)
    debug{"GAMMA " + alpha + " " + r.mkString(" ")}
    val myMat = Matrix(m.rows*numClasses,m.columns * numClasses)
    val numAlpha = m.rows 
    (0 until numClasses).foreach{i=>
      myMat.viewPart(i * numAlpha,i * numAlpha,numAlpha,numAlpha).assign(m) * r(i)
    }
    println(myMat)
    myMat
  }

  def matReq[A <: BioEnum](m:MatReq[A])=m match{
    case MatReq(n,Some(m),Some(p)) =>
      if (mat.isEmpty){
        mat = Some(applyMat(m,alpha,p))
      }
      if (processedPi.isEmpty){
        processedPi = Some(applyPi(p))
      }
      MatReq(n,mat,processedPi)
    case m:MatReq[A]=>
      problem("Unfilled " + m)
      m
  }

  def unclean{
    mat = None
    processedPi = None
  }
  def updateParam(p:ParamChanged[_]){
    p match {
      case ParamChanged(`shapeName`,d:Double) =>
        alpha=d
        unclean
    }
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
        case Unclean(a:ActorParamComponent) =>
          rec.foreach{v=>
            v !? Unclean(self)
          }
          reply(Unclean(a))
        case QMatReq(n,m,p)=>
          if (rec1Map contains n.id){
            rec1Map(n.id) forward QMatReq(n,m,p)
          }else {
            warning{"Map does not contain " + n.id + " " + n}
            sender ! QMatReq(n,m,p) 
          }
        case MatReq(n,m,p)=>
          if(rec1Map contains n.id){
            rec1Map(n.id) forward MatReq(n,m,p)
          }else {
            if (! n.isRoot){warning{"Map does not contain " + n.id + " " + n}}
            sender ! MatReq(n,m,p) 
          }
        case m:NMessage=>
          if (rec1Map contains m.n){
            rec1Map(m.n) forward m
          }else {
            println(self + " received unexpected message " + m)
          }
        case a:Any=> 
          println(self + " received unexpected message " + a)
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
/*
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
  */

}

class BasicSingleExpActorModel[A <: BioEnum](tree:Tree[A],branchLengthParams:ActorTreeComponent[A],rec1:Option[Actor]) extends ActorModelComponent{
  val blName = branchLengthParams.name
  val nodes = tree.descendentNodes.toArray //root has no length

  override def start={
    branchLengthParams addActor this
    if (rec1.isDefined){rec1.get.start}
    super.start
  }
  var eigen:MatrixExponential=null 

  def act{
    main(None,nodes.foldLeft(Map[Int,Double]()){(m,n)=>m+((n.id,n.lengthTo))})
  }
  def main(eigen:Option[MatrixExponential],lengths:Map[Int,Double]){
    case class ExpReq[A <: BioEnum](n:Node[A],e:MatrixExponential,lengthTo:Double,pi:Vector)
      react{
        case q:QMatReq[_] =>
          sender ! q
          main(eigen,lengths)
        case MatReq(n,Some(m),Some(pi)) => 
          if (n.isRoot){//don't need exp(qt)
            if (rec1.isDefined){rec1.get forward MatReq(n,None,Some(pi))}
            else { sender ! MatReq(n,None,Some(pi))}
            main(eigen,lengths)
          }
          val myEigen = if (eigen.isEmpty){
            case object Get
            val actorX = new Actor{
              def act{
                react{
                  case Get => reply(new MatExpNormal(m,pi,None)) // send this calculation out onto a normal actor so Param doesn't do too much computation if user has lowered thread count
                }
              }
            }
            actorX.start
            (actorX !? Get).asInstanceOf[MatExpNormal]//match { case a:MatExpNormal => a}
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
            reply(Unclean(a))
          }
          main(None,lengths)
        case ParamChanged(`blName`,a:Array[Double])=>
          if (rec1.isDefined){rec1.get forward Unclean(branchLengthParams)} else {reply(Unclean(branchLengthParams))}
          val nodeMap = nodes.zip(a).foldLeft(Map[Int,Double]()){(m,t)=>m + ((t._1.id,t._2))}
          main(eigen,nodeMap)
      }
  }
}

trait ActorParamComponent extends ParamActor{
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
      reply(Array(lower))
      main(param)
    case Upper =>
      reply(Array(upper))
      main(param)
    }
  }
}


object JoinedParamWrapper{
  def apply(pList:List[ActorParamComponent])={
    val x = new JoinedParamWrapper(pList)
    x.start
    x
  }
}
class JoinedParamWrapper private(pList:List[ActorParamComponent]) extends ActorParamComponent{
  override def start ={
    super.start
  }
  val name = JoinedParam(pList.map{_.name})
  lazy val paramArray=getParam2D.toArray
  lazy val numParams =paramArray.map{_.size}
    def getParam2D= pList.map{p=>(p !? RequestOpt).asInstanceOf[Array[Double]]}
    def flatten(a:Array[Array[Double]]):Array[Double]={
      flatten(a.toList)
    }
    def flatten(a:List[Array[Double]]):Array[Double]={
      //jeez this is lazy - but it probably isn't going to slow things down too much
      a.map{_.toList}.flatten[Double].toArray 
    }
    def setParam(from:Array[Double],to:Array[Array[Double]])={
      val iter = from.elements
      for (i<-0 until to.length){
        for (j<-0 until to(i).length){
          to(i)(j)=iter.next
        }
      }
    }
  lazy val lower = pList.map{p=>(p !? Lower).asInstanceOf[Array[Double]].toList}.flatten[Double].toArray
  lazy val upper = pList.map{p=>(p !? Upper).asInstanceOf[Array[Double]].toList}.flatten[Double].toArray
  def act{
    loop{
    react{
    case RequestParam=>
      reply(ParamChanged(name,flatten(getParam2D)))
    case ParamUpdate(x:Array[Double])=>
      val to = paramArray.map{i=>new Array[Double](i.size)}.toArray
      setParam(x,to)
      pList.elements.zip(to.elements).foreach{t=>
        val (pActor,pArray)=t
        pActor forward OptUpdate(pArray)
      }
    case OptUpdate(x:Array[Double]) =>
      val to = paramArray.map{i=>new Array[Double](i.size)}
      setParam(x,to)
      pList.elements.zip(to.elements).foreach{t=>
        val (pActor,pArray)=t
        pActor forward OptUpdate(pArray)
      }
    case RequestOpt =>
      reply(flatten(getParam2D))
    case Lower =>
      reply(lower)
    case Upper =>
      reply(upper)
    }}
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

case class ReadSerial(x:Seq[Double])
case object WriteSerial

abstract class AbstractActorParam[A] extends ActorParamComponent with Logging{
  override def !?(msg:Any)={
    debug{this + " !? "+ msg}
    val ans = super.!?(msg)
    debug{this + " finished" }
    ans 
  }

  type Param = {
    def getParams:Array[Double]
    def setParams(a:Array[Double]):Unit
  }
  def internal:Param
  def myParam:A
  override def toString={
    name.toString + " = " + internal.getParams.toList.mkString(",")
  }

  def readSerial(x:Seq[Double])=internal setParams x.toArray
  def writeSerial:Seq[Double]=internal.getParams

  def lower:Double
  def upper:Double
  def lowerArray = Array.make(internal.getParams.length,lower)
  def upperArray = Array.make(internal.getParams.length,upper)
  //override this if raw params are different!
  def setRaw(p:Array[Double]){internal setParams p}
  /*def waitOn(i:Int,o:OutputChannel[Any]){
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
  }*/
  private val myHandler:PartialFunction[Any,Unit] = {
        case RequestParam=>
         reply(ParamChanged(name,myParam))
        case ParamUpdate(x:Array[Double])=> 
          setRaw(x)
          modelComp.foreach{c=>
            c !? ParamChanged(name,myParam)
          }
          reply('ok)
        //  waitOn(modelComp.length,sender)
        case ParamUpdate(v:Vector)=>
          debug{"ParamUpdate 1"}
          setRaw(v.toArray)
          debug{"ParamUpdate 2"}
          modelComp.foreach{c=>
            debug{"ParamUpdate c"}
            c !? ParamChanged(name,myParam)
          }
          debug{"ParamUpdate 3"}
          reply('ok)
        case OptUpdate(x)=>
          internal setParams x
          modelComp.foreach{c=>
            c !? ParamChanged(name,myParam)
          }
          reply('ok)
        case RequestOpt=>
          reply(internal.getParams)
        case Lower=>
          reply(lowerArray)
        case Upper=>
          reply(upperArray)
        case ReadSerial(s) =>
          readSerial(s)
          modelComp.foreach{c=>
            c !? ParamChanged(name,myParam)
          }
          reply('ok)
        case WriteSerial =>
          reply(writeSerial)
  }
  def handler:PartialFunction[Any,Unit]= {
     myHandler 
  }

  def act{
    loop{
      react{handler}
    }
  }
}
class ActorPiComponent(pi:Vector,val name:ParamName) extends AbstractActorParam[Vector]{

  override def handler:PartialFunction[Any,Unit]={
    super.handler
  }


  class PiParam(pi:Vector){

    var medianIndex=0
    def fromFit(array:Array[Double],medianIndex:Int)={
      val exponentiated =  array.map{i=>Math.exp(i)}
      val total = (0.0D /: exponentiated){_+_} + Math.exp(0.0D)
      ((0 to medianIndex-1 ).map{i=> exponentiated(i)/total}.toList ++ List(Math.exp(0.0D)/total) ++ (medianIndex to array.length-1).map{i=> exponentiated(i)/total}).toArray
    }
  
    def toFit(pi:Vector):(List[Double],Int)={
      val t  = (pi.toList.zipWithIndex.toList.sort{_._1<_._1})(pi.size/2)
      medianIndex = t._2
      def toFitness={i:Int=>Math.log(pi(i)/pi(medianIndex))}
      ((0 to medianIndex-1).map{toFitness}.toList ++ (medianIndex+1 to pi.size-1).map{toFitness}.toList,
       medianIndex)
    }
    
    def setParams(a:Array[Double]){pi assign fromFit(a,medianIndex)}
    def getParams={
      val (l,m)=toFit(pi)
      medianIndex=m
      l.toArray
    }
    def view=pi.copy
    def setPi(p:Array[Double]){pi assign p}

  
  }
  val internal = new PiParam(pi)

  def lower = -10 // this is log-odds scaled prob
  def upper = 10
  def myParam = pi
  override def setRaw(p:Array[Double]){
    internal setPi p
  }
  override def toString=name.toString + " = " + internal.view.toArray.mkString(",")
  override def writeSerial = internal.view.toList
  override def readSerial(x:Seq[Double])=internal.setPi(x.toArray)
}
class ActorSComponent(s:Matrix,val name:ParamName) extends AbstractActorParam[Matrix] with SMatUtil{
  class SMatParam(s:Matrix) extends SMatUtil{
    var cache:Option[Array[Double]]=None
    def getParams={ if (cache.isEmpty){cache = Some(linearSMat(s).toArray)}; cache.get}
    def setParams(a:Array[Double]) = setSMat(a,s)
  }
  val internal = new SMatParam(s)
  def myParam = s
  def lower = 0.0
  def upper = 200
  private val myHandler:PartialFunction[Any,Unit]={
      case ParamUpdate(m:Matrix)=>
        setRaw(linearSMatFull(m).toArray)
        modelComp.foreach{c=>
          c !? ParamChanged(name,myParam)
        }
        //waitOn(modelComp.length,sender)
        reply('ok)
  }
  override def handler={
    myHandler orElse (super.handler)
  }
}

class ActorFullSComponent(s:Matrix,val name:ParamName) extends AbstractActorParam[Matrix] with SMatUtil{
  class FullSMatParam(s:Matrix) extends SMatUtil{
    var cache:Option[Array[Double]]=None
    def getParams={ if (cache.isEmpty){cache = Some(linearSMatFull(s).toArray)}; cache.get}
    def setParams(a:Array[Double]) = {setSMat(a,s);cache=None}
  }
  val internal = new FullSMatParam(s)
  def lower = 0.0D
  def upper = 200
  def myParam = s
  private val myHandler:PartialFunction[Any,Unit]={
      case RequestParam=>
          reply(ParamChanged(name,myParam))
      case ParamUpdate(m:Matrix)=>
        setRaw(linearSMatFull(m).toArray)
        modelComp.foreach{c=>
          c !? ParamChanged(name,myParam)
        }
        //waitOn(modelComp.length,sender)
        reply('ok)
  }
  override def handler={
    myHandler orElse (super.handler)
  }
}


class ActorTreeComponent[B <: BioEnum](tree:Tree[B],name:ParamName) extends ActorArrayComponent(tree.getBranchLengths.toArray,name,0.0,100.0)
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
        c !? ParamChanged(name,myParam)
      }
     // waitOn(modelComp.length,sender)
      reply('ok)
  }
  override def handler={
    myHandler orElse (super.handler)
  }
}
class ActorArrayComponent(param:Array[Double],val name:ParamName,val lower:Double, val upper:Double) extends AbstractActorParam[Array[Double]]{
  class ArrayParam{
    var myP=param
    def getParams = myP.clone
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
    val components = new BasicActorModel(pi,s, new NormaliserActorModel(new BasicSingleExpActorModel(tree,branchLength,None)))
    val pList = List(pi,s,branchLength)
    val pMap = pList.map{p => (p.name,p)}.foldLeft(Map[ParamName,ActorParamComponent]()){_+_}
    new ActorModel(tree,components,pMap)
  }
}
object GammaModel{
  def apply[A <: BioEnum](tree:Tree[A]):ActorModel[A]=apply(tree,WAG.pi,WAG.S)
  def apply[A <: BioEnum](tree:Tree[A],piVect:Vector,sMat:Matrix):ActorModel[A]={
    val pi = new ActorPiComponent(piVect,Pi(0))
    val s = new ActorSComponent(sMat,S(0))
    val branchLength = new ActorTreeComponent(tree,BranchLengths(0))
    val alpha = new ActorGammaComponent(0.5D,Alpha(0))
    val components = new BasicActorModel(pi,s, 
      new GammaActorModel(alpha,tree.alphabet.numClasses,
       new NormaliserActorModel(
       new BasicSingleExpActorModel(tree,branchLength,None))))
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
         new NormaliserActorModel(
          new BasicSingleExpActorModel(tree,branchLength,None)))))
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
         new NormaliserActorModel(
          new THMMActorModel(sigma,numClasses,
            new BasicSingleExpActorModel(tree,branchLength,None))))))
    val pList = List(pi,s,branchLength,alpha,invarPrior,sigma)
    val pMap = pList.map{p => (p.name,p)}.foldLeft(Map[ParamName,ActorParamComponent]()){_+_}
     new ActorModel(tree,components,pMap)
  }
}

object BranchSpecificThmmModel{
  def apply[A <: BioEnum](tree:Tree[A]):ActorModel[A]={
    val nodeMap=tree.descendentNodes.foldLeft[Map[Int,Int]](IntMap[Int]()){(m,n)=>
        m + ((n.id,n.id))
    }
    apply(tree,nodeMap)
  }
  def apply[A <: BioEnum](tree:Tree[A],map:Map[Int,Int]):ActorModel[A]={
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
          new NormaliserActorModel(
          new ForkActor(tree,modelMap)))))
          //modelMap(1))))


    new ActorModel(tree,components,pMap)

   
  }
  def apply[A <: BioEnum](tree:Tree[A],l:List[Int]):ActorModel[A]={
    apply(tree,(tree::tree.descendentNodes).map{n => (n.id,if (l contains n.id){1}else{0})}.foldLeft(Map[Int,Int]()){_+_})
  }
}
class BadParameterException(s:String) extends Exception(s)

trait PSetter[A]{
  def internal:ActorParamComponent
  def update(x:Seq[Double]){
    x.elements.zipWithIndex.foreach{t=>
      update(t._2,t._1)
    }
  }

  def readSerial(x:Seq[Double]){internal !? ReadSerial(x)}
  def writeSerial = {internal !? WriteSerial}.asInstanceOf[Seq[Double]]
  def update(i:Int,x:Double) { throw new BadParameterException("Can't accept 1D parameter location")}
  def update(i:Int,j:Int,x:Double) { throw new BadParameterException("Can't accept 2D parameter location")}
  def apply(i:Int):Double = throw new BadParameterException("Can't accept 1D parameter location")
  def apply(i:Int,j:Int):Double = throw new BadParameterException("Can't accept 2D parameter location")
  def apply():Double = throw new BadParameterException("Can't accept 0D parameter location")
  def toList:List[Double]
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
  def length = numArguments
  def apply(d:Array[Double])
  def currentArgs:Array[Double]
  def latestArgs:Array[Double]
  lazy val random = new org.apache.commons.math.random.MersenneTwister
  def randomise{
    apply(
      lower.zip(upper).map{t=>
        random.nextDouble * (t._2-t._1) + t._1
      }.toArray
    )
  }
  def randomise(lowerX:Double,upperX:Double){
    apply(
      lower.zip(upper).map{t=> ((t._1 max lowerX),(t._2 min upperX))}.map{t=>
        random.nextDouble * (t._2-t._1) + t._1
      }.toArray
    )
  }

  def shuffle{
    val prevArgs = currentArgs
    val newArgs = prevArgs.map{i=>
      prevArgs(random nextInt currentArgs.length)
    }
    apply(newArgs)
  }
}



class ActorModel[A <: BioEnum](t:Tree[A],components:ActorModelComponent,val paramMap:Map[ParamName,ActorParamComponent]) extends Logging{

  val tree = t.splitAln(6)
  tree.foreach{_.start}
  val params = paramMap.keys.toList
  paramMap.values.foreach{_.start}
  components.start


  def rawReq(a:Any)={
    Actor.actor{
      receive{
        case a:Any=>
        components ! a
        var a2:Any=null
        receive{
          case b:Any=>
            a2=b
          }
        reply(a2)
        exit
        }
      } !? a
  }

  def qMat(i:Int)={
    case class GetQ(i:Int)


    val actor = Actor.actor{
      receive{
        case GetQ(i)=>
        val baseTree = tree(0)
        components ! QMatReq(baseTree(i),None,None)
        var matrix:Matrix=null
        receive{
          case QMatReq(node,Some(m),Some(pi))=>
            matrix = m
          }
        reply(matrix)
        exit
        }
      }
    (actor !? GetQ(i)).asInstanceOf[Matrix]
    }

  def eMat(i:Int)={
    case class GetQ(i:Int)


    val actor = Actor.actor{
      receive{
        case GetQ(i)=>
        val baseTree = tree(0)
        components ! MatReq(baseTree(i),None,None)
        var matrix:Matrix=null
        receive{
          case MatReq(node,Some(m),Some(pi))=>
            matrix = m
          }
        reply(matrix)
        exit
        }
      }
    (actor !? GetQ(i)).asInstanceOf[Matrix]
    }

  debug{"Param Map " + paramMap}
  def getParam(p:ParamName):ActorParamComponent={
    val ans = p match {
      case SingleParam(p2)=>SingleParamWrapper(getParam(p2))
      case JoinedParam(p2)=>JoinedParamWrapper(p2.map{getParam})
      case m:AllSingle => {
       JoinedParamWrapper( paramMap.filter{m matches _._1}.map{t=>SingleParamWrapper(t._2)}.toList)
     }
      case m:ParamMatcher => {
        JoinedParamWrapper( paramMap.filter{m matches _._1}.map{_._2}.toList)
      }
      case p  => paramMap(p)
    }
    ans
  }
  
  def setParamsFrom(other:ActorModel[A]){
    other.paramMap.foreach{t=> val (k,v)=t
      if (paramMap.contains(k) && !k.isInstanceOf[BranchLengths]){
        debug{"Setting " + k + " " + v}
        val p = (v !? RequestParam).asInstanceOf[ParamChanged[_]]
        debug{"Param " + p}
        this(k)=p
      }
    }
  }

  def from(other:ActorModel[A])={
    this << other
    apply(BranchLengths(0)) << other(BranchLengths(0))
    this
  }


  def from(other:String)={
    val pLookup = paramMap.map{t=>(t._1.toString,t._1)}.foldLeft(Map[String,ParamName]()){_+_}
    other.lines.foreach{s=>
      val p = s.split("=").map{_.trim}
      val name=pLookup.get(p(0))
      if (name.isDefined){
        apply(name.get)  readSerial p(1).split(",").map(_.toDouble)
      }else {
        warning{"Ignoring apparent Parameter " + p(0)}
      }
    }
    this
  }

  def <<(other:ActorModel[A]){setParamsFrom(other)}

  val paramLengthMap = paramMap.keys.map{t=>(t,optGet(t).length)}.foldLeft(Map[ParamName,Int]()){(m,p)=>m+((p._1,p._2))}

  def update(p:ParamName,x:Array[Double]) {
    getParam(p) !? ParamUpdate[Array[Double]](x)
  }
  def update(p:ParamName,x:Double) {
    getParam(p) !? ParamUpdate(Array(x))
  }
  def update(p:ParamName,x:Vector) {
    getParam(p) !? ParamUpdate(x.toArray)
  }
  def update(p:ParamName,x:Matrix) {
    getParam(p) !? ParamUpdate(x)
  }
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
      lazy val lower = (param !? Lower).asInstanceOf[Array[Double]]
      lazy val upper = (param !? Upper).asInstanceOf[Array[Double]]
      lazy val numArguments=lower.length
      lazy val currentArgs=(param !? RequestOpt).asInstanceOf[Array[Double]]
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
          s(d.slice(pointer,pointer + len))
          pointer = pointer + len
        }
      }
      override def toString = sub.map{_.toString}.mkString(" ")
    }
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
      def toList = getP.toList
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
      def toList = getP.toList
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
        val ans = getP(i,j)
        ans
      }
      def toList = getP.toList
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
        val ans = getP
        ans
      }
      override def update(i:Int,x:Double)={
        assert(i==0)
        update(x)
      }
      def toList = getP :: Nil
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
        for (t<-tree){
          t ! LogLikelihoodCalc(components,self)
        }
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

  def paramString=paramMap.map{t=> t._2.name + " = " + (t._2 !? WriteSerial).asInstanceOf[Seq[Double]].mkString(",")}.mkString("\n")

  override def toString = paramString + "\nlog-likelihood: " + logLikelihood

  def optimise(params:ParamName*):Double={
    optimise(params.toList)
  }
  def optimise(params:List[ParamName],tolfx:Double,tolx:Double):Double={
    import ModelOptimiser._
    ModelOptimiser.optimise(getConjugateDirection,params,this,tolfx,tolx)
  }
  def optimiseQuick(params:List[ParamName]):Double=optimise(params,1E-1,1E-1)
  def optimise(params:List[ParamName]):Double={
    import ModelOptimiser._
    ModelOptimiser.optimise(getConjugateDirection,params,this)
  }
  def optimise(params:List[List[ParamName]])(implicit m:Manifest[List[List[ParamName]]]):Double={
    var end = logLikelihood
    var start = end
    do {
      start=end
      params.foreach{pset=>
        optimise(pset) 
      }
      end = logLikelihood
    } while (end - start > 0.001)
    end
  }
  def optimiseCustom(params:List[List[ParamName]],tolfx:Double,tolx:Double)(implicit m:Manifest[List[List[ParamName]]]):Double={
  var end = logLikelihood
    var start = end
    do {
      start=end
      params.foreach{pset=>
        optimise(pset,tolfx,tolx) 
      }
      end = logLikelihood
    } while (end - start > 0.001)
    end
  }
  def optimiseQuick(params:List[List[ParamName]])(implicit m:Manifest[List[List[ParamName]]]):Double={
    var end = logLikelihood
    var start = end
    do {
      start=end
      params.foreach{pset=>
        optimiseQuick(pset) 
      }
      end = logLikelihood
    } while (end - start > 0.001)
    end
  }


//  def switchingRate(nodeID:Int)={
//    logLikelihood //make sure parameters are set
//    (components !? SwitchingRate(nodeID)).asInstanceOf[Option[Double]].get
//  }

  def switchingRate(nodeID:Int)={
    val myTree=tree(0)
    val ans = rawReq(QMatReq(myTree(nodeID),None,None)).asInstanceOf[QMatReq[_]]
    val qMat = ans.m.get
    val pi = ans.pi.get
    Math.max(qMat.rate(pi),1.0)-1.0
  }

  def switchingBranchLengths={
    val myTree = currentTree
    val branchLengths=apply(BranchLengths(0)).toList
    myTree.descendentNodes.map{i:Node[_]=>i.id}.map{i=> (i,switchingRate(i) * myTree(i).lengthTo)}.foldLeft(Map[Int,Double]()){_+_}
  }
  def switchingTree={
    t setBranchLengths switchingBranchLengths
  }
  def currentTree={
    t setBranchLengths (apply(BranchLengths(0)).toList)
  }
  
}
case class MatReq[A <: BioEnum](n:Node[A],m:Option[Matrix],pi:Option[Vector]){
  def toQMatReq=QMatReq(n,m,pi)
}
case class QMatReq[A <: BioEnum](n:Node[A],m:Option[Matrix],pi:Option[Vector]){// extends MatReq(n2,m2,pi2) //used for skipping exponentiation
  def toMatReq=MatReq(n,m,pi)
}
object NewMatReq{
  def apply[A <: BioEnum](n:Node[A])=MatReq(n,None,None)
}



