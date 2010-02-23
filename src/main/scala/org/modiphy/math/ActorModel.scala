package org.modiphy.math
import EnhancedMatrix._
import org.modiphy.tree._
import org.modiphy.tree.DataParse.Tree
import org.modiphy.sequence._
import scala.collection.immutable.{Map,IntMap}
import scala.actors.Actor
import scala.actors.Actor._
import scala.actors.OutputChannel

abstract class ActorModelComponent extends Actor{
}
case class Params(list:List[ActorParamComponent])
object GetParams extends Params(Nil)
class MatBuilder(m:Option[Map[Node[_],Matrix]])
case class MatReq(n:Node[_],m:Option[Matrix],pi:Option[Vector])
object NewMatReq{
  def apply(n:Node[_])=MatReq(n,None,None)
}
case object Unclean
case object RequestOpt
object RequestParam 
abstract class ParamName

case class Pi(i:Int) extends ParamName
object Pi extends Pi(0)
case class S(i:Int) extends ParamName
object S  extends S(0)
case class BranchLengths(i:Int) extends ParamName
object BranchLengths extends BranchLengths(0)
case class Alpha(i:Int) extends ParamName
object Alpha extends Alpha(0)
case class InvarPrior(i:Int) extends ParamName
object InvarPrior extends InvarPrior(0)
case class Sigma(i:Int) extends ParamName
object Sigma extends Sigma(0)

case class SingleParam(p:ParamName) extends ParamName

case class ParamUpdate[A](d:A)
case class ParamChanged[A](p:ParamName,param:A)

object Lower
object Upper


/**
 Optimiser may have different view of params (log of real params, for example)
*/
case class OptUpdate(d:Array[Double])

class BasicActorModel(piParam:ActorPiComponent,sParam:ActorSComponent,reciever:Actor) extends ActorModelComponent{
  val piName = piParam.name
  val sName = sParam.name
  override def start={
    piParam.start
    sParam.start
    piParam addActor this
    sParam addActor this
    reciever.start
    super.start
  }
  def act{
    piParam ! RequestParam
    sParam ! RequestParam
    initialise(None,None)
  }
  def initialise(s:Option[Matrix],pi:Option[Vector]){
    if (s.isDefined && pi.isDefined){
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
        reciever forward Params(piParam :: sParam :: l)
        main(s,pi,mat)
      case MatReq(n,oldM,oldPi)=>
        val myMat = if (mat.isDefined){
          mat.get
        }else {
          s.sToQ(pi)
        }
        reciever forward MatReq(n,Some(myMat),Some(pi))
        main(s,pi,Some(myMat))
        case ParamChanged(piName,v:Vector)=>
          reciever forward Unclean
          main(s,v.copy,None)
        case ParamChanged(sName,m:Matrix)=>
          reciever forward Unclean
          main(m.copy,pi,None)
     }
  }
}

class THMMActorModel(sigmaParam:ActorFullSComponent,numClasses:Int,reciever:Actor) extends ActorModelComponent{
  val sigmaName = sigmaParam.name
  override def start={
    sigmaParam addActor this
    sigmaParam.start
    reciever.start
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
    sigmaParam ! RequestParam
    initialise(None)
  }
  def initialise(sigma:Option[Matrix]){
    react{
      case ParamChanged(sigmaName,a:Matrix)=>main(a,None)
    }
  }
  def main(sigma:Matrix,mat:Option[Matrix]){
    react{
      case Params(l)=> 
        reciever forward Params(sigmaParam :: l)
        main(sigma,mat)
      case ParamChanged(sigmaName,array:Matrix)=>
        reciever forward Unclean
        main(array,None)
      case MatReq(n,Some(m),Some(p))=>
        val myMat = if (mat.isDefined){
          mat.get
        }else {
          applyMat(m,p,sigma)
        }
        reciever forward MatReq(n,Some(myMat),Some(p))
        main(sigma,Some(myMat))
      case Unclean => 
        reciever forward Unclean 
        main(sigma,None)
    }
  }
}




class InvarActorModel(priorParam:ActorProbComponent,piParam:ActorPiComponent,numClasses:Int,reciever:Actor) extends ActorModelComponent{
  val priorName=priorParam.name
  val piName = piParam.name
  override def start={
    priorParam.start
    priorParam addActor this
    piParam.start
    piParam addActor this
    reciever.start
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
    piParam ! RequestParam
    priorParam ! RequestParam
    initialise(None,None)
  }
  def initialise(prior:Option[Double],pi:Option[Vector]){
    if (prior.isDefined && pi.isDefined){
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
        reciever forward Params(priorParam :: piParam :: l)
        main(prior,rawpi,processedPi,mat)
      case ParamChanged(priorName,d:Double)=>
        reciever forward Unclean
        main(d,rawpi,None,None)
      case ParamChanged(piName,v:Vector)=>
        reciever forward Unclean
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
        reciever forward MatReq(n,Some(myMat),Some(myPi))
        main(prior,rawpi,Some(myPi),Some(myMat))
      case Unclean => 
        reciever forward Unclean 
        main(prior,rawpi,None,None)
    }
  }
}

class GammaActorModel(shape:ActorGammaComponent,numClasses:Int,reciever:Actor) extends ActorModelComponent{
  val shapeName = shape.name
  override def start={
    shape.start
    shape addActor this
    reciever.start
    super.start
  }
  val gamma = new Gamma(numClasses)
  
  def applyMat(m:Matrix,alpha:Double,pi:Vector):Matrix={
    val r = gamma(alpha)
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
    shape ! RequestParam
    react{
      case ParamChanged(shapeName,alpha:Double)=>main(alpha,None,None)
    }
  }
  def main(alpha:Double,mat:Option[Matrix],pi:Option[Vector]){
    react{
      case Params(l)=>
        reciever forward Params(shape :: l)
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
        reciever forward MatReq(n,Some(myMat),Some(myPi))
        main(alpha,Some(myMat),Some(myPi))
      case ParamChanged(shapeName,d:Double)=>
        reciever forward Unclean
        main(d,None,pi)
      case Unclean =>
        reciever forward Unclean
        main(alpha,None,None)
    }
  }
}

class ForkActor[A <: BioEnum](tree:Tree[A],receiverMap:Map[Node[A],Actor]) extends ActorModelComponent{
  val numRec = receiverMap.values.toList.length
  override def start={
    receiverMap.values.foreach{v=>
      v.start
    }
    super.start
  }

  def act{
    loop{
      react{
        case Params(l)=>
          receiverMap.values.foreach{v=>
            v ! Params(l)             
          }
          getListReplies(sender,numRec,Nil)
        case Unclean =>
          receiverMap.values.foreach{v=>
            v ! Unclean
          }
          getCleanReplies(sender,numRec)
        case MatReq(n,m,p)=>
          if(receiverMap contains n){
            receiverMap(n) forward MatReq(n,m,p)
          }else {
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

  def getCleanReplies(output:OutputChannel[Any],toDo:Int){
    if (toDo==0){
      output ! Unclean
      act
    }
    react{
      case Unclean => getCleanReplies(output,toDo-1)
    }
  }

}

class BasicSingleExpActorModel[A <: BioEnum](tree:Tree[A],branchLengthParams:ActorTreeComponent[A],reciever:Option[Actor]) extends ActorModelComponent{
  val blName = branchLengthParams.name
  val nodes = tree.descendentNodes.toArray //root has no length

  override def start={
    branchLengthParams addActor this
    branchLengthParams.start
    if (reciever.isDefined){reciever.get.start}
    super.start
  }
  var eigen:MatrixExponential=null 

  def act{
    main(None,nodes.foldLeft(Map[Node[A],Double]()){(m,n)=>m+((n,n.lengthTo))})
  }
  def main(eigen:Option[MatrixExponential],lengths:Map[Node[A],Double]){
    case class ExpReq(n:Node[_],e:MatrixExponential,lengthTo:Double,pi:Vector)
      react{
        case Params(l)=>
          if (reciever.isDefined){reciever.get forward Params(branchLengthParams :: l)}
          else {sender ! Params(branchLengthParams :: l)}
          main(eigen,lengths)
        case MatReq(n,Some(m),Some(pi)) => 
          if (n.isRoot){//don't need exp(qt)
            sender ! MatReq(n,None,Some(pi))
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
                  sender ! MatReq(n,Some(ans),Some(pi)) 
                  exit
              }
            }
          }.start forward ExpReq(n,myEigen,lengths(n),pi)
          main(Some(myEigen),lengths)
        case Unclean=>
          if (reciever.isDefined){reciever.get forward Unclean} else { reply(Unclean)}
          main(None,lengths)
        case ParamChanged(blName,a:Array[Double])=>
          if (reciever.isDefined){reciever.get forward Unclean} else {reply(Unclean)}
          val nodeMap = nodes.zip(a).foldLeft(Map[Node[A],Double]()){(m,t)=>m + ((t._1,t._2))}
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
      sender ! ParamChanged(name,param)
    case ParamUpdate(x:Array[Double])=>
      val newParam = Array.make(param.length,x(0))
      p forward OptUpdate(newParam)
      main(newParam)
    case OptUpdate(x:Array[Double]) =>
      val newParam = Array.make(param.length,x(0))
      p forward OptUpdate(newParam)
      main(newParam)
    case RequestOpt =>
      sender ! Array(param(0))
    }
  }
}

abstract class AbstractActorParam[A] extends ActorParamComponent{
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
  private val myHandler:PartialFunction[Any,Unit] = {
        case RequestParam=>
          sender ! ParamChanged[A](name,myParam)
        case ParamUpdate(x:Array[Double])=> 
          setRaw(x)
          modelComp.foreach{c=>
            c !? ParamChanged[A](name,myParam)
          }
          reply('ok)
        case OptUpdate(x)=>
          internal setParams x
          modelComp.foreach{c=>
            c !? ParamChanged[A](name,myParam)
          }
          reply('ok)
        case RequestOpt=>
          sender ! internal.getParams
        case Lower=>
          sender ! lowerArray
        case Upper=>
          sender ! upperArray
  }
  def handler:PartialFunction[Any,Unit]= {
     myHandler 
  }

  def act{
    loop{
      react(handler)
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
        setRaw(linearSMat(m).toArray)
        modelComp.foreach{c=>
          c !? ParamChanged(name,myParam)
        }
        reply('ok)
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
      case ParamUpdate(m:Matrix)=>
        setRaw(linearSMatFull(m).toArray)
        modelComp.foreach{c=>
          c !? ParamChanged(name,myParam)
        }
        reply('ok)
  }
  override def handler={
    myHandler orElse (super.handler)
  }
}


class ActorTreeComponent[B <: BioEnum](tree:Tree[B],val name:ParamName) extends AbstractActorParam[Array[Double]]{
  def lower =0.0
  def upper = 100.0
  class TreeParam{
    val nodes = tree.descendentNodes.toArray
    var myP=nodes.map{n:Node[B]=>n.lengthTo}
    def getParams=myP
    def setParams(a:Array[Double]){Array.copy(a,0,myP,0,Math.min(a.length,myP.length))}
  }
  val internal = new TreeParam
  def myParam = internal.getParams
  private val myHandler:PartialFunction[Any,Unit]={
        case RequestParam=>
          sender ! ParamChanged(name,myParam)
        case ParamUpdate(x:Array[Double])=> 
          setRaw(x)
          modelComp.foreach{c=>
            c !? ParamChanged(name,myParam)
          }
          reply('ok)
  }
  override def handler={
    myHandler orElse (super.handler)
  }
}
class ActorDoubleComponent(param:Double,val name:ParamName,val lower:Double,val upper:Double) extends AbstractActorParam[Double]{
  class DoubleParam{
    var myP:Double=param
    def getParams=Array(myP)
    def setParams(a:Array[Double]){myP=a(0)}
  }
  val internal = new DoubleParam
  def myParam = internal.getParams(0)
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
class ActorProbComponent(prob:Double,name:ParamName) extends ActorDoubleComponent(prob,name,0,1)
class ActorGammaComponent(alpha:Double,name:ParamName) extends ActorDoubleComponent(alpha,name,0.01,1000)
object SimpleModel{
  def apply[A <: BioEnum](tree:Tree[A])={
    val pi = new ActorPiComponent(WAG.pi,Pi)
    val s = new ActorSComponent(WAG.S,S)
    val branchLength = new ActorTreeComponent(tree,BranchLengths)
    val components = new BasicActorModel(pi,s, new BasicSingleExpActorModel(tree,branchLength,None))
    new ActorModel(tree,components)
  }
}
object GammaModel{
  def apply[A <: BioEnum](tree:Tree[A])={
    val pi = new ActorPiComponent(WAG.pi,Pi)
    val s = new ActorSComponent(WAG.S,S)
    val branchLength = new ActorTreeComponent(tree,BranchLengths)
    val alpha = new ActorGammaComponent(0.5D,Alpha)
    val components = new BasicActorModel(pi,s, 
      new GammaActorModel(alpha,tree.alphabet.numClasses,
       new BasicSingleExpActorModel(tree,branchLength,None)))
    new ActorModel(tree,components)
  }
}
object InvarGammaModel{
  def apply[A <: BioEnum](tree:Tree[A])={
    val pi = new ActorPiComponent(WAG.pi,Pi)
    val s = new ActorSComponent(WAG.S,S)
    val branchLength = new ActorTreeComponent(tree,BranchLengths)
    val alpha = new ActorGammaComponent(0.5D,Alpha)
    val invarPrior = new ActorProbComponent(0.2D,InvarPrior)
    val components = new BasicActorModel(pi,s,
      new GammaActorModel(alpha,tree.alphabet.numClasses-1,
        new InvarActorModel(invarPrior,pi,tree.alphabet.numClasses,
          new BasicSingleExpActorModel(tree,branchLength,None))))
     new ActorModel(tree,components)
  }
}

object InvarThmmModel{
def apply[A <: BioEnum](tree:Tree[A])={
  val numClasses = tree.alphabet.numClasses
    val pi = new ActorPiComponent(WAG.pi,Pi)
    val s = new ActorSComponent(WAG.S,S)
    val branchLength = new ActorTreeComponent(tree,BranchLengths)
    val alpha = new ActorGammaComponent(0.5D,Alpha)
    val invarPrior = new ActorProbComponent(0.2D,InvarPrior)
    val sigma = new ActorFullSComponent(Matrix(numClasses,numClasses),Sigma)
    val components = new BasicActorModel(pi,s,
      new GammaActorModel(alpha,numClasses-1,
        new InvarActorModel(invarPrior,pi,numClasses,
          new THMMActorModel(sigma,numClasses,
            new BasicSingleExpActorModel(tree,branchLength,None)))))
     new ActorModel(tree,components)
  }
}

object BranchSpecificThmmModel{
  def apply[A <: BioEnum](tree:Tree[A]):ActorModel={
    var i = -1
    val nodeMap=tree.descendentNodes.foldLeft(Map[Node[A],Int]()){(m,n)=>
        i=i+1
        m + ((n,i))
    }
    apply(tree,nodeMap)
  }
  def apply[A <: BioEnum](tree:Tree[A],map:Map[Node[A],Int]):ActorModel={
    val numClasses = tree.alphabet.numClasses
    val pi = new ActorPiComponent(WAG.pi,Pi)
    val s = new ActorSComponent(WAG.S,S)
    val branchLength = new ActorTreeComponent(tree,BranchLengths)
    val alpha = new ActorGammaComponent(0.5D,Alpha)
    val invarPrior = new ActorProbComponent(0.2D,InvarPrior)
    val modelNums = map.values.toList.removeDuplicates
    val modelMapTmp = modelNums.map{a=>
      val sigma = new ActorFullSComponent(Matrix(numClasses,numClasses),Sigma(a))
      val ans = new THMMActorModel(sigma,numClasses, new BasicSingleExpActorModel(tree,branchLength,None))
      (a,ans)
    }.foldLeft[Map[Int,Actor]](IntMap[Actor]()){_+_}
      
    val modelMap = (map.keys).map{n=> (n,modelMapTmp(map(n)))}.foldLeft(Map[Node[A],Actor]()){_+_}
    val components = new BasicActorModel(pi,s,
      new GammaActorModel(alpha,numClasses-1,
        new InvarActorModel(invarPrior,pi,numClasses,
          new ForkActor(tree,modelMap))))


     new ActorModel(tree,components)

   
  }
}
class BadParameterException(s:String) extends Exception(s)

trait PSetter{
  def update(i:Int,x:Double) { throw new BadParameterException("Can't accept 1D parameter location")}
  def update(i:Int,j:Int,x:Double) { throw new BadParameterException("Can't accept 2D parameter location")}
  def apply(i:Int):Double = throw new BadParameterException("Can't accept 1D parameter location")
  def apply(i:Int,j:Int):Double = throw new BadParameterException("Can't accept 2D parameter location")
  def apply():Double = throw new BadParameterException("Can't accept 0D parameter location")
  def update(x:Double) { throw new BadParameterException("Can't accept 0D parameter location")}
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

class ActorModel(tree:Tree[_],components:ActorModelComponent){
  components.start
  tree.start
  val params = (components !? GetParams).asInstanceOf[Params].list.removeDuplicates
  val paramMap = params.foldLeft(Map[ParamName,ActorParamComponent]()){(m,p)=>m+((p.name,p))}
  def getParam(p:ParamName)={
    p match {
      case SingleParam(p2)=>SingleParamWrapper(paramMap(p2))
      case p => paramMap(p)
    }
  }
  val paramLengthMap = paramMap.keys.map{t=>(t,optGet(t).length)}.foldLeft(Map[ParamName,Int]()){(m,p)=>m+((p._1,p._2))}
  def update(p:ParamName,x:Array[Double]) {getParam(p) !? ParamUpdate[Array[Double]](x)}
  def update(p:ParamName,x:Double) {getParam(p) !? ParamUpdate(Array(x))}
  def update(p:ParamName,x:Vector) {getParam(p) !? ParamUpdate(x.toArray)}
  def update(p:ParamName,x:Matrix) {getParam(p) !? ParamUpdate(x)}
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
          param !? ParamUpdate(d)
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
          println(i + " " + len + " " + s)
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
    class APSetter extends PSetter{
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
    def getP={(p !? RequestParam).asInstanceOf[ParamChanged[Array[Double]]].param}
    class APSetter extends PSetter{
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
    class APSetter extends PSetter{
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
    class APSetter extends PSetter{
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


  def logLikelihood = {
    object Send
    Actor.actor{
      receive{
      case Send => 
        (tree ! LogLikelihoodCalc(components,self))
        var ans:Option[Double]=None
        receive{
          case LogLikelihood(d) => 
            ans=Some(d)
        }
        reply{ans.get}
      }
      exit
    } !? Send
  }.asInstanceOf[Double]

  def paramString=paramMap.map{t=> t._2.toString}.mkString("\n")

  def optimise(params:ParamName*)={
    import ModelOptimiser._
    ModelOptimiser.optimise(getConjugateDirection,params.toList,this)
  }

  
}

