package org.modiphy.math
import EnhancedMatrix._
import org.modiphy.tree._
import org.modiphy.tree.DataParse.Tree
import org.modiphy.sequence._
import scala.collection.immutable.Map
import scala.actors.Actor
import scala.actors.Actor._

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
case object RequestParam
case object RequestOpt

abstract class ParamName

case object Pi extends ParamName
case object S  extends ParamName
case object BranchLengths  extends ParamName
case object Alpha extends ParamName
case object InvarPrior extends ParamName
case object Sigma extends ParamName

case class ParamUpdate[A](d:A)
case class ParamChanged[A](p:ParamName,param:A)


/**
 Optimiser may have different view of params (log of real params, for example)
*/
case class OptUpdate(d:Array[Double])

class BasicActorModel(piParam:ActorPiComponent,sParam:ActorSComponent,reciever:Actor) extends ActorModelComponent{
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
      case ParamChanged(Pi,v:Vector)=>initialise(s,Some(v))
      case ParamChanged(S,m:Matrix)=>initialise(Some(m),pi)
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
        case ParamChanged(Pi,v:Vector)=>
          reciever forward Unclean
          main(s,v.copy,None)
        case ParamChanged(S,m:Matrix)=>
          reciever forward Unclean
          main(m.copy,pi,None)
     }
  }
}

class THMMActorModel(sigmaParam:ActorFullSComponent,numClasses:Int,reciever:Actor) extends ActorModelComponent{
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
      case ParamChanged(Sigma,a:Matrix)=>main(a,None)
    }
  }
  def main(sigma:Matrix,mat:Option[Matrix]){
    react{
      case Params(l)=> 
        reciever forward Params(sigmaParam :: l)
        main(sigma,mat)
      case ParamChanged(Sigma,array:Matrix)=>
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
      case ParamChanged(InvarPrior,d:Double)=>initialise(Some(d),pi)
      case ParamChanged(Pi,v:Vector)=>initialise(prior,Some(v))
    }
  }
  def main(prior:Double,rawpi:Vector,processedPi:Option[Vector],mat:Option[Matrix]){
    react{
      case Params(l)=> 
        reciever forward Params(priorParam :: piParam :: l)
        main(prior,rawpi,processedPi,mat)
      case ParamChanged(InvarPrior,d:Double)=>
        reciever forward Unclean
        main(d,rawpi,None,None)
      case ParamChanged(Pi,v:Vector)=>
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
      case ParamChanged(Alpha,alpha:Double)=>main(alpha,None,None)
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
      case ParamChanged(Alpha,d:Double)=>
        reciever forward Unclean
        main(d,None,pi)
      case Unclean =>
        reciever forward Unclean
        main(alpha,None,None)
    }
  }
}

class BasicSingleExpActorModel[A <: BioEnum](tree:Tree[A],branchLengthParams:ActorTreeComponent[A],reciever:Option[Actor]) extends ActorModelComponent{
  override def start={
    branchLengthParams.start
    if (reciever.isDefined){reciever.get.start}
    super.start
  }
  var eigen:MatrixExponential=null 

  def act{
    main(None)
  }
  def main(eigen:Option[MatrixExponential]){
    case class ExpReq(n:Node[_],e:MatrixExponential,lengthTo:Double,pi:Vector)
      react{
        case Params(l)=>
          if (reciever.isDefined){reciever.get forward Params(branchLengthParams :: l)}
          else {sender ! Params(branchLengthParams :: l)}
          main(eigen)
        case MatReq(n,Some(m),Some(pi)) => 
          if (n.isRoot){//don't need exp(qt)
            sender ! MatReq(n,None,Some(pi))
            main(eigen)
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
                //  println("Node " + n.id + " " + ans(0,0))
                //  println("Send Mat " + ans)
                  sender ! MatReq(n,Some(ans),Some(pi)) 
                  exit
              }
            }
          }.start forward ExpReq(n,myEigen,n.lengthTo,pi)
          main(Some(myEigen))
        case Unclean=>
          if (reciever.isDefined){reciever.get forward Unclean} else { reply(Unclean)}
          main(None)
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
  def myParam = pi
  override def setRaw(p:Array[Double]){
    internal setPi p
  }
  override def toString=name.toString + " " + internal.view
}
class ActorSComponent(s:Matrix,val name:ParamName) extends AbstractActorParam[Matrix] with SMatUtil{
  val internal = new SMatParam(s,name.toString)
  def myParam = s
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


class ActorTreeComponent[A <: BioEnum](tree:Tree[A],val name:ParamName) extends AbstractActorParam[Array[Double]]{
  val internal = new TreeParamControl[A](tree)
  def myParam = internal.getParams
}
class ActorDoubleComponent(param:Double,val name:ParamName) extends AbstractActorParam[Double]{
  class DoubleParam{
    var myP:Double=param
    def getParams=Array(myP)
    def setParams(a:Array[Double]){myP=a(0)}
  }
  val internal = new DoubleParam
  def myParam = internal.getParams(0)
}
class ActorArrayComponent(param:Array[Double],val name:ParamName) extends AbstractActorParam[Array[Double]]{
  class ArrayParam{
    var myP=param
    def getParams = myP.toArray
    def setParams(a:Array[Double]){Array.copy(a,0,myP,0,Math.min(a.length,myP.length))}
  }
  val internal = new ArrayParam
  def myParam = internal.getParams
}
class ActorProbComponent(prob:Double,name:ParamName) extends ActorDoubleComponent(prob,name)
class ActorGammaComponent(alpha:Double,name:ParamName) extends ActorDoubleComponent(alpha,name)
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

class ActorModel(tree:Tree[_],components:ActorModelComponent){
  components.start
  val params = (components !? GetParams).asInstanceOf[Params].list.removeDuplicates
  val paramMap = params.foldLeft(Map[ParamName,ActorParamComponent]()){(m,p)=>m+((p.name,p))}
  def update(p:ParamName,x:Array[Double]) {paramMap(p) !? ParamUpdate[Array[Double]](x)}
  def update(p:ParamName,x:Double) {paramMap(p) !? ParamUpdate(Array(x))}
  def update(p:ParamName,x:Vector) {paramMap(p) !? ParamUpdate(x.toArray)}
  def update(p:ParamName,x:Matrix) {paramMap(p) !? ParamUpdate(x)}
  def optUpdate(p:ParamName,x:Array[Double]){paramMap(p) !? OptUpdate(x)}

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

  //model(Pi)=Array(0.1,0.2)
}

