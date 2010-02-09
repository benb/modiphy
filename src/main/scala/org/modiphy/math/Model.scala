package org.modiphy.math

import cern.colt.matrix._
import org.modiphy.tree._
import org.modiphy.sequence._
import org.modiphy.math.EnhancedMatrix._
import org.modiphy.tree.DataParse._
import org.modiphy.util._
import scala.collection.immutable.IntMap
import tlf.Logging
import org.modiphy.util.PList._

/**
 Observer pattern. Anything that implements receiveUpdate(subject) can observe a Subject
*/
trait Subject {
  type Observer = { def receiveUpdate(subject:Subject) }
  private var observers = List[Observer]()
  def addObserver(o:Observer){observers=o::observers}
  def notifyObservers{observers.foreach{_.receiveUpdate(this)}}
}

class InvalidMatrixException(m:String) extends RuntimeException(m)

abstract class ParamControl extends Subject{
  def getParams:Array[Double]
  def setParams(a:Array[Double]):Unit
  def numParams=getParams.length
  def name:String
  def view:Vector
  /**
   check to see whether there is a difference before recaluclating
   useful if we combine parameters into a single large array
  */
  def softSetParams(a:Array[Double]):Unit
  def softSetParams(from:Array[Double],to:Array[Double])={
    if (
      !(from deepEquals to)
    ){
      setParams(to)
    }
  }


  def lower = Math.MIN_DOUBLE //default
  def upper = Math.MAX_DOUBLE //default

  def softLower = lower
  def softUpper = upper

  def lowerBound(i:Int)=lower
  def upperBound(i:Int)=upper

  lazy val lowerBounds = (0 until numParams).map{i=> lowerBound(i)}
  lazy val upperBounds = (0 until numParams).map{i=> upperBound(i)}

  def softLowerBound(i:Int)=softLower
  def softUpperBound(i:Int)=softUpper

  private def contained(i:Int)={val p = getParams(i); p >= lowerBound(i) && p <= upperBound(i)}
  def cromulent = (0 until numParams).foldLeft(true){(b,p)=> b && contained(p)}

  def makeCromulent{
    if (! cromulent){
      setParams(
        getParams.zipWithIndex.map{t=> if (t._1.isNaN || t._1 < softLowerBound(t._2)){softLowerBound(t._2)}else if (t._1 > softUpperBound(t._2)){softUpperBound(t._2)} else {t._1}}
      )
    }
  }
  def distFromCromulence={
    getParams.zipWithIndex.map{t=>
      if (t._1 < lowerBound(t._2)){t._1 - lowerBound(t._2)}else if (t._1 > upperBound(t._2)){upperBound(t._2) - t._1 } else {0.0D}
    }
  }
}

class BasicParamControl(a:Array[Double],val name:String) extends ParamControl{
  def getParams=a.toArray//copy
  def setParams(x:Array[Double]){
    x.copyToArray(a,0)
    notifyObservers
  }
  def view=a.toVector
  def softSetParams(x:Array[Double]){
    softSetParams(a,x)
  }
}

trait LogParamControl extends ParamControl{
  def sensibleMin = -5E2 // exp(sensibleMin) should be seems to be a small but non-zero number so opt is not flat
  def sensibleMax = 400D // again, exp(sensibleMax) should be < Inf
  abstract override def getParams=super.getParams.map{Math.log}
  abstract override def setParams(a:Array[Double])=super.setParams(a.map{Math.exp})
  override def upper = {val ans = Math.log(super.upper);if (ans.isNaN || ans > sensibleMax){sensibleMax}else{ans}}
  override def lower = {val ans = Math.log(super.lower);if (ans.isNaN || ans < sensibleMin){sensibleMin}else{ans}}
}

class TreeParamControl[A <: BioEnum](t:Tree[A]) extends ParamControl{
  var params = t.getBranchLengths
  def getParams = params.toArray
  def setParams(a:Array[Double]){
    params = a.toList
    notifyObservers
  }
  def getLatestTree:Tree[A] = {t.setBranchLengths(params)}
  val name="Branch Lengths"
  def view = getLatestTree.getBranchLengths.toVector
  def softSetParams(a:Array[Double])=softSetParams(params.toArray,a)
  override def lower = 0.0D
  override def upper = 10.0D
}

class LogTreeParamControl[A <: BioEnum](t:Tree[A]) extends TreeParamControl[A](t) with LogParamControl
/**
 The basic model, composed of an sMatrix, a pi matrix, and the MathComponent that converts the S into a Q matrix.
*/

class MixtureModel[A <: BioEnum](val models:List[Model[A]],val priors:PiComponent) extends BasicModel[A]{
  val nodeDependent = models.foldLeft(false){(a,b)=>a||b.nodeDependent} 
  assert(nodeDependent==false,"Components must be node independent to be used with MixtureModel")
  models.foreach{_.addObserver(this)}
  // tree assumed to be identical for all models
  def tree = models.head.tree
  def logLikelihood= {
    var subLikelihoods = models.map{a=>a.realLikelihoods}
    var patternCounts = models.head.tree.aln.pCount
    var total=0.0D
    val pri = priors(tree).toList
    while (subLikelihoods.head != Nil){
      total+=(patternCounts.head * Math.log(subLikelihoods.map{_.head}.zip(pri).foldLeft(0.0D){(a,b)=>a+(b._1*b._2)}))
      subLikelihoods = subLikelihoods.map{_.tail}
      patternCounts=patternCounts.tail
    }
    total
  }
  val params = (priors.getParams ++ models.map{_.params}.flatten[ParamControl].removeDuplicates.toList).toArray
  def cromulent = models.foldLeft(true){_ && _.cromulent}

  def setPi(a:Array[Double]){
    models.foreach{_.setPi(a)}
  }

}
class ComposeModel[A <: BioEnum](pi:PiComponent,s:SComponent,maths:MathComponent,var tree:Tree[A]) extends Model[A]{
  val alphabet = tree.alphabet
  val params = (s.getParams ++ pi.getParams ++ maths.getParams ++ List(treeParamControl)).toArray
  List(pi,s,maths).foreach{m =>
    m.addObserver(this)
  }

  val nodeDependent = pi.nodeDependent || s.nodeDependent || maths.nodeDependent
  def careAboutNodeDependence=true
  assert(!(careAboutNodeDependence) || !(nodeDependent),"Components must be node independent to be used with ComposeModel")
  override def qMat(node:Node[A])={
    maths.sToQ(node)(s(node),pi(node)) 
  }

  override def cromulent={
    println(params.map{c => c.cromulent.toString + " " + c.name})
    params.foldLeft(true){(a,b) => a && b.cromulent}
  }

  def setPi(a:Array[Double]){pi.setPi(a)}

  def sMat=s(tree)
  def sMat(node:Node[A])=s(node)
  def getParamName=""
  def getParamName(i:Int)=""

  def piVals(node:Node[A])=pi(node)
}

class NodeDependentComposeModel[A <: BioEnum](pi:PiComponent,s:SComponent,maths:MathComponent,tree:Tree[A]) extends ComposeModel[A](pi,s,maths,tree){
  override def careAboutNodeDependence=false
}


/**
 A building block for a model
*/
abstract class MComponent extends Subject{
  var clean=false
  def getParams:List[ParamControl]
  def receiveUpdate(s:Subject)={changedParam(s)}
  def changedParam(s:Subject)={notifyObservers; clean=false}//by default fire up the chain
  def cromulent=true
  val nodeDependent=true
}

abstract class SComponent extends MComponent{
  def paramName:String
  def apply(node:Node[_]):Matrix
}

class BasicSMatComponent(sMat:Matrix) extends SComponent with SMatUtil{
  def apply(node:Node[_])=sMat
  val paramName="S Matrix"
  val param = new SMatParam(sMat,"SMat")
  param.addObserver(this)
  def getParams=List(param)
  override val nodeDependent=false
}
abstract class MathComponent extends MComponent{
  def sToQ(node:Node[_])(start:Matrix,pi:Vector):Matrix 
}
class ScaledMathComponent(scale:SingleComponent) extends DefaultMathComponent{
  scale.addObserver(this)
  override def normalValue(node:Node[_])=scale(node)
}

class DefaultMathComponent extends MathComponent{
  var cacheQ:Matrix=null
  def normalValue(node:Node[_])=1.0D
  def sToQ(node:Node[_])(sMat:Matrix,pi:Vector)={
    if (!clean){
      cacheQ=sMat.sToQ(pi).normalize(pi,normalValue(node))
    }
    cacheQ
  }
  def getParams=List()
  override val nodeDependent=false
}
object PiComponent{
  def makeSensible(startPi:Vector)={
    val smallestPi=0.001D
    val initial = startPi.copy
    (0 until initial.size).foreach{i=> 
      if (initial(i) < smallestPi){
        initial(i)=smallestPi
      }
    }
    initial.normalize(1)

  }
}

abstract class VectorComponent extends MComponent{
  def apply(node:Node[_]):Vector
}
abstract class SingleComponent extends MComponent{
  def apply(node:Node[_]):Double
}
abstract class PiComponent extends VectorComponent{
  def setPi(a:Array[Double])
}

class PiParam(pi:Vector,val name:String) extends ParamControl{
  def view=pi.copy
  def this(pi:Vector)=this(pi,"Pi")
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

  var medianIndex=0
  var lastGotParams:Array[Double]=null
  def getParams:Array[Double]={
    val (l,m)=toFit(pi)
    medianIndex=m
    lastGotParams=l.toArray
    l.toArray
  }

  def setParams(a:Array[Double])={
    pi assign fromFit(a,medianIndex)
    notifyObservers
  }
  
  def setPi(a:Array[Double])={
    pi assign a
    notifyObservers
  }
  def softSetParams(a:Array[Double])={
   softSetParams(lastGotParams,a) 
  }

  override def softLower = -10.0D
  override def softUpper = 10.0D
}

class BasicPiComponent(startPi:Vector) extends PiComponent{
  val pi = PiComponent.makeSensible(startPi)
  val param = new PiParam(pi)
  param.addObserver(this)
  val getParams=List(param)
  def apply(node:Node[_])=pi
  override val cromulent=true // should not be possible to get completely awful pis
  override val nodeDependent=false

  def setPi(a:Array[Double]){
    param.setPi(a)
  }
}

class PriorPiComponent(startPi:PiComponent,numClasses:Int,numAlpha:Int) extends PiComponent{
  def this(startPi:PiComponent,alphabet:BioEnum)=this(startPi,alphabet.numClasses,alphabet.numAlpha)
  assert(startPi.nodeDependent==false)//not implemented node-dependent pis here
  override val nodeDependent=false

  val pi=Vector(numAlpha * numClasses)
  startPi.addObserver(this)
  var prior = Vector((0 until numClasses).map{i=> 1.0D/numClasses}.toArray)
  val param:ParamControl = new PiParam(prior,"Prior Pi")
  param.addObserver(this)
  def apply(node:Node[_]) = {
    if (!clean){
      recaluclatePi
      clean=true
    }
    pi
  }
  def recaluclatePi{
    (0 until numClasses).foreach{i=>
      pi.viewPart(i * numAlpha ,numAlpha).assign(startPi(null))*prior(i)
      //println("Added class " + i + "\n" + pi)
    }
  }
  def getParams=param :: startPi.getParams

  def getView(start:Int):PiComponent=getView(start,numAlpha*numClasses)
  def getView(start:Int,stop:Int):PiComponent={
    val outer=this
    new PiComponent{
      def apply(node:Node[_])=outer.apply(node).viewPart(start,stop-start).copy
      def getParams=List()
      def setParams(array:Array[Double])=throw new IllegalArgumentException("Can't optimise PriorPiView")
      override val nodeDependent = outer.nodeDependent
      outer.addObserver(this)
      def setPi(a:Array[Double]){}
    }
  }
  def setPi(a:Array[Double])=startPi.setPi(a)
}

class FirstOnlyPiParam(pi:Vector,val name:String) extends ParamControl{
  def getParams=Array(pi(0))
  def setParams(a:Array[Double]){
    pi(0)=a(0)
    for (i<-1 until pi.size){
      pi(i)=(1.0D-pi(0))/(pi.size-1)
    }
    notifyObservers
  }
  def view=pi.copy
  def softSetParams(a:Array[Double]){
    if (a(0)!=pi(0)){
      setParams(a)
    }
  }
  override def lower=0.0D
  override def softUpper=1.0D-1E-5
  override def upper=1.0D
}


class FlatPriorPiComponent(startPi:PiComponent,numClasses:Int,numAlpha:Int) extends PriorPiComponent(startPi,numClasses,numAlpha){
  def this(startPi:PiComponent,alphabet:BioEnum)=this(startPi,alphabet.numClasses,alphabet.numAlpha)
  override def getParams = startPi.getParams
}
class FlatPrior(numClasses:Int) extends PiComponent{
  val priors = Vector(numClasses)
  priors.assign(1.0D/numClasses)
  def apply(n:Node[_])=priors
  def getParams=Nil
  def setParams(a:Array[Double]){}
  def setPi(a:Array[Double]){

  }
}


class FirstPriorPiComponent(startPi:PiComponent,numClasses:Int,numAlpha:Int) extends PriorPiComponent(startPi,numClasses,numAlpha){
  def this(startPi:PiComponent,alphabet:BioEnum)=this(startPi,alphabet.numClasses,alphabet.numAlpha)
  val myParam = new FirstOnlyPiParam(prior,"First Prior") //with LogParamControl
  myParam.addObserver(this)
  override def getParams = myParam :: startPi.getParams

  override def cromulent = {
    if (!clean){
      recaluclatePi
      clean=true
    }
    val sum = pi.zSum
    pi.toArray.foldLeft(true){_ && _ > -Math.EPS_DOUBLE} && pi.zSum > 0.99999 && pi.zSum < 1.00001
  }
}

abstract class SMatComponent extends MComponent with SMatUtil{
  def sMat:Matrix
  override def cromulent = !(sMat exists {i=> i < 0.0D}) 
  def linearSMat:Seq[Double]=linearSMat(sMat)
  def apply(node:Node[_]):Matrix
}

class GammaMathComponent(a:Double,numClasses:Int,numAlpha:Int,pi:PiComponent,s:SComponent) extends MathComponent{
  def this(a:Double,alphabet:BioEnum,pi:PiComponent,s:SComponent)=this(a,alphabet.numClasses,alphabet.numAlpha,pi,s)
  pi.addObserver(this)
  s.addObserver(this)
  override def cromulent = alpha(0)>0.0D && alpha(0)<1E5 && !(alpha(0).isNaN) &&  super.cromulent
  def matLength = numClasses * numAlpha
  val alpha = Array(a)
  def gMath = new Gamma(numClasses)
  val param=new BasicParamControl(alpha,"Alpha Shape") {override def lower=1E-6; override def upper=100.0D}//with LogParamControl {override def lower=1.01D;override def upper=11.5}
  def getParams=List(param)
  param.addObserver(this)
  val internalS=Matrix(matLength,matLength)
  override val nodeDependent = false
  var cachedQ:Matrix=null

  def rates=gMath(alpha(0))
  
  def sToQ(node:Node[_])(sMat:Matrix,pi:Vector)={
      if (nodeDependent || !(clean)){
        val r = rates
        (0 until numClasses).foreach{i=>
          internalS.viewPart(i * numAlpha,i * numAlpha,numAlpha,numAlpha).assign(sMat)*r(i)
        }
        cachedQ= internalS.sToQ(pi).normalize(pi)
        clean=true
      }
    
    cachedQ
  }
  
}

class InvariantMathComponent(numAlpha:Int,pi:PiComponent,s:SComponent,base:GammaMathComponent) extends  MathComponent{
  var cachedQ:Matrix=null
  pi.addObserver(this)
  s.addObserver(this)
  base.addObserver(this)
  def sToQ(node:Node[_])(sMat:Matrix,pi:Vector)={
    if (!(clean)){
      val startQ = base.sToQ(node)(sMat,pi.viewPart(numAlpha,pi.size-numAlpha))
      if (cachedQ==null){
        cachedQ = Matrix(startQ.rows + numAlpha, startQ.rows + numAlpha)
      }
      cachedQ assign 0.0D
      cachedQ.viewPart(numAlpha,numAlpha,startQ.rows,startQ.rows).assign(startQ)
      cachedQ = cachedQ.normalize(pi)
      clean=true
    }
    cachedQ
  }
  def getParams=base.getParams
  override val nodeDependent=false
  override def cromulent = base.cromulent && super.cromulent
}

class THMMGammaMathComponent(gMath:MathComponent,cMat:Matrix,alphabet:BioEnum) extends MathComponent{
  def scale(node:Node[_])=1.0D
  override val nodeDependent = false
  var cachedQ:Matrix=null
  val param = new FullSMatParam(cMat,"CMat")
  param.addObserver(this)
  gMath.addObserver(this)
  override def sToQ(node:Node[_])(sMat:Matrix,pi:Vector)={
    import alphabet._
    if (!clean || nodeDependent){
      val qStart = gMath.sToQ(node)(sMat,pi)
        for (i <- 0 until numClasses){
        for (j <- i+1 until numClasses){
          for (x <- 0 until numAlpha){
            qStart(i * numAlpha + x, j * numAlpha + x) = cMat(i,j) * scale(node) *  pi(j * numAlpha + x) // i->j transition

            qStart(j * numAlpha + x, i * numAlpha + x) = cMat(i,j) * scale(node) * pi(i * numAlpha + x) // j->i transition
          }
        }
      }
      cachedQ=qStart.fixDiag
  //    println("CMAT " + cMat)
  //    println("QMAT " + cachedQ)
      clean=true
    }
    cachedQ
    
  }
  override def cromulent = gMath.cromulent && super.cromulent
  def getParams=param::gMath.getParams
}

class THMMGammaMathComponentBranch[A <: BioEnum](gMath:MathComponent,cMat:Matrix,alphabet:BioEnum,scales:Array[Double],tree:Tree[A]) extends THMMGammaMathComponent(gMath,cMat,alphabet){
  assert(scales.length==tree.getBranchLengths.length)
  override val nodeDependent = true
  val nodeOrder:Map[Int,Int]=tree.nodes.filter{n:Node[A]=> !(n.isRoot)}.map{_.id}.zipWithIndex.foldLeft[Map[Int,Int]](IntMap[Int]()){_+_}
  override def scale(node:Node[_])={
    assert(nodeOrder isDefinedAt node.id)
    scales(nodeOrder(node.id)) / node.lengthTo // now appears independent of branch length to optimiser
  }

  val scaleParam=new BasicParamControl(scales,"THMM Scale") {override def lower=0.0D; override def softUpper=1E100} 
  override def getParams=scaleParam::super.getParams
  override def cromulent = gMath.cromulent && super.cromulent
}

object ModelFact{
  def basic[A <: BioEnum](pi:Vector,s:Matrix,tree:Tree[A])={
    val sC=new BasicSMatComponent(s)
    val piC = new BasicPiComponent(pi)
    new ComposeModel(piC,sC,new DefaultMathComponent,tree)
  }
  def gamma[A <: BioEnum](pi:Vector,s:Matrix,alpha:Double,tree:Tree[A])={
    val sC=new BasicSMatComponent(s)
    //val cC=new SMatComponent(cC)
    val piC = new FlatPriorPiComponent(new BasicPiComponent(pi),tree.alphabet)
    val gammaC=new GammaMathComponent(alpha,tree.alphabet,piC,sC)
    new ComposeModel(piC,sC,gammaC,tree)
  }
  def gammaMixture[A <: BioEnum](pi:Vector,s:Matrix,alpha:Double,tree:Tree[A],numSC:Int)={
    val sC = new BasicSMatComponent(s)
    val piC = new BasicPiComponent(pi)
    val gammaC = new GammaMixtureComponent(alpha,numSC)
    import tree.alphabet._
    val modelList = (for (i <- 0 until numSC) yield new ComposeModel(piC,sC,new ScaledMathComponent(gammaC.getView(i)),tree)).toList
    new MixtureModel(modelList,new FlatPrior(numSC))
  }
  def thmm[A <: BioEnum](pi:Vector,s:Matrix,alpha:Double,cMat:Matrix,tree:Tree[A])={
    val sC = new BasicSMatComponent(s)
    val piC = new BasicPiComponent(pi)
    val piC2 = new FlatPriorPiComponent(piC,tree.alphabet)
    val gammaC = new GammaMathComponent(alpha,tree.alphabet,piC2,sC)
    val thmmMath = new THMMGammaMathComponent(gammaC,cMat,tree.alphabet)
    new ComposeModel(piC2,sC,thmmMath,tree)
  }
  def invarThmm[A <: BioEnum](pi:Vector,s:Matrix,alpha:Double,cMat:Matrix,tree:Tree[A])={
    val sC = new BasicSMatComponent(s)
    val piC = new BasicPiComponent(pi)
    val piC2 = new FirstPriorPiComponent(piC,tree.alphabet)
    val view = piC2.getView(tree.alphabet.numAlpha)
    val gammaC = new GammaMathComponent(alpha,tree.alphabet.numClasses-1,tree.alphabet.numAlpha,view,sC)
    val invariantMath = new InvariantMathComponent(tree.alphabet.numAlpha,piC2,sC,gammaC)
    val thmmMath = new THMMGammaMathComponent(invariantMath,cMat,tree.alphabet)
    new ComposeModel(piC2,sC,thmmMath,tree)

  }
  def invarThmmBS[A <: BioEnum](pi:Vector,s:Matrix,alpha:Double,cMat:Matrix,tree:Tree[A])={
    val sC = new BasicSMatComponent(s)
    val piC = new BasicPiComponent(pi)
    val piC2 = new FirstPriorPiComponent(piC,tree.alphabet)
    val view = piC2.getView(tree.alphabet.numAlpha)
    val gammaC = new GammaMathComponent(alpha,tree.alphabet.numClasses-1,tree.alphabet.numAlpha,view,sC)
    val invariantMath = new InvariantMathComponent(tree.alphabet.numAlpha,piC2,sC,gammaC)
    val blScale = tree.getBranchLengths.toArray
    val thmmMath = new THMMGammaMathComponentBranch(invariantMath,cMat,tree.alphabet,blScale,tree)
    new NodeDependentComposeModel(piC2,sC,thmmMath,tree)
  }
}

/**
 A basic model, but not necessarily one that can be passed to a tree to do likelihood calculations

*/
abstract class BasicModel[A <: BioEnum] extends Logging with Subject{

  def setPi(a:Array[Double])
  var clean=false
 
  def params:Array[ParamControl]

  def nodeDependent:Boolean
  def cromulent:Boolean
  def getParams=params.map{_.getParams}.toList
  def getParam(s:String)=params.filter{_.name.toLowerCase==s.toLowerCase}.toList
  def setParams(i:Int)(a:Array[Double]){params(i).setParams(a)}
  def getParams(i:Int)=params(i).getParams

  def receiveUpdate(s:Subject){
    clean=false
    notifyObservers
  }
  def logLikelihood:Double
}
abstract class Model[A <: BioEnum] extends BasicModel[A]{

  val treeParamControl = new LogTreeParamControl(tree)
  treeParamControl.addObserver(this)

  override def receiveUpdate(s:Subject){
    if (s==treeParamControl){
      tree = treeParamControl.getLatestTree
    }
    super.receiveUpdate(s)
  }
 
 
  var tree:Tree[A]
  def piVals(node:Node[A]):Vector
  def piVals:Vector=piVals(tree) //assume we want pi at root unless specified
  def qMat(node:Node[A]):Matrix
  def sMat:Matrix
  def getMat(node:Node[A])={ debug{"e^Qt: t==" + node.lengthTo}; qMat(node).exp(node.lengthTo) }

  def getParamName(i:Int):String

  def sMat(node:Node[A]):Matrix

  def likelihoods(node:Node[A]):List[Vector]={
    val childLkl = node.children.map{i:Node[A]=>(i.likelihoods(this),i.lengthTo)}.toList

    val intermediates= childLkl.zip(node.children.toList).map{t=>
    val ((siteVectorList,length),c)=t // list of vectors  1 for each site
    //println ("siteVectorList " + siteVectorList)
    //println("Q = " + qMat)
      debug{"Q = " + qMat(c)}
      val matrix = getMat(c) //e^Qt
      siteVectorList.map{siteVector=>
        val alphabet = c.alphabet

        debug{"e^Qt=" + matrix}
        val ret = new Array[Double](alphabet.matLength) 

        (0 to alphabet.matLength - 1).foreach{i=>
          ret(i)=siteVector.zDotProduct(matrix.viewRow(i))
        }
        Vector(ret)
      }
    }.toList
    val ans = intermediates.head
    val func = new cern.colt.function.DoubleDoubleFunction{def apply(x:Double,y:Double)=x*y}
    intermediates.tail.foreach{list2=>
      ans.zip(list2).map{t=> // not really a map but used for parallel reasons
        val (vec,vec2)=t
        vec.assign(vec2,func)
      }
    }
    ans
  }

  def likelihoods:List[Vector]=tree.likelihoods(this)
  def realLikelihoods=tree.realLikelihoods(this)
  def logLikelihood:Double={
    /*
    val (p:Option[List[List[Double]]],diff:Double)=
  if (cromulent){
    (None,0.0D) 
  }else {
    val ans = (Some(params.map{_.getParams}),params.map{_.distFromCromulence.toList}.toList.flatten[Double].foldLeft(0.0D){(a,b)=>a+b.abs})
      println(params.map{a => (a.distFromCromulence.toList, a.name)}.toList)
    params.foreach{_.makeCromulent}
    println("Penalty " + ans._2)
    ans
  }*/
  if (cromulent){
    val lnL = try{
      tree.logLikelihood(this)
    }catch {
      case mex:org.apache.commons.math.MathException => Math.NaN_DOUBLE
      case iex:InvalidMatrixException =>{
        println("Trying invalid params " + iex)
        Math.NaN_DOUBLE
      }

    }/*finally{
      if (p.isDefined){
        params.toList.zip(p.get).foreach{t:(ParamControl,List[Double])=> t._1.setParams(t._2.toArray)}
      }
    }*/
    lnL
  }else {
    -1E100
  }
  }
  def makeCromulent=params.foreach{_.makeCromulent}
}

class SMatParam(sMat:Matrix,val name:String) extends ParamControl with SMatUtil{
  def getParams={
    lastParams = linearSMat(sMat).toArray
    lastParams.toArray
  }
  var lastParams:Array[Double]=linearSMat(sMat).toArray
  def setParams(a:Array[Double])={
    setSMat(a,sMat)
    notifyObservers
  }
  def view=getParams.toVector
  def softSetParams(a:Array[Double]){softSetParams(lastParams,a)}
  override def lower=0.0D
  override def upper=1E20
}
class FullSMatParam(sMat:Matrix,name:String) extends SMatParam(sMat,name){
  override def getParams=linearSMatFull(sMat).toArray
}
trait SMatUtil{
  def linearSMatFull(sMatInst:Matrix):Seq[Double]={
    val list = (0 to sMatInst.rows -2).map{ i=> // don't skip last elment...
       val list1 = sMatInst(i).viewPart(i+1,sMatInst.columns - i - 1).toList
       list1
    }.toList.flatten[Double]
    list
  }
  def linearSMat(sMatInst:Matrix):Seq[Double]={
    val list = (0 to sMatInst.rows -3).map{ i=> // skip last elment...
       val list1 = sMatInst(i).viewPart(i+1,sMatInst.columns - i - 1).toList
       list1
    }.toList.flatten[Double]
    list
  }
def setSMat(array:Array[Double],sMatInst:Matrix){
    val iter:Iterator[Double]=array.elements
    (0 to sMatInst.rows-1).foreach{i=>
      (i+1 to sMatInst.columns-1).foreach{j=>
        if (iter.hasNext){
          sMatInst(i,j)=iter.next
        }else{ //by passing in a shorter array, remaining elements set to 1
          sMatInst(i,j)=1
        }
      }
    }
    sMatInst
  }
}

object Gamma{
  val cache = new CacheMap[(Int,Double),Array[Double]](100)
}

class Gamma(numCat:Int){
  import org.apache.commons.math.distribution.ChiSquaredDistributionImpl
  import cern.jet.stat.Gamma.incompleteGamma
  import Gamma._
  val chi2=new ChiSquaredDistributionImpl(1.0D)
  def chiSquareInverseCDF(prob:Double,df:Double)={
     chi2.setDegreesOfFreedom(df)
     val ans = chi2.inverseCumulativeProbability(prob) 
     ans
  }
  def gammaInverseCDF(prob:Double,alpha:Double,beta:Double)=chiSquareInverseCDF(prob,2.0*(alpha))/(2.0*(beta))


  def apply(shape:Double):Array[Double]={
    cache.getOrElseUpdate((numCat,shape),gamma(shape))
  }
  def gamma(shape:Double):Array[Double]={
    if (shape==Math.POS_INF_DOUBLE){
      (0 until numCat).map{a=>1.0}.toArray
    }else {
      val alpha = shape
      val beta = shape
      val factor=alpha/beta*numCat
      val freqK=new Array[Double](numCat)
        val rK=new Array[Double](numCat)

        (0 until numCat-1).foreach{i=>
          freqK(i)=gammaInverseCDF((i+1.0)/numCat, alpha, beta);
          freqK(i) = incompleteGamma(alpha+1,freqK(i)*beta)
        }

        rK(0) = freqK(0)*factor;
        rK(numCat-1) = (1-freqK(numCat-2))*factor;
        (1 until numCat-1).foreach{i=>
        rK(i) = (freqK(i)-freqK(i-1))*factor;
      }
      // println("RATES " + rK.toList)
      rK
    }
  }
}

class GammaMixtureComponent(initAlpha:Double,numCat:Int) extends VectorComponent{
  val param = new BasicParamControl(Array(initAlpha),"Alpha Shape")
  def getParams=List(param)
  val gammaMath = new Gamma(numCat)
  var gammaVect:Vector=Vector(numCat)
  def apply(n:Node[_])={
    if (!clean){
      gammaVect.assign(gammaMath(param.getParams(0)))
    }
    gammaVect
  }

  def getView(i:Int)={
    assert(i < numCat)
    val out=this
    val ans = new SingleComponent{
      def apply(n:Node[_])=out(n)(i)
      def getParams=List()
      override val nodeDependent = out.nodeDependent
    }
    ans.addObserver(this)
    ans
  }
  override val nodeDependent=false
}


