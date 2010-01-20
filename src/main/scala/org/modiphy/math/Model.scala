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

trait Subject {
  type Observer = { def receiveUpdate(subject:Subject) }
  private var observers = List[Observer]()
  def addObserver(o:Observer){observers=o::observers}
  def notifyObservers{observers.foreach{_.receiveUpdate(this)}}
}

abstract class ParamControl extends Subject{
  def getParams:Array[Double]
  def setParams(a:Array[Double]):Unit
}
class BasicParamControl(a:Array[Double]) extends ParamControl{
  def getParams=a.toArray//copy
  def setParams(x:Array[Double]){
    x.copyToArray(a,0)
    notifyObservers
  }
}
/**
 The basic model, composed of an sMatrix, a pi matrix, and the MathComponent that converts the S into a Q matrix.
*/
class ComposeModel[A <: BioEnum](pi:PiComponent,s:SComponent,maths:MathComponent,var tree:Tree[A]) extends Model[A]{
  val alphabet = tree.alphabet
  val params = (s.getParams ++ pi.getParams ++ maths.getParams).toArray
  List(pi,s,maths).foreach{m =>
    m.addObserver(this)
  }
  var clean=false
  def receiveUpdate(s:Subject){clean=false}
  val nodeDependent = pi.nodeDependent || s.nodeDependent || maths.nodeDependent
  assert(nodeDependent==false,"Components must be node independent to be used with ComposeModel")
  var qMatCache:Matrix=null
  override def qMat(node:Node[A])={
    maths.sToQ(node)(s(node),pi(node)) 
  }

  override def cromulent=pi.cromulent && s.cromulent && maths.cromulent
  def getParams=params.map{_.getParams}.toList
  def setParams(i:Int)(a:Array[Double]){params(i).setParams(a)}
  def getParams(i:Int)=params(i).getParams



  def sMat=s(tree)
  def sMat(node:Node[A])=s(node)
  def getParamName=""
  def getParamName(i:Int)=""

  def piVals(node:Node[A])=pi(node)
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
  val param = new SMatParam(sMat)
  param.addObserver(this)
  def getParams=List(param)
  override val nodeDependent=false
}

abstract class MathComponent extends MComponent{
  def sToQ(node:Node[_])(start:Matrix,pi:Vector):Matrix 
}

class DefaultMathComponent extends MathComponent{
  var cacheQ:Matrix=null
  def sToQ(node:Node[_])(sMat:Matrix,pi:Vector)={
    if (!clean){
      cacheQ=sMat.sToQ(pi).normalize(pi)
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

abstract class PiComponent extends MComponent{
  def apply(node:Node[_]):Vector
}

class PiParam(pi:Vector) extends ParamControl{
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
  def getParams:Array[Double]={
    val (l,m)=toFit(pi)
    medianIndex=m
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
}

class BasicPiComponent(startPi:Vector) extends PiComponent{
  val pi = PiComponent.makeSensible(startPi)
  val param = new PiParam(pi)
  param.addObserver(this)
  val getParams=List(param)
  def apply(node:Node[_])=pi
  override val cromulent=true // should not be possible to get completely awful pis
  val paramName="Pi Values"
  override val nodeDependent=false
}

class PriorPiComponent(startPi:PiComponent,alphabet:BioEnum) extends PiComponent{
  assert(startPi.nodeDependent==false)//not implemented node-dependent pis here
  override val nodeDependent=false

  val pi=Vector(alphabet.matLength)
  startPi.addObserver(this)
  val numClasses = alphabet.numClasses
  var prior = Vector((0 until numClasses).map{i=> 1.0D/numClasses}.toArray)
  val param = new PiParam(prior)
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
      pi.viewPart(i * alphabet.numAlpha ,alphabet.numAlpha).assign(startPi(null))*prior(i)
      //println("Added class " + i + "\n" + pi)
    }
  }
  def getParams=param :: startPi.getParams
}
class FlatPriorPiComponent(startPi:PiComponent,alphabet:BioEnum) extends PriorPiComponent(startPi,alphabet){
  override def getParams = startPi.getParams
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
  def matLength = numClasses * numAlpha
  val alpha = Array(a)
  def gMath = new Gamma(numClasses)
  val param=new BasicParamControl(alpha)
  def getParams=List(param)
  param.addObserver(this)
  val internalS=Matrix(matLength,matLength)
  override val nodeDependent = false
  var cachedQ:Matrix=null
  
  def sToQ(node:Node[_])(sMat:Matrix,pi:Vector)={
    if (!(clean)){
    val rates = gMath(alpha(0))
    (0 until numClasses).foreach{i=>
      internalS.viewPart(i * numAlpha,i * numAlpha,numAlpha,numAlpha).assign(sMat)*rates(i)
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
      val startQ = base.sToQ(node)(sMat,pi)
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
  def getParams=List()
  override val nodeDependent=false
}

class THMMGammaMathComponent(gMath:GammaMathComponent,cMat:Matrix,alphabet:BioEnum) extends MathComponent{
  override val nodeDependent = false
  var cachedQ:Matrix=null
  val param = new FullSMatParam(cMat)
  param.addObserver(this)
  gMath.addObserver(this)
  override def sToQ(node:Node[_])(sMat:Matrix,pi:Vector)={
    import alphabet._
    if (!clean){
      val qStart = gMath.sToQ(node)(sMat,pi)
        for (i <- 0 until numClasses){
        for (j <- i+1 until numClasses){
          for (x <- 0 until numAlpha){
            qStart(i * numAlpha + x, j * numAlpha + x) = cMat(i,j) * pi(j * numAlpha + x) // i->j transition

            qStart(j * numAlpha + x, i * numAlpha + x) = cMat(i,j) * pi(i * numAlpha + x) // j->i transition
          }
        }
      }
      cachedQ=qStart.fixDiag
      clean=true
      //println("Q:\n" + cachedQ)
    }
    cachedQ
    
  }
  def getParams=param::gMath.getParams
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
  def thmm[A <: BioEnum](pi:Vector,s:Matrix,alpha:Double,cMat:Matrix,tree:Tree[A])={
    val sC = new BasicSMatComponent(s)
    val piC = new BasicPiComponent(pi)
    val piC2 = new FlatPriorPiComponent(piC,tree.alphabet)
    val gammaC = new GammaMathComponent(alpha,tree.alphabet,piC2,sC)
    val thmmMath = new THMMGammaMathComponent(gammaC,cMat,tree.alphabet)
    new ComposeModel(piC2,sC,thmmMath,tree)
  }
}



class UnknownParamException(i:Int) extends java.lang.RuntimeException("Unknown parameter " + i.toString)
trait Model[A <: BioEnum] extends Logging{
  
  var tree:Tree[A]
  def cromulent:Boolean

  def setParams(i:Int)(a:Array[Double]):Unit

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
  def logLikelihood=if (cromulent){tree.logLikelihood(this)}else{Math.NEG_INF_DOUBLE}
  def getParams(i:Int):Array[Double]
  def getParams:List[Array[Double]]
}

trait TraitModel[A <: BioEnum] extends Model[A] with SMatUtil{
  def getParams(i:Int):Array[Double]={Array()}
  def sMat(node:Node[A]):Matrix=sMat
  val careful=false
  def likelihoodTree = tree
  def getParams:List[Array[Double]] = List()
  def setParams(paramSet:Int)(params:Array[Double]){}
  def getParamName(i:Int):String=""

  def cromulent = true
  def piVals(node:Node[A])=piVals
  def qMat(node:Node[A]):Matrix={sMat(node).sToQ(pi).normalize(pi)}
  def qMat:Matrix=qMat(tree) //identical for all nodes in basic model, so can cache

  def pi:Vector
  def getTree = likelihoodTree
  override def toString ="Model:\n:"
 
}

class SMatParam(sMat:Matrix) extends ParamControl with SMatUtil{
  def getParams=linearSMat(sMat).toArray
  def setParams(a:Array[Double])={
    setSMat(a,sMat)
    notifyObservers
  }
}
class FullSMatParam(sMat:Matrix) extends SMatParam(sMat){
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

trait OptBranchLengths[A <: BioEnum] extends TraitModel[A]{
  var tree:Tree[A]
  private val paramNum = super.getParams.length
  override def getParamName(i:Int)={
    if (i != paramNum){super.getParamName(i)}
    else{"Branch Lengths"}
  }
  override def setParams(i:Int)(a:Array[Double])={
    if (i != paramNum){super.setParams(i)(a)}
    else {tree = tree.setBranchLengths(a.toList).setRoot}
  }
  override def getParams:List[Array[Double]]={
    super.getParams ++ List(tree.getBranchLengths.toArray)
  }
  override def toString={
    super.toString + "\n" + tree.toString
  }
  override def cromulent=super.cromulent && !(tree.getBranchLengths.exists{_<0.0D})
}

trait OptBranchScale[A <: BioEnum] extends TraitModel[A]{
  var tree:Tree[A]
  var treeScale:Double

  override def getMat(node:Node[A])={ qMat(node).exp(node.lengthTo * treeScale) }
  private val paramNum = super.getParams.length

  override def getParamName(i:Int)={
    if (i != paramNum){super.getParamName(i)}
    else{"Branch Scale"}
  }
  override def setParams(i:Int)(a:Array[Double])={
    if (i != paramNum){super.setParams(i)(a)}
    else {treeScale = a(0)}
  }
  override def getParams:List[Array[Double]]={
    super.getParams ++ List(Array(treeScale))
  }
  override def toString={
    super.toString + "\n" + tree.toString + " (scale: " + treeScale.toString + " )"
  }
  override def cromulent=super.cromulent && treeScale >= 0.0D
}

trait CombineOpt[A <: BioEnum] extends TraitModel[A]{
  def toCombine:List[Int]
  //following params need to be lazy so they are not evaluated until the rest of the model is set up
  // a mapping from new indicies to old
  lazy val paramMap=super.getParams.zipWithIndex.filter{t=>val (p,i)=t;!(toCombine.exists{j=>j==i})}.map(_._2).zipWithIndex.foldLeft(IntMap[Int]():Map[Int,Int]){(m,t)=>m+((t._2,t._1))}
  private lazy val paramNum = super.getParams.length - toCombine.length
  lazy val paramLengths = super.getParams.zipWithIndex.filter{t=>toCombine.exists{j=>j==t._2}}.map{_._1.length}
  override def getParams:List[Array[Double]]={
    val start = super.getParams
    //not super efficient but this should be rarely called
    start.zipWithIndex.filter{t=>val (array,i)=t; !toCombine.exists{j=>j==i}}.map{_._1} ++ List(start.zipWithIndex.filter{t=>toCombine.exists{j=>j==t._2}}.map{_._1.toList}.flatten[Double].toArray)
  }

  override def getParamName(i:Int)={
    if (i != paramNum){super.getParamName(i)}
    else{"Combined: " + toCombine.map{super.getParamName(_)}}
  }

  override def setParams(i:Int)(a:Array[Double]){
    if (paramMap contains i){
      super.setParams(paramMap(i))(a)
    }else {
      if (i==paramNum){
        var listPtr = a.toList
        paramLengths.zip(toCombine).foreach{t=>
          val (length,index)=t
          super.setParams(index)(listPtr.take(length).toArray)
          listPtr = listPtr drop length
        }
        
      }
      else {
        throw new Exception("Unsupported parameter value " + i)
      }
    }
  }

}
/*
trait ExposeOpt[A <: BioEnum] extends TraitModel[A]{
  val exposed:List[Int]
  override def getParams:List[Array[Double]]={
    super.getParams.zipWithIndex.filter{t=>
        val (p,i)=t
        exposed contains i
      }.sort{(i,j)=>exposed.find(_==i._2).get < exposed.find(_==j._2).get}.map{_._1}
  }

  override def getParamName(i:Int)={
    if (i != paramNum){super.getParamName(i)}
    else{"Expose: " + toCombine.map{getParamName}}
  }

  override def setParams(i:Int)(a:Array[Double]){
    super.setParams(exposed(i))(a)
  }
}

*/


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
     chi2.inverseCumulativeProbability(prob)
  }
  def gammaInverseCDF(prob:Double,alpha:Double,beta:Double)=chiSquareInverseCDF(prob,2.0*(alpha))/(2.0*(beta))


  def apply(shape:Double):Array[Double]={
    cache.getOrElseUpdate((numCat,shape),gamma(shape))
  }
  def gamma(shape:Double):Array[Double]={
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


trait AlternateModel[A <: BioEnum] extends TraitModel[A]{
  val nodeIDs:Set[Int]
  val altModel:Model[A]
  private lazy val startParamNum=super.getParams.length
  private lazy val endParamNum=startParamNum + altModel.getParams.length
  override def qMat(n:Node[A])={
    if (nodeIDs contains n.id){
      debug{"Using altmodel as " + nodeIDs + " contains " + n.id}
      altModel.qMat(n)
    }else {
      debug{"Using regular model as " + nodeIDs + " does not contain " + n.id}
      super.qMat(n)
    }
  }
  override def getParams:List[Array[Double]]={
    //default to last parameter of altModel
    super.getParams ++ altModel.getParams
  }

  override def getParamName(i:Int)={
    if (i>=startParamNum && i < endParamNum){"Alternate Model " + altModel.getParamName(i-startParamNum)}
    else {super.getParamName(i)}

  }
  override def setParams(i:Int)(a:Array[Double])={
    if (i>=startParamNum && i < endParamNum){
      altModel.setParams(i-startParamNum)(a)
    }else {
      super.setParams(i)(a)
    }
  }
  override def cromulent = super.cromulent && altModel.cromulent
  override def toString = {
    super.toString + "\n" + 
    "alt: {\n " + altModel.toString + "\n"+
    "} "
  }
}

class GammaModel[A <: BioEnum](val pi:Vector,val sMat:Matrix,var tree:Tree[A],val alpha:Array[Double]) extends SiteClassPiModel[A] with GammaSMat[A]
trait SiteClassSubstitutionsScaled[A <: BioEnum] extends SiteClassSubstitutions[A]{
  private val paramNum=super.getParams.length

  override def getParams:List[Array[Double]]={
    super.getParams ++ List(Array(norm))
  }

  override def getParamName(i:Int)={
    if (i==paramNum){"SC Substitution Scale "}
    else {super.getParamName(i)}
  }

  
  override def setParams(i:Int)(a:Array[Double]){
    if (i==paramNum){norm=a(0)}
    else{super.setParams(i)(a)}
  }

  override def cromulent=super.cromulent && norm>0.0D
}

trait SiteClassSubstitutions[A <: BioEnum] extends GammaSMat[A]{ //do not use directly, use a descenent
  val rateChangeS:Matrix
  val priors:Vector

  var norm=1.0D
  
  def getSCQMat(node:Node[A])={
    val scQ = rateChangeS.sToQ(priors)
    scQ.normalize(priors,norm)
  }

  override def qMat(node:Node[A])={
   val qStart=super.qMat(node) // we assume a nice, scaled to 1 qMat
   val scQ = getSCQMat(node)
   val numClasses = tree.alphabet.numClasses
    val numAlpha = tree.alphabet.numAlpha
    (0 until numClasses).foreach{i:Int=>
      (0 until numClasses).foreach{j:Int =>
        (0 until numAlpha).foreach{k:Int =>
          if (i!=j){
            val realI = i * numAlpha + k
            val realJ = j * numAlpha + k

            qStart(realI,realJ)=scQ(i,j) * pi(k)
          }
        }
      }
    }
    qStart.fixDiag
  }

  private val paramNum=super.getParams.length

  override def getParams:List[Array[Double]]={
    super.getParams ++ List(linearSMat(rateChangeS).toArray)
  }
  
  override def getParamName(i:Int)={
    if (i==paramNum){"SC Substitution SMat "}
    else {super.getParamName(i)}
  }
  override def setParams(i:Int)(a:Array[Double]){
    if (i==paramNum){setSMat(a,rateChangeS)}
    else{super.setParams(i)(a)}
  }

  override def cromulent=super.cromulent && !(rateChangeS.exists{i=>i<0})
}
trait SiteClassSubstitutionsSeparateBranchLength[A <: BioEnum] extends SiteClassSubstitutions[A]{
  private val paramNum=super.getParams.length
  
  var bl2:Array[Double] = tree.descendentNodes.map{_.lengthTo}.toArray // initialize to branch lengths
  var nodeMap = recalculateNodeMap
  def recalculateNodeMap = tree.descendentNodes.zip(bl2.toList).foldLeft(Map[Node[A],Double]()){_+_}
  override def getParams:List[Array[Double]]={
    super.getParams ++ List(bl2)
  }
  
  override def getParamName(i:Int)={
    if (i==paramNum){"SC Substitution Per Branch Rate "}
    else {super.getParamName(i)}
  }

  override def setParams(i:Int)(a:Array[Double]){
    if (i==paramNum){bl2=a;nodeMap = recalculateNodeMap}
    else{super.setParams(i)(a)}
  }

  override def getSCQMat(node:Node[A])={
    val scQ = rateChangeS.sToQ(priors)
    scQ.normalize(priors,nodeMap(node)/node.lengthTo) // so if the parameter is set to double the 'real' branch length then the sc changes take place at twice the rate
  }

  override def cromulent=super.cromulent && bl2.foldLeft(true){_ && _>0}
}

class THMM[A <: BioEnum](pi:Vector, sMat:Matrix,tree:Tree[A],alpha:Array[Double], val rateChangeS:Matrix) extends GammaModel[A](pi,sMat,tree,alpha) with SiteClassSubstitutionsScaled[A]{
  val priors = Vector((for (i <- 0 until tree.alphabet.numClasses) yield 1.0D/tree.alphabet.numClasses).toList)
}

class THMM2[A <: BioEnum](pi:Vector, sMat:Matrix,tree:Tree[A],alpha:Array[Double],val rateChangeS:Matrix)  extends GammaModel[A](pi,sMat,tree,alpha) with SiteClassSubstitutionsSeparateBranchLength[A]{
  val priors = Vector((for (i <- 0 until tree.alphabet.numClasses) yield 1.0D/tree.alphabet.numClasses).toList)
  }

trait GammaSMat [A <: BioEnum] extends SMat[A]{

  val alpha:Array[Double]
  val alphabet = tree.alphabet
  val numClasses = alphabet.numClasses
  private val paramNum = super.getParams.length
  val gamma = new Gamma(numClasses)

  override def qMat(node:Node[A]) = {
    if (cachedQMat.isEmpty){
      cachedQMat=Some(gammaQMat)
    }
    cachedQMat.get
  }

  var cachedQMat:Option[Matrix]=None

  def gammaQMat={
   // println("HELLO")
   // println("PI: " + pi)
    debug{"Gamma matrix with gamma value " + alpha(0)}
    val gammaVals = gamma(alpha(0)).toList
    val piCat = pi.toArray
    (0 to pi.size-1).foreach{i=> piCat(i)/=numClasses}
    val qStart = Matrix(pi.size * numClasses, pi.size * numClasses)
    (0 to numClasses-1).toList.zip(gammaVals).foreach{t=>
      val (i,rate)=t
      val subQ=sMat.sToQ(piCat).normalize(pi,rate)
      val part = qStart.viewPart(i * pi.size,i*pi.size,pi.size,pi.size)
      part.assign(subQ)
    }
    if (qStart exists {d=>d.isNaN}){
      println("ERROR ")
      println(getParams.map{_.toList}.mkString(","))
      println("produces sMat " + sMat)
      println("and pi " + piVals(tree))
      println("gamma Rates " + gammaVals)
      println("qMat " + qStart)
    }
    qStart
  }
  
  override def getParams:List[Array[Double]]={
    super.getParams ++ List(alpha)
  }

  override def cromulent  = {
    super.cromulent && alpha(0) > 0.0D && alpha(0) < 200.0D && gamma(alpha(0)).toList.find{_.isNaN}.isEmpty
    //no point having alpha go higher - count as infinity
  }

  override def getParamName(i:Int)={
    if (i==paramNum){"GammaSMat"}
    else {super.getParamName(i)}
  }
  override def setParams(i:Int)(a:Array[Double])={
    cachedQMat=None//wipe cache
    debug{"Setting main alpha => " + a(0)}
    if (i != paramNum){super.setParams(i)(a)}
    else {alpha(0) = a(0)}
  }

  override def toString = {
    super.toString + "\n"+
    "shape: " + alpha(0)  
  }
}

trait SiteClassPiModel[A <: BioEnum] extends FullPiModel[A]{
  override def piVals = {
    val numClasses = tree.alphabet.numClasses
    val piStart = Vector(pi.size * numClasses)
    val piCat = pi.toArray
  //  println(numClasses + " " + piStart + " " + pi)
    (0 to piCat.length-1).foreach{i=> piCat(i)/=numClasses}
    (0 to numClasses-1).foreach{i=> 
    //  println("piStart.viewPart("+i+" * "+pi.size+","+pi.size+").assign("+piCat+")")
      piStart.viewPart(i * pi.size,pi.size).assign(piCat)}
    piStart
  }
}

trait FullPiModel[A <: BioEnum] extends TraitModel[A]{
  val pi:Vector
  final val piParam = super.getParams.length

  override def toString = {
    super.toString + "\n"+
    "Pi: " + pi.toString
  }
  var medianIndex:Int=0
  override def cromulent = super.cromulent &&  (pi.zSum < 1.00000001) && (pi.zSum > 0.99999999)

  override def setParams(i:Int)(a:Array[Double])={
    if (i==piParam){setPi(a)}
    else{throw new UnknownParamException(i)}
  }

  override def getParams:List[Array[Double]]={
    super.getParams ++ List(toFit.toArray)
  }

  override def getParamName(i:Int)={
    if (i==piParam){"Pi Values"}
    else {super.getParamName(i)}
  }

  def setPi(array:Array[Double]){
    val exponentiated =  array.map{i=>Math.exp(i)}
    val total = (0.0D /: exponentiated){_+_} + Math.exp(0.0D)
    val values = ((0 to medianIndex-1 ).map{i=> exponentiated(i)/total}.toList ++ List(Math.exp(0.0D)/total) ++ (medianIndex to array.length-1).map{i=> exponentiated(i)/total}).toArray
    pi assign values
  }


  def toFit:List[Double]={                                                                   
    val t  = (pi.toList.zipWithIndex.toList.sort{_._1<_._1})(pi.size/2)
    medianIndex = t._2
    def toFitness={i:Int=>Math.log(pi(i)/pi(medianIndex))}
    (0 to medianIndex-1).map{toFitness}.toList ++ (medianIndex+1 to pi.size-1).map{toFitness}.toList
  }
}

class EnhancedModel[A <: BioEnum](val pi:Vector,val sMat:Matrix, var tree:Tree[A]) extends TraitModel[A] with FullPiModel[A] with SMat[A]

trait SMat[A <: BioEnum] extends TraitModel[A]{
  val sMat:Matrix
  final val sMatParam=super.getParams.length

  override def toString = {
    super.toString + "\n" + 
    "S: " + sMat.toString 
  }

  import cern.colt.function.DoubleFunction

  override def cromulent = super.cromulent && !(sMat exists {i=> i < 0.0D}) 

  def linearSMat:Seq[Double]=linearSMat(sMat)
 
  override def getParamName(i:Int)={
    if (i==sMatParam){"S Values"}
    else {super.getParamName(i)}
  }

  override def getParams:List[Array[Double]]={
    super.getParams ++ List(linearSMat.toArray)
  }
  
  override def setParams(i:Int)(a:Array[Double]){
    if (i==sMatParam){setSMat(a)}
    else{super.setParams(i)(a)}
  }
  def setSMat(array:Array[Double]){setSMat(array,sMat)}

}



