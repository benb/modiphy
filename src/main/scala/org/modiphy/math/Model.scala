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

abstract class SingleComponent extends MComponent{
  def apply(node:Node[_]):Double
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


abstract class SMatComponent extends MComponent with SMatUtil{
  def sMat:Matrix
  override def cromulent = !(sMat exists {i=> i < 0.0D}) 
  def linearSMat:Seq[Double]=linearSMat(sMat)
  def apply(node:Node[_]):Matrix
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
        }//else{ //by passing in a shorter array, remaining elements unset 
          //sMatInst(i,j)=1
       // }
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

