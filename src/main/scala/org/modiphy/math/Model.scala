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


class InvalidMatrixException(m:String) extends RuntimeException(m)


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

