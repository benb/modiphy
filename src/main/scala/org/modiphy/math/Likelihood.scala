package org.modiphy.math
import EnhancedMatrix._
import org.modiphy.tree._
import org.modiphy.sequence._
import cern.colt.matrix.DoubleMatrix1D

class LikelihoodCalc[A <: BioEnum](tree:Node[A],model:Matrix,alphabet:A){
  def count = alphabet.filter{i=>alphabet.isReal(i)}.foldLeft(0){(i,j)=>i+1}

  def likelihood(i:Int)={
        1.0  
  }
}

object BasicLikelihoodCalc{
  val func = new cern.colt.function.DoubleDoubleFunction{def apply(x:Double,y:Double)=x*y}

  def partialLikelihoodCalc(end:List[Vector],matrix:Matrix)={
    val width = matrix.rows
    val rows = new Array[DoubleMatrix1D](width)
    (0 until width).foreach{i=> rows(i)=matrix viewRow i}
    end.map{siteVector=>
        val ret = Vector(width)

        for (i<-0 until width){
          ret set (i,siteVector.zDotProduct(rows(i)))
        }
        ret
      }
  }
  def combinePartialLikelihoods(intermediates:List[List[Vector]])={
    val ans = intermediates.head
    intermediates.tail.foreach{list2=>
      ans.zip(list2).map{t=> // not really a map but used for parallel reasons
        val (vec,vec2)=t
        vec.assign(vec2,func)
      }
    }
    ans

  }
}
