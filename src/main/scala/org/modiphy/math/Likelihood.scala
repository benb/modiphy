package org.modiphy.math
import EnhancedMatrix._
import org.modiphy.tree._
import org.modiphy.sequence._

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
    end.map{siteVector=>
        val ret = Vector(width)

        for (i<-0 until width){
          ret(i)=siteVector.zDotProduct(matrix.viewRow(i))
        }
        ret
      }
  }
  def combinePartialLikelihoods(intermediates:List[List[Vector]])={
    val ans = intermediates.head
    intermediates.tail.foreach{list2=>
      ans.zip(list2).foreach{t=> 
        val (vec,vec2)=t
        vec.assign(vec2,func)
      }
    }
    ans

  }
  def likelihoods(pl:List[Vector],pi:Vector)={
    pl.map{vec=>
     pi.elements.zipWithIndex.map{t=>
        val(p,i)=t
        vec(i)*p
     }.foldLeft(0.0D){_+_}
   }
 }

  def logLikelihood(pl:List[Vector],pi:Vector)={
    likelihoods(pl,pi).foldLeft(0.0D){_+Math.log(_)}
  }
}
