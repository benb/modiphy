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
  def apply[A<:BioEnum](matMap:MatrixPi[A],myNode:Node[A]):List[Vector]={
   val partialLikelihoods:List[List[Vector]]=myNode.children.map{_.likelihoods(matMap)}
    val intermediates=partialLikelihoods.zip(myNode.children.toList).map{t=>
      val (siteVectorList,c)=t // list of vectors  1 for each site
      val matrix = matMap.getMat(c) //e^Qt
      partialLikelihoodCalc(siteVectorList,matrix)
    }.toList
    combinePartialLikelihoods(intermediates)
  }

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
      ans.zip(list2).map{t=> // not really a map but used for parallel reasons
        val (vec,vec2)=t
        vec.assign(vec2,func)
      }
    }
    ans

  }
}
