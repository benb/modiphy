package org.modphy.math
import org.modphy.tree._
import org.modphy.sequence._
import EnhancedMatrix._
import org.apache.commons.math.analysis._
import org.apache.commons.math.optimization.GoalType._
import org.apache.commons.math.optimization._
import org.apache.commons.math.optimization.direct._

class MultFunction(f:Array[Double]=>Double) extends MultivariateRealFunction{
  def value(point:Array[Double])=f(point)
}

object Optimiser{
  def sMatMapper(original:Matrix)(array:Array[Double]):Matrix={
    val iter:Iterator[Double]=array.elements
    val newmat=original.like
    (0 to newmat.rows-1).foreach{i=>
      (i+1 to newmat.columns-1).foreach{j=>
        if (iter.hasNext){
          newmat(i,j)=iter.next
        }else{ //by passing in a shorter array, remaining elements set to 1
          newmat(i,j)=1
        }
      }
    }
    newmat
  }

  def opt(opt:MultivariateRealOptimizer)(start:Seq[Double],f:Array[Double]=>Double)={
    val result = opt.optimize(new MultFunction(f),MAXIMIZE,start.toArray)
    (result.getPoint,result.getValue)
  }
  def nelderMead=opt(new NelderMead)_

  def optMat[A <: BioEnum](opt:MultivariateRealOptimizer)(start:Seq[Double],pi:Vector,mapper:Array[Double]=>Matrix,tree:INode[A])={
    def func(point:Array[Double]) = {
      if (point.findIndexOf{i=> i < 0.0D} > -1){Math.NEG_INF_DOUBLE}
      else{
        val ans = tree.mkLkl(mapper(point).sToQ(pi).normalize,pi).logLikelihood
        //println("f("+ point.mkString(",") + ") = " + ans)
        //if (ans.isNaN){throw new org.apache.commons.math.FunctionEvaluationException(point)}
        if (ans.isNaN){Math.NEG_INF_DOUBLE}else{ans}
      }
    }
    val result = opt.optimize(new MultFunction(func),MAXIMIZE,start.toArray)
    (mapper(result.getPoint),result.getValue)
  }

  def optMatNelderMead[A <: BioEnum] = optMat[A](new NelderMead)_

}
