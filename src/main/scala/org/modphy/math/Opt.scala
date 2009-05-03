package org.modphy.math
import org.modphy.tree._
import org.modphy.sequence._
import EnhancedMatrix._
import org.apache.commons.math.analysis._
import org.apache.commons.math.optimization.GoalType._
import org.apache.commons.math.optimization._
import org.apache.commons.math.optimization.direct._
import cern.colt.function.DoubleFunction


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

  def piMapper(original:Vector)(array:Array[Double]):Vector={
    val newVec = original.like 
    println("oldvec " + newVec)
    newVec.assign(array)
    println("newvec " + newVec)
    val total = newVec.zSum
    newVec.assign(new DoubleFunction(){def apply(x:Double)=x/total})
    newVec
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
    (mapper(result.getPoint),result.getValue,result.getPoint)
  }
  def optPi[A <: BioEnum](opt:MultivariateRealOptimizer)(start:Seq[Double],sMat:Matrix,mapper:Array[Double]=>Vector,tree:INode[A])={
    def func(point:Array[Double]) = {
      if (point.findIndexOf{i=> i <= 0.0D} > -1){Math.NEG_INF_DOUBLE}
      else{
        val pi = mapper(point)
        val ans = tree.mkLkl(sMat.sToQ(pi).normalize,pi).logLikelihood
        //println("f("+ point.mkString(",") + ") = " + ans)
        //if (ans.isNaN){throw new org.apache.commons.math.FunctionEvaluationException(point)}
        if (ans.isNaN){Math.NEG_INF_DOUBLE}else{ans}
      }
    }
    val result = opt.optimize(new MultFunction(func),MAXIMIZE,start.toArray)
    (mapper(result.getPoint),result.getValue,result.getPoint)
  }

  def optModel[A <: BioEnum](optFactory: => MultivariateRealOptimizer)(start:(Seq[Double],Seq[Double]),model:(Vector,Matrix),mapper:(Array[Double]=>Vector,Array[Double]=>Matrix),tree:INode[A])={
    var startPi = start._1.toArray
    var startMat=start._2.toArray
    def getPi=mapper._1(startPi)
    def getQ=mapper._2(startMat).sToQ(getPi).normalize
    println("Start Q " + getQ)
    var newLkl=tree.mkLkl((getQ,getPi)).logLikelihood
    var startLkl=Math.NEG_INF_DOUBLE
    val cutoff=0.01
   while (newLkl - startLkl > cutoff){
      startLkl = newLkl
      println("Starting " + getPi + "\n" + getQ + "\n" + newLkl)
      val ans = optPi(optFactory)(startPi,getQ,mapper._1,tree)
      startPi = ans._3
      println("Opt " + getPi + "\n" + getQ + "\n" + ans._2)
      val ans2 = optMat(optFactory)(startMat,getPi,mapper._2,tree)
      startMat=ans._3
      newLkl=ans2._2
      println("Opt " + getPi + "\n" + getQ + "\n" + ans2._2)
    }
    (mapper._1(startPi.toArray),mapper._2(startMat.toArray),newLkl,startPi,startMat)
  }

  


  def optMatNelderMead[A <: BioEnum] = optMat[A](new NelderMead)_

}
