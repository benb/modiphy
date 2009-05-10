package org.modiphy.math
import org.modiphy.tree._
import org.modiphy.sequence._
import EnhancedMatrix._
import org.apache.commons.math.analysis._
import org.apache.commons.math.optimization.GoalType._
import org.apache.commons.math.optimization._
import org.apache.commons.math.optimization.direct._
import cern.colt.function.DoubleFunction
import DataParse._
import tlf.Logging


class MultFunction(f:Array[Double]=>Double) extends MultivariateRealFunction{
  def value(point:Array[Double])=f(point)
}

object ModelOptimiser extends Logging{

  def optimise[A <: BioEnum](optFactory: => MultivariateRealOptimizer)(model:Model[A])={
    var startLkl = model.logLikelihood
    var newLkl=startLkl
    do {
      startLkl = newLkl
      val startModels = model.getParams
      info{"START OPT" + model + "\n" + model.logLikelihood}
      startModels.zipWithIndex.foreach{t=> 
        val (start,index)=t
        val result = optFactory.optimize(new MultFunction({ d:Array[Double]=>
            model.setParams(index)(d)
            val lkl = model.logLikelihood
            debug{"f" + index +": " + d.toList.mkString(",") + " => " + lkl}
            if (lkl.isNaN){Math.NEG_INF_DOUBLE}else{ lkl}
        }
        ),MAXIMIZE,start.toArray)
        model.setParams(index)(result.getPoint)
        info{"OPT " + index + "\n" +  model + "\n" + model.logLikelihood}
        finest{model.likelihoods.mkString(" ")}
        finest{model.realLikelihoods.mkString(" ")}
        newLkl = result.getValue
        true
      }
    } while (newLkl - startLkl > 0.03)
    model
  }
  def nelderMead[A <: BioEnum](model:Model[A])=optimise[A](getNelderMead)(model)
  def multiDirectional[A <: BioEnum](model:Model[A])=optimise[A](getMultiDirectional)(model)
  object DefaultConvergence extends RealConvergenceChecker{
    import   org.apache.commons.math.optimization.RealPointValuePair
    def converged(iter:Int,oldPt:RealPointValuePair,newPt:RealPointValuePair)={iter > 1000 || newPt.getValue - oldPt.getValue < 0.01} 
  }
  def getNelderMead={
    val opt = new NelderMead
    opt.setConvergenceChecker(DefaultConvergence)
    opt
    }

  def getMultiDirectional={
    val opt = new MultiDirectional
    opt.setConvergenceChecker(DefaultConvergence)
    opt

    }

}


