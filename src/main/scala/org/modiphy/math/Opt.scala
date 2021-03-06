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

import dr.math.{MultivariateFunction,MultivariateMinimum,UnivariateMinimum,UnivariateFunction}


abstract class Gradient extends MultivariateVectorialFunction{
  def apply(point:Array[Double]):Array[Double]
  def value(point:Array[Double])=apply(point)
}
class FuncWrapper(model:ActorModel[_],p:OptPSetter) extends MultivariateFunction with MultivariateRealFunction with Logging{
  def setBest=apply(best)
  var best:Array[Double]=p.latestArgs
  var bestLnL = -1E100
  var count = 0
  val length = latestArgs.length
  val lower = (0 until length).map{i=>p.lower(i)}.toArray
  val upper = (0 until length).map{i=>p.upper(i)}.toArray
  var logT = System.currentTimeMillis()
  var logT1 = logT
  def getUpperBound(i:Int)=upper(i)
  def getLowerBound(i:Int)=lower(i)
  //def this(model:ActorModel,p:OptPSetter)=this(model,p,{s:String=>})
  def apply(point:Array[Double])={
    val t1 = System.currentTimeMillis
    finest{"out " + (t1 - logT1)}
    p(point)
    val ans = model.logLikelihood
    val t2 = System.currentTimeMillis
    if (ans > bestLnL){
      best = point.clone
      finest{"NEW BEST" + best + " " + ans}
      bestLnL = ans
    }
    if (count % 100==0){
      val newT = System.currentTimeMillis()
      info{count + "t: " + (newT - logT) + " f: " + p + " " + ans}
      info{"tSingle: " + (t2 - t1)}
      logT = newT
    }else {
      extra{count + " f: " + p + " " + ans}
    }
    count=count+1
    logT1 = System.currentTimeMillis
    finest{"in " +  (logT1 - t1)}
    ans
  }
  def evaluate(point:Array[Double])= -apply(point)
  def inBounds(point:Array[Double])={
    point.zip(lower).foldLeft(true){(bool,t) => bool && t._1 >= t._2} && 
    point.zip(upper).foldLeft(true){(bool,t) => bool && t._1 <= t._2} 

  }
  def value(point:Array[Double]) = {
    if (inBounds(point)){
      apply(point)
    }else {
      if (inBounds(point.map{_.abs})){
        apply(point.map{_.abs})
      }else {
      extra{"f: " + point.mkString(" ") + " out of bounds"}
      extra{"lower " + lower.mkString(" ")}
      extra{"upper " + upper.mkString(" ")}

      -1e100
    }
  }
}
  def getNumArguments=p.numArguments
  def latestArgs=p.latestArgs
}

object ModelOptimiser extends Logging{

  
  def optimise[A <: BioEnum](optFactory: => MultivariateMinimum,p:ParamName,model:ActorModel[_]):Double={
    optimise(optFactory,p::Nil,model)
  }
  def optimise[A <: BioEnum](optFactory: => MultivariateMinimum,pList:List[ParamName],model:ActorModel[_]):Double={
    optimise(optFactory,pList,model,1E-4,1E-3)
  }
  def optimise[A <: BioEnum](optFactory: => MultivariateMinimum,pList:List[ParamName],model:ActorModel[_],tolfx:Double,tolx:Double):Double={
    val startParams = model.optSetter(pList)
      val func = new FuncWrapper(model,startParams)
      if (func.length==1){
        getNelderMead.optimize(func,MAXIMIZE,func.latestArgs)
      }else {
        optFactory.optimize(func,func.latestArgs,tolfx,tolx)
      }
      func.setBest 
  }
  def optimiseSequential[A <: BioEnum](optFactory: => MultivariateMinimum,pListList:List[List[ParamName]],model:ActorModel[_]):Double={
     var start = model.logLikelihood
     var end = start
     do {
        start=end
        pListList.foreach{l:List[ParamName]=>
          end=optimise(optFactory,l,model)
        }
     } while (end - start > 0.02)
     end
  }

  object DefaultConvergence extends RealConvergenceChecker{
    import   org.apache.commons.math.optimization.RealPointValuePair
    def converged(iter:Int,oldPt:RealPointValuePair,newPt:RealPointValuePair)={iter > 30000 || newPt.getValue - oldPt.getValue < 0.01} 
  }
  def getNelderMead={
    val opt = new NelderMead
    opt.setConvergenceChecker(DefaultConvergence)
    opt
  }
  def getConjugateDirection={
    new dr.math.ConjugateDirectionSearch
  }
  def getConjugateGradient={
    new dr.math.ConjugateGradientSearch
  }

  def getMultiDirectional={
    val opt = new MultiDirectional
    opt.setConvergenceChecker(DefaultConvergence)
    opt

  }
}
