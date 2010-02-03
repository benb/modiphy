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

class JoinedParam(p:Array[ParamControl]) extends ParamControl{
  def this(p2:List[ParamControl])=this(p2.toArray)
  val lengths = p.map{_.getParams.length}
  val pointers = lengths.foldLeft(List(0)){(p,a)=> (a+p.head)::p}.reverse

  val subArrays = p.map{_.getParams.toArray}


  val backing = new Array[Double](lengths.foldLeft(0){_+_})
  def update{
    var pointer = 0
    for (t<-lengths.zip(p)){
      println(t._2.getParams.length + " " + 0 + " " + pointer + " " + t._1 + " " + backing.length)
      Array.copy(t._2.getParams,0,backing,pointer,t._1)
      pointer+=t._1
    }
  }
  def getParams={
    update
    backing.toArray //copy
  }
  def setParams(a:Array[Double])={
    for (t<-pointers.zip(pointers.tail).zip(p.toList)){
      t._2.softSetParams(a.subArray(t._1._1,t._1._2))
    }
  }
  def softSetParams(a:Array[Double])=setParams(a)
  def view=null
  val name = "Joined param: " + p.map{_.name}.mkString(" ")
}
class RestrictedParam(p:ParamControl,f:Array[Double]=>Array[Double],initialParam:Array[Double]) extends ParamControl{
  val currParam = initialParam.toArray
  def getParams = currParam.toArray
  
  private def copyToInternal(a:Array[Double]){
    Array.copy(a,0,currParam,0,a.length)
  }
  
  def setParams(a:Array[Double]){copyToInternal(a); p.setParams(f(a))}
  def softSetParams(a:Array[Double]){copyToInternal(a); p.softSetParams(f(a))}
  def view = null
  def name="Restricted " + p.name
}
  
object RestrictedParamFunctions{
  def oneToMany(num:Int):Array[Double]=>Array[Double]=
    {d:Array[Double] =>
      val x = new Array[Double](num)
      for (i<-0 until x.length){x(i)=d(0)}
      x
    }
  /**
    First array value is 0<->[1,2,3,4..k]
    Second array value is for [1,2,3...k] <->[1,2,3...k]
    Intended for having separate invar<->gamma and between gamma switching rates
  */
  def dualThmm(numSC:Int)={
    val size = numSC*(numSC-1)/2
    val ans = {d:Array[Double] => 
      val x = new Array[Double](size)
      for (i<-0 until numSC -1){x(i)=d(0)}
      for (i<-numSC until size){x(i)=d(1)}
      x
    }
    ans
   }
}

object ModelOptimiser extends Logging{

  
  def optimise[A <: BioEnum](paramList:Seq[ParamControl],optFactory: => MultivariateRealOptimizer)(model:ComposeModel[A])={
    var startLkl = model.logLikelihood
    var newLkl=startLkl

    do {
      startLkl = newLkl
      val startModels = paramList
      info{"START OPT " + model + " \n" + model.logLikelihood}

      startModels.foreach{start=> 
        val startP=start.getParams
        val result = optFactory.optimize(new MultFunction({ d:Array[Double]=>
            start.setParams(d)
            val lkl = model.logLikelihood
            extra{"f(" + start.name +"): " + d.toList.mkString(",") + " => " + lkl}
            if (lkl.isNaN){Math.NEG_INF_DOUBLE}else{ lkl}
        }
        ),MAXIMIZE,startP.toArray)
        start.setParams(result.getPoint)
        println("OPT " + start.name + " " + model.logLikelihood)
        println(model)
        finest{model.likelihoods.mkString(" ")}
        finest{model.realLikelihoods.mkString(" ")}
        newLkl = result.getValue
        true
      }

    } while (newLkl - startLkl > 0.03)
    model
  }
  



  def nelderMead[A <: BioEnum](paramSet:Seq[ParamControl],model:ComposeModel[A]):ComposeModel[A]=optimise[A](paramSet,getNelderMead)(model)
  def multiDirectional[A <: BioEnum](paramSet:Seq[ParamControl],model:ComposeModel[A]):ComposeModel[A]=optimise[A](paramSet,getMultiDirectional)(model)

  def nelderMead[A <: BioEnum](model:ComposeModel[A]):ComposeModel[A]=nelderMead(model.params,model)
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


