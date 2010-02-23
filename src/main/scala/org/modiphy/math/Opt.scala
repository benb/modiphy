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
class FuncWrapper(model:ActorModel,p:OptPSetter) extends MultivariateFunction with MultivariateRealFunction with Logging{
  var count = 0
  val length = latestArgs.length
  val lower = (0 until length).map{i=>p.lower(i)}.toArray
  val upper = (0 until length).map{i=>p.upper(i)}.toArray
  def getUpperBound(i:Int)=upper(i)
  def getLowerBound(i:Int)=lower(i)
  //def this(model:ActorModel,p:OptPSetter)=this(model,p,{s:String=>})
  def apply(point:Array[Double])={
    p(point)
    val ans = model.logLikelihood
    if (count % 100==0){
      info{count + " f: " + p + " " + ans}
    }else {
      extra{count + " f: " + p + " " + ans}
    }
    count=count+1
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
      extra{"f: " + point.mkString(" ") + " out of bounds"}
      extra{"lower " + lower.mkString(" ")}
      extra{"upper " + upper.mkString(" ")}

      -1e100
    }
  }
  def getNumArguments=p.numArguments
  def latestArgs=p.latestArgs
}

class MultFunction(f:Array[Double]=>Double,val numArg:Int,p:ParamControl,badVal:Double) extends DifferentiableMultivariateRealFunction with MultivariateFunction{
 def this(f:Array[Double]=>Double,numArg:Int,p:ParamControl)=this(f,numArg,p,-1E100)

  def makeCromulent(point:Array[Double]){
    p.setParams(point)
    p.makeCromulent
  }
  def apply(point:Array[Double])=value(point)
  def value(point:Array[Double])={val ans = f(point); if (ans.isNaN){badVal}else{ans}}
  def evaluate(point:Array[Double])=value(point)
  def partialDerivative(k:Int)={
    val outer = this
    new MultivariateRealFunction{
      def value(point:Array[Double])={
        val h = 1E-4
        val hplus = point.toArray
        val hminus = point.toArray
        hplus(k)=hplus(k)+h
        hminus(k)=hminus(k)-h
        (outer.value(hplus)-outer.value(hminus))/h/2.0
      }
    }
  }
  def gradient():Gradient={
    new Gradient{
      def apply(point:Array[Double])={
        for (k<-0 until point.length) yield partialDerivative(k).value(point)
      }.toArray
    }
  }
  def gradient(point:Array[Double]):Array[Double]=gradient()(point)
  def getNumArguments=numArg
  def getLowerBound(i:Int)=p.lowerBound(i)
  def getUpperBound(i:Int)=p.upperBound(i)
  def negative = new MultFunction({a:Array[Double]=> -f(a)},numArg,p,-badVal)
}
class SingleFunction(m:MultFunction) extends UnivariateFunction{
  def apply(x:Double)=evaluate(x)
  def evaluate(x:Double)={val ans = m.evaluate(Array(x)); if (ans < 1E-100 || ans==Math.NEG_INF_DOUBLE){1E-100}else{ans}}//hack to make BEAST behave
  def getLowerBound=m.getLowerBound(0)
  def getUpperBound=m.getUpperBound(0)
  
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
  val lowerB = p.toList.map{i=>i.lowerBounds.toList}.flatten[Double].toArray
  val upperB = p.toList.map{i=>i.upperBounds.toList}.flatten[Double].toArray
  override def lowerBound(i:Int)=lowerB(i)
  override def upperBound(i:Int)=upperB(i)
}



class RestrictedParam(p:ParamControl,f:Array[Double]=>Array[Double],initialParam:Array[Double]) extends ParamControl{
  var currParam = initialParam.toArray
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

  type Optimizer ={
    def optimize(f:MultFunction,initial:Array[Double]):RealPointValuePair
  }
  abstract class OptimWrapper{
    def optimize(f:MultFunction,initial:Array[Double]):RealPointValuePair={
      f.numArg match {
        case 1=>uniOpt.optimize(f,initial)
        case _=>optimizeImpl(f,initial)
      }
    }
    lazy val uniOpt = new OptimWrapperApache(new NelderMead)
    def optimizeImpl(f:MultFunction,initial:Array[Double]):RealPointValuePair
/*    def optimizeImpl(f:SingleFunction,initial:Double):RealPointValuePair = {
      println("BOUND " + f.getLowerBound + " " + f.getUpperBound)
      println("1.0D " + f(1.0D))
      println("-1.0D " + f(-1.0D))
      val ans = uniOpt.optimize(f,1E-12,f.getLowerBound,f.getUpperBound)
      new RealPointValuePair(Array(ans),f(ans))
    }*/
  }
  class OptimWrapperApache(o:MultivariateRealOptimizer) {
    def optimize(f:MultFunction,initial:Array[Double]):RealPointValuePair = 
      o.optimize(f,MAXIMIZE,initial)
  }
  class OptimWrapperBeast(o:MultivariateMinimum) extends OptimWrapper{
    def optimizeImpl(f:MultFunction,initial:Array[Double]):RealPointValuePair = {
      val point = Array(initial :_*)
      f.makeCromulent(point)

      import sbt.TrapExit
      val ans = TrapExit({ o.optimize(f.negative,point,0E-2,1E-6)},new sbt.ConsoleLogger)
      if (ans == 0){
        println("SUCCESS! " + ans)
        new RealPointValuePair(point,f.value(point))
      }else {
        val point2 = Array(initial :_*)
        println("Trying " + point2.mkString(","))

        (new OptimWrapperApache(new NelderMead)).optimize(f,point2)
        f(point2)
        new RealPointValuePair(point2,f.value(point2))
      }
    }
  }
  implicit def RealPointValuePairToTuple(r:RealPointValuePair)=(r.getPoint,r.getValue)
  implicit def WrapOpt(o:MultivariateRealOptimizer)=new OptimWrapperApache(o)
  implicit def WrapOpt(o:MultivariateMinimum)=new OptimWrapperBeast(o)

  
  def optimise[A <: BioEnum](paramList:Seq[ParamControl],optFactory: => Optimizer)(model:ComposeModel[A])={
    var startLkl = model.logLikelihood
    var newLkl=startLkl
    val startModels = paramList

    do {
      startLkl = newLkl
      info{"START OPT " + model + " \n" + model.logLikelihood}

      startModels.foreach{start=> 
        val startP=start.getParams
        val result = optFactory.optimize(new MultFunction({ d:Array[Double]=>
            start.setParams(d)
            val lkl = model.logLikelihood
            extra{"f(" + start.name +"): " + d.toList.mkString(",") + " => " + lkl}
            if (lkl.isNaN){Math.NEG_INF_DOUBLE}else{ lkl}
        },start.getParams.length,start
        ),startP.toArray)
        start.setParams(result._1)
        println("OPT " + start.name + " " + model.logLikelihood)
        println(model)
        finest{model.likelihoods.mkString(" ")}
        finest{model.realLikelihoods.mkString(" ")}
        newLkl = result._2
        true
      }

    } while (newLkl - startLkl > 0.03 && !(startModels.drop(1).isEmpty))//only loop if there is >1 model
    model
  }
  def optimise[A <: BioEnum](optFactory: => MultivariateMinimum,p:ParamName,model:ActorModel):Double={
    optimise(optFactory,p::Nil,model)
  }
  def optimise[A <: BioEnum](optFactory: => MultivariateMinimum,pList:List[ParamName],model:ActorModel):Double={
      val startParams = model.optSetter(pList)
      val func = new FuncWrapper(model,startParams)
      if (func.length==1){
        getNelderMead.optimize(func,MAXIMIZE,func.latestArgs)
      }else {
        optFactory.optimize(func,func.latestArgs,1E-4,1E-3)
      }
      model.logLikelihood
  }
  def optimiseSequential[A <: BioEnum](optFactory: => MultivariateMinimum,pListList:List[List[ParamName]],model:ActorModel):Double={
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




  def nelderMead[A <: BioEnum](paramSet:Seq[ParamControl],model:ComposeModel[A]):ComposeModel[A]=optimise[A](paramSet,getNelderMead)(model)
  def multiDirectional[A <: BioEnum](paramSet:Seq[ParamControl],model:ComposeModel[A]):ComposeModel[A]=optimise[A](paramSet,getMultiDirectional)(model)
  def conjugateGradient[A <: BioEnum](paramSet:Seq[ParamControl],model:ComposeModel[A]):ComposeModel[A]=optimise[A](paramSet,getConjugateGradient)(model)
  def conjugateDirection[A <: BioEnum](paramSet:Seq[ParamControl],model:ComposeModel[A]):ComposeModel[A]=optimise[A](paramSet,getConjugateDirection)(model)



  def nelderMead[A <: BioEnum](model:ComposeModel[A]):ComposeModel[A]=nelderMead(model.params,model)
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

/*import ModelOptimiser._
object NewtonRaphson{
  def iterate(f:DifferentiableMultivariateRealFunction,x:Array[Double])={
    val fx = f.value(x)
    x.zip(f.gradient.value(x)).map{t=> 
      t._1-fx/t._2
    }
  }
  
  def optimize(f:DifferentiableMultivariateRealFunction,direction:GoalType,initial:Array[Double],iter:Int,tol:Double,maxIter:Int):RealPointValuePair={
    if (iter>maxIter){
      new RealPointValuePair(initial,f.value(initial)) 
    }else{
      val newLoc = iterate(f,initial)
      if (
        newLoc.zip(initial).foldLeft(true){(ans,t)=>
          ans && (t._1-t._2).abs < tol
        }
      ){
        new RealPointValuePair(initial,f.value(initial))
      }else {
        optimize(f,direction,newLoc,iter+1,tol,maxIter)
      }
    }
  }
  def optimize(f:DifferentiableMultivariateRealFunction,direction:GoalType,initial:Array[Double]):RealPointValuePair=optimize(f,direction,initial,0,1E-7,30000)

}
*/





