package org.modphy.math

import cern.colt.matrix._
import org.modphy.tree._
import org.modphy.sequence._
import org.modphy.math.EnhancedMatrix._
import org.modphy.tree.DataParse._
import scala.collection.immutable.IntMap

trait Model[A <: BioEnum]{
  def logLikelihood:Double
  def getParams:List[Array[Double]] = List()
  def setParams(paramSet:Int)(params:Array[Double]){}
  def cromulent = true
  def tree:Tree[A]
  def likelihoods=tree.mkLkl(this).likelihoods
  def likelihoods(node:CalcLikelihoodNode[A]):List[Vector]
  def realLikelihoods=tree.mkLkl(this).realLikelihoods
  def piVals:Vector
}

trait OptBranchLengths[A <: BioEnum] extends Model[A]{
  var tree:Tree[A]
  val treeParamNum = super.getParams.length
  override def setParams(i:Int)(a:Array[Double])={
    if (i != treeParamNum){super.setParams(i)(a)}
    else {tree = tree.setBranchLengths(a.toList).setRoot}
  }
  override def getParams:List[Array[Double]]={
    super.getParams ++ List(tree.getBranchLengths.toArray)
  }
  override def toString={
    super.toString + "\n" + tree.toString
  }
  override def cromulent=super.cromulent && !(tree.getBranchLengths.exists{_<0.0D})
}
trait OptBranchScale[A <: BioEnum] extends Model[A]{
  var tree:Tree[A]
  val treeScaleParamNum = super.getParams.length
  val originalBranchLengths = tree.getBranchLengths
  override def setParams(i:Int)(a:Array[Double])={
    if (i != treeScaleParamNum){super.setParams(i)(a)}
    else {tree = tree.setBranchLengths(originalBranchLengths.map{i=>i*a(0)}).setRoot}
  }
  def scale = tree.getBranchLengths.head/originalBranchLengths.head 
  override def getParams:List[Array[Double]]={
    super.getParams ++ List(Array(scale))
  }
  override def toString={
    super.toString + "\n" + tree.toString + " (scale: " + scale.toString + " )"
  }
  override def cromulent=super.cromulent && !(tree.getBranchLengths.exists{_<0.0D})
}

trait CombineOpt[A <: BioEnum] extends Model[A]{
  def toCombine:List[Int]
  // a mapping from new indicies to old
  val paramMap=super.getParams.zipWithIndex.filter{t=>val (p,i)=t;!(toCombine.exists{j=>j==i})}.map(_._2).zipWithIndex.foldLeft(IntMap[Int]():Map[Int,Int]){(m,t)=>m+((t._2,t._1))}
  val thisParam = getParams.length-1
  val paramLengths = super.getParams.zipWithIndex.filter{t=>toCombine.exists{j=>j==t._2}}.map{_._1.length}
  override def getParams:List[Array[Double]]={
    val start = super.getParams
    //not super efficient but this should be rarely called
    start.zipWithIndex.filter{t=>val (array,i)=t; !toCombine.exists{j=>j==i}}.map{_._1} ++ List(start.zipWithIndex.filter{t=>toCombine.exists{j=>j==t._2}}.map{_._1.toList}.flatten[Double].toArray)
  }
  override def setParams(i:Int)(a:Array[Double]){
    if (paramMap contains i){
      super.setParams(paramMap(i))(a)
    }else {
      if (i==thisParam){
        var listPtr = a.toList
        paramLengths.zip(toCombine).foreach{t=>
          val (length,index)=t
          super.setParams(index)(listPtr.take(length).toArray)
          listPtr = listPtr drop length
        }
        
      }
      else {
        throw new Exception("Unsupported parameter value " + i)
      }
    }
  }

}



object Gamma{
  val cache = new scala.collection.jcl.WeakHashMap[(Int,Double),Seq[Double]]()
}
class Gamma(numCat:Int){
  import net.sf.doodleproject.numerics4j.statistics.distribution.GammaDistribution

  import Gamma._
  val gamma = new GammaDistribution(1.0,1.0)

    def apply(shape:Double)={
    if (shape < 150.0D){
      cache.getOrElseUpdate((numCat,shape),{
        
         
          gamma setAlpha (shape + 1E-10)
          gamma setBeta (shape + 1E-10)
          //+ 1E-10 seems to avoid convergence exceptions!

          val inter = (0 to numCat-1).map{i => gamma inverseCumulativeProbability (i * 2.0 + 1.0)/(2.0 * numCat)}
          val mean = inter.foldLeft(0.0){_+_}/numCat
          inter.map{_/mean}}
        )
    }else {
      (1 to numCat).map(i=> 1.0D).toList // assume alpha => inf
    }

  }
}


class GammaModel[A <: BioEnum](pi:Vector,sMat:Matrix,tree:Tree[A],shape:Double) extends EnhancedModel(pi,sMat,tree){
  var alpha = shape
  val alphabet = tree.alphabet
  val numClasses = alphabet.numClasses
  val gamma = new Gamma(numClasses)
  val gammaParamNum = super.getParams.length


  override def piVals = {
    val piStart = Vector(pi.size * numClasses)
    val piCat = pi.toArray
    (0 to piCat.length-1).foreach{i=> piCat(i)/=numClasses}
    (0 to piCat.length-1).foreach{i=> piStart.viewPart(i * pi.size,pi.size).assign(piCat)}
    piStart
  }

  override def qMat(node:CalcLikelihoodNode[A]) = {
    val gammaVals = gamma(alpha).toList
    val piCat = pi.toArray
    (0 to pi.size-1).foreach{i=> piCat(i)/=numClasses}
    val qStart = Matrix(pi.size * numClasses, pi.size * numClasses)
    (0 to numClasses-1).toList.zip(gammaVals).foreach{t=>
      val (i,rate)=t
      val part = qStart.viewPart(i * numClasses,i*numClasses,pi.size,pi.size)
      part.assign(sMat.sToQ(piCat).normalize(rate))
    }
    //println(qStart)
    qStart
  }
  
  override def getParams:List[Array[Double]]={
    super.getParams ++ List(Array(alpha))
  }

  override def cromulent  = {
    import net.sf.doodleproject.numerics4j.exception.ConvergenceException 
    val ok = super.cromulent && alpha > 0.0D && alpha < 200.0D 
    val converges:Boolean = try{gamma(alpha);true}catch{case e: Exception => false}
    ok && converges
    //no point having alpha go higher - count as infinity
  }
  override def setParams(i:Int)(a:Array[Double])={
    if (i != gammaParamNum){super.setParams(i)(a)}
    else {alpha = a(0)}
  }

  override def toString = {
    "Pi: " +  pi.toString + "\n"+
    "S: " + sMat.toString + "\n"+
    "shape: " + alpha  
  }
}


class EnhancedModel[A <: BioEnum](val pi:Vector,val sMat:Matrix, var tree:Tree[A]) extends Model[A]{


  var medianIndex:Int=0

  /**
   The values that should be used in the likelihood calc - i.e. numClasses * numChars long
   */
  def piVals = pi
  override def toString = {
    "Pi: " +  pi.toString + "\n"+
    "S: " + sMat.toString 
  }
  def getTree = tree.mkLkl(this)

  import cern.colt.function.DoubleFunction

  override def cromulent = !(sMat exists {i=> i < 0.0D}) && (pi.zSum < 1.000001) && (pi.zSum > 0.9999999)

  def logLikelihood=if (cromulent){tree.mkLkl(this).logLikelihood}else{Math.NEG_INF_DOUBLE}

  def qMat(node:CalcLikelihoodNode[A])=sMat.sToQ(pi).normalize
  def getMat(node:CalcLikelihoodNode[A])={
    qMat(node).exp(node.lengthTo)
  }

  def likelihoods(node:CalcLikelihoodNode[A])={
    val childLkl = node.childElements.map{i:LikelihoodNode[A]=>(i.likelihoods,i.lengthTo)}
    val intermediates= childLkl.map{t=>
    val (siteVectorList,length)=t // list of vectors - 1 for each site
    //println ("siteVectorList " + siteVectorList)
    //println("Q = " + qMat)
    val matrix = getMat(node) //e^Qt
    //println("e^Qt=" + matrix)
    siteVectorList.map{siteVector=>
     val alphabet = node.alphabet
    
     val ret = DoubleFactory1D.dense.make(alphabet.matLength) 
        (0 to alphabet.matLength-1).foreach{i=>
          (0 to alphabet.matLength-1).foreach{j=>
            ret(j)=ret(j) + siteVector(i) * matrix(j,i)
          }
        }
        //println("MAP => " +siteVector + " " + ret)
        ret
      }
    }.toList
    val ans = intermediates.head
    intermediates.tail.foreach{list2=>
    ans.zip(list2).foreach{t=>
          val (vec,vec2)=t
          (0 to vec.size-1).foreach{base=>vec(base)=vec(base)*vec2(base)}
      }
    }
    ans

  }

  def linearSMat:Seq[Double]={
    val list = (0 to sMat.rows -3).map{ i=> // skip last elment...
       val list1 = sMat(i).viewPart(i+1,sMat.columns - i - 1).toList
       list1
    }.toList.flatten[Double]

    list
  }

  override def getParams:List[Array[Double]]={
    toFit.toArray :: linearSMat.toArray :: List()
  }
  
  override def setParams(i:Int)(a:Array[Double]){
   i match {
    case 0 => setPi(a)
    case 1 => setSMat(a)
    case _ => throw new Exception("Param index out of bounds")
   }

  }

 def setSMat(array:Array[Double]){
    val iter:Iterator[Double]=array.elements
    (0 to sMat.rows-1).foreach{i=>
      (i+1 to sMat.columns-1).foreach{j=>
        if (iter.hasNext){
          sMat(i,j)=iter.next
        }else{ //by passing in a shorter array, remaining elements set to 1
          sMat(i,j)=1
        }
      }
    }
    sMat
  }

  def setPi(array:Array[Double]){
    val exponentiated =  array.map{i=>Math.exp(i)}
    val total = (0.0D /: exponentiated){_+_} + Math.exp(0.0D)
    val values = ((0 to medianIndex-1 ).map{i=> exponentiated(i)/total}.toList ++ List(Math.exp(0.0D)/total) ++ (medianIndex to array.length-1).map{i=> exponentiated(i)/total}).toArray
    pi assign values
  }


  def toFit={                                                                   
    val t  = (pi.toList.zipWithIndex.toList.sort{_._1<_._1})(pi.size/2)
    medianIndex = t._2
    def toFitness={i:Int=>Math.log(pi(i)/pi(medianIndex))}
    (0 to medianIndex-1).map{toFitness}.toList ++ (medianIndex+1 to pi.size-1).map{toFitness}.toList
  }



}


