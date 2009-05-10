package org.modphy.math

import cern.colt.matrix._
import org.modphy.tree._
import org.modphy.sequence._
import org.modphy.math.EnhancedMatrix._
import org.modphy.tree.DataParse._
import scala.collection.immutable.IntMap
import tlf.Logging

class UnknownParamException(i:Int) extends java.lang.RuntimeException("Unknown parameter " + i.toString)
trait Model[A <: BioEnum] extends Logging{
  def logLikelihood=if (cromulent){tree.mkLkl(this).logLikelihood}else{Math.NEG_INF_DOUBLE}
  def getParams:List[Array[Double]] = List()
  def setParams(paramSet:Int)(params:Array[Double]){}
  def cromulent = true
  var tree:Tree[A]
  def qMat(node:CalcLikelihoodNode[A])={
    sMat.sToQ(pi).normalize

  }

  def getMat(node:CalcLikelihoodNode[A])={ qMat(node).exp(node.lengthTo) }
  def likelihoods(node:CalcLikelihoodNode[A])={
    val childLkl = node.childElements.map{i:LikelihoodNode[A]=>(i.likelihoods,i.lengthTo)}
    val intermediates= childLkl.map{t=>
    val (siteVectorList,length)=t // list of vectors  1 for each site
    //println ("siteVectorList " + siteVectorList)
    //println("Q = " + qMat)
    val matrix = getMat(node) //e^Qt
    //println("e^Qt=" + matrix)
    siteVectorList.map{siteVector=>
      val alphabet = node.alphabet

      val ret = DoubleFactory1D.dense.make(alphabet.matLength) 
        (0 to alphabet.matLength - 1).foreach{i=>
          (0 to alphabet.matLength - 1).foreach{j=>
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
        (0 to vec.size - 1).foreach{base=>vec(base)=vec(base)*vec2(base)}
      }
    }
    ans

  }

  def likelihoods=tree.mkLkl(this).likelihoods
  def realLikelihoods=tree.mkLkl(this).realLikelihoods

  def pi:Vector
  def sMat:Matrix
  def piVals=pi
  def getTree = tree.mkLkl(this)
  override def toString ="Model:\n:"
  
}

trait OptBranchLengths[A <: BioEnum] extends Model[A]{
  var tree:Tree[A]
  private val paramNum = super.getParams.length
  override def setParams(i:Int)(a:Array[Double])={
    if (i != paramNum){super.setParams(i)(a)}
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
  private val paramNum = super.getParams.length
  val originalBranchLengths = tree.getBranchLengths
  override def setParams(i:Int)(a:Array[Double])={
    if (i != paramNum){super.setParams(i)(a)}
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
  //following params need to be lazy so they are not evaluated until the rest of the model is set up
  // a mapping from new indicies to old
  lazy val paramMap=super.getParams.zipWithIndex.filter{t=>val (p,i)=t;!(toCombine.exists{j=>j==i})}.map(_._2).zipWithIndex.foldLeft(IntMap[Int]():Map[Int,Int]){(m,t)=>m+((t._2,t._1))}
  private lazy val paramNum = super.getParams.length - toCombine.length
  lazy val paramLengths = super.getParams.zipWithIndex.filter{t=>toCombine.exists{j=>j==t._2}}.map{_._1.length}
  override def getParams:List[Array[Double]]={
    val start = super.getParams
    //not super efficient but this should be rarely called
    start.zipWithIndex.filter{t=>val (array,i)=t; !toCombine.exists{j=>j==i}}.map{_._1} ++ List(start.zipWithIndex.filter{t=>toCombine.exists{j=>j==t._2}}.map{_._1.toList}.flatten[Double].toArray)
  }
  override def setParams(i:Int)(a:Array[Double]){
    if (paramMap contains i){
      super.setParams(paramMap(i))(a)
    }else {
      if (i==paramNum){
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
  import pal.substmodel.GammaRates

  import Gamma._

    val gamma = new GammaRates(numCat,1.0)
    def apply(shape:Double):Seq[Double]={
    if (shape < 150.0D){
      cache.getOrElseUpdate((numCat,shape),{
          gamma setParameter (shape,0)
          gamma getRates
        }
        )
    }else {
      (1 to numCat).map(i=> 1.0D).toList // assume alpha => inf
    }

  }
}

trait AlternateModel[A <: BioEnum] extends Model[A]{
  val nodeIDs:Set[Int]
  val altModel:Model[A]
  def altModelParam:Int = altModel.getParams.length-1
  private val paramNum=super.getParams.length
  override def qMat(n:CalcLikelihoodNode[A])={
    if (nodeIDs contains n.id){
      debug{"Using altmodel as " + nodeIDs + " contains " + n.id}
      altModel.qMat(n)
    }else {
      debug{"Using regular model as " + nodeIDs + " does not contain " + n.id}
      super.qMat(n)
    }
  }
  override def getParams:List[Array[Double]]={
    //default to last parameter of altModel
    super.getParams ++ List(altModel.getParams.last)
  }
  override def setParams(i:Int)(a:Array[Double])={
    if (i==paramNum){
      altModel.setParams(altModelParam)(a)
    }else {
      super.setParams(altModelParam)(a)
    }
  }
  override def cromulent = super.cromulent && altModel.cromulent
  override def toString = {
    super.toString + "\n" + 
    "alt: {\n " + altModel.toString + "\n"+
    "} "
  }
}

class GammaModel[A <: BioEnum](val pi:Vector,val sMat:Matrix,var tree:Tree[A]) extends SiteClassPiModel[A] with GammaSMat[A]

trait GammaSMat [A <: BioEnum] extends SMat[A]{

  var alpha=1.0
  val alphabet = tree.alphabet
  val numClasses = alphabet.numClasses
  private val paramNum = super.getParams.length
  val gamma = new Gamma(numClasses)

  override def qMat(node:CalcLikelihoodNode[A]) = gammaQMat
  def gammaQMat={
    debug{"Gamma matrix with gamma value " + alpha}
    val gammaVals = gamma(alpha).toList
    val piCat = pi.toArray
    (0 to pi.size-1).foreach{i=> piCat(i)/=numClasses}
    val qStart = Matrix(pi.size * numClasses, pi.size * numClasses)
    (0 to numClasses-1).toList.zip(gammaVals).foreach{t=>
      val (i,rate)=t
      val part = qStart.viewPart(i * numClasses,i*numClasses,pi.size,pi.size)
      part.assign(sMat.sToQ(piCat).normalize(rate))
    }
    if (qStart exists {d=>d.isNaN}){
      println("ERROR ")
      println(getParams.map{_.toList}.mkString(","))
      println("produces sMat " + sMat)
      println("and pi " + piVals)
      println("gamma Rates " + gammaVals)
      println("qMat " + qStart)
    }
    qStart
  }
  
  override def getParams:List[Array[Double]]={
    super.getParams ++ List(Array(alpha))
  }

  override def cromulent  = {
    super.cromulent && alpha > 0.0D && alpha < 200.0D && gamma(alpha).toList.find{_.isNaN}.isEmpty
    //no point having alpha go higher - count as infinity
  }

  override def setParams(i:Int)(a:Array[Double])={
    debug{"Setting main alpha => " + a(0)}
    if (i != paramNum){super.setParams(i)(a)}
    else {alpha = a(0)}
  }

  override def toString = {
    super.toString + "\n"+
    "shape: " + alpha  
  }
}

trait SiteClassPiModel[A <: BioEnum] extends FullPiModel[A]{
  override def piVals = {
    val numClasses = tree.alphabet.numClasses
    val piStart = Vector(pi.size * numClasses)
    val piCat = pi.toArray
    (0 to piCat.length-1).foreach{i=> piCat(i)/=numClasses}
    (0 to piCat.length-1).foreach{i=> piStart.viewPart(i * pi.size,pi.size).assign(piCat)}
    piStart
  }
}

trait FullPiModel[A <: BioEnum] extends Model[A]{
  val pi:Vector
  final val piParam = super.getParams.length

  override def toString = {
    super.toString + "\n"+
    "Pi: " + pi.toString
  }
  var medianIndex:Int=0
  override def cromulent = super.cromulent &&  (pi.zSum < 1.00000001) && (pi.zSum > 0.99999999)

  override def setParams(i:Int)(a:Array[Double])={
    if (i==piParam){setPi(a)}
    else{throw new UnknownParamException(i)}
  }

  override def getParams:List[Array[Double]]={
    super.getParams ++ List(toFit.toArray)
  }

  def setPi(array:Array[Double]){
    val exponentiated =  array.map{i=>Math.exp(i)}
    val total = (0.0D /: exponentiated){_+_} + Math.exp(0.0D)
    val values = ((0 to medianIndex-1 ).map{i=> exponentiated(i)/total}.toList ++ List(Math.exp(0.0D)/total) ++ (medianIndex to array.length-1).map{i=> exponentiated(i)/total}).toArray
    pi assign values
  }


  def toFit:List[Double]={                                                                   
    val t  = (pi.toList.zipWithIndex.toList.sort{_._1<_._1})(pi.size/2)
    medianIndex = t._2
    def toFitness={i:Int=>Math.log(pi(i)/pi(medianIndex))}
    (0 to medianIndex-1).map{toFitness}.toList ++ (medianIndex+1 to pi.size-1).map{toFitness}.toList
  }
}

class EnhancedModel[A <: BioEnum](val pi:Vector,val sMat:Matrix, var tree:Tree[A]) extends Model[A] with FullPiModel[A] with SMat[A]

trait SMat[A <: BioEnum] extends Model[A]{
  val sMat:Matrix
  final val sMatParam=super.getParams.length

  override def toString = {
    super.toString + "\n" + 
    "S: " + sMat.toString 
  }

  import cern.colt.function.DoubleFunction

  override def cromulent = super.cromulent && !(sMat exists {i=> i < 0.0D}) 

  def linearSMat:Seq[Double]={
    val list = (0 to sMat.rows -3).map{ i=> // skip last elment...
       val list1 = sMat(i).viewPart(i+1,sMat.columns - i - 1).toList
       list1
    }.toList.flatten[Double]
    list
  }

  override def getParams:List[Array[Double]]={
    super.getParams ++ List(linearSMat.toArray)
  }
  
  override def setParams(i:Int)(a:Array[Double]){
    if (i==sMatParam){setSMat(a)}
    else{super.setParams(i)(a)}
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
}


