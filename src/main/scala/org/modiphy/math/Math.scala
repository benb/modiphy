package org.modiphy.math
import cern.colt.matrix._
import cern.colt.matrix.DoubleFactory2D._
import cern.colt.function.DoubleFunction
import cern.colt.matrix.linalg.{Algebra,EigenvalueDecomposition}
import org.modiphy.util._


case class Foo(d:Double)

object EnhancedMatrix{
  type Matrix=DoubleMatrix2D
  type Matrix1D=DoubleMatrix1D
  implicit def enhanceMatrix(m:Matrix)=new EnhancedMatrix(m)
  implicit def enhanceMatrix1D(m:Matrix1D)=new EnhancedMatrix1D(m)
  implicit def enhanceArray(m:Array[Double])=new EnhancedArray(m)
  implicit def enhanceList(m:List[Double])=new EnhancedList(m)
  private lazy val algebra = new Algebra

  class EnhancedList(l:List[Double]){
    def toMatrix1D=Matrix1D(l)
  }
  class EnhancedArray(a:Array[Double]){
    def toMatrix1D=Matrix1D(a)
  }

}

import EnhancedMatrix._

//Doesn't play well with other actors 
//(react from channel belonging to other actor error)
//needs to be recoded as pure actor
//TODO
trait CachedMatrixExponential extends MatrixExponential{
  import scala.actors.Actor
  import scala.actors.OutputChannel
  import scala.actors.Actor._
  case class CacheReq(t:Double)
  def realExp(t:Double) = super.exp(t)
  case class Calc(o:OutputChannel[Any],r:CacheReq)
  case class Answer(t:Double,m:Matrix)
  case object ExitNow
  class CacheActor extends Actor{
    val cache = new SoftCacheMap[Double,Matrix](500)
    def act{
        loop{
          react{
            case ExitNow => Actor.exit
            case Answer(t,m)=>
            println("GOT ANSWER")
              cache+((t,m)) 
            case CacheReq(t) =>
            println("GOT CACHEREQ")
              val cacheLookup = cache.get(t)
              if (cacheLookup.isDefined){
                reply(Answer(t,cacheLookup.get))
              }
              else {
                Actor.actor{
                  react{
                    case Calc(o,CacheReq(t))=>
                      val ans = realExp(t)
                      o ! Answer(t,ans)
                      reply(Answer(t,ans))
                  }
                  exit
                }
              } ! Calc(sender,CacheReq(t))
          }
        }
    }
    
  }
  val actor = new CacheActor
  actor.start
  override def exp(t:Double)={
    println("SENDING")
    (actor !? CacheReq(t)).asInstanceOf[Answer].m
  }
  def exit { actor ! ExitNow }
}



object Matrix1D{
  def apply(i:Int):Matrix1D=cern.colt.matrix.DoubleFactory1D.dense.make(i)
  def apply(a:Double*):Matrix1D=Matrix1D(a.toArray)
  def apply(a:Array[Double]):Matrix1D=cern.colt.matrix.DoubleFactory1D.dense.make(a)
  def apply(l:List[Double]):Matrix1D=Matrix1D(l.toArray)
}
object Matrix{
  def apply(i:Int,j:Int)=dense.make(i,j)
}

trait MatrixExponential{
  def pi:Matrix1D
  def u:Matrix
  def v:Matrix
  def q:Matrix
  lazy val qNorm = u * d * v
  def rawD:Matrix
  def scale:Option[Double]
  lazy val d = {
    if (scale.isDefined){
      def normFact = - dense.diagonal(q).zDotProduct(pi) / scale.get
      def normFunc = new DoubleFunction(){def apply(v:Double)=v/normFact}
      rawD.copy.assign(normFunc)
    }else {
      rawD
    }
  }
  def exp(t:Double)={
    (u * d.expVals(t)) * v
  }
}


class MatExpNormal(val q:Matrix,val pi:Matrix1D,val scale:Option[Double]) extends MatrixExponential{
  val eigen = new EigenvalueDecomposition(q)  
  val algebra = new Algebra
  val u = eigen.getV
  val rawD = eigen.getD
  val v = algebra.inverse(u)
}
class MatExpScale(m:MatrixExponential,val s:Double) extends MatrixExponential{
  def pi = m.pi
  def u = m.u
  def v = m.v
  def q = m.q
  def rawD = m.rawD
  val scale = Some(s)
}
class BasicMatExpScale(val u:Matrix,val rawD:Matrix,val v:Matrix,val pi:Matrix1D, s:Double) extends MatrixExponential{
  def this(u:Matrix,rawD:Matrix1D,v:Matrix,pi:Matrix1D,scale:Double)= this(u,sparse.diagonal(rawD),v,pi,scale)
  val scale = Some(s)
  def q = u * rawD * v
}


class MatExpYang(val q:Matrix,val pi:Matrix1D,val scale:Option[Double]) extends MatrixExponential{
  import cern.colt.matrix.DoubleFactory2D.sparse

    val funcRoot = new DoubleFunction(){def apply(v:Double)=Math.sqrt(v)}
  val funcRec = new DoubleFunction(){def apply(v:Double)=1.0D/v}
  val piRootVec =  pi.copy.assign(funcRoot)
    val piRoot = sparse.diagonal(piRootVec)
    val inversePiRoot:Matrix = sparse.diagonal(piRootVec.copy.assign(funcRec))

    val a = ((piRoot * q) * inversePiRoot)//.symmetrise


  
  val eigen = new EigenvalueDecomposition(a)
    val u:Matrix = inversePiRoot * eigen.getV
  val v:Matrix = eigen.getV.viewDice * piRoot
  val rawD = eigen.getD
//  println("MatExpYang" + d(0,0) + " " + q(0,0))
}



object MatExp{
  private lazy val algebra = new Algebra
  private lazy val cache:SoftCacheMap[String,(EigenvalueDecomposition,Matrix)] = new SoftCacheMap(10)
  def decomp(m:Matrix)={
    if (m.exists{d=> d==Math.NEG_INF_DOUBLE || d==Math.POS_INF_DOUBLE || d.isNaN}) throw new InvalidMatrixException("Invalid matrix")
    cache.getOrElseUpdate(m.toString,{
        val e = new EigenvalueDecomposition(m); (e,algebra.inverse(e.getV))
     })
  }
  def exp(m:Matrix,t:Double) = {
    val (eigen,vprime)=decomp(m)
    val v = eigen.getV
    v * (eigen.getD expVals t) * vprime
  }
}
class EnhancedMatrix(d:DoubleMatrix2D){
  type ID={def id:Int}
  def apply(i:Int)=d.viewRow(i)
  def exp(t:Double)=MatExp.exp(d,t)
  def expVals(t:Double)=sparse.diagonal(dense.diagonal(d).assign( new DoubleFunction(){def apply(arg:Double)={Math.exp(t * arg)}}))
  def *(m:Matrix)=algebra.mult(d,m)
  def *(n:Double)=d.assign(new DoubleFunction(){def apply(v:Double)=v * n}) 
  def update(i:Int,j:Int,v:Double):Unit=d.set(i,j,v)
  def update(i:ID,j:ID,v:Double):Unit=d.set(i.id,j.id,v)
  def apply(i:Int,j:Int):Double=d.get(i,j)
  def apply(i:ID,j:ID):Double=apply(i.id,j.id)
  def symmetrise(s:Matrix)={
    val s2 = s.copy
    (0 to d.columns-1).foreach{i=>
      (i+1 to d.columns-1).foreach{j=>
        s2(j,i)=s2(i,j)
      }
    }
    s2
  }
  def sToQ(pi:{def apply(i:Int):Double})={
    val q = symmetrise(d)
    (0 to d.columns-1).foreach{i=>
      q.viewColumn(i).assign(new DoubleFunction(){
          def apply(j:Double)=j*pi(i)
      })
    }
    (0 to q.rows-1).foreach{i=> 
      q(i,i)=0
      q(i,i)= -(q.viewRow(i).zSum)
    }
    q
  }
  def normalize(v:Matrix1D):Matrix=normalize(v,1.0D)
  def normalize(v:Matrix1D,overall:Double):Matrix={
    val sum = rate(v)
    d.copy.assign(new DoubleFunction(){def apply(d:Double)= overall * d/sum})
  }
  def rate(v:Matrix1D):Double={
    -dense.diagonal(d).zDotProduct(v)
  }

  def exists(f: Double=>Boolean):Boolean={
   this.elements.exists{f}
  }
  def toList=elements.toList
  

  def diagonal=dense.diagonal(d)

  def toUpper={
    val s = d.like
    (0 to d.rows-1).foreach{i=>
      (i to d.columns-1).foreach{j=>
        s(i,j)=d(j,i) 
      }
    }
  s
  }

  def fixDiag={
    rowElements.toList.zipWithIndex.foreach{t=>
      val (v,i)=t
      v(i)=0
      v(i)= -v.zSum
    }
    d 
  }
  def rowElements={
    for (i <- 0 until d.rows) yield d.viewRow(i)
  }
  def elements={
    val outer = this
    new Iterator[Double]{
      var i=0
      var j=0
      def hasNext={i < d.rows && j < d.columns}
      private def inc = {
        j+=1
        if (j >= d.columns){
          j=0
          i+=1
        }
      }
      def next={val ans = d.getQuick(i,j);inc;ans}
    }
  }

}
class EnhancedMatrix1D(d:DoubleMatrix1D){
  type ID={def id:Int}
  def toList = (0 to d.size -1).map{i=>d.get(i)}.toList
  def elements=d.toArray.elements
  def apply(i:Int):Double=d.get(i)
  def apply(i:ID):Double=apply(i.id)
  def update(i:Int,v:Double):Unit=d.set(i,v)
  def update(i:ID,v:Double):Unit=update(i.id,v)
  def *(n:Double)=d.copy.assign(new DoubleFunction(){def apply(v:Double)=v * n}) 
  def /(n:Double)={*(1/n)} 
  def normalize(a:Double):Matrix1D={
    val sum = d.zSum
    d.copy.assign(new DoubleFunction(){def apply(o:Double)=a*o/sum})
  }

  def normalize:Matrix1D=normalize(1.0D)
}
