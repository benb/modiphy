package org.modiphy.math
import cern.colt.matrix._
import cern.colt.matrix.DoubleFactory2D._
import cern.colt.function.DoubleFunction
import cern.colt.matrix.linalg.{Algebra,EigenvalueDecomposition}
import org.modiphy.util._

object EnhancedMatrix{
  type Matrix=DoubleMatrix2D
  type Vector=DoubleMatrix1D
  implicit def enhanceMatrix(m:Matrix)=new EnhancedMatrix(m)
  implicit def enhanceVector(m:Vector)=new EnhancedVector(m)
  private lazy val algebra = new Algebra
}
import EnhancedMatrix._
object Vector{
  def apply(i:Int):Vector=cern.colt.matrix.DoubleFactory1D.dense.make(i)
  def apply(a:Double*):Vector=Vector(a.toArray)
  def apply(a:Array[Double]):Vector=cern.colt.matrix.DoubleFactory1D.dense.make(a)
  def apply(l:List[Double]):Vector=Vector(l.toArray)
}
object Matrix{
  def apply(i:Int,j:Int)=dense.make(i,j)
}

object MatExp{
  private lazy val algebra = new Algebra
  private lazy val cache:SoftCacheMap[String,(EigenvalueDecomposition,Matrix)] = new SoftCacheMap(10)
  def decomp(m:Matrix)={
    cache.getOrElseUpdate(m.toString,{val e = new EigenvalueDecomposition(m); (e,algebra.inverse(e.getV))})
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
  def normalize(v:Vector):Matrix=normalize(v,1.0D)
  def normalize(v:Vector,overall:Double):Matrix={
    val sum = -dense.diagonal(d).zDotProduct(v)
    d.copy.assign(new DoubleFunction(){def apply(d:Double)= overall * d/sum})
  }

  def exists(f: Double=>Boolean):Boolean={
   d.toArray.toList.exists{row:Array[Double] => row.toList.exists(f)}
  }

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

}
class EnhancedVector(d:DoubleMatrix1D){
  type ID={def id:Int}
  def toList = (0 to d.size -1).map{i=>d.get(i)}.toList
  def elements=d.toArray.elements
  def apply(i:Int):Double=d.get(i)
  def apply(i:ID):Double=apply(i.id)
  def update(i:Int,v:Double):Unit=d.set(i,v)
  def update(i:ID,v:Double):Unit=update(i.id,v)
  def *(n:Double)=d.assign(new DoubleFunction(){def apply(v:Double)=v * n}) 
  def /(n:Double)={*(1/n)} 
  def normalize(a:Double):Vector={
    val sum = d.zSum
    d.copy.assign(new DoubleFunction(){def apply(o:Double)=a*o/sum})
  }

  def normalize:Vector=normalize(1.0D)
}
