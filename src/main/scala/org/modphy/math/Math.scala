package org.modphy.math
import cern.colt.matrix._
import cern.colt.matrix.DoubleFactory2D._
import cern.colt.function.DoubleFunction
import cern.colt.matrix.linalg.{Algebra,EigenvalueDecomposition}

object Model{
  def GTR={
    val vector = Vector(4)
    vector assign 0.25
    val matrix = Matrix(4,4)
    matrix assign 1.0
    (vector,matrix)
  }
}

object EnhancedMatrix{
  type Matrix=DoubleMatrix2D
  type Vector=DoubleMatrix1D
  implicit def enhanceMatrix(m:Matrix)=new EnhancedMatrix(m)
  implicit def enhanceVector(m:Vector)=new EnhancedVector(m)
  private lazy val algebra = new Algebra
}
object Vector{
  def apply(i:Int)=cern.colt.matrix.DoubleFactory1D.dense.make(i)
}
object Matrix{
  def apply(i:Int,j:Int)=dense.make(i,j)
}

import EnhancedMatrix._
object MatExp{
  private lazy val algebra = new Algebra
  import scala.collection.jcl.WeakHashMap
  private lazy val cache:WeakHashMap[Matrix,(EigenvalueDecomposition,Matrix)] = new WeakHashMap()
  def decomp(m:Matrix)={
    cache.getOrElseUpdate(m,{val e = new EigenvalueDecomposition(m); (e,algebra.inverse(e.getV))})
  }
  def exp(m:Matrix,t:Double) = {
    val (eigen,vprime)=decomp(m)
    //println(eigen)
    val v = eigen.getV
    //println("M " + m)
    v * (eigen.getD expVals t) * vprime
  }
}
class EnhancedMatrix(d:DoubleMatrix2D){
  type ID={def id:Int}
  def apply(i:Int)=d.viewRow(i)
  def exp(t:Double)=MatExp.exp(d,t)
  def expVals(t:Double)=sparse.diagonal(dense.diagonal(d).assign( new DoubleFunction(){def apply(arg:Double)={Math.exp(t * arg)}}))
  def *(m:Matrix)=algebra.mult(d,m)
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
  def normalize:Matrix=normalize(1.0D)
  
  def normalize(i:Double):Matrix={
    val normFact = -(dense.diagonal(d).zSum)
    d.copy.assign(new DoubleFunction(){def apply(d:Double)= i * d/normFact})
  }

  def exists(f: Double=>Boolean):Boolean={
   d.toArray.toList.exists{row:Array[Double] => row.toList.exists(f)}
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
}
