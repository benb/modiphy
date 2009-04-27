package org.modphy.math
import cern.colt.matrix._
import cern.colt.matrix.DoubleFactory2D._
import cern.colt.function.DoubleFunction
import cern.colt.matrix.linalg.{Algebra,EigenvalueDecomposition}


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

import EnhancedMatrix._
object MatExp{
  private lazy val algebra = new Algebra
  import scala.collection.jcl.WeakHashMap
  private lazy val cache:WeakHashMap[Matrix,(EigenvalueDecomposition,Matrix)] = new WeakHashMap()
  def decomp(m:Matrix)=cache.getOrElseUpdate(m,{val e = new EigenvalueDecomposition(m); (e,algebra.inverse(e.getV))})
  def exp(m:Matrix,t:Double) = {
    val (eigen,vprime)=decomp(m)
    println(eigen)
    val v = eigen.getV
    println("M " + m)
    v * (eigen.getD expVals t) * vprime
  }
}
class EnhancedMatrix(d:DoubleMatrix2D){
  def exp(t:Double)=MatExp.exp(d,t)
  def expVals(t:Double)=sparse.diagonal(dense.diagonal(d).assign( new DoubleFunction(){def apply(arg:Double)={Math.exp(t * arg)}}))
  def *(m:Matrix)=algebra.mult(d,m)
  def update(i:Int,j:Int,v:Double)=d.set(i,j,v)
  def apply(i:Int,j:Int)=d.get(i,j)
  def symmetrise(s:Matrix)={
    val s2 = s.copy
    (0 to d.columns-1).foreach{i=>
      (i+1 to d.columns-1).foreach{j=>
        s2(j,i)=s2(i,j)
      }
    }
    s2
  }
  def sToQ(pi:Vector)={
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

}
class EnhancedVector(d:DoubleMatrix1D){
  def elements=d.toArray.elements
  def apply(i:Int)=d.get(i)
  def update(i:Int,v:Double)=d.set(i,v)
}
