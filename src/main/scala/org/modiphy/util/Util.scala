package org.modiphy.util
import scala.actors.Futures._
import org.apache.commons.collections.LRUMap


class PList[A](l:List[A]){                          
  def pMap[B](f:A=>B)={l.map{i=>future(f(i))}.map{_()}}
}
class PStream[A](s:Stream[A]){
  val chunk =20
  def pMap[B](f:A=>B)={
    s.map{i=>future(f(i))}.toList.map{_()}.toStream 
  }

}
object PList{
  implicit def ListToPList[A](l:List[A])=new PList(l)
  implicit def StreamToPStream[A](s:Stream[A])=new PStream(s)
}

class CacheMap[A,B](defaultSize:Int) extends scala.collection.mutable.Map[A,B]{
  val map = new org.apache.commons.collections.LRUMap(defaultSize)
  def get(a:A)={
    val b = map.get(a)
    if (b==null){
      None
    }else{
      Some(b.asInstanceOf[B])
    }
  }
  def -=(a:A)={
    map remove a
    this
  }
  def +=(kv:(A,B))={
    map put (kv._1,kv._2)
    this
  }
  def iterator={new JclIterator(map.iterator).map{_.asInstanceOf[A]}.map{a=>(a,get(a).get)}}
}
import scala.ref.SoftReference
class SoftCacheMap[A,B <: AnyRef](defaultSize:Int) extends scala.collection.mutable.Map[A,B]{
  val map = new CacheMap[A,SoftReference[B]](defaultSize)
  def +=(kv:(A,B))={
    map put (kv._1,new SoftReference[B](kv._2))
    this
  }
  def get(a:A)={
    val b = map.get(a)
    if (b.isDefined){
      b.get.get
    }else{
      None
    }
  }
  def -=(a:A)={
    map -= a
    this
  }
  def iterator=map.iterator.map{t=>(t._1,t._2.get)}.filter{_._2.isDefined}.map{t=>(t._1,t._2.get)}
}

class JclIterator[A](i:java.util.Iterator[A]) extends Iterator[A]{
  def next=i.next 
  def hasNext = i.hasNext
}

object Test{
  def lrt(lkl1:Double,lkl2:Double,df:Int)={
    import org.apache.commons.math.distribution._
    val d = - 2.0 * (lkl1 - lkl2)
    val chi = new ChiSquaredDistributionImpl(df)
    1.0D - chi.cumulativeProbability(d)
  }
  def aic(lkl:Double,numParams:Int)={
    2 * (numParams.toDouble - lkl)
  }
}
