package org.modphy.sequence

abstract class BioEnum(names:String*) extends Enumeration(names: _*){
  def isReal(a:Value)=true
  def getNums(a:Value):List[Int]
  val length:Int=filter{isReal}.toList.length
  def fromString(s:String):BioSeq[Value]
}

object DNA extends BioEnum("A","G","C","T","N","-"){
  val A,G,C,T,N,GAP=Value
  override def isReal(a:Base)=(a!=N && a!=GAP)
  type Base = Value
  def getNums(a:Base)=if (isReal(a)){List(a.id)}else{(A.id to T.id).toList}
  def fromString(s:String)=new BioSeq(s.toList.map{i=> DNA valueOf i.toString getOrElse(DNA.N)})
  override val length=4
}

class BioSeq[A](seq:Seq[A]) extends Seq[A]{
  def length = seq.length
  def elements = seq.elements
  def apply(i:Int)=seq(i)
}

object BioSeq{
  def string2DNA(s:String)=new BioSeq(s.toList.map{i=> DNA valueOf i.toString getOrElse(DNA.N)})
}

object EnhancedIterator{
  implicit def MakeEnhancedIterator[A](i: Iterator[A])=new EnhancedIterator(i.buffered)
}
class EnhancedIterator[A](i:BufferedIterator[A]){
  def takeWhileSafe(f:A=>Boolean):List[A]={
    takeWhileSafe(List(),f)
  }
  def takeWhileSafe(l:List[A],f:A=>Boolean):List[A] ={
    if (i.hasNext && f(i.head)){
      takeWhileSafe(i.next::l,f)
    }else {
      l.reverse
    }
  }
}

class Fasta(source:Iterator[String]) extends Iterator[(String,String)]{
  import EnhancedIterator._
  val iter=source.map{i=>i.trim}.buffered
  def next = {
    val name = iter.next.trim.drop(1)
      val seq = iter.takeWhileSafe{s:String => !(s startsWith (">"))}.mkString("")
      (name,seq)
    }
    def hasNext = iter.hasNext
  }
