package org.modiphy.sequence
import org.modiphy.math._

abstract class BioEnum(names:String*) extends Enumeration(names: _*){
  def isReal(a:Value)=true
  def getNums(a:Value):List[Int]
  val matLength:Int
  def parseString(s:String):BioSeq[Value]
  def matElements:List[Value]
  val numClasses=1
  def unknown:Value
}

class SiteClassDNA(override val numClasses:Int) extends BioEnum("A","G","C","T","N","-"){
  type Base = Value
  val A,G,C,T,N,GAP=Value
  override val matLength=numClasses*4
  override def isReal(a:Base)=((a!=N) && (a!=GAP))
  def getNums(a:Base)=if (isReal(a)){(a.id to (matLength-1) by 4).toList}else{(0 to matLength-1).toList}
  def parseString(s:String)=new BioSeq(s.toList.map{i=> valueOf(i.toString).getOrElse(N)})
  def matElements=List(A,G,C,T)
  val unknown=N
}

object DNA extends SiteClassDNA(1)
class SiteClassAA(override val numClasses:Int) extends BioEnum("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","X","-"){
  type AminoAcid = Value
  val A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,X,GAP=Value
  override val matLength=numClasses*20
  override def isReal(a:AminoAcid)=((a!=X) && (a!=GAP))
  def getNums(a:AminoAcid)=if (isReal(a)){(a.id to (matLength-1) by 20).toList}else{(0 to matLength-1).toList}
  def parseString(s:String)=new BioSeq(s.toList.map{i=> valueOf(i.toString).getOrElse(N)})
  def matElements=List(A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V)
  val unknown=X
}
object AA extends SiteClassAA(1)

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

class Maf(source:Iterator[String]) extends Iterator[MafAln]{
  val iter = source.filter{x=> !(x startsWith "#")}.buffered
  def next = {
    new MafAln(iter)
  }
  def hasNext = {while (iter.hasNext && iter.head.matches("\\s+")){iter.next};iter.hasNext}
}

/**
 Represents a single entry of a MAF format alignment
 http://genome.ucsc.edu/FAQ/FAQformat#format5
*/
class MafAln(source:BufferedIterator[String]){
  /**
   Key->Value pairs from the Alignment Block Line
  */
  val aLine = {
    val line = source.next
    assert(line startsWith "a ")
    line.split("\\s+").toList.tail.foldLeft(Map[String,String]()){(m,i)=>val s = i.split("="); m+((s(0),s(1)))}
  }
  /**
   Map of (Name->aligned sequence) 
  */
  var seqs=Map[String,String]()

  /**
   iterator over sequences
  */
  def elements = seqs.elements
  
  while (source.hasNext && !(source.head.startsWith("a"))){
    val line = source.next
    if (line startsWith "s"){
      val s = line.split("\\s+")
      seqs=seqs+((s(1),s(6)))
    }else if (line startsWith "i"){
      true //ignore for now
    }else if (line startsWith "q"){
      true //ignore for now
    }else {
      true //ignore this line
    }
  }

}

class Alignment[A<:BioEnum](val map:Map[String,String],val alphabet:A){

  val stateMap=map.map{t=>(t._1,t._2.split("").drop(1).map{alphabet.valueOf(_).getOrElse(alphabet.unknown)}.toList)}.foldLeft(Map[String,List[alphabet.Value]]()){_+_}
  def get(a:String):List[alphabet.Value]=stateMap.get(a).getOrElse(Nil)

  def apply(a:String)=map(a)
  lazy val columns:List[List[String]] = new FlippedIterator(map.values.map{s=>s.split("").drop(1).elements}.toList).toList
  def getF={
    val alphaList:List[List[alphabet.Value]] = columns.map(_.map{alphabet.valueOf(_).getOrElse(alphabet.unknown)}) //mapped to alphabet values
    val countMaps:List[Map[alphabet.Value,Int]] = alphaList.map(_.foldLeft(Map[alphabet.Value,Int]()){(m,l)=>m.update(l,1+m.getOrElse(l,0))}) //mapped to count of each letter in column

    val filteredCountMaps:List[Map[alphabet.Value,Int]] = countMaps.map{_.filter{t=>alphabet.isReal(t._1)}}//remove gaps/unknown
    val freqMaps:List[Map[alphabet.Value,Double]] = filteredCountMaps.map{a=>(a.values.foldLeft(0.0D){_+_},a)}.map{a=>a._2.map{t=>(t._1,t._2.toDouble/a._1)}.foldLeft(Map[alphabet.Value,Double]()){_+_}} //convert to frequencies at each column
    val (total,totalFreqMap) = freqMaps.foldLeft((0,Map[alphabet.Value,Double]())){(a,b)=>
      (a._1+1,b.foldLeft(a._2){(c,d)=>
        c+((d._1,c.getOrElse(d._1,0.0D)+d._2))
      })
    }//total frequencies
    
    val ans = totalFreqMap.map(t=>(t._1,t._2/total)).foldLeft(Map[alphabet.Value,Double]()){_+_}
    ans
    }

  def getFPi={
    val f = getF
    Vector(alphabet.matElements.map{f(_)})
  }
    
}

class FlippedIterator[A](l:List[Iterator[A]]) extends Iterator[List[A]]{
  def next = l.map{_.next}
  def hasNext = l.foldLeft(true){(a,b)=>a && b.hasNext}
}
