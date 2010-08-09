package org.modiphy.sequence
import org.modiphy.math._
import org.modiphy.math.EnhancedMatrix._
import org.modiphy.math.Vector

abstract class BioEnum(names:String*) extends Enumeration(names: _*){
  def isReal(a:Value)=true
  def getNums(a:Value):List[Int]
  val matLength:Int
  def parseString(s:String):BioSeq[Value]
  def matElements:List[Value]
  val numClasses=1
  val numAlpha=matElements.length
  def unknown:Value
}

case class SiteClassDNA(override val numClasses:Int) extends BioEnum("A","G","C","T","N","-"){
  type Base = Value
  val A,G,C,T,N,GAP=Value
  override val matLength=numClasses*4
  override def isReal(a:Base)=((a!=N) && (a!=GAP))
  def getNums(a:Base)=if (isReal(a)){(a.id to (matLength-1) by numAlpha).toList}else{(0 to matLength-1).toList}
  def parseString(s:String)=new BioSeq(s.toList.map{i=> valueOf(i.toString).getOrElse(N)})
  def matElements=List(A,G,C,T)
  val unknown=N
}

object DNA extends SiteClassDNA(1)
case class SiteClassAA(override val numClasses:Int) extends BioEnum("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","X","-"){
  type AminoAcid = Value
  val A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,X,GAP=Value
  override val matLength=numClasses*20
  override def isReal(a:AminoAcid)=((a!=X) && (a!=GAP))
  def getNums(a:AminoAcid)=if (isReal(a)){(a.id to (matLength-1) by numAlpha).toList}else{(0 to matLength-1).toList}
  def parseString(s:String)=new BioSeq(s.toList.map{i=> valueOf(i.toString).getOrElse(N)})
  def matElements=List(A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V)
  val unknown=X
}
object AA extends SiteClassAA(1)

class BioSeq[A](seq:Seq[A]) extends Seq[A]{
  def length = seq.length
  def iterator = seq.iterator
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

class PamlAlignment(source:Iterator[String]) extends Iterator[(String,String)]{
  source.next // skip first line
  def next = {
    val l = source.next.split("\\s+")
    (l(0),l(1))
  }
  def hasNext = source.hasNext
}
class UnsupportedAlignmentFormatException(s:String) extends Exception(s)

object GenAlnParser{
  def apply(source:Iterator[String]):Iterator[(String,String)]={
    val iter = source.buffered
    if (iter.head.trim matches "[0-9]+\\s+[0-9]+"){//paml first line
      new PamlAlignment(iter)
    }else if(iter.head matches ">.*") {//fasta
      new Fasta(iter)
    }else {
      throw new UnsupportedAlignmentFormatException("Can't read alignment with header:\n" + iter.head)
    }

  }
}


class Fasta(source:Iterator[String]) extends Iterator[(String,String)]{
  import EnhancedIterator._
  val iter=source.map{i=>i.trim}.buffered
  def next = {
    val name = iter.next.trim.drop(1).trim.split("\\s+")(0)
      val seq = iter.takeWhileSafe{s:String => !(s startsWith (">"))}.mkString("").toUpperCase
      (name,seq)
    }
    def hasNext = iter.hasNext
    def toAlignment[A <: BioEnum](alphabet:A)={new Alignment(this.foldLeft(Map[String,String]()){_+_},alphabet)}
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


abstract class AbsAlignment[A<:BioEnum]{
  import scala.collection.immutable.{Map,TreeMap}
  val alphabet:A
  type Letter = alphabet.Value
  def get(a:String):List[Letter] 
  def patterns:(TreeMap[String,List[Letter]],List[Int])
  def split(i:Int):List[Alignment[A]]
}


class Alignment[A<:BioEnum](m:Map[String,String],val alphabet:A){
  type Letter=alphabet.Value
  import scala.collection.immutable.{Map,TreeMap}
  val map = m.foldLeft(TreeMap[String,String]()){_+_} // use TreeMap to keep sorted


  val stateMap=map.map{t=>(t._1,t._2.split("").drop(1).map{alphabet.valueOf(_).getOrElse(alphabet.unknown)}.toList)}.foldLeft(Map[String,List[Letter]]()){_+_}
  private def columns:List[List[String]] = new FlippedIterator(map.values.map{s=>s.split("").drop(1).elements}.toList).toList
  val alphaList:List[List[Letter]] = columns.map(_.map{alphabet.valueOf(_).getOrElse(alphabet.unknown)}) //mapped to alphabet values
  
  def get(a:String):List[Letter]=stateMap.get(a).getOrElse(Nil)
  def colPatterns=alphaList.foldLeft(Map[List[Letter],Int]()){(m,l)=>
    m(l)=m.getOrElse(l,0)+1  
  }.toList

  val patterns={
    val l=colPatterns
    val i=new FlippedIterator(l.map{_._1.elements}).toList
    (
      map.keySet.toList.zip(i).foldLeft(TreeMap[String,List[Letter]]()){(m,t)=>m+((t._1,t._2.toList))}, 
      l.map{_._2}
    )
  }
  val pCount={
    //println("PATTERNS") 
    //println(patterns._1)
    //println(patterns._2)
    patterns._2
  }

  def length = patterns._2.foldLeft(0){_+_}
  /**
   return (List[patternCount],Map(Name,Patterns)) taking npatterns after dropping skip
  */
  def sub(skip:Int,nPatterns:Int)={
    (patterns._2.drop(skip).take(nPatterns),
      patterns._1.map{t=> val(name,patterns)=t
        (name,patterns.drop(skip).take(nPatterns))
      }.foldLeft(Map[String,List[Letter]]()){_+_}
    )
  }

  def split(i:Int):List[Alignment[A]]={
    val size = patterns._1.values.iterator.next.length
    for (j <- 0 until (size - (size % i)) by (size/i)) yield {
      val t = if (j < (size - (size % i)) - (size/i)){
        sub(j,(size/i))
      }else {
        sub(j,(size/i)+(size%i))
      }
        val subMap = t._2.map{tx => 
          (tx._1, tx._2.zip(t._1).map{ty=> ty._1.toString * ty._2}.mkString(""))
        }.foldLeft(Map[String,String]()){_+_}
      new Alignment(subMap,alphabet)
    }
  }.toList

  def getPatterns(a:String)=patterns._1(a)


  def apply(a:String)=map(a)
  def getF={

    val letterCounts:List[Map[Letter,Int]] = map.values.map{
      alphabet.parseString(_).foldLeft(Map[Letter,Int]()){
        (m,l)=>m(l)=m.getOrElse(l,0)+1
      }
    }.toList//map each sequence to a map of counts of each letter

    val pMap = letterCounts.map{m=>
      val propMap = Map[Letter,Double]()
        val total = m.values.foldLeft(0){_+_}
        m.foldLeft(propMap){(m2,t)=> 
          m2(t._1)=t._2.toDouble/total 
        }
    }
    val num = map.size
    val ans = pMap.foldLeft(Map[Letter,Double]()){(m,m2)=>
      m2.foldLeft(m){(m3,t)=>m3(t._1)=m3.getOrElse(t._1,0.0D)+t._2/num}
    }
    ans


    /*
    val countTreeMaps:List[Map[Letter,Int]] = alphaList.map(_.foldLeft(Map[Letter,Int]()){(m,l)=>m.update(l,1+m.getOrElse(l,0))}) //mapped to count of each letter in column

    val filteredCountTreeMaps:List[Map[Letter,Int]] = countTreeMaps.map{_.filter{t=>alphabet.isReal(t._1)}}//remove gaps/unknown
    val freqMaps:List[Map[Letter,Double]] = filteredCountTreeMaps.map{a=>(a.values.foldLeft(0.0D){_+_},a)}.map{a=>a._2.map{t=>(t._1,t._2.toDouble/a._1)}.foldLeft(Map[Letter,Double]()){_+_}} //convert to frequencies at each column
    val (total,totalFreqMap) = freqMaps.foldLeft((0,Map[Letter,Double]())){(a,b)=>
      (a._1+1,b.foldLeft(a._2){(c,d)=>
        c+((d._1,c.getOrElse(d._1,0.0D)+d._2))
      })
    }//total frequencies
    
    val ans = totalFreqMap.map(t=>(t._1,t._2/total)).foldLeft(Map[Letter,Double]()){_+_}
    ans
    */
    }

  def getFPi={
    val f = getF
    Vector(alphabet.matElements.map{f.getOrElse(_,0.0)}).normalize(1.0D)
  }

  def toFasta={
    map.foldLeft(""){(s,t)=>s + ">"+t._1+"\n"+t._2+"\n"}
  }
    
}

class FlippedIterator[A](l:List[Iterator[A]]) extends Iterator[List[A]]{
  def next = l.map{_.next}
  def hasNext = l.foldLeft(true){(a,b)=>a && b.hasNext}
}
