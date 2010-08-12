package org.modiphy.sequence
import org.modiphy.math._
import org.modiphy.math.EnhancedMatrix._
import org.modiphy.math.Vector

object FlippedIterator{
  def apply[A](l:Iterable[Iterator[A]])(implicit m:Manifest[Iterable[Iterator[A]]])=new FlippedIterator(l)
  implicit def MakeIterator[A](l:Iterable[Iterable[A]])=l.map{_.iterator}
//  def apply[A](l:Iterable[Iterable[A]])(implicit m:Manifest[Iterable[Iterable[A]]])=new FlippedIterator(l.map{_.iterator})
}
class FlippedIterator[A](l:Iterable[Iterator[A]]) extends Iterator[Iterable[A]]{
  def next = l.map{_.next}
  def hasNext = l.foldLeft(true){(a,b)=>a && b.hasNext}
}

import FlippedIterator.MakeIterator
abstract class BioEnum(names:String*) extends Enumeration(names: _*){
  def isReal(a:Value)=true
  def getNums(a:Value):List[Int]
  val matLength:Int
  def parseString(s:String):IndexedSeq[Value]
  def matElements:List[Value]
  val numClasses=1
  val numAlpha=matElements.length
  def unknown:Value
}

case class SiteClassDNA(override val numClasses:Int) extends BioEnum("A","G","C","T","N","-"){
  import scala.collection.immutable.{Vector,VectorBuilder}
  type Base = Value
  val A,G,C,T,N,GAP=Value
  override val matLength=numClasses*4
  override def isReal(a:Base)=((a!=N) && (a!=GAP))
  def getNums(a:Base)=if (isReal(a)){(a.id to (matLength-1) by numAlpha).toList}else{(0 to matLength-1).toList}
  def parseString(s:String)=s.map{i=> valueOf(i.toString).getOrElse(N)}.foldLeft(new VectorBuilder[Base]){_ += _}.result
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
  def parseString(s:String)=s.map{i=> valueOf(i.toString).getOrElse(X)}.toIndexedSeq
  def matElements=List(A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V)
  val unknown=X
}
object AA extends SiteClassAA(1)


object BioSeq{
  def string2DNA(s:String)=s.toList.map{i=> DNA valueOf i.toString getOrElse(DNA.N)}.toIndexedSeq
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
    def toAlignment(alphabet:BioEnum)={SimpleAlignment(this.foldLeft(Map[String,String]()){_+_},alphabet)}
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


abstract class Alignment{
  import scala.collection.immutable.{Map,TreeMap}
  val alphabet:BioEnum
  type Letter = BioEnum#Value 
  def patterns:(Map[String,Seq[BioEnum#Value]],List[Int])
  def split(i:Int):List[Alignment]
  def pCount=patterns._2
  def getPatterns(a:String)=patterns._1(a)
  
  def length = patterns._2.foldLeft(0){_+_}

  def getF={

    val letterCounts:List[Map[Letter,Int]] = patterns._1.values.map{list=> list.zip(patterns._2).foldLeft(Map[Letter,Int]()){(m,t)=>
      val (letter,count)=t
      m.updated(letter,m.getOrElse(letter,0)+count)
    }}.toList//map each sequence to a map of counts of each letter

    val pMap = letterCounts.map{m=>
      val propMap = Map[Letter,Double]()
        val total = m.values.foldLeft(0){_+_}
        m.foldLeft(propMap){(m2,t)=> 
          m2.updated(t._1,t._2.toDouble/total)
        }
    }
    val num = patterns._1.size
    val ans = pMap.foldLeft(Map[Letter,Double]()){(m,m2)=>
      m2.foldLeft(m){(m3,t)=>m3.updated(t._1,m3.getOrElse(t._1,0.0D)+t._2/num)}
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

  def sequenceNames = patterns._1.keys

  def getFPi={
    val f = getF
    Vector(alphabet.matElements.map{f.getOrElse(_,0.0)}).normalize(1.0D)
  }
}

abstract class SimpleAlignment(alphabet:BioEnum) extends SplitAlignment(alphabet){
  val pointers:List[Int]
  def get(seqName:String)={
    val vec = sitePatterns(seqName)
    pointers.map{vec}
  }
  def toFasta={
    sequenceNames.map{s=> (s,get(s).mkString)}.foldLeft(""){(s,t)=>s + ">"+t._1+"\n"+t._2+"\n"}
  }
  override def canEqual(other:Any)=other.isInstanceOf[SimpleAlignment]
  override def equals(other:Any)={
    other match {
      case that:SimpleAlignment=> (that canEqual this) && (that.pointers == this.pointers) && super.equals(that)
      case _ => false
    }
  }
 
}

abstract class SplitAlignment(val alphabet:BioEnum) extends Alignment{
  val sitePatterns:Map[String,scala.IndexedSeq[BioEnum#Value]]
  val counts:List[Int]

  def canEqual(other:Any)=other.isInstanceOf[SplitAlignment]
  override lazy val hashCode = 41 + alphabet.hashCode + sitePatterns.hashCode + counts.hashCode
  override def equals(other:Any)={
    other match {
      case that:SplitAlignment=> (that canEqual this) && (that.sitePatterns==this.sitePatterns) && (that.counts == this.counts)
      case _ => false
    }
  }
  


  def patterns=(sitePatterns,counts)
  def split(i:Int):List[Alignment]={
    val originalLength = sitePatterns.values.head.length
    val size = (originalLength/i) + {if (originalLength%i>0) {1} else {0}}
    val p = sitePatterns.map{t=> (t._1,t._2.grouped(size))}.map{t=>(t._1,t._2.zip(counts.grouped(size)))}.map{ t=>//produce Iterator[(String,Iterable[List[Vector[Letter]],List[List[Int]]])]
      t._2.map{i=>(t._1,i)} // produce Seq[(String,(Vector[Letter],List[Int]))]
    }
    val x = p.map{_.toList}.toList
    val pIter = FlippedIterator(x)
    pIter.map{slice:Iterable[(String,(IndexedSeq[BioEnum#Value],List[Int]))]=>
      new SplitAlignment(alphabet){
        val sitePatterns:Map[String,scala.IndexedSeq[this.alphabet.Value]] = slice.foldLeft(Map[String,scala.IndexedSeq[this.alphabet.Value]]()){(m,t)=>
          m updated (t._1,t._2._1.asInstanceOf[IndexedSeq[this.alphabet.Value]])
        }      
        val counts = slice.head._2._2
      }
    }.toList
  }
}

object SimpleAlignment{
  // using this pattern:
  // http://www.scala-lang.org/node/6353
  trait pack {val alphabet:BioEnum;val patterns:Map[String,scala.IndexedSeq[BioEnum#Value]]}
  trait raw {val alphabet:BioEnum;val alignment:Map[String,scala.IndexedSeq[BioEnum#Value]]}

  def apply(p:pack,countList:List[Int],pointerList:List[Int]):SimpleAlignment={
    import p._
    new SimpleAlignment(alphabet){
      val counts = countList
      val pointers = pointerList
      val sitePatterns = patterns.asInstanceOf[Map[String,scala.IndexedSeq[alphabet.Value]]]
    }
  }

  def apply(p:raw):SimpleAlignment = {
    import p._
    val patterns = FlippedIterator(alignment.values.toList).foldLeft(Map[List[BioEnum#Value],Int]()){(m,s)=> val l = s.toList; m.updated(l,m.getOrElse(l,0)+1)}
    val mySitePatterns = alignment.keys.iterator.zip(FlippedIterator(patterns.keys.toList)).foldLeft(Map[String,scala.Vector[BioEnum#Value]]()){(m,t)=>
      m.updated(t._1,new scala.collection.immutable.VectorBuilder().++=(t._2).result)
    }
    val hash = patterns.keys.zipWithIndex.foldLeft(Map[List[BioEnum#Value],Int]()){_+_}
    val myCounts = patterns.values.toList
    val myPointers = FlippedIterator(alignment.values).map{a=>hash(a.toList)}.toList
    new SimpleAlignment(alphabet){
      val counts = myCounts
      val pointers = myPointers
      val sitePatterns = mySitePatterns.asInstanceOf[Map[String,scala.Vector[alphabet.Value]]]
    }
  }
  def apply(aln:Map[String,String],alpha:BioEnum):SimpleAlignment={
    apply(new raw{
        val alphabet=alpha
        val alignment = aln.map{t=>(t._1,alphabet.parseString(t._2))}.foldLeft[Map[String,IndexedSeq[alphabet.Value]]](Map[String,scala.IndexedSeq[alphabet.Value]]()){_+_}
      }
    )
  }

  /*
  TODO
  def apply(m:Map[String,String],alphabet:A):SimpleAlignment={
  }*/
}

