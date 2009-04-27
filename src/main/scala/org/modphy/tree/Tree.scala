package org.modphy.tree
import org.modphy.sequence._
import org.modphy.math.EnhancedMatrix._
import cern.colt.matrix.DoubleFactory1D


import scala.util.parsing.combinator._
class TreeParser[A <: BioEnum](m:Map[String,String],alphabet:A) extends JavaTokenParsers{
  def root: Parser[INode[A,Node[A]]] = "("~node~rep(","~>node)~");" ^^ {
    case "("~node1~nodeList~");" => new INode[A,Node[A]](node1::nodeList,alphabet,0.0)
  }
  def node: Parser[Node[A]] = leaf | "("~node~rep(","~>node)~"):"~floatingPointNumber ^^ {
    case "("~node1~nodeList~"):"~length => new INode[A,Node[A]](node1::nodeList,alphabet,length.toDouble)
  } 
  def seqName: Parser[String] = regex(new scala.util.matching.Regex("[a-zA-Z0-9_.+-]+"))
  def leaf: Parser[Leaf[A]] = seqName~":"~floatingPointNumber ^^ {case name~":"~length => new Leaf(name,m(name),alphabet,length.toDouble)}
}


object DataParse{
  def apply[A <: BioEnum](tree:String,alignment:Iterator[String],alphabet:A)={
    val aln = new Fasta(alignment)
    val seqMap = aln.foldLeft(Map[String,String]()){_+_}
    val root = new TreeParser[A](seqMap,alphabet){def parseAll=parse(root,tree)}.parseAll.get
    (root,aln)
  }
}



object Likelihood{
  type Model=(Matrix,Vector)//Q matrix and pi Values
}
import Likelihood._


abstract class Node[A <: BioEnum]{
  def mkLkl(mod:Model):Node[A] with LikelihoodNode[A]
  def lengthTo:Double
  def child(i:Int):Option[Node[A]]
  def name:String
}

trait LikelihoodNode[A <: BioEnum]{
  val alphabet:A
  def likelihoods:List[Vector]
}

trait CalcLikelihoodNode[A <: BioEnum] extends LikelihoodNode[A]{
  //return a list of the lists of probabilties for each site
  def childElements:Iterator[Node[A] with LikelihoodNode[A]]
  def length(i:Int):Double
  val model:Model
  val qMatrix=model._1
  val pi=model._2
  
  override lazy val likelihoods:List[Vector]={



    val childLkl = childElements.map{i=>(i.likelihoods,i.lengthTo)}
    val intermediates= childLkl.map{t=>
    val (siteVectorList,length)=t // list of vectors - 1 for each site
    println("Q = " + qMatrix)
    val matrix = qMatrix exp length //e^Qt
    println("e^Qt=" + matrix)
    siteVectorList.map{siteVector=>
      val ret = DoubleFactory1D.dense.make(alphabet.length) 
        (0 to alphabet.length-1).foreach{i=>
          (0 to alphabet.length-1).foreach{j=>
            ret(j)=ret(j) + siteVector(i) * matrix(i,j)
          }
        }
        ret
      }
    }
    intermediates.map{vecList=>
      vecList.tail.foldLeft(vecList.head){(a,b)=>(0 to a.size-1).foreach{i=>a(i)*=b(i)};a}
    }.toList
  }
}

trait LeafLikelihoodNode[A <: BioEnum] extends LikelihoodNode[A]{
  val sequence:BioSeq[alphabet.Value]
  override lazy val likelihoods:List[Vector]={
    val elements = alphabet.length
    sequence.map{a=>
      val vec = DoubleFactory1D.dense.make(alphabet.length)
        println("Vector " + vec)
        alphabet.getNums(a).foreach{i=>vec(i)=1.0D}
        println(vec)
        vec
    }.toList
  }
}


class INode[A <: BioEnum,B <:Node[A]](children:List[B],val alphabet:A,val lengthTo:Double) extends Node[A]{
  val name=""
  def childElements:Iterator[B] = children.elements
  def child(i:Int)=if (i < children.length){Some(children(i))}else{None}
  def length(i:Int):Double=children(i).lengthTo
  def length(n:Node[A]):Double=children.find{node=>node==n}.get.lengthTo
  def mkLkl(mod:Model):INode[A,Node[A] with LikelihoodNode[A]] with LikelihoodNode[A]=new INode[A,Node[A] with LikelihoodNode[A]](children.map{t=>t.mkLkl(mod)}.toList,alphabet,lengthTo) with CalcLikelihoodNode[A]{val model=mod}
}

class Leaf[A <: BioEnum](val name:String,val seq:String,val alphabet:A,val lengthTo:Double) extends Node[A]{
  def mkLkl(mod:Model):Leaf[A] with LikelihoodNode[A]=new Leaf[A](name,seq,alphabet,lengthTo) with LeafLikelihoodNode[A]{val sequence=alphabet.fromString(seq)}
  def child(i:Int)=None
  override def toString=name+":"+lengthTo
}


