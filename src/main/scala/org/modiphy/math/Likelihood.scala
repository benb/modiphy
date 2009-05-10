package org.modiphy.math
import EnhancedMatrix._
import org.modiphy.tree._
import org.modiphy.sequence._

class LikelihoodCalc[A <: BioEnum](tree:Node[A],model:Matrix,alphabet:A){
  def count = alphabet.filter{i=>alphabet.isReal(i)}.foldLeft(0){(i,j)=>i+1}

  def likelihood(i:Int)={
        1.0  
  }
}
