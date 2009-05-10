package org.modiphy.tree
import org.modiphy.sequence._
import org.scalatest.Suite
import org.modiphy.tree._
import org.modiphy.math._
import org.modiphy.math.EnhancedMatrix._
import org.apache.commons.math.optimization.direct._


class ModelSuite extends Suite{
  def testMode{
   val testMat = Matrix(20,20)
     testMat(3,7)= -2
   assert(testMat.exists{i=>i < 0.0})

   println(DNA.getNums(DNA.A))
   println(DNA.getNums(DNA.G))
   println(DNA.getNums(DNA.C))
   println(DNA.getNums(DNA.T))
   assert(DNA.getNums(DNA.A)==List(0))
   assert(DNA.getNums(DNA.G)==List(1))
   assert(DNA.getNums(DNA.C)==List(2))
   assert(DNA.getNums(DNA.T)==List(3))


    val tree = "((a:0.1,b:0.1):0.2,d:0.3);"
    val data = DataParse(tree,DNA)
    val model = new EnhancedModel(Vector(4), Matrix(4,4),data)
    model.setSMat(Array(1.2,0.1,1.1,1.2,0.8))
    model.pi.assign(Array(0.25,0.25,0.25,0.25))
    println("Model " + model)
    assert(model.cromulent)
    model.setSMat(Array(1.2,-0.1,1.1,1.2,0.8))
    println("Model " + model)
    assert(!model.cromulent)
    println("Cromulent " + model.cromulent)


  }
}


