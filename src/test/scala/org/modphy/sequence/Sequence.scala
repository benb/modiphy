package org.modphy.sequence
import org.scalatest.Suite
import org.modphy.sequence.DNA._

class SeqSuite extends Suite{
    def testSequence{
      assert(DNA.matlength==4) 
    }
}


