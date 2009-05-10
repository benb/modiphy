package org.modiphy.sequence
import org.scalatest.Suite
import org.modiphy.sequence.DNA._

class SeqSuite extends Suite{
    def testSequence{
      assert(DNA.matLength==4) 
    }
}


