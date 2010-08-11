import org.scalatest.FunSuite
import org.scalatest.matchers.ShouldMatchers._
import org.modiphy.math._
import org.modiphy.tree._
import ModelData._
import org.modiphy.math.EnhancedMatrix._


class AlignmentSuite extends FunSuite {
   val (tree,aln) = DataParse(treeStr,alnStr.lines,new org.modiphy.sequence.SiteClassAA(1))

   test("Split"){
       val split = aln.split(10)
       println(split)
       println(aln.patterns._2.length)
       println(split.head.patterns._2.length)
       split.length should be (10)
       split.zip(aln.split(10)).foreach{t=> t._1 should equal (t._2)}

       println(split.map{_.length})
       split.foldLeft(0){_+_.length} should be (aln.length)
   }
 }
