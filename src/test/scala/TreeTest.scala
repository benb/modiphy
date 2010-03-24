import org.scalatest.FunSuite
import org.scalatest.matchers.ShouldMatchers._
import org.modiphy.math._
import org.modiphy.tree._
import ModelData._
import org.modiphy.math.EnhancedMatrix._


class TreeSuite extends FunSuite {
   val (tree,aln) = DataParse(treeStr,alnStr.lines,new org.modiphy.sequence.SiteClassAA(1))

   test("BranchLengths"){
     val tree2 = tree.copy.setBranchLengths(tree.getBranchLengths)
     tree2.getBranchLengths should be (tree.getBranchLengths)
   }
 }
