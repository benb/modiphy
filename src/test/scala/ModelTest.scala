import org.scalatest.FunSuite
import org.scalatest.matchers.ShouldMatchers._
import org.modiphy.math._
import org.modiphy.tree._
import ModelData._


class ModelSuite extends FunSuite {
   
  val (tree,aln) = DataParse(treeStr,alnStr.lines,new org.modiphy.sequence.SiteClassAA(1))
  val model = ModelFact.basic(WAG.pi,WAG.S,tree)//,Array(0.53))
  val (tree4,aln4) = DataParse(treeStr,alnStr.lines,new org.modiphy.sequence.SiteClassAA(4))
  val model2 = ModelFact.gamma(WAG.pi,WAG.S,0.5,tree4)



  test("log likelihood of gamma model should match PAML") (model2.logLikelihood should be (-5808.929978 plusOrMinus 0.01))//from PAML
  test("log likelihood of basic model should match PAML") (model.logLikelihood should be (-6057.892394 plusOrMinus 0.01))//from PAML

}
