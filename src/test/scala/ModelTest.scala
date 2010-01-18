import org.scalatest.FunSuite
import org.scalatest.matchers.ShouldMatchers._
import org.modiphy.math._
import org.modiphy.tree._
import ModelData._


class ModelSuite extends FunSuite {
   
  val (tree,aln) = DataParse(treeStr,alnStr.lines,new org.modiphy.sequence.SiteClassAA(1))
  val model = ModelFact.basic(WAG.pi,WAG.S,tree)


  val (tree4,aln4) = DataParse(treeStr,alnStr.lines,new org.modiphy.sequence.SiteClassAA(4))
  val sC=new BasicSMatComponent(WAG.S)
  val gammaC=new GammaMathComponent(0.5,tree4.alphabet)
  val piC = new BasicPiComponent(WAG.pi)
  val piCG = new FlatPriorPiComponent(piC,tree4.alphabet)
  val model2 = new ComposeModel(piCG,sC,gammaC,tree4)

  test("log likelihood of basic model should match PAML") (model.logLikelihood should be (-6057.892394 plusOrMinus 0.001))//from PAML
  test("log likelihood of gamma model should match PAML") (model2.logLikelihood should be (-5808.929978 plusOrMinus 0.001))//from PAML
  test("updating gamma model to gamma+F model should match PAML"){
    piC.getParams(0).setPi(Array(0.038195,0.070238,0.054858,0.072802,0.037939,0.046398,0.080749,0.048962,0.017175,0.043066,0.085106,0.069726,0.015124,0.046142,0.028198,0.073571,0.044604,0.024096,0.049474,0.053576))
    model2.logLikelihood should be (-5810.399586 plusOrMinus 0.001)
  }

}
