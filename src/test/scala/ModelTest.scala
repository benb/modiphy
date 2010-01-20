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
  val piC = new BasicPiComponent(WAG.pi)
  val piCG = new FlatPriorPiComponent(piC,tree4.alphabet)
  val gammaC=new GammaMathComponent(0.5,tree4.alphabet,piCG,sC)
  val model2 = new ComposeModel(piCG,sC,gammaC,tree4)

  test("log likelihood of basic model should match PAML") (model.logLikelihood should be (-6057.892394 plusOrMinus 0.001))//from PAML
  test("log likelihood of gamma model should match PAML") (model2.logLikelihood should be (-5808.929978 plusOrMinus 0.001))//from PAML
  test("updating gamma model to gamma+F model should match PAML"){
    piC.getParams(0).setPi(Array(0.038195,0.070238,0.054858,0.072802,0.037939,0.046398,0.080749,0.048962,0.017175,0.043066,0.085106,0.069726,0.015124,0.046142,0.028198,0.073571,0.044604,0.024096,0.049474,0.053576))
    model2.logLikelihood should be (-5810.399586 plusOrMinus 0.001)
  }

    val mat = Matrix(4,4)
    val thmmMath = new THMMGammaMathComponent(gammaC,mat,tree4.alphabet)
    val model3 = new ComposeModel(piCG,sC,thmmMath,tree4)

  test("thmm should give same answers with 0 c matrix") {
    model3.logLikelihood should be (model2.logLikelihood plusOrMinus 0.001)
  }
  test ("thmm matrix should have 6 parameters for 4 classes"){// as it is not normalised
    thmmMath.getParams.head.getParams.length should be (6)
  }
  test("thmm should give improvment with slow rate of change") {
    //assume a conservative rate of change between fast site classes is probably better than none at all
    //this is not a definite but seems to be true for this dataset
    //and possible that slow<->fast might be close to 0!
    thmmMath.getParams.head.setParams(Array(0,0,0,0,0,0.05))
    val ans1 = model2.logLikelihood
    model3.logLikelihood should be > (ans1)
  }
  test("gamma model should have the right param controls"){
    val params = model2.params
    params.length should be (4)//pi + sMat + gamma(shape) + branches
    val pSorted = params.toList.sort{(a,b)=>a.getParams.length<b.getParams.length}
    pSorted(0).getParams.length should be (1)
    pSorted(1).getParams.length should be (19)
    pSorted(2).getParams.length should be (19) // branch lengths
    pSorted(3).getParams.length should be (189)
  }
  test("thmm model should have the right param controls"){
    val params = model3.params
    params.length should be (5)//pi + sMat + gamma(shape) + cmat + branches
    val pSorted = params.toList.sort{(a,b)=>a.getParams.length<b.getParams.length}
    pSorted(0).getParams.length should be (1)
    pSorted(1).getParams.length should be (6)
    pSorted(2).getParams.length should be (19)
    pSorted(3).getParams.length should be (19) //branch lengths
    pSorted(4).getParams.length should be (189)
  }
  test("Gamma + Invariant model should match phyml"){
    val (tree5,aln5) = DataParse(treeStr,alnStr.lines,new org.modiphy.sequence.SiteClassAA(5))
    val piMe = new BasicPiComponent(WAG.pi)
    val pi = new PriorPiComponent(piMe,tree5.alphabet)
    val gammaMe=new GammaMathComponent(0.5,4,20,pi.getView(20,100),sC)
    val model4 = new ComposeModel(pi,sC,new InvariantMathComponent(20,pi,sC,gammaMe),tree5)
    model4.logLikelihood should be (-5824.746968 plusOrMinus 0.01) // phyml -d aa -o n -i png1-aln.phy -u png1.tre -a 0.5 -m WAG -v 0.2
    var priorInv = 0.4D
    def priorGamma = (1.0D - priorInv)/4.0D
    model4.params(1).asInstanceOf[PiParam].setPi(Array(priorInv,priorGamma,priorGamma,priorGamma,priorGamma))
    model4.logLikelihood should be (-5864.879865 plusOrMinus 0.01) // phyml -d aa -o n -i png1-aln.phy -u png1.tre -a 0.5 -m WAG -v 0.4
    priorInv = 0.3D
    model4.params(1).asInstanceOf[PiParam].setPi(Array(priorInv,priorGamma,priorGamma,priorGamma,priorGamma))
    model4.logLikelihood should be (-5841.318438 plusOrMinus 0.01) // phyml -d aa -o n -i png1-aln.phy -u png1.tre -a 0.5 -m WAG -v 0.4
  }
  test("Branch Length changes"){
    val param = model3.params(4)
    val bls = param.getParams
    val bl = bls(0)
    bls(0)=1.0
    val start=model3.logLikelihood
    param.setParams(bls)
    model3.logLikelihood should be < (start) // artificially setting a branch length to exp(1.0) should be bad for the likelihood
    bls(0)=bl
    param.setParams(bls)
    model3.logLikelihood should be (start)
  }
}
