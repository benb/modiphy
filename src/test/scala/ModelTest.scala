import org.scalatest.FunSuite
import org.scalatest.matchers.ShouldMatchers._
import org.modiphy.math._
import org.modiphy.tree._
import ModelData._
import org.modiphy.math.EnhancedMatrix._
import tlf.Logging
import org.modiphy.math.LazyP._


class ModelSuite extends FunSuite {
   
  val (tree,aln) = DataParse(treeStr,alnStr.lines,new org.modiphy.sequence.SiteClassAA(1))
  val model = org.modiphy.math.SimpleModel(tree)

//  Logging.toStdout 
//  Logging setLevel "finest"


  val (tree4,aln4) = DataParse(treeStr,alnStr.lines,new org.modiphy.sequence.SiteClassAA(4))
  val (tree5,aln5) = DataParse(treeStr,alnStr.lines,new org.modiphy.sequence.SiteClassAA(5))

  test("BS complex model with restricted params should match simpler Sigma model"){
    val simpleModel = InvarThmmModel(tree5)
    simpleModel(InvarPrior)=0.123
    simpleModel(BranchLengths)(1)=0.5462
    simpleModel(Alpha)=0.3
    simpleModel(Sigma)=Array.make(10,0.4232)

    val simpleModel2=InvarThmmModel(tree5)
    simpleModel2 << simpleModel
    simpleModel2(BranchLengths) << simpleModel(BranchLengths)
    val lkl1 = simpleModel.logLikelihood

    simpleModel2.logLikelihood should be (lkl1 plusOrMinus 1E-7)

    val complexModel = BranchSpecificThmmModel(tree5,2::3::4::Nil)
    complexModel.logLikelihood should not (be (simpleModel.logLikelihood plusOrMinus 1E-5))
    complexModel << simpleModel
    complexModel(BranchLengths) << simpleModel(BranchLengths)
    complexModel.logLikelihood should not (be (simpleModel.logLikelihood plusOrMinus 1E-5))
    complexModel(Sigma(1)) << simpleModel(Sigma(0))
    complexModel.logLikelihood should be (simpleModel.logLikelihood plusOrMinus 1E-5)

  }

  test("OptUpdate"){
    val simpleModel = InvarThmmModel(tree5)
    val start = simpleModel.logLikelihood
    val opt = simpleModel.optSetter(InvarPrior(0)::BranchLengths(0)::Alpha(0)::Nil)
    val startP = opt.latestArgs
    startP(startP.length-1)=0.2
    opt(startP)
    simpleModel.logLikelihood should not (be (start plusOrMinus 1E-4))
    startP(startP.length-1)=0.5
    opt(startP)
    simpleModel.logLikelihood should be (start plusOrMinus 1E-4)
  }
 
  test("log likelihood of basic model should match PAML") (model.logLikelihood should be (-6057.892394 plusOrMinus 0.001))//from PAML
  test("TufA Gamma Model"){
    val (tree,aln)=DataParse(tufaTree,tufaAln.lines,new org.modiphy.sequence.SiteClassAA(4))
    val pi = Vector(0.066487,0.050798,0.021011,0.069415,0.001862,0.033511,0.095479,0.089628,0.012500,0.061436,0.093617,0.053457,0.032979,0.038830,0.051064,0.027926,0.067287,0.000266,0.025532,0.106915)

    //val gammaModel = ModelFact.gamma(pi,WAG.S,0.496466,tree)
    val gammaModel = GammaModel(tree)//
    gammaModel(Alpha(0))=0.496466D
    gammaModel(Pi(0))=pi
    gammaModel.logLikelihood should be (-2900.328678 plusOrMinus 1E-6)
  }
  val model2 = GammaModel(tree4)
  test("log likelihood of gamma model should match PAML") {
    model2(Alpha(0))=0.5//piCG,sC,gammaC,tree4)
    model2.logLikelihood should be (-5808.929978 plusOrMinus 0.001)
  }//from PAML
    val plusF=aln.getFPi.toArray
  test("updating gamma model to gamma+F model should match PAML"){
    //val plusF=Array(0.038195,0.070238,0.054858,0.072802,0.037939,0.046398,0.080749,0.048962,0.017175,0.043066,0.085106,0.069726,0.015124,0.046142,0.028198,0.073571,0.044604,0.024096,0.049474,0.053576)
    model2(Pi(0))=plusF
    model2.logLikelihood should be (-5810.399586 plusOrMinus 0.001)
  }
  
  test("complex model with restricted params should match simpler model"){
    val model2B = InvarGammaModel(tree5)
    model2B(InvarPrior(0))=0.0D
    model2B(Pi(0))=plusF
    //val model2B = ModelFact.invarThmm(Vector(plusF),WAG.S,0.5,Matrix(5,5),tree5)
   // model2B.getParam("First Prior").head.setParams(Array(0))
    model2B.logLikelihood should be (-5810.399586 plusOrMinus 0.001)
  }

  test("Gamma + Invariant model should match phyml"){
    val (tree5,aln5) = DataParse(treeStr,alnStr.lines,new org.modiphy.sequence.SiteClassAA(5))
    val model4 = InvarGammaModel(tree5)
    model4.params.length should be (5) // pi firstPrior sMat bl shape

    model4.logLikelihood should be (-5824.746968 plusOrMinus 0.01) // phyml -d aa -o n -i png1-aln.phy -u png1.tre -a 0.5 -m WAG -v 0.2

    model4(InvarPrior)=0.4D
    model4.logLikelihood should be (-5864.879865 plusOrMinus 0.01) // phyml -d aa -o n -i png1-aln.phy -u png1.tre -a 0.5 -m WAG -v 0.4
    model4(InvarPrior)=0.3D
    model4.logLikelihood should be (-5841.318438 plusOrMinus 0.01) // phyml -d aa -o n -i png1-aln.phy -u png1.tre -a 0.5 -m WAG -v 0.4
  }

  test("Single Param Wrapper works"){
    val (tree5,aln5) = DataParse(treeStr,alnStr.lines,new org.modiphy.sequence.SiteClassAA(5))
    val model4 = BranchSpecificThmmModel(tree5)
    model4(Sigma(2))(2,3) should be (0.0)
    model4(SingleParam(Sigma(2)))=1.0
    model4(Sigma(2))(2,3) should be (1.0)
    model4(SingleParam(Sigma(2)))=5.0
    model4(Sigma(2))(1,3) should be (5.0)
  }

  test("Copying Params"){

    val (tree5,aln5) = DataParse(treeStr,alnStr.lines,new org.modiphy.sequence.SiteClassAA(5))
    val (tree5a,aln5a) = DataParse(treeStr,alnStr.lines,new org.modiphy.sequence.SiteClassAA(5))
    val model4 = InvarThmmModel(tree5)
    model4(Pi(0))=tree5.aln.getFPi
    val model4b = BranchSpecificThmmModel(tree5a)
    model4b(Pi(0)) << model4(Pi(0))
    model4b.logLikelihood should be (model4.logLikelihood plusOrMinus 1E-5)
    model4b(Pi(0))=WAG.pi
    model4b.logLikelihood should not (be (model4.logLikelihood plusOrMinus 1E-5))

    model4b << model4

    model4b(BranchLengths)(0)=3.0

    model4b.logLikelihood should not (be (model4.logLikelihood plusOrMinus 1E-5))
    model4(BranchLengths) << model4b(BranchLengths)
                                    
    model4b.logLikelihood should be (model4.logLikelihood plusOrMinus 1E-5)
  }                                 



  test("THMM.SI"){
    val (tree5,aln5)=DataParse(pfTree,pfAln.lines,new org.modiphy.sequence.SiteClassAA(5))
    val thmmsi=InvarThmmModel(tree5)

    val sigma = Matrix(5,5)
    sigma(0,4)=2.415327
    sigma(1,4)=2.415327
    sigma(2,4)=2.415327
    sigma(3,4)=2.415327

    val pi = Vector(Array(0.024191,0.002492,0.002932,0.002492,0.001906,0.002492,0.006304,0.023018,0.002346,0.026683,0.034307,0.008943,0.007037,0.014808,0.005278,0.018326,0.013928,0.007477,0.007917,0.020379)).normalize(1)
    println("PI " + pi)
    thmmsi(Pi(0))=pi
    thmmsi(Sigma(0))=sigma
    thmmsi(Alpha(0))=3.270690
    thmmsi(InvarPrior(0))=0.066963
    thmmsi.logLikelihood should be (-2972.196109 plusOrMinus 1e-1)
    //can't hope for much better as we only know params/bls to a certain number of decimal places
    //this log-likelihood is taken from an optimised model of Leaphy's
  }

  test("BS complex model with restricted params should match simpler model"){
    model2(Pi(0))=plusF
    model2.logLikelihood should be (-5810.399586 plusOrMinus 0.001)
  //  val model2B = ModelFact.invarThmmBS(Vector(plusF),WAG.S,0.5,Matrix(5,5),tree5)
    val typ = new org.modiphy.sequence.SiteClassAA(5)
    val (tree5a,aln5a) = DataParse(treeStr,alnStr.lines,typ)
    var i= -1
    val model2C=BranchSpecificThmmModel(tree5a,(tree5a::tree5a.descendentNodes).foldLeft(Map[Int,Int]()){(m,n)=>i=i+1;m+((n.id,i))})
    model2C(InvarPrior(0))=0.0
    model2C(Pi(0))=plusF
    model2C.logLikelihood should be (-5810.399586 plusOrMinus 0.001)
  }

 test("Branch Length changes"){
    val (tree5,aln5)=DataParse(pfTree,pfAln.lines,new org.modiphy.sequence.SiteClassAA(5))
    val model3=InvarThmmModel(tree5)
    val start=model3.logLikelihood
    val startBL = model3(BranchLengths(0))
    val bl = startBL(0)
    startBL(0)=10.0D
  //  model3(BranchLengths)=startBL
    model3.logLikelihood should be < (start) // artificially setting a branch length to exp(1.0) should be bad for the likelihood
    startBL(0)=bl
    println(startBL)
  //  model3(BranchLengths)=startBL
    model3.logLikelihood should be (start plusOrMinus 1e-4)
  }

/*
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
  test("gamma model should give same answer whether mixture or bigmat"){
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

*/
}
