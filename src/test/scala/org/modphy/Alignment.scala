package org.modphy.math
import  org.modphy.sequence._
import org.scalatest.Suite
import org.modphy.tree._
import org.modphy.math._
import org.modphy.math.EnhancedMatrix._
import org.apache.commons.math.optimization.direct._



class AlignmentSuite extends Suite{


  val maffile = """##maf version=1 scoring=autoMZ.v1
a score=297559.000000
s hg18.chr6                  5000 364 + 170899992 GATCTTATATAACTGTGAGATTAATCTCAGATAATGACACAAAATATAGTGAAGTTGGTAAGTTATTTAGT-AAAGCTCATGAAAATTGTGCCCTCCATTCCCATATAATTTATTAATTGTCTAGGAACTTCCACATACATTGCCTCAATTTATCTTTCAACAAC-TTGTGTGTTATATTTTGGAATACAGATACAAAGTTA---TTATGCTTTCAAAATATTCTTTTG-CTAATTCTTAGAACAAAGAAAGGCATA--------AATATATTAGTATTTGTGTACATCTGTTCCTTCCTGTGTGACC---CTAAGTTTAGTAGAAGAAAGGAGAGAAAATATAGC---CTAGCttataaatttaaaaaaaaatttatttGGTC
s rheMac2.chr7           87625102 364 - 169801366 -ATCTTATATCACTGTGCGATTAATCTCAGATAATGACATAAAATATAGTGAAGTTGGTAAGTTATTTAGTAAAAGCTCATGAAAATTGTGCCTTTCATTCCCATATAATTTAGTAATCGTCTAGGAACTTTCACATACACTGCCTCAATTTATCTTTCAACAAC-TTGTGTGTTATATTTTGGAATATAGATACAAACTTA---TCATGCTTTCAAAATATTATTTTG-TTAATTCTTAAAACAAAGAAAGGCATA--------AACA---TAGTATTTGTGTACACCTGTACCTTCCTGTGTGTTCAGTAGAAGTTTAGTAGAAGAAAGGAGAGAAAATATAGC---CTAGCTTATAAATTTTTTACTATTTTTATTTGACC
q rheMac2.chr7                                    -99999999999999999999999999999999999999999999999999889999999999999999999999999999999999999999999999999999999999999999999999999999999999999996999999999999999999999999-999999999999999999999999999999999999---999999999999999999999999-999999999999999999999999999--------9999---99999999999999999999999999999999999999999999999999999999999999999999999999---99999999999999999999999999999999999
i rheMac2.chr7           N 0 C 0
s panTro2.chrUn           9727052 363 +  58616431 GATCTTATATAACTGTGAGATTAATCTCAGATAATGACACAAAATATAGTGAAGTTGGTAAGTTATTTAGT-AAAGCTCATGAAAATTGTGCCTTCCATTCCCATATAATTTAGTAATTGTCTAGGAACTTCCACATACATTGCCTCAATTTATCTTTCAACAAC-TTGTGTGTTATATTTTGGAATACAGATACAAAGTTA---TTATGCTTTCAAAATATTCTTTTG-CTAATTCTTGGAACAAAGAAAGGCATA--------AATATATTAGTATTTGTGTACACCTGTTCCTTCCTGTGTGATC---CTAAGTTTAGTAGAAGAAAGGAGAGAAAATATAGC---CTAGCttataaatt-aaaaaaaaatttatttGGTC
q panTro2.chrUn                                   99999999999999999999999999999999999999999999999999999999999999999999999-999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999-999999999999999999999999999999999999---999999999999999999999999-999999999999999999999999999--------9999999999999999999999999999999999999999999---99999999999999999999999999999999999---99999999999999-99999999999999999999
i panTro2.chrUn          N 0 C 0

a score=101230.000000
s hg18.chr6                  5364 153 + 170899992 CATTTTGTGAAAAA----CATAAAAAAAGAACTGTCACAT-CTTAATTTAAAAAATATATGCTTAGTGGTAAG--GAGATATATGTCAACTTT
TAAGAGG-TTGAAAAACAAACGCCTCCCATTATA--AGTT---TATACTTCA-CCTCCCACCACTATAACAACC
s rheMac2.chr7           87625466 158 - 169801366 CATTTTGTGAAAAA----CATAAAAAAAGAACTGTCACATGCTTAATTTAAAAAATATATGCTTATTGGTAAG--GAGATATATGTCAACTTT
TAAGAGGTTTGAAAAACAAATGCCTCCCATTATAAGCGTT---CATACTTCACCCTCCCACCACTGTAACAACC
q rheMac2.chr7                                    89999999999999----9999999999999999999999999999999999999999999999999999999--999999999999999999
9999999999999999999999999999999999999999---9999999999999999999999999999999
i rheMac2.chr7           C 0 C 0
s panTro2.chrUn           9727415 153 +  58616431 CATTTTGTGAAAAA----CATAAAAAAAGAACTGTCACAT-CTTAATTTAAAAAATACATGCTTAGTGGTAAG--GAGATATATGTCAACTTT
TAAGAGG-TTGAAAAACAAACGCCTCCTATTATA--AGTT---TATACTTCA-CCTCCCACCACTATAACAACC
q panTro2.chrUn                                   99999999999999----9999999999999999999999-99999999999999999999999999999999--999999999999999999
9999999-99999999999999999999999999--9999---999999999-999999999999999999999
i panTro2.chrUn          C 0 C 0
s ponAbe2.chr15             15026 155 -  99152023 CATTTTGTGAAAAA----CATAAAAAAAGAACTGTCACAT-CTTAATTTAAAAAATATATGCTTAGTGGTAAG--GAGATATATGTCAACTTT
TAAGAGG-TTGAAAAACAAACGCCTCCCATTATAAGAGTT---TATACTTCA-cctcccaccactataacaacc
q ponAbe2.chr15                                   99999999999999----9999999999999999999999-99999999999999999999999999999999--999999999999999999
9999999-99999999999999999999999999999999---999999999-999999999999999999999
i ponAbe2.chr15          C 0 C 0
s micMur1.scaffold_35565     1629 159 +     10343 CATTTTCTTTAAAAACTTCAAGGAAAAAAAACTGCCATATGCCTAACTCAAAA---ATATATTTAGTGATAAGGTGATATATATGTCTACTTT
TAGGAGG-CTTAAAAAATAACTCTTCTCATTGCAATAGTG---CATATTTCA-TGTCCTACCACCATAGCACCC
q micMur1.scaffold_35565                          99999999999999999998999999999999996979999999999999999---9999999999799556655894446766646565699
6613734-55649999765435524555444233546335---453449442-463554432322244245272
i micMur1.scaffold_35565 I 1154 N 0
s choHof1.scaffold_37150     5108 123 -     13819 catTTGACCTAAAAACCTCAAGTAAGAAACTCCTTCCTGT-CTTGCTC---AAAATTTTTGCTTAATGTTATG--------------------
-----------AAGAAAAAATTCCTTTCATTGTAACATTTTGCTATAGTTCA-CTTCCTATCATTA--------
q choHof1.scaffold_37150                          9999999999999999999999999999999999999999-9999999---9999999999999999999999--------------------
-----------99999999999999999999999999999999999999999-9999999999999--------
i choHof1.scaffold_37150 C 0 C 0

"""
      def testAlign{

        val mafReader = new Maf(maffile.split("\n").elements.buffered)

        val alignment1=mafReader.next
        val alignment2=mafReader.next
        val tree = """(((((((((((((((((hg18:0.003731,panTro2:0.005501):0.013010,gorGor1:0.01):0.01,
          ponAbe2:0.020000):0.020000,rheMac2:0.031571):0.01,calJac1:0.01):0.020000,
        tarSyr1:0.100000):0.050000,(micMur1:0.084110,
        otoGar1:0.145437):0.033956):0.020000,tupBel1:0.203975):0.020000,
      (((((mm9:0.104920,rn4:0.109421):0.020000,dipOrd1:0.200000):0.050000,
        cavPor3:0.150000):0.050000,speTri1:0.150000):0.020000,(oryCun1:0.100000,
      ochPri2:0.200000):0.050000):0.020000):0.020000,(((vicPac1:0.100000,
      (turTru1:0.120000,bosTau4:0.162368):0.050000):0.050000,((equCab2:0.150000,
        (felCat3:0.098674,canFam2:0.114682):0.002783):0.050000,(myoLuc1:0.142600,
        pteVam1:0.142460):0.003381):0.007170):0.030000,(eriEur1:0.279121,
      sorAra1:0.309867):0.023929):0.040000):0.01,(((loxAfr2:0.110021,
      proCap1:0.150000):0.003000,echTel1:0.265218):0.066339,(dasNov2:0.178799,
    choHof1:0.200000):0.090000):0.053170):0.213469,monDom4:0.320721):0.088647,
  ornAna1:0.488110):0.118797,((galGal3:0.23,taeGut1:0.16):0.16,
  anoCar1:0.513962):0.093688):0.151358,xenTro2:0.778272):0.174596,
(((tetNig1:0.203933,fr2:0.239587):0.203949,(gasAcu1:0.314162,
  oryLat2:0.501915):0.055354):0.346008,danRer5:0.730028):0.174596):0.1,
petMar1:0.1);"""

       import DNA._

       val aln = alignment1.seqs.foldLeft(Map[String,String]()){(i,j)=>i+((j._1.replaceAll("\\..*",""),j._2))} //get rid of chromosome info
       val data = DataParse.dropNodes[DNA.type](tree,aln,DNA)

         println("TREE: " + data._1)

       val pi=Vector(4)
       for (i<-DNA.matElements){
         pi(i)=0.25
       }

             val sMat = Matrix(4,4)
       sMat assign 1
       val model = new EnhancedModel[DNA.type](pi,sMat,data._1) with OptBranchLengths[DNA.type] 

       val testModel1 = model.logLikelihood
       println("__________________MODEL1")
       println(model.pi)
       println(model.sMat)
       model.pi(0)=0.001
       model.pi(1)=0.001
       model.pi(2)=0.001
       model.pi(3)=1.0-model.pi(0)-model.pi(1)-model.pi(2)

       println("__________________MODEL2")
       println(model.pi)
       println(model.sMat)

       val testModel2 = model.logLikelihood
       println("Likelihoods " + testModel1 + " " + testModel2)

       assert(testModel2<testModel1)



       println("START "  + model)

       println("FIT: " + model.toFit.toList.mkString(","))
       assert (!model.toFit.contains{i:Double => i < -0.00001} && !model.toFit.contains{i:Double=> i < 0.00001})
       val test = model.toFit.toArray
       test(2)=1.0D
       println("FIT2: " + test.mkString(","))
       model setPi test
       println("PI2: " + model.pi)
       val pis = model.pi.toList.sort{_>_}
       println("PI2B: " + pis.mkString(","))
       assert(pis(0)>0.25)
       assert(!(pis.tail.contains{i:Double=>i>0.25}))
       

       ModelOptimiser.nelderMead[DNA.type](model)

       val lkl1=model.logLikelihood
       val branchLengths = model.getParams(2)
       model.setParams(2)(branchLengths)
       assert((model.logLikelihood-lkl1).abs < 0.0000001)

       println("OPIMISED " + model)

        //alignment is AT biased
       assert(model.pi(A) > 0.33)
       assert(model.pi(T) > 0.33)

       //GAMMA
       object GammaDNA extends SiteClassDNA(4)
       val gammaModel = new GammaModel[GammaDNA.type](model.pi.copy.assign(0.25),model.sMat.copy.assign(1),data._1.setNewDataType[GammaDNA.type](GammaDNA),1.0) with OptBranchLengths[GammaDNA.type]

       println("BaseLkl " + model.logLikelihood)

       ModelOptimiser.nelderMead[GammaDNA.type](gammaModel)

       println("GammaLkl " + gammaModel.logLikelihood)

       //gammaModel.alpha=0.1


       println(model)
       println(model.logLikelihood)
       println(gammaModel)
       println(gammaModel.logLikelihood)


       assert(gammaModel.logLikelihood > model.logLikelihood)

       
       assert(mafReader.hasNext==false)

      }
    }

