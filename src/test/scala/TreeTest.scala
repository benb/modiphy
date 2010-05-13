import org.scalatest.FunSuite
import org.scalatest.matchers.ShouldMatchers._
import org.modiphy.math._
import org.modiphy.tree._
import ModelData._
import org.modiphy.math.EnhancedMatrix._


class TreeSuite extends FunSuite {
   val (tree,aln) = DataParse(
    "(Chlamydomona:0.082441#0,Chlorella:0.050529#1,((Euglena:0.099344#5,Nephroselmis:0.088455#10):0.010341#27,(((Mycobacteriu:0.527327#8,(Mycoplasma:0.338658#9,(Ecoli:0.214706#4,Rickettsia:0.253023#13):0.050486#16):0.061350#17):0.143567#18,(Synechocysti:0.079101#14,(Thermosynech:0.057044#15,(Gloeobacter:0.179811#6,Prochlorococ:0.109110#12):0.038575#19):0.009400#20):0.043594#21):0.033610#22,(Guillardia:0.067446#7,(Porphyra:0.064926#11,(Cyanidioschy:0.076726#2,Cyanidium:0.107181#3):0.028576#23):0.013723#24):0.010431#25):0.060348#28):0.000000#26);",
    """>Chlamydomona
GRIVQIIGPVLDIVF
>Chlorella
GRITQIIGPVLDVSF
>Cyanidioschy
GKVVQIIGPVLDIEF
>Cyanidium     
GKVIQIIGPVLDIVF
>Ecoli        
GKIVQVIGAVVDVEF
>Euglena
GIILQIIGPVMDISF
>Gloeobacter
GYITQIIGPVVDAEF
>Guillardia
GYITQIIGPVVDVEF
>Mycobacteriu
GRVVRVTGPVVDVEF
>Mycoplasma 
GKVHQVIGPVVDVIF
>Nephroselmis 
GKIVQIVGPVMDVAF
>Porphyra  
GSVTQIIGPVLDIAF
>Prochlorococ
GVIRQVIGPVLDVEF
>Rickettsia 
GKITQVISAVVDVKF
>Synechocysti 
GKITQVIGPVIDAQF
>Thermosynech
GFITQVIGPVVDIEF
""".lines
,new org.modiphy.sequence.SiteClassAA(1))

   test("Node IDs"){
     tree.children.map{_.label.get}.sort{_ < _} should be (List(0,1,26))
     tree.children.sort{_.label.get > _.label.get}.head.children.map{_.label.get}.sort{_ < _} should be (List(27,28))
   }
  test("Branch Lengths"){
     val t2 = tree.setBranchLengths(tree.getBranchLengths)
     println(tree.getBranchLengths)
     println(t2.getBranchLengths)
     t2.getBranchLengths should be (tree.getBranchLengths)
     t2.toString should be (tree.toString)
   }

 }
