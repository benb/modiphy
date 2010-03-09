import org.scalatest.FunSuite
import org.scalatest.matchers.ShouldMatchers._
import org.modiphy.math._
import org.modiphy.tree._
import org.modiphy.math.EnhancedMatrix._
import ModelData._
import scala.actors.Actor
import scala.actors.Actor._


class ActorModelSuite extends FunSuite {
  val (tree,data)=DataParse(treeStr,alnStr.lines,new org.modiphy.sequence.SiteClassAA(1))

  test ("Actor Model should give correct matrix exp"){
    case class Test(b:Branch[_])
    class ActorTest extends Actor{
      tree.startTree
      val actorPi = new ActorPiComponent(WAG.pi,Pi(0))
      val actorS = new ActorSComponent(WAG.S,S(0))
      val branchLengthComp =  new ActorTreeComponent(tree,BranchLengths(0))
      val part1 = new BasicActorModel(actorPi,actorS,new NormaliserActorModel(new BasicSingleExpActorModel(tree,branchLengthComp,None)))
      actorPi.start
      actorS.start
      part1.start
      def act{
        loop{
          react{
            case Test(branch) => part1 forward NewMatReq(branch)
            }
          }
        }
      }

      val actor = new ActorTest
      actor.start
      val ans = (actor !? Test(tree.bList.head.myBranch)).asInstanceOf[MatReq]
      val mat1 = ans.m.get
      val me = new MatExpYang(WAG.S.sToQ(WAG.pi),WAG.pi,Some(1.0))
      val mat2 = me.exp(tree.bList.head.dist)

      mat1.elements.zip(mat2.elements).foreach{t=>
        t._1 should be (t._2 plusOrMinus 1E-7) //could be MatExpNormal or other impl that might not give exactly the same answer
      }

      val b1= tree.bList.filter{!_.down.isInstanceOf[Leaf[_]]}.head
      val b2 = b1.down.bList.filter{_.myBranch != b1.myBranch}.head

       val mat3 = (actor !? Test(b2.myBranch)).asInstanceOf[MatReq].m.get
           mat3.elements.zip(me.exp(b2.dist).elements).foreach{t=>
           t._1 should be (t._2 plusOrMinus 1E-7)
         }
    }
    

  }
