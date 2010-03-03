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
    case class Test(n:Node[_])
    class ActorTest extends Actor{
      val actorPi = new ActorPiComponent(WAG.pi,Pi(0))
      val actorS = new ActorSComponent(WAG.S,S(0))
      val branchLengthComp =  new ActorTreeComponent(tree,BranchLengths(0))
      val part1 = new BasicActorModel(actorPi,actorS,new BasicSingleExpActorModel(tree,branchLengthComp,None))
      actorPi.start
      actorS.start
      part1.start
      def act{
        loop{
          react{
            case Test(node) => part1 forward NewMatReq(node)
            }
          }
        }
      }

      val actor = new ActorTest
      actor.start
      val ans = (actor !? Test(tree.children(0))).asInstanceOf[MatReq]
      val mat1 = ans.m.get
      val me = new MatExpYang(WAG.S.sToQ(WAG.pi),WAG.pi,Some(1.0))
      val mat2 = me.exp(tree.children(0).lengthTo)

      mat1.elements.zip(mat2.elements).foreach{t=>
        t._1 should be (t._2 plusOrMinus 1E-7) //could be MatExpNormal or other impl that might not give exactly the same answer
      }

       val mat3 = (actor !? Test(tree.children(0).children(1))).asInstanceOf[MatReq].m.get
           mat3.elements.zip(me.exp(tree.children(0).children(1).lengthTo).elements).foreach{t=>
           t._1 should be (t._2 plusOrMinus 1E-7)
         }
    }
    

  }
