import org.scalatest.FunSuite
import org.scalatest.matchers.ShouldMatchers._
import org.modiphy.math._
import org.modiphy.tree._
import ModelData._


class MathSuite extends FunSuite {
  val gamma = new Gamma(4)
  test ("Gamma calculation should match PAML"){
    gamma(0.5).toList.zip(List(0.033388, 0.251916, 0.820268, 2.894428)).foreach{t=>
      t._1 should be (t._2 plusOrMinus 0.001)
    }
  }
}
