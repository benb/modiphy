import org.scalatest.FunSuite
import org.scalatest.matchers.ShouldMatchers._
import org.modiphy.math._
import org.modiphy.tree._
import org.modiphy.math.EnhancedMatrix._
import ModelData._


class ParamSuite extends FunSuite {
  class Obs {
    var count = 0
    def receiveUpdate(s:Subject){
      count=count+1
    }
  }
  import Math._
  test("FirstOnlyPiParam should be sane") {
    val pi = Vector(10)
    for (i<-0 until pi.size){pi(i)=0.1}
    val obs = new Obs
    val param = new FirstOnlyPiParam(pi,"first param") with LogParamControl
    param.addObserver(obs)

    param.name should be ("first param")
    param.setParams(Array(log(0.5D)))
    obs.count should be (1)
    for (i <- 1 until pi.size){
      pi(i) should be (0.5D/9.0D plusOrMinus 1e-7)
    }
    param.setParams(Array(log(0.01D)))
    obs.count should be (2)
    for (i <- 1 until pi.size){
      pi(i) should be (0.99D/9.0D plusOrMinus 1e-7D)
    }
  }
    
}
