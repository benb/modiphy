package org.modipht.newtree
import org.modiphy.math._
import scala.actors.Actor._ 
import scala.actors._ 
import org.modiphy.math.EnhancedMatrix._

case object BranchLength
case class UpdateDist(d:Double)
case class Unclean[A <: BioEnum](direction:Node[A])
case object UpdateMat
case class LikelihoodCalcDir[A <: BioEnum](model:ActorModelComponend,b:DirBranch[A])
case class CalculatedPartialLikelihoods(pl:List[Vector])

class DirBranch[A <: BioEnum](val to:Node[A],var dist:Double,val from:Node[A]) extends Actor{
  def act{
    main(None)
  }
  //partial likelihoods go up the tree - so from -> to
  def main(pl:Option[List[Vector]]){
    react{
      case BranchLength=>
        reply(dist)
        main(pl)
      case UpdateDist(d)=>
        dist=d
        reply('ok)
        main(None)
      case UpdateMat=>
        to !? Unclean(this)
        reply('ok)
        main(None)
      case LikelihoodCalc(model)=>
        if (!pl.isDefined){
          //get PL up to 'descendent' node
          from ! LikelihoodCalc(model)
          //get Matrix
          model ! NewMatReq(this)
          getAns(sender)
        }else {
          sender ! CalculatedPartialLikelihoods(pl.get)
          main(pl)
        }
      }
    }
  def getAns1(replyTo:OutputChannel[Any]){
    react{
      case MatReq(_,mat,pi)=>
        getAns2(replyTo,mat,pi)
    }
  }
  def getAns2(replyTo:OutputChannel[Any],mat:Matrix,pi:Vector){
    react{
      case CalculatedPartialLikelihoods(pl)=>
        val ans = BasicLikelihoodCalc.partialLikelihoodCalc(pl2,eMat)
        replyTo ! CalculatedPartialLikelihoods(ans)
        main(Some(ans))
    }
  }
}

class Branch[A <: BioEnum](a:Node[A],b:Node[A],var dist:Double) extends Actor{
  val endA = new DirBranch(a,dist,b)
  val endB = new DirBranch(b,dist,a)
  def act{
    loop{  
      react{
        case BranchLength=>reply(dist)
        case UpdateDist(d)=>
          a !? UpdateDist(d)
          b !? UpdateDist(d)
          reply('ok)
      }
    }
  }
}


class Node[A]{
}
class INode[A](val id:Int) extends Node[A]{
  private var branchEnds:List[DirBranch[A]]=Nil

  def act{
    main(new scala.collection.immutable.HashMap[DirBranch[A],List[Vector]]())
  }
  type Cache=Map[Option[DirBranch[A]],List[Vector]] 
  def main(m:Map[Option[DirBranch[A]],List[Vector]]){
    react{
      case Unclean(branch)=>
        branchEnds.filter{_ != branch}.foreach{c=> c !? Unclean(this)}
        reply('ok)
        main(m)
      case LogLikelihoodCalc(model,actor) =>
        if (m.contains(None)){
          sender ! BasicLikelihoodCalc.logLikelihood(m(None),model.pi)
        }else {
          var i = 0
          branchEnds.filter(_ != branch).foreach{c => c ! LikelihoodCalc(model);i=i+1}
          logLikelihood(m,Nil,i,sender,model.pi)
        }
      case LikelihoodCalcDir(model,branch)=>
        if (m.contains(branch)){
          sender ! CalculatedPartialLikelihoods(m(Some(branch)))
        }else {
        //request partial likelihoods along each branch
          var i=0
          branchEnds.filter(_ != branch).foreach{c => c ! LikelihoodCalc(model);i=i+1}
          partialLikelihoods(branch,m,Nil,i,sender)
        }
      }
    }
    
    def logLikelihood(cache:Cache,plList:List[List[Vector]],toDo:Int,replyTo:OutputChannel[Any],pi:Vector){
      if (toDo==0){
        val ans = BasicLikelihoodCalc.combinePartialLikelihoods(plList)
        replyTo ! BasicLikelihoodCalc.logLikelihood(plList,pi)
        main(cache + ((None,ans)))
      }else {
        react{
          case CalculatedPartialLikelihoods(pl)=>
            logLikelihood(cache,pl::plList,toDo-1,replyTo,pi)
        }
      }
    }

    def partialLikelihoods(branch:DirBranch[A],cache:Cache,plList:List[List[Vector]],toDo:Int,replyTo:OutputChannel[Any]){
      if (toDo==0){
        val ans = BasicLikelihoodCalc.combinePartialLikelihoods(plList)
        replyTo ! CalculatedPartialLikelihoods(ans)
        main(cache + ((Some(branch),ans)))
      }else {
        react{
          case CalculatedPartialLikelihoods(pl) =>
            partialLikelihoods(branch,cache,pl::plList,toDo-1,replyTo)
        }
      }
    }

}
class Leaf[A](val id:Int) extends Node[A]{
  var branchEnd:Option[DirBranch[A]]=None
}
