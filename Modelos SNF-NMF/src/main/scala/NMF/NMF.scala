import breeze.linalg.{DenseMatrix, sum}
import breeze.numerics._
import org.apache.spark.{SparkConf, SparkContext}
import breeze.linalg.DenseVector

object Trainer {
  def jnmf(xs: Array[DenseMatrix[Double]], r: Int = 2, n: Int = 50000, eps: Double = 0.0001, epsEval: Int = 1): JNMFModel = {
    val hs0 = xs.map(x => DenseMatrix.rand[Double](r, x.cols))
    val w0 = DenseMatrix.rand[Double](xs(0).rows, r)
    val cost0 = JNMFModel(w0, hs0).cost(xs)

    val (_, w, hs, _) = (1 to n).foldOrStop((xs, w0, hs0, cost0)) { case ((xs, w, hs, cost), i) =>
      val nhs = hs.zip(xs).map { case (h, x) =>
        h *:* ((w.t * x) /:/ ((w.t * w) * h))
      }
      val (num, den) =
        nhs.zip(xs).foldLeft((DenseMatrix.zeros[Double](w.rows, w.cols), DenseMatrix.zeros[Double](w.cols, w.cols))) {
          case ((num, den), (h, x)) => (num :+= x * h.t, den :+= h * h.t)
        }
      val nw = w *:* (num /:/ (w * den))
      if (i % epsEval == 0) {
        var currCost = JNMFModel(w, hs).cost(xs)
        val cmp = (cost-currCost)/currCost
        println(s"$i: $currCost $cost $cmp")
        if (cmp <= eps && i > 1)
          None
        else
          Some((xs, nw, nhs, currCost))
      } else {
        Some((xs, nw, nhs, cost))
      }
    }
    JNMFModel(w, hs)
  }

  implicit class FancyIterable[A](it: Iterable[A]) {
    def foldOrStop[B](zero: B)(f: (B, A) => Option[B]): B = {
      val ii = it.iterator
      var b = zero
      while (ii.hasNext) {
        val x = ii.next
        val nb = f(b, x)
        if (nb == None) return b
        b = nb.getOrElse(b)
      }
      b
    }
  }

  def clustering(jnmf: JNMFModel): (DenseMatrix[Double], Seq[Int]) = {
    val w = jnmf.w
    val labels = DenseVector.zeros[Int](w.rows)
    for(i <- 0 until w.rows) {
      val max = w(i, ::).inner.toArray.max
      for(j <- 0 until w.cols) {
        if(w(i, j) == max)
          labels(i) = j+1
      }
    }
    (w, labels.toArray.toSeq)
  }
}

case class JNMFModel(w: DenseMatrix[Double], hs: Array[DenseMatrix[Double]]) {
  def inverse(): Array[DenseMatrix[Double]] = {
    hs.map(w * _)
  }

  def cost(xs: Array[DenseMatrix[Double]]): Double = {
    xs.zip(hs).foldLeft(0d) { case (cost, (x, h)) => cost + sum(pow(x - w * h, 2)) }
  }
}