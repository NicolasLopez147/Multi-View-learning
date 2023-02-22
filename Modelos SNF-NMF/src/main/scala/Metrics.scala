import scala.math.{pow, sqrt}
import breeze.linalg._
import breeze.numerics._
import breeze.stats._

object Metrics {
  def euclidean_distance(x: Seq[Double], y: Seq[Double]): Double = {
    sqrt(x.zip(y).map { case (a, b) => pow(a - b, 2) }.sum)
  }

  private def centroid(points: Seq[Seq[Double]]): Seq[Double] = {
    points.transpose.map(_.sum / points.size)
  }

  def davies_bouldin_index(dataDense: DenseMatrix[Double], labels: Seq[Int]): Double = {
    val k = labels.toSet.size
    var data: List[List[Double]] = List()
    for (i <- 0 until dataDense.rows) {
      data = data :+ dataDense(i,::).t.toArray.toList
    }
    val clusters = (0 until k).map { i =>
      val points = data.zip(labels).collect { case (x, j) if j == (i+1) => x }
      (i, points, centroid(points))
    }

    val s = clusters.map { case (_, points, c) =>
      sqrt(points.map(x => pow(euclidean_distance(x, c), 2)).sum / points.size)
    }

    val m = for {
      (i, pi, ci) <- clusters
      (j, pj, cj) <- clusters
      if i != j
    } yield (i, j, (s(i) + s(j) / euclidean_distance(ci, cj)))

    m.groupBy(_._1).mapValues(_.map(_._3).max).values.sum/k
  }

  def silhouette(data: DenseMatrix[Double], labels: Seq[Int]): Double = {
    // number of points
    val n = labels.length

    val distanceMatrix = DenseMatrix.zeros[Double](n, n)

    for (i <- 0 until n) {
      for (j <- 0 until n) {
        distanceMatrix(i, j) = euclidean_distance(data(i, ::).t.toArray.toSeq, data(j, ::).t.toArray.toSeq)
      }
    }

    // calculate the average distance to all other points in the same cluster for each point
    val a = DenseVector.zeros[Double](n)
    for (i <- 0 until n) {
      val clusterPoints = labels.zipWithIndex.filter(_._1 == labels(i)).map(_._2)
      if (clusterPoints.length == 1) {
        a(i) = 0.0
      } else {
        a(i) = clusterPoints.map(j => distanceMatrix(i, j)).sum / (clusterPoints.length - 1)
      }
    }

    // calculate the average distance to all points in the closest other cluster for each point
    val b = DenseVector.zeros[Double](n)
    for (i <- 0 until n) {
      val otherClusters = labels.distinct.filter(_ != labels(i))
      if (otherClusters.isEmpty) {
        b(i) = 0.0
      } else {
        val closestCluster = otherClusters.map(c => (c, labels.zipWithIndex.filter(_._1 == c).map(_._2).map(j => distanceMatrix(i, j)).sum / labels.count( _ == c))).minBy(_._2)._1
        b(i) = labels.zipWithIndex.filter(_._1 == closestCluster).map(j => distanceMatrix(i, j._2)).sum / labels.count( _ == closestCluster)
      }
    }

    // calculate the silhouette score for each point
    val s = DenseVector.zeros[Double](n)
    for (i <- 0 until n) {
      s(i) = if (a(i) < b(i)) {
        1 - a(i) / b(i)
      } else if (a(i) > b(i)) {
        b(i) / a(i) - 1
      } else {
        0.0
      }
    }

    // return the average silhouette score
    mean(s)
  }

  def logRank(data: Seq[Seq[Double]], labels: Seq[Int]) = {
  }

  def centroidMSE(centroid: Array[Double], data: Array[Array[Double]]): Double = {
    val n = data.length
    val distances = data.map(datum => euclidean_distance(centroid, datum))
    val sumOfSquaredDistances = distances.map(d => math.pow(d, 2)).sum
    sumOfSquaredDistances / n
  }

  def psnr(dataDense: DenseMatrix[Double], labels: Array[Int]): Double = {
    val k = labels.toSet.size

    var data: List[List[Double]] = List()
    for (i <- 0 until dataDense.rows) {
      data = data :+ dataDense(i,::).t.toArray.toList
    }

    val centroids = (0 until k).map { i =>
      val points = data.zip(labels).collect { case (x, j) if j == (i+1) => x }
      (centroid(points))
    }
    val mseValues = centroids.indices.map { i =>
      val centroid = centroids(i)
      val dataForCentroid = data.zip(labels).filter { case (_, label) => label == i }.map(_._1)
      centroidMSE(centroid.toArray, dataForCentroid.toArray.map(_.toArray))
    }
    val maxSignal = data.map(_.max).max
    val mseAverage = mseValues.sum / mseValues.length
    10 * math.log10(math.pow(maxSignal, 2) / mseAverage)
  }


  def runtime() = {}
}




