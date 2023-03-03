import breeze.linalg._
import breeze.numerics._

object KKMeans {
  var kernel_function = (x: DenseVector[Double], y: DenseVector[Double]) => (x dot y);

  def set_kernel_function(kernel_function_type: String) : Unit = {
    kernel_function_type match {
      case "linear" => this.kernel_function = (x: DenseVector[Double], y: DenseVector[Double]) => (x dot y)
      case "polynomial" => this.kernel_function = (x: DenseVector[Double], y: DenseVector[Double]) => pow((x dot y) + 1, 3)
      case "rbf" => this.kernel_function = (x: DenseVector[Double], y: DenseVector[Double]) => exp(-pow(norm(x - y), 2))
      case "sigmoid" => this.kernel_function = (x: DenseVector[Double], y: DenseVector[Double]) => tanh((x.t * y) + 1)
      case _ => throw new IllegalArgumentException("Invalid kernel function: " + kernel_function)
    }
  }

  def kernel_matrix(data: DenseMatrix[Double]) : DenseMatrix[Double] = {
    val n = data.rows
    val kernel_matrix = DenseMatrix.zeros[Double](n, n)
    for (i <- 0 until n) {
      for (j <- 0 until n) {
        kernel_matrix(i, j) = kernel_function(data(i, ::).t, data(j, ::).t)
      }
    }
    kernel_matrix
  }

  def kernel_k_means(data: DenseMatrix[Double], k: Int, max_iterations: Int = 100): List[Int] = {
    val kernel = kernel_matrix(data)
    val n = data.rows
    val m = data.cols
    val centroids = DenseMatrix.zeros[Double](k, m)
    val labels = DenseVector.zeros[Int](n)
    val old_labels = DenseVector.zeros[Int](n)

    // Initialize centroids
    for (i <- 0 until k) {
      centroids(i, ::) := data(i, ::)
    }

    // Initialize labels
    for (i <- 0 until n) {
      labels(i) = i % k
    }

    var iterations = 0
    while (iterations < max_iterations && labels != old_labels) {
      // Save old labels
      old_labels := labels

      // Update centroids
      for (i <- 0 until k) {
        val cluster = (0 until n).filter(labels(_) == i)
        if (cluster.nonEmpty) {
          centroids(i, ::) := mean(kernel(::, cluster), Axis._1).t
        }
      }

      // Update labels
      for (i <- 0 until n) {
        val distances = DenseVector.zeros[Double](k)
        for (j <- 0 until k) {
          distances(j) = kernel(i, j) - 2 * (centroids(j, ::) dot data(i, ::).t) + (centroids(j, ::) dot centroids(j, ::).t)
        }
        labels(i) = argmin(distances)
      }

      iterations += 1
    }

    labels.toArray.toList
  }
}
