import breeze.linalg.{DenseMatrix, csvread}
import org.scalatest.funsuite.AnyFunSuite

class KKMeansTest extends AnyFunSuite {

  test("KKMeans debería calcular el número correcto de centroides para los datos de entrada") {
    val entrada: DenseMatrix[Double] = csvread(new java.io.File("src/test/resources/clusters.csv"))
    val k: Int = 2
    val kkMeans = new KKMeans()
    
    // Set the kernel function and compute the kernel matrix
    kkMeans.set_kernel_function("rbf")
    val kernelMatrix = kkMeans.kernel_matrix(entrada)

    // Compute the initial centroids and run the kernel k-means algorithm
    val initialCentroids = kkMeans.k_means_pp(kernelMatrix, k)
    val labels = kkMeans.kernel_k_means(kernelMatrix, k)

    // Count the number of unique labels
    val numCentroids = labels.toSet.size

    // Assert that the number of centroids is equal to k
    assert(numCentroids == k)
  }

}