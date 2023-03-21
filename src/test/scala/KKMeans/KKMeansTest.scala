import breeze.linalg.{DenseMatrix, csvread}
import org.scalatest.funsuite.AnyFunSuite
import KKMeans.KKMeans

class KKMeansTest extends AnyFunSuite {

  test("KKMeans debería calcular el número correcto de centroides para los datos de entrada") {
    val entrada: DenseMatrix[Double] = csvread(new java.io.File("src/test/resources/clusters.csv"))
    val k: Int = 2
    
    // Set the kernel function and compute the kernel matrix
    KKMeans.set_kernel_function("rbf")
    val kernelMatrix = KKMeans.kernel_matrix(entrada)

    // Compute the initial centroids and run the kernel k-means algorithm
    val initialCentroids = KKMeans.k_means_pp(kernelMatrix, k)
    val labels = KKMeans.kernel_k_means(kernelMatrix, k)

    // Count the number of unique labels
    val numCentroids = labels.toSet.size

    // Assert that the number of centroids is equal to k
    assert(numCentroids == k)
  }

}