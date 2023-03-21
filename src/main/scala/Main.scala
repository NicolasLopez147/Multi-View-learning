import SNF.SNF
import NMF.Trainer
import KKMeans.KKMeans
import breeze.linalg.{DenseMatrix, sum}
import breeze.numerics._
import breeze.linalg.DenseVector

/**
  * This is an example of how to use the SNF model and NMF model. To compare the results of the SNF model with the NMF model,
  * the Kernel K-Means model is used to cluster the data. The results are compared using the metrics: Silhouette, Davies-Bouldin,
  * PSNR, Ball-Hall, Calinski-Harabasz, Hartigan and Xu.
  */
object Main {
  // Number of iterations for the Kernel K-Means model
  val iterations = 100
  // Number of clusters for the Kernel K-Means model
  val k = 30

  /**
    * main method
    *
    * @param args
    */
  def main(args: Array[String]): Unit = {
    //creaate 3 random matrices to use nmf
    val X: Array[DenseMatrix[Double]] = Array(
    DenseMatrix.rand[Double](50, 50),
    DenseMatrix.rand[Double](50, 50)
  )

    // Implement JNMF model
   val modeloEntrenado = Trainer.jnmf(X,r=6, eps = 0.0001 )
    val X_reconstructed = modeloEntrenado.inverse()

    //subtract each matrix from the reconstructed matrix
    val X1_diff = X(0) - X_reconstructed(0)
    val X2_diff = X(1) - X_reconstructed(1)


    // Apply JNMF
    val suma = sum(abs(X1_diff))
    val suma2 = sum(abs(X2_diff))
    println("Suma: " + suma)
    println("Suma2: " + suma2)
  }

  def mainTemp():Unit = {
 // Load data
    val (expMatriz,methyMatriz,mirnaMatriz) = LoadData.tratarDatos(Seq("resources"))

    // Implement SNF model
    val modeloSNF = new SNF(exp = expMatriz, methy = methyMatriz, mirna = mirnaMatriz)

    // Apply SNF
    val (
      matrizEstatusPromedio,
      expMatrizEstatus,
      methyMatrizEstatus,
      mirnaMatrizEstatus,
    ) = modeloSNF.aplicarSNF(m = 0.5,porcentaje = 5,iteraciones = 20)

    // Apply Kernel K-Means
    val est_prom_labels = KKMeans.kernel_k_means(matrizEstatusPromedio, k, iterations)

    // Build a table with the results based on the metrics
    val table:Seq[Seq[Any]] = Seq(
      Seq("Metrica", "Kernel"),
      Seq("Silhouette", Metrics.silhouette(matrizEstatusPromedio, est_prom_labels)),
      Seq("Davies-Bouldin", Metrics.davies_bouldin_index(matrizEstatusPromedio, est_prom_labels)),
      Seq("PSNR", Metrics.psnr(matrizEstatusPromedio, est_prom_labels)),
      Seq("Ball-Hall", Metrics.ball_hall(matrizEstatusPromedio, est_prom_labels)),
      Seq("Calinski-Harabasz", Metrics.calinski_harabasz(matrizEstatusPromedio, est_prom_labels)),
      Seq("Hartigan", Metrics.hartigan(matrizEstatusPromedio, est_prom_labels)),
      Seq("Xu", Metrics.xu(matrizEstatusPromedio, est_prom_labels)),
    )

    // Print the table for the SNF model
    printTable(table);

    // Implement NMF model
    val NMFmodel = Trainer.jnmf(Array(expMatriz, methyMatriz, mirnaMatriz), r = k)

    // Apply NMF
    val W_matrix = NMFmodel.w

    // Apply Kernel K-Means
    val w_labels = KKMeans.kernel_k_means(W_matrix, W_matrix.cols, iterations)

    // Build a table with the results based on the metrics
    val table2:Seq[Seq[Any]] = Seq(
      Seq("Metrica", "Kernel"),
      Seq("Silhouette", Metrics.silhouette(W_matrix, w_labels)),
      Seq("Davies-Bouldin", Metrics.davies_bouldin_index(W_matrix, w_labels)),
      Seq("PSNR", Metrics.psnr(W_matrix, w_labels)),
      Seq("Ball-Hall", Metrics.ball_hall(W_matrix, w_labels)),
      Seq("Calinski-Harabasz", Metrics.calinski_harabasz(W_matrix, w_labels)),
      Seq("Hartigan", Metrics.hartigan(W_matrix, w_labels)),
      Seq("Xu", Metrics.xu(W_matrix, w_labels)),
    )

    // Print the table for the NMF model
    printTable(table2);

  }

  /**
    * Print a table with the results of the metrics
    *
    * @param table
    */
  def printTable[T](table: Seq[Seq[T]]): Unit = {
    val colWidths = table.transpose.map(_.map(_.toString.length).max)
    table.foreach { row =>
      println(row.zip(colWidths).map { case (elem, width) => elem.toString.padTo(width, ' ') }.mkString(" "))
    }
  }
}