import LoadData.DataMa
import SNF.SNF
import breeze.linalg.DenseMatrix
import breeze.numerics.round

object Main {
  def main(args: Array[String]): Unit = {
    // Cargar datos
    val (expMatriz,methyMatriz,mirnaMatriz) = LoadData.tratarDatos(Seq("aml"))

    KKMeans.set_kernel_function("polynomial")
    // Implementar model SNF
    val modeloSNF = new SNF(exp = expMatriz, methy = methyMatriz, mirna = mirnaMatriz)
    val (matrizEstatusPromedio,expMatrizEstatus,methyMatrizEstatus,mirnaMatrizEstatus, expKernelDispersa,methyKernelDispersa,mirnaKernelDispersa) = modeloSNF.aplicarSNF(m = 0.5,porcentaje = 5,iteraciones = 20)
    
    val exp_est_labels = KKMeans.kernel_k_means(expMatrizEstatus, 20, 100)
    val exp_kernel_labels = KKMeans.kernel_k_means(expKernelDispersa, 20, 100)
    println("Metrica \t Kernel \t Estatus")
    print("Silhouette \t")
    print(Metrics.silhouette(expMatrizEstatus, exp_est_labels))
    print("\t")
    println(Metrics.silhouette(expKernelDispersa, exp_kernel_labels))
    print("Davies-Bouldin \t")
    print(Metrics.davies_bouldin_index(expMatrizEstatus, exp_est_labels))
    print("\t")
    println(Metrics.davies_bouldin_index(expKernelDispersa, exp_kernel_labels))
    print("PSNR \t")
    print(Metrics.psnr(expMatrizEstatus, exp_est_labels))
    print("\t")
    println(Metrics.psnr(expKernelDispersa, exp_kernel_labels))
    print("Ball-Hall \t")
    print(Metrics.ball_hall(expMatrizEstatus, exp_est_labels))
    print("\t")
    println(Metrics.ball_hall(expKernelDispersa, exp_kernel_labels))
    print("Calinski-Harabasz \t")
    print(Metrics.calinski_harabasz(expMatrizEstatus, exp_est_labels))
    print("\t")
    println(Metrics.calinski_harabasz(expKernelDispersa, exp_kernel_labels))
    print("Hartigan \t")
    print(Metrics.hartigan(expMatrizEstatus, exp_est_labels))
    print("\t")
    println(Metrics.hartigan(expKernelDispersa, exp_kernel_labels))
    print("Xu \t")
    print(Metrics.xu(expMatrizEstatus, exp_est_labels))
    print("\t")
    println(Metrics.xu(expKernelDispersa, exp_kernel_labels))

    //Implementar modelo NMF
    /*val modelo = Trainer.jnmf(Array(expMatriz.t,methyMatriz.t,mirnaMatriz.t),r = 50)
    println(modelo.cost(Array(expMatriz.t,methyMatriz.t,mirnaMatriz.t)))*/


    /*val points = DenseMatrix(
      (1.0, 1.0),
      (3.0, 3.0),
      (9.0, 9.0),
      (7.0, 7.0)
    )

    val labels = Seq(1,1,2,2)

    println(Metrics.silhouette(points, labels))*/
  }
}