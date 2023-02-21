import LoadData.DataMa
import SNF.SNF
import breeze.linalg.DenseMatrix
import breeze.numerics.round

object Main {
  def main(args: Array[String]): Unit = {
    // Cargar datos
    val (expMatriz,methyMatriz,mirnaMatriz) = LoadData.tratarDatos(Seq("aml"))


    // Implementar model SNF
    val modeloSNF = new SNF(exp = expMatriz, methy = methyMatriz, mirna = mirnaMatriz)
    val (matrizEstatusPromedio,expMatrizEstatus,methyMatrizEstatus,mirnaMatrizEstatus) = modeloSNF.aplicarSNF(m = 0.5,porcentaje = 5,iteraciones = 20)

    for (elemento <- matrizEstatusPromedio.activeIterator) println(elemento)

    //Implementar modelo NMF
    /*val modelo = Trainer.jnmf(Array(expMatriz.t,methyMatriz.t,mirnaMatriz.t),r = 50)
    println(modelo.cost(Array(expMatriz.t,methyMatriz.t,mirnaMatriz.t)))*/


  }
}