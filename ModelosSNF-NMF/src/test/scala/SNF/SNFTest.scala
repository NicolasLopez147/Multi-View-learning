import org.scalatest.funsuite.AnyFunSuite
import breeze.linalg._
import breeze.numerics._
import breeze.stats.distributions._
import SNF.SNF

class SNFTest extends AnyFunSuite {

  test("SNF calcula la matriz de similitud entre los datos de entrada") {
    // Lee los datos de entrada desde resources, son archivos csv
    val entradas = Array(
      csvread(new java.io.File("src/test/resources/data1.csv")),
      csvread(new java.io.File("src/test/resources/data2.csv")),
      csvread(new java.io.File("src/test/resources/data3.csv"))
    )

    val salidasEsperadas = Array(
      csvread(new java.io.File("src/test/resources/W.csv")),
      csvread(new java.io.File("src/test/resources/W1.csv")),
      csvread(new java.io.File("src/test/resources/W2.csv")),
      csvread(new java.io.File("src/test/resources/W3.csv"))
    )

    val data1: DenseMatrix[Double] = entradas(0)
    val data2: DenseMatrix[Double] = entradas(1)
    val data3: DenseMatrix[Double] = entradas(2)

    val W: DenseMatrix[Double] = salidasEsperadas(0)
    val W1: DenseMatrix[Double] = salidasEsperadas(1)
    val W2: DenseMatrix[Double] = salidasEsperadas(2)
    val W3: DenseMatrix[Double] = salidasEsperadas(3)

    val mu = 0.5
    val percent = 5
    val maxIter = 20
    val modeloSNF = new SNF(exp = data1, methy = data2, mirna = data3)
    val (
      matrizEstatusPromedio,
      data1MatrizEstatus,
      data2MatrizEstatus,
      data3MatrizEstatus
    ) =
      modeloSNF.aplicarSNF(m = mu, porcentaje = percent, iteraciones = maxIter)

    // Define the percentage of error epsilon
    val epsilon = 0.05

    // Calculate the maximum absolute difference between the matrices
    val diffW = max(abs(W - matrizEstatusPromedio))
    val diffW1 = max(abs(W1 - data1MatrizEstatus))
    val diffW2 = max(abs(W2 - data2MatrizEstatus))
    val diffW3 = max(abs(W3 - data3MatrizEstatus))

    // Check that the maximum absolute difference is less than epsilon
    assert(diffW <= epsilon * max(abs(W)))
    assert(diffW1 <= epsilon * max(abs(W1)))
    assert(diffW2 <= epsilon * max(abs(W2)))
    assert(diffW3 <= epsilon * max(abs(W3)))

  }

}

