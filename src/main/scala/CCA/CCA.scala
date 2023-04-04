package CCA

import breeze.linalg.{DenseMatrix, eig, inv, sum}
import breeze.numerics.{abs, pow}
import breeze.stats.{mean, stddev}

class CCA (matriz:DenseMatrix[Double]){
  private def calcularMatrizMedia(matriz: DenseMatrix[Double]): DenseMatrix[Double] = {
    val matrizMedia = matriz.copy
    for (j <- 0 until matriz.cols) {
      val media = mean(matriz(::, j))
      for (i <- 0 until matriz.rows) {
        matrizMedia(i, j) = media
      }
    }
    matrizMedia
  }

  private def calcularMatrizCovarianza(x: DenseMatrix[Double], y: DenseMatrix[Double], mediaMatrizX: DenseMatrix[Double], mediaMatrizY: DenseMatrix[Double]): DenseMatrix[Double] = {
    ((x - mediaMatrizX).t * (y - mediaMatrizY)).map(l => l/x.rows)
  }
  private def normalizarDatos(m:DenseMatrix[Double]): DenseMatrix[Double] = {
    (m - mean(m)) / stddev(m)
  }

  private def calcularMatrizCorrelacion(u:DenseMatrix[Double],v:DenseMatrix[Double]): DenseMatrix[Double] = {
    val uNorm = normalizarDatos(u)
    val vNorm = normalizarDatos(v)
    val mediaMatrizU = calcularMatrizMedia(uNorm)
    val mediaMatrizV = calcularMatrizMedia(vNorm)
    val covarianzaU = abs(pow(sum(calcularMatrizCovarianza(uNorm,uNorm,mediaMatrizU,mediaMatrizU)),0.5))
    val covarianzaV = abs(pow(sum(calcularMatrizCovarianza(vNorm,vNorm,mediaMatrizV,mediaMatrizV)),0.5))
    val covarianzaUV = abs(sum(calcularMatrizCovarianza(uNorm,vNorm,mediaMatrizU,mediaMatrizV)))
    val covarianzaVU = abs(sum(calcularMatrizCovarianza(vNorm,uNorm,mediaMatrizV,mediaMatrizU)))

    DenseMatrix((covarianzaU,covarianzaUV),(covarianzaVU,covarianzaV))

  }

  private def calcularPesos(x: DenseMatrix[Double], y: DenseMatrix[Double]): (DenseMatrix[Double], DenseMatrix[Double]) = {
    val mediaMatrizX = calcularMatrizMedia(x)
    val mediaMatrizY = calcularMatrizMedia(y)

    val sigmaS = List(calcularMatrizCovarianza(x, x, mediaMatrizX, mediaMatrizX),
      calcularMatrizCovarianza(x, y, mediaMatrizX, mediaMatrizY),
      calcularMatrizCovarianza(y, y, mediaMatrizY, mediaMatrizY))

    val matrizSigmaX: DenseMatrix[Double] = inv(sigmaS(0)) * sigmaS(1) * inv(sigmaS(2)) * sigmaS(1).t
    val matrizSigmaY: DenseMatrix[Double] = inv(sigmaS(2)) * sigmaS(1).t * inv(sigmaS(0)) * sigmaS(1)

    (calcularEiginvalores(matrizSigmaX, sigmaS(0)), calcularEiginvalores(matrizSigmaY, sigmaS(2)))
  }

  private def calcularEiginvalores(matrizSigma: DenseMatrix[Double], matrizCovarianza: DenseMatrix[Double]): DenseMatrix[Double] = {
    val eigenValores = eig(matrizSigma).eigenvalues.asDenseMatrix
    val restriccion = eigenValores * matrizCovarianza * eigenValores.t
    eigenValores.map(x => x / pow(restriccion(0, 0), 0.5))
  }

  private def calcularVariablesCanonicas(datos: DenseMatrix[Double], pesos: DenseMatrix[Double]): DenseMatrix[Double] = {
    if (datos.cols != pesos.rows) {
      datos * pesos.t
    } else
      datos * pesos
  }

  def aplicarAnalisisCorrelacionCanonica(): DenseMatrix[Double] = {
    val matrizT = matriz.t
    var max:Double = 0
    var correlationMatrixMax = matriz

    for (i <- 5 until matrizT.cols){
      val x = matrizT(0 until matrizT.rows - 1, 0 until i)
      val y = matrizT(0 until matrizT.rows - 1, i until matrizT.cols)


      val (weightX, weightY) = calcularPesos(x, y)

      val u = calcularVariablesCanonicas(x, weightX)
      val v = calcularVariablesCanonicas(y, weightY)

      val correlationMatrix = calcularMatrizCorrelacion(u, v)
      if (correlationMatrix(0,1) > max){
        max = correlationMatrix(0,1)
        correlationMatrixMax = correlationMatrix
      }
    }
    correlationMatrixMax
  }

}
