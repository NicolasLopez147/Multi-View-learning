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
    val filas = x.rows
    val s = (x - mediaMatrizX).t * (y - mediaMatrizY)
    for (i <- 0 until s.rows) {
      for (j <- 0 until s.cols) {
        s(i, j) = s(i, j) / filas
      }
    }
    s
  }
  private def normalizarDatos(m:DenseMatrix[Double]): DenseMatrix[Double] = {
    val avg = mean(m)
    val dev = stddev(m)
    (m - avg) / dev
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

    val matrizCorrelacion = DenseMatrix((covarianzaU,covarianzaUV),(covarianzaVU,covarianzaV))
    matrizCorrelacion
  }

  private def calcularPesos(x: DenseMatrix[Double], y: DenseMatrix[Double]): (DenseMatrix[Double], DenseMatrix[Double]) = {
    val mediaMatrizX = calcularMatrizMedia(x)
    val mediaMatrizY = calcularMatrizMedia(y)
    val matrizCovarianzaXX = calcularMatrizCovarianza(x, x, mediaMatrizX, mediaMatrizX)
    val matrizCovarianzaYY = calcularMatrizCovarianza(y, y, mediaMatrizY, mediaMatrizY)
    val sigmaS = List(matrizCovarianzaXX, calcularMatrizCovarianza(x, y, mediaMatrizX, mediaMatrizY), matrizCovarianzaYY)
    val matrizSigmaX: DenseMatrix[Double] = inv(sigmaS(0)) * sigmaS(1) * inv(sigmaS(2)) * sigmaS(1).t
    val matrizSigmaY: DenseMatrix[Double] = inv(sigmaS(2)) * sigmaS(1).t * inv(sigmaS(0)) * sigmaS(1)


    val eigenValoresX = calcularEiginvalores(matrizSigmaX, matrizCovarianzaXX)
    val eigenValoresY = calcularEiginvalores(matrizSigmaY, matrizCovarianzaYY)

    (eigenValoresX, eigenValoresY)
  }

  private def calcularEiginvalores(matrizSigma: DenseMatrix[Double], matrizCovarianza: DenseMatrix[Double]): DenseMatrix[Double] = {
    val eigenValores = eig(matrizSigma).eigenvalues.asDenseMatrix
    val restriccion = eigenValores * matrizCovarianza * eigenValores.t
    for (i <- 0 until eigenValores.rows) {
      eigenValores(0, i) = eigenValores(0, i) / pow(restriccion(0, 0), 0.5)
    }
    eigenValores
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
    var indice = 0

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
        indice = i
      }
    }
    println(indice)
    correlationMatrixMax
  }

}
