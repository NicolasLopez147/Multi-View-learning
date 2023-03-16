package SNF

import breeze.linalg.{DenseMatrix, InjectNumericOps}

class SNF (val exp:DenseMatrix[Double],val methy:DenseMatrix[Double],val mirna:DenseMatrix[Double]){

  /*private def normalizar(m: DenseMatrix[Double]): DenseMatrix[Double] = {
    val avg = mean(m)
    val dev = stddev(m)
    (m - avg) / dev
  }*/
  private def normalizar(matriz: DenseMatrix[Double]): DenseMatrix[Double] = {

    for (i <- 0 until matriz.rows) {
      var media: Double = 0
      var contador = 0
      for (j <- 0 until matriz.cols) {
        media += matriz(i, j)
        contador += 1
      }
      media = media / contador
      var varianza: Double = 0
      contador = 0
      for (j <- 0 until matriz.cols) {
        varianza += math.pow(matriz(i, j) - media, 2)
        contador += 1
      }
      varianza = math.pow(varianza / contador, 0.5)
      if (varianza == 0){
        varianza = 0.00000000000000001
      }
      for (j <- 0 until matriz.cols) {
        matriz(i, j) = (matriz(i, j) - media) / varianza
      }
    }
    matriz
  }

  private def calcularDistancia(matriz: DenseMatrix[Double]): DenseMatrix[Double] = {
    val p = DenseMatrix.zeros[Double](matriz.rows, matriz.rows)
    for (i <- 0 until matriz.rows) {
      for (j <- 0 until matriz.rows) {
        var aux: Double = 0
        for (k <- 0 until matriz.cols) {
          aux += math.pow(matriz(i, k) - matriz(j, k), 2)
        }
        p(i, j) = math.pow(aux, 0.5)
      }
    }
    p
  }

  private def calcularVecinos(filas: Integer, porcentaje: Double): Integer = {
    var vecinos = (filas * porcentaje / 100).toInt
    if (vecinos == 0) {
      vecinos = 1
    }
    vecinos
  }

  private def calcularPromedioDistancia(p: DenseMatrix[Double], porcentaje: Double): DenseMatrix[Double] = {
    val pMedia = DenseMatrix.zeros[Double](p.rows, 1)
    val vecinos: Integer = calcularVecinos(p.rows, porcentaje)


    for (i <- 0 until  p.rows ) {
      val fila: Array[Double] = new Array[Double](p.cols)
      for (j <- 0 until  p.cols ) {
        fila(j) = p(i, j)
      }
      val ordFila = fila.sorted
      var media: Double = 0
      for (k <- 1 to vecinos) {
        media += ordFila(k)
      }
      pMedia(i, 0) = media / vecinos
    }
    pMedia
  }

  private def calcularE(i: Integer, j: Integer, pMedia: DenseMatrix[Double], p: DenseMatrix[Double]): Double = {
    val e: Double = (pMedia(i.toInt, 0) + pMedia(j.toInt, 0) + p(i.toInt, j.toInt)) / 3
    e
  }

  private def calcularKernel(p: DenseMatrix[Double], m: Double, porcentaje: Double): DenseMatrix[Double] = {
    val kernel = DenseMatrix.zeros[Double](p.rows, p.rows)
    val pMedia: DenseMatrix[Double] = calcularPromedioDistancia(p, porcentaje)

    for (i <- 0 until p.rows) {
      for (j <- 0 until p.cols) {
        kernel(i, j) = math.exp(-math.pow(p(i, j), 2) / (calcularE(i, j, pMedia, p) * m))
      }
    }
    kernel
  }

  private def calcularSumatoria(w: DenseMatrix[Double]): DenseMatrix[Double] = {
    val matrizSumatoria = DenseMatrix.zeros[Double](w.rows, 1)

    for (i <- 0 until w.rows) {
      var sumatoria: Double = 0
      for (j <- 0 until w.cols) {
        if (j != i) sumatoria += w(i, j)
      }
      matrizSumatoria(i, 0) = sumatoria
    }
    matrizSumatoria
  }

  private def calcularKernelCompleto(matrizKernel: DenseMatrix[Double]): DenseMatrix[Double] = {
    val kernelCompleto = DenseMatrix.zeros[Double](matrizKernel.rows, matrizKernel.cols)

    val matrizSumatoria = calcularSumatoria(matrizKernel)
    for (i <- 0 until kernelCompleto.rows) {
      for (j <- 0 until kernelCompleto.cols) {
        if (i != j) {
          kernelCompleto(i, j) = matrizKernel(i, j) /:/ (2 * matrizSumatoria(i, 0))
        } else {
          kernelCompleto(i, j) = 0.5
        }
      }
    }

    kernelCompleto
  }

  private def calcularKernelDispersa(matrizKernel: DenseMatrix[Double], porcentaje: Double): DenseMatrix[Double] = {
    val kernelDisperso = DenseMatrix.zeros[Double](matrizKernel.rows, matrizKernel.cols)
    val vecinos = calcularVecinos(matrizKernel.rows, porcentaje)

    for (i <- 0 until matrizKernel.rows) {
      val fila: Array[Double] = new Array[Double](matrizKernel.cols)
      for (j <- 0 until matrizKernel.cols) {
        fila(j) = matrizKernel(i, j)
      }
      val ordFila = fila.sorted.dropRight(fila.length - vecinos-1)
      var sumatoria: Double = 0
      for (i <- 1 until ordFila.length) {
        sumatoria += ordFila(i)
      }
      if (sumatoria == 0){
        sumatoria = 0.0000000000000000000001
      }
      for (j <- 0 until matrizKernel.cols) {
        if (ordFila.contains(matrizKernel(i, j))) {
          kernelDisperso(i, j) = matrizKernel(i, j) / sumatoria
        }
      }
    }

    kernelDisperso
  }

  private def calcularMatrizEstatus(kernelCompletoVista1: DenseMatrix[Double], kernelCompletoVista2: DenseMatrix[Double], kernelDisperso: DenseMatrix[Double]): DenseMatrix[Double] = {
    val sumatoriaMatrices = (kernelCompletoVista1 + kernelCompletoVista2) * 0.5

    val matrizEstatus = kernelDisperso * sumatoriaMatrices * kernelDisperso.t

    matrizEstatus
  }
  def aplicarSNF(m:Double,porcentaje:Double,iteraciones:Integer): (DenseMatrix[Double],DenseMatrix[Double],DenseMatrix[Double],DenseMatrix[Double],DenseMatrix[Double],DenseMatrix[Double],DenseMatrix[Double]) = {

    val expKernel = calcularKernel(calcularDistancia(normalizar(exp).t),m,porcentaje)
    val methyKernel = calcularKernel(calcularDistancia(normalizar(methy).t),m,porcentaje)
    val mirnaKernel = calcularKernel(calcularDistancia(normalizar(mirna).t),m,porcentaje)

    // Calcular matriz kernel completa
    var expKernelCompleto = calcularKernelCompleto(expKernel)
    var methyKernelCompleto = calcularKernelCompleto(methyKernel)
    var mirnaKernelCompleto = calcularKernelCompleto(mirnaKernel)

    // Calcular matriz kernel dispersa
    val expKernelDispersa = calcularKernelDispersa(expKernel, porcentaje)
    val methyKernelDispersa = calcularKernelDispersa(methyKernel, porcentaje)
    val mirnaKernelDispersa = calcularKernelDispersa(mirnaKernel, porcentaje)

    // Calcular matriz estatus
    var expMatrizEstatus = DenseMatrix.zeros[Double](expKernelCompleto.rows, expKernelCompleto.cols)
    var methyMatrizEstatus = DenseMatrix.zeros[Double](methyKernelCompleto.rows, methyKernelCompleto.cols)
    var mirnaMatrizEstatus = DenseMatrix.zeros[Double](mirnaKernelCompleto.rows, mirnaKernelCompleto.cols)

    for (_ <- 1 to iteraciones) {
      expMatrizEstatus = calcularMatrizEstatus(methyKernelCompleto, mirnaKernelCompleto, expKernelDispersa)
      methyMatrizEstatus = calcularMatrizEstatus(expKernelCompleto, mirnaKernelCompleto, methyKernelDispersa)
      mirnaMatrizEstatus = calcularMatrizEstatus(methyKernelCompleto, expKernelCompleto, mirnaKernelDispersa)

      expKernelCompleto = expMatrizEstatus
      methyKernelCompleto = methyMatrizEstatus
      mirnaKernelCompleto = mirnaMatrizEstatus
    }
    val matrizEstatusPromedio = (expMatrizEstatus+methyMatrizEstatus+mirnaMatrizEstatus)*0.3333
    (matrizEstatusPromedio,expMatrizEstatus,methyMatrizEstatus,mirnaMatrizEstatus)
  }
}
