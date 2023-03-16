import breeze.io.CSVReader
import breeze.linalg._
import breeze.numerics.NaN

import java.io.{FileInputStream, InputStreamReader}

object LoadData {

  private val MISSING_VALUE = 0

  def tratarDatos(types: Seq[String]): (DenseMatrix[Double],DenseMatrix[Double],DenseMatrix[Double]) = {

    println("CARGANDO LOS DATOS")

    val cnts = loadRappoport(types)
    val (vista1 ,vista2, vista3 ) = igualadorEstandarizador(cnts(0)(0),cnts(0)(1),cnts(0)(2))

    val vista1Final = outlier_removal(vista1)
    val vista2Final = outlier_removal(vista2)
    val vista3Final = outlier_removal(vista3)

    println("DATOS CARGADOS Y TRATADOS")
    (vista1Final,vista2Final,vista3Final)
  }

  private def igualadorEstandarizador(A: DataMa, B: DataMa, C: DataMa): (DenseMatrix[Double], DenseMatrix[Double],DenseMatrix[Double]) = {
    val commonColNames = A.colNames.intersect(B.colNames).intersect(C.colNames)
    val filteredAData = DenseMatrix.zeros[Double](A.data.rows, commonColNames.length)
    val filteredBData = DenseMatrix.zeros[Double](B.data.rows, commonColNames.length)
    val filteredCData = DenseMatrix.zeros[Double](C.data.rows, commonColNames.length)

    commonColNames.zipWithIndex.map { case (colName, i) =>
      val aCol = A.colNames.indexOf(colName)
      val bCol = B.colNames.indexOf(colName)
      val cCol = C.colNames.indexOf(colName)
      filteredAData(::, i) := A.data(::, aCol)
      filteredBData(::, i) := B.data(::, bCol)
      filteredCData(::, i) := C.data(::, cCol)
    }

    (filteredAData, filteredBData, filteredCData)
  }

  private def outlier_removal(m: DenseMatrix[Double]): DenseMatrix[Double] = {
    var rows_removed: List[Int] = List[Int]()
    var count_missing_values = 0
    var m_res = m.copy
    for (i <- 0 until m.rows) {
      for (j <- 0 until m.cols) {
        if (m(i, j) == MISSING_VALUE) {
          count_missing_values = count_missing_values + 1
        }
      }
      if (count_missing_values > m.cols * 0.2) {
        rows_removed = rows_removed.appended(i)
      }
      count_missing_values = 0
    }
    for (i <- Range(rows_removed.length-1,0,-1)){
      val aux = m_res.delete(rows_removed(i), Axis._0)
      m_res = aux
    }
    m_res
  }

  private def loadRappoport(types: Seq[String]): Seq[Seq[DataMa]] = {

    val omics = Seq("exp", "methy", "mirna")

    val cnts = types.map(c => {
      omics.map(o =>
        csvreader(
          new InputStreamReader(
            new FileInputStream(
              //s"C:/Users/NicolasL/IdeaProjects/SNFDefinito/src/main/$c/$o"
              
              s"C:/Users/Jhonatan/Documents/Mis documentos/Multiview learning/GITHUB_Multi-View-learning/ModelosSNF-NMF/src/main/$c/$o"
              //s"./SNF"
            )
          )
        )
      )
    })
    cnts
  }

  case class DataMa(
                     rowNames: Seq[String],
                     colNames: Seq[String],
                     data: DenseMatrix[Double]
                   )

  private def csvreader(
                 reader: java.io.Reader,
                 separator: Char = ' ',
                 quote: Char = '"',
                 escape: Char = '\\',
                 skipLines: Int = 1,
                 skipCols: Int = 1,
                 rowNameIndex: Int = 0
               ): DataMa = {
    var mat = CSVReader.read(reader, separator, quote, escape)
    mat = mat.takeWhile(line =>
      line.nonEmpty && line.head.nonEmpty
    ) // empty lines at the end
    reader.close()
    val ma = if (mat.isEmpty) {
      DenseMatrix.zeros[Double](0, 0)
    } else {
      DenseMatrix.tabulate(mat.length-2 , mat.head.length ) {
        case (i, j) =>
          val v = mat(i +2)(j + 1)
          if (v.length > 0 && v != "NA")
            v.toDouble
          else
            NaN
      }
    }
    val rowNames = (skipLines until mat.length).map(i => mat(i)(rowNameIndex))
    val colNames = (skipCols until mat.head.length).map(i => mat(0)(i))
    DataMa(rowNames, colNames, ma)
  }

  /*def csvreader(
                 reader: java.io.Reader,
                 separator: Char = ' ',
                 quote: Char = '"',
                 escape: Char = '\\',
                 skipLines: Int = 1,
                 skipCols: Int = 1,
                 rowNameIndex: Int = 0): DataMa = {
    var mat = CSVReader.read(reader, separator, quote, escape)
    mat = mat.takeWhile(line => line.nonEmpty && line.head.nonEmpty) // empty lines at the end
    reader.close()
    val ma = if (mat.isEmpty) {
      DenseMatrix.zeros[Double](0, 0)
    } else {
      DenseMatrix.tabulate(mat.length - skipLines, mat.head.length - skipCols) { case (i, j) =>
        val v = mat(i + skipLines)(j + skipCols)
        if (v.length > 0 && v != "NA")
          v.toDouble
        else
          NaN
      }
    }
    val rowNames = (skipLines until mat.length).map(i => mat(i)(rowNameIndex))
    val colNames = (skipCols until mat.head.length).map(i => mat(0)(i))
    DataMa(rowNames, colNames, ma)
  }*/
}
