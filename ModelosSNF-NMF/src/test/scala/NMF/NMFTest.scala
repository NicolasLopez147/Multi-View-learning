
import org.scalatest.funsuite.AnyFunSuite
import breeze.linalg._
import breeze.numerics._
import breeze.stats.distributions._

//configura la semilla aleatoria para que los resultados sean reproducibles


class NMFTest extends AnyFunSuite {
  
  //rand es una instancia de RandBasis que se usa para generar números aleatorios dado un valor de semilla 
  implicit val rand: RandBasis = RandBasis.withSeed(0)

  //X es un arreglo de dos matrices densas de 4x3 y 8x3 respectivamente , que se usará para probar el método nmf
  val X: Array[DenseMatrix[Double]] = Array(
    DenseMatrix.rand[Double](4, 3),
    DenseMatrix.rand[Double](8, 3),
    DenseMatrix.rand[Double](7, 3)
  )

  val epsilon: Double = 0.001

  test("NMF devuelve matrices W y las Hs que reconstruyen X con un error epsilon") {
    val modelo = Trainer.jnmf(X, eps = epsilon)
    val errorDeReconstruccion = (X, modelo.inverse()).zipped.map {
      case (x, x_hat) => sum(pow(x - x_hat, 2)) / sum(pow(x, 2))
    }.sum / X.length * 100
    assert(errorDeReconstruccion <= epsilon)
  }
}

/* Primero, creamos una matriz aleatoria X con dos elementos: 
    una matriz de 4x3 y otra de 8x3. Luego, definimos el límite 
    de error epsilon en 0.001. Llamamos al método jnmf con estos
     valores y obtenemos el objeto modelo que contiene las matrices W y H.*/

/* A continuación, calculamos el error de reconstrucción entre X 
y el producto de W y cada una de las matrices H. Iteramos a través del arreglo X
 y calculamos el error cuadrático medio entre X y X_hat, donde X_hat 
 es el producto de W y la matriz H para ese elemento de X. Finalmente,
  calculamos la media de estos valores y lo multiplicamos
   por 100 para obtener el porcentaje de error.*/

/* Finalmente, comprobamos que este porcentaje de error
 sea menor o igual al valor de epsilon. Si se cumple
  esta condición, nuestro test pasa.*/
