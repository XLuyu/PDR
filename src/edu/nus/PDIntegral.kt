package edu.nus
import java.lang.Math.*
import org.apache.commons.math3.complex.*

class PDIntegeral(val l1:Double, val l2:Double, val d1:Double, val d2:Double, val I:Int){
    constructor(A:Mapping, B:Mapping): this(
            l1 = A.binCount*PDMEGA.binSize.toDouble(),
            l2 = B.binCount*PDMEGA.binSize.toDouble(),
            d1 = B.pos-A.pos.toDouble(),
            d2 = (if (B.forward()==A.forward()) 1 else -1)*(B.alignmentStart-A.alignmentStart).toDouble(),
            I = if (B.forward()==A.forward()) 1 else -1
            ) { }
    operator fun Double.plus(c:Complex) = Complex(this).add(c)!!
    operator fun Double.minus(c:Complex) = Complex(this).subtract(c)!!
    operator fun Double.times(c:Complex) = Complex(this).multiply(c)!!
    operator fun Complex.plus(c:Complex) = this.add(c)!!
    operator fun Complex.minus(c:Complex) = this.subtract(c)!!
    operator fun Complex.plus(c:Double) = this.add(c)!!
    operator fun Complex.minus(c:Double) = this.subtract(c)!!
    private fun log(x:Double) = Complex(x).log()
    private fun integeralRect(y:Double, Y:Double, x:Double, X:Double, d:Double, D:Double, Id:Double, ID:Double, sign:Double):Double {
      // integrate((T*(a+L*y-x)/(b+K*y-x)), (x,q,p),(y,w,v))
      val (q,p,w,v)=listOf(x,X,y,Y)
      val (a,b,L,K)=listOf(d,D,Id,ID)
      val s = (
        - v*(v*(K - L) - 2*(K*v - L*v - 2*a + 2*b)*log(-K*v - b + p)) 
        + v*(v*(K - L) - 2*(K*v - L*v - 2*a + 2*b)*log(-K*v - b + q)) 
        + w*(w*(K - L) - 2*(K*w - L*w - 2*a + 2*b)*log(-K*w - b + p)) 
        - w*(w*(K - L) - 2*(K*w - L*w - 2*a + 2*b)*log(-K*w - b + q)) 
        + 4*(p*v - p*w - q*v + q*w) 
        - 2*K*v*(-2*K*a + K*b + K*p + L*b - L*p) 
        + 2*K*v*(-2*K*a + K*b + K*q + L*b - L*q) 
        + 2*K*w*(-2*K*a + K*b + K*p + L*b - L*p) 
        - 2*K*w*(-2*K*a + K*b + K*q + L*b - L*q) 
        + 2*(b - p)*(-2*K*a + K*b + K*p + L*b - L*p)*log(2*K*(K*v + b - p)) 
        - 2*(b - p)*(-2*K*a + K*b + K*p + L*b - L*p)*log(2*K*(K*w + b - p)) 
        - 2*(b - q)*(-2*K*a + K*b + K*q + L*b - L*q)*log(2*K*(K*v + b - q)) 
        + 2*(b - q)*(-2*K*a + K*b + K*q + L*b - L*q)*log(2*K*(K*w + b - q))).divide(4*sign)
        if (abs(s.imaginary)>abs(s.real)/1000000+0.1) throw Exception("Integral error")
        return s.real
    }
    private fun integeral11(y:Double, Y:Double, x:Double, X:Double, d:Double, D:Double, sign:Double):Double {
        // integrate((a+y-x)/(b+y-x), (x,q,p+y),(y,w,v))
        val (q,p,w,v)=listOf(x,X,y,Y)
        val (a,b)=listOf(d,D)
        val s= (  v*v/2
                - w*w/2
                + v*(a - b)*log(-b + q - v)
                - w*(a - b)*log(-b + q - w)
                - v*(a*log(-b + p) + a - b*log(-b + p) - b - p + q)
                + w*(a*log(-b + p) + a - b*log(-b + p) - b - p + q)
                + (a - b)*(b - q)*log(b - q + v)
                - (a - b)*(b - q)*log(b - q + w))
        if (abs(s.imaginary)>abs(s.real)/1000000+0.1) throw Exception("Integral error")
        return sign*s.real
    }
    private fun integeral1_1(y:Double, Y:Double, x:Double, X:Double, d:Double, D:Double, sign:Double):Double {
        // integrate((a+y-x)/(b-y-x), (x,q,p+y),(y,w,v))
        val (q,p,w,v)=listOf(x,X,y,Y)
        val (a,b)=listOf(d,D)
        val s= (p*v - p*w - q*v + q*w
                - v*v/2 - v*(a - q) - v*(-2*a + b + p)/2
                + w*w/2 + w*(a - q) + w*(-2*a + b + p)/2
                + v*(-a + b + v)*log(-b + p - 2*v)
                - w*(-a + b + w)*log(-b + p - 2*w)
                - v*(-a + b + v)*log(-b + q - v)
                + w*(-a + b + w)*log(-b + q - w)
                + (a - q)*(b - q)*log(b - q + v)
                - (a - q)*(b - q)*log(b - q + w)
                + (b - p)*(-2*a + b + p)*log(2*b - 2*p + 4*v).divide(4.0)
                - (b - p)*(-2*a + b + p)*log(2*b - 2*p + 4*w).divide(4.0))
        if (abs(s.imaginary)>abs(s.real)/1000000+0.1) throw Exception("Integral error")
        return sign*s.real
    }
    private fun integeral_1_1(y:Double, Y:Double, x:Double, X:Double, d:Double, D:Double, sign:Double):Double {
        // integrate((a-y-x)/(b-y-x), (x,q,p-y),(y,w,v))
        val (q,p,w,v)=listOf(x,X,y,Y)
        val (a,b)=listOf(d,D)
        val s= (- v*v/2
                + w*w/2
                + v*(a - b)*log(-b + q + v)
                - w*(a - b)*log(-b + q + w)
                - v*(a*log(-b + p) + a - b*log(-b + p) - b - p + q)
                + w*(a*log(-b + p) + a - b*log(-b + p) - b - p + q)
                - (a - b)*(b - q)*log(-b + q + v)
                + (a - b)*(b - q)*log(-b + q + w))
        if (abs(s.imaginary)>abs(s.real)/1000000+0.1) throw Exception("Integral error")
        return sign*s.real
    }
    private fun integeral_11(y:Double, Y:Double, x:Double, X:Double, d:Double, D:Double, sign:Double):Double {
        // integrate((a-y-x)/(b+y-x), (x,q,p-y),(y,w,v))
        val (q,p,w,v)=listOf(x,X,y,Y)
        val (a,b)=listOf(d,D)
        val s= (p*v - p*w - q*v + q*w
                - v*v/2 - v*(a - q) - v*(-2*a + b + p)/2
                + w*w/2 + w*(a - q) + w*(-2*a + b + p)/2
                + v*(-a + b + v)*log(-b + p - 2*v)
                - w*(-a + b + w)*log(-b + p - 2*w)
                - v*(-a + b + v)*log(-b + q - v)
                + w*(-a + b + w)*log(-b + q - w)
                + (a - q)*(b - q)*log(b - q + v)
                - (a - q)*(b - q)*log(b - q + w)
                + (b - p)*(-2*a + b + p)*log(2*b - 2*p + 4*v).divide(4.0)
                - (b - p)*(-2*a + b + p)*log(2*b - 2*p + 4*w).divide(4.0))
        if (abs(s.imaginary)>abs(s.real)/1000000+0.1) throw Exception("Integral error")
        return sign*s.real
    }
    private fun integral11doubleY(y:Double, Y:Double, x:Double, X:Double, d:Double, D:Double, sign:Double):Double {
        // integrate((a+y-x)/(b+y-x), (x,q+y,p+y),(y,w,v))
        val (q,p,w,v)=listOf(x,X,y,Y)
        val (a,b)=listOf(d,D)
        val s = (v - w)*(p - q - (a - b)*log(-b + p) + (a - b)*log(-b + q))
        if (abs(s.imaginary)>abs(s.real)/1000000+0.1) throw Exception("Integral error")
        return sign*s.real
    }
    fun getPD():Double{
        var sum = 0.0
        try {
            if (I == 1) {
                val bp1 = -d2
                val bp2 = bp1 + l1
                val bpm1 = -(d1 + d2) / 2
                val bpm2 = bpm1 + l1
                if (d2 > d1) return integeralRect(y = 0.0, Y = l2, x = 0.0, X = l1, d = d1, D = d2, Id = 1.0, ID = 1.0, sign = 1.0) //#1
                if (bp2 < l2) sum += integeralRect(y = max(0.0, bp2), Y = l2, x = 0.0, X = l1, d = min(d1, d2), D = max(d1, d2), Id = 1.0, ID = 1.0, sign = 1.0) //#1
                if (0 < bp2 && bp1 < l2) sum += integeral11(y = max(0.0, bp1), Y = min(bp2, l2), x = 0.0, X = d2, d = d2, D = d1, sign = 1.0) //#2
                if (bpm2 > bp1) {
                    if (0 < bp2 && bpm2 < l2) sum += -integeral11(y = max(0.0, bpm2), Y = min(bp2, l2), x = l1, X = d2, d = d2, D = d1, sign = -1.0) //#3
                    if (0 < bpm2 && bp1 < l2) sum += integral11doubleY(y = max(0.0, bp1), Y = min(bpm2, l2), x = d2, X = (d2 + d1) / 2, d = d2, D = d1, sign = -1.0) // #3.5
                    if (0 < bp1 && bpm1 < l2) sum += integeral11(y = max(0.0, bpm1), Y = min(bp1, l2), x = 0.0, X = (d2 + d1) / 2, d = d2, D = d1, sign = -1.0) //#4
                } else {
                    if (0 < bp2 && bp1 < l2) sum += -integeral11(y = max(0.0, bp1), Y = min(bp2, l2), x = l1, X = d2, d = d2, D = d1, sign = -1.0) //#3
                    if (0 < bp1 && bpm2 < l2) sum += integeralRect(y = max(0.0, bpm2), Y = min(bp1, l2), x = 0.0, X = l1, d = d2, D = d1, Id = 1.0, ID = 1.0, sign = -1.0) // #3.5
                    if (0 < bpm2 && bpm1 < l2) sum += integeral11(y = max(0.0, bpm1), Y = min(bpm2, l2), x = 0.0, X = (d2 + d1) / 2, d = d2, D = d1, sign = -1.0) //#4
                }
                if (0 < bpm2 && bpm1 < l2) sum += -integeral11(y = max(0.0, bpm1), Y = min(bpm2, l2), x = l1, X = (d1 + d2) / 2, d = d1, D = d2, sign = -1.0) // #5
                if (0 < bpm1) sum += integeralRect(y = 0.0, Y = min(bpm1, l2), x = 0.0, X = l1, d = d1, D = d2, Id = 1.0, ID = 1.0, sign = -1.0) //#6
            } else {
                val bp2 = d2
                val bp1 = bp2 - l1
                val bpm1 = (d2 - d1) / 2
                val bpm2 = (d1 + d2) / 2
                val rightBound = min(l1, bpm2)
                if (bp2 < l2 && 0 < rightBound) sum += integeralRect(y = max(0.0, bp2), Y = l2, x = 0.0, X = rightBound, d = d2, D = d1, Id = -1.0, ID = 1.0, sign = -1.0) // #1
                if (0 < bp2 && bp1 < l2) {
                    sum += -integeral_11(y = max(0.0, bp1), Y = min(bp2, l2), x = rightBound, X = d2, d = d2, D = d1, sign = -1.0) // #2
                    sum += integeral_11(y = max(0.0, bp1), Y = min(bp2, l2), x = 0.0, X = d2, d = d2, D = d1, sign = 1.0) // #4
                }
                if (bpm2 < l1) sum += integeralRect(y = 0.0, Y = l2, x = max(bpm2, 0.0), X = l1, d = d1, D = d2, Id = 1.0, ID = -1.0, sign = -1.0) // #3
                if (0 < bp1 && bpm1 < l2) sum += integeralRect(y = max(0.0, bpm1), Y = min(bp1, l2), x = 0.0, X = l1, d = d2, D = d1, Id = -1.0, ID = 1.0, sign = 1.0) // #5
                if (bpm1 > 0) sum += integeralRect(y = 0.0, Y = min(bpm1, l2), x = 0.0, X = l1, d = d1, D = d2, Id = 1.0, ID = -1.0, sign = 1.0) // #6
            }
        } catch (e:Exception){
            println("$l1 $l2 $d1 $d2 $I")
            throw e
        }
        return sum
    }

}
fun main(args: Array<String>) {
    val a = PDIntegeral(1000.0,7548000.0,7.7471E7,-7547499.0,1)
    println(a.getPD())
}