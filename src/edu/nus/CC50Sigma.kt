package edu.nus

import java.lang.Math.*
import kotlin.system.exitProcess

class CC50Sigma(private val A:Mapping, private val B:Mapping, private val threshold:Int){
    fun getCC50(): Int {
        if (A.forward()==B.forward()) {
            val d1 = B.start-A.start
            val d2 = (B.alignmentStart-A.alignmentStart)*(if (A.forward()) 1 else -1)
            val cut = d1-threshold
            if (d2<=cut) return 0
            if (0<=cut && cut<A.getLength()) return min(A.getLength()-cut, B.getLength())
            if (-B.getLength()<cut && cut<0) return min(A.getLength(), B.getLength()+cut)
            return 0
        } else {
            val d1 = B.start-A.start
            val d2 = (B.alignmentStart-A.alignmentStart)*(if (A.forward()) 1 else -1)
            val cut = d1-threshold
            val x = (d2+cut+1)/2 // open right bound, first no-solution x
            val y = (d2-cut+1)/2 // open right bound, first no-solution ys
            if (0<=cut && cut<A.getLength() && y>0) return min(y,min(A.getLength()-cut, B.getLength()))
            if (-B.getLength()<cut && cut<0 && x>0) return min(x,min(A.getLength(), B.getLength()+cut))
            return 0
        }
    }
    fun bruteforce(): Int{
        var ans = 0
        val dire1 = if (A.forward()) 1 else -1
        val dire2 = if (B.forward()) 1 else -1
        for (i in 0 until A.getLength()) {
            val j = A.start + i + threshold - B.start
            if (j < 0 || B.getLength() <= j) continue
            if (((B.alignmentStart + j * dire2)-(A.alignmentStart + i * dire1))*dire1>0) ans += 1
        }
        return ans
    }
//            println("A: [${A.start},${A.end}] -> [${A.alignmentStart},${A.alignmentEnd}]")
//            println("B: [${B.start},${B.end}] -> [${B.alignmentStart},${B.alignmentEnd}]")
//            println("T: $threshold")
//            println("Formula:    $a")
//            println("BruteForce: $b")
//            exitProcess(1)
}