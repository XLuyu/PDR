
package edu.nus

import com.github.ajalt.clikt.core.*
import com.github.ajalt.clikt.parameters.arguments.argument
import com.github.ajalt.clikt.parameters.options.*
import com.github.ajalt.clikt.parameters.types.*
import java.io.File

class PDMEGA : CliktCommand(help =
"""FILES is a list of sequencing files (fasta or fastq, may be with .gz suffix) given in following format:
    + t_1.fq t_2.fq - bg1_1.fq.gz bg1_2.fq.gz - bg_s2_1.fa bg_s2_2.fa
    where target sample is given with leading '+' and each background sample is given with leading '-'
""") {
    companion object {
        var binSize = 0
    }
    val threads by option(help = "Threads to use, default by the CPU core number").int().default(Runtime.getRuntime().availableProcessors())
    val K by option("-k", help = "bin size").int().default(1000)
    val output by option("-o", help = "Output file path").file().default(File("DesirableOut.fasta"))
    val tmpDir: File by option("-d", help = "Temporary folder for intermediate results").file(exists = false).default(File("PDMEGATmp"))
    val reference by argument(help = "Reference genome").file(exists = true)
    val assembly by argument(help = "Assembly to evaluate").file(exists = true)
    val bwa by option(help = "BWA executable path if not in environment variable PATH").default("bwa")

    override fun run() {
        tmpDir.mkdirs()
        PDMEGA.binSize = K
        println("[Info] ====== ChopReference ======")
        val refSeg = File(tmpDir,"ReferenceFragment.fasta")
        val refInfo = FastaSplitter(reference).chopAndWriteFasta(K,refSeg)
        println("[Info] finished successfully")
        println("[Info] ====== BWA ======")
        val bwaSam = File(tmpDir,"bwa.sam")
        val runtime = kotlin.system.measureTimeMillis {
            if (File(assembly.path+".sa").exists()){
                println("[Info] index file for ${assembly.path} exists, use it")
            } else {
                println("[CommandLine] $bwa index ${assembly.path}")
                Runtime.getRuntime().exec("$bwa index ${assembly.path}").errorStream.bufferedReader().forEachLine { println(it) }
            }
            println("[CommandLine] $bwa mem -t $threads ${assembly.path} ${refSeg.path} -o ${bwaSam.path}")
            val sm = Runtime.getRuntime().exec("$bwa mem -t $threads ${assembly.path} ${refSeg.path} -o ${bwaSam.path}")
            sm.errorStream.bufferedReader().forEachLine { println(it) }
            sm.inputStream.bufferedReader().forEachLine { println(it) }
            sm.waitFor()
        }
        println("[Info] successfully finished in ${runtime/1000} seconds")
        println("[Info] ====== Mapping Analysis ======")
        val mappingAnalyzer = MappingAnalyzer(bwaSam, refInfo)
        mappingAnalyzer.run()
//        tmpDir.deleteRecursively()
    }
}

fun main(args: Array<String>) = PDMEGA().main(if (args.isEmpty()) arrayOf("--help") else args)