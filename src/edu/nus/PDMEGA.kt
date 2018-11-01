
package edu.nus

import com.github.ajalt.clikt.core.*
import com.github.ajalt.clikt.parameters.arguments.argument
import com.github.ajalt.clikt.parameters.options.*
import com.github.ajalt.clikt.parameters.types.*
import me.tongfei.progressbar.ProgressBar
import java.io.File

class PDMEGA : CliktCommand() {
    companion object {
        var binSize = 0
    }
    val threads by option(help = "Threads to use, default by the CPU core number").int().default(Runtime.getRuntime().availableProcessors())
    val K by option("-k", help = "bin size").int().default(1000)
//    val output by option("-o", help = "Output file path").file().default(File("DesirableOut.fasta"))
    val tmpDir: File by option("-d", help = "Temporary folder for intermediate results").file(exists = false).default(File("PDMEGATmp"))
    val reference by argument(help = "Reference genome").file(exists = true)
    val assembly by argument(help = "Assembly to evaluate").file(exists = true)
    val aligner by option("-a", help = "executable path of aligner (BWA or minimap2), default: bwa").default("bwa")

    private fun bwa(refBlock: File, totalBin: Int):File {
        val sam = File(tmpDir,"RefToAsm.sam")
        if (File(assembly.path+".sa").exists()){
            println("[Info] index file for ${assembly.path} exists, use it")
        } else {
            println("[CommandLine] $aligner index ${assembly.path}")
            Runtime.getRuntime().exec("$aligner index ${assembly.path}").errorStream.bufferedReader().forEachLine { println(it) }
        }
        println("[CommandLine] $aligner mem -t $threads ${assembly.path} ${refBlock.path} -o ${sam.path}")
        val sm = Runtime.getRuntime().exec("$aligner mem -t $threads ${assembly.path} ${refBlock.path} -o ${sam.path}")
        val pb = ProgressBar("BWA", totalBin.toLong())
        sm.errorStream.bufferedReader().forEachLine {
            val match = "Processed (\\d+) reads".toRegex().find(it)
            if (match!=null) pb.stepBy(match.groupValues[1].toLong())
        }
        sm.waitFor()
        pb.close()
        return sam
    }
    private fun minimap2(refBlock: File, totalBin: Int):File {
        val sam = File(tmpDir,"RefToAsm.sam")
        val cmd = "$aligner -t $threads -a -x sr ${assembly.path} ${refBlock.path}"
        println("[CommandLine] $cmd")
        val sm = ProcessBuilder(cmd.split(" "))
        sm.redirectOutput(sam)
        val proc = sm.start()
        val pb = ProgressBar("Minimap2", totalBin.toLong())
        proc.errorStream.bufferedReader().forEachLine {
            val match = "mapped (\\d+) sequences".toRegex().find(it)
            if (match!=null) pb.stepBy(match.groupValues[1].toLong())
        }
        proc.waitFor()
        pb.close()
        return sam
    }
    override fun run() {
        val runtime = kotlin.system.measureTimeMillis {
            tmpDir.mkdirs()
            PDMEGA.binSize = K
            println("[Info] ====== Reference Breaking ======")
            val refBlock = File(tmpDir,"ReferenceBlock.fasta")
            val refInfo = FastaSplitter(reference).chopAndWriteFasta(K,refBlock)
            val totalBin = refInfo.sumBy { it.binCount }
            println("[Info] Successfully finished")
            println("[Info] ====== Alignment ======")
            val sam = if (aligner.endsWith("bwa")) bwa(refBlock, totalBin) else minimap2(refBlock, totalBin)
            println("[Info] successfully finished")
            println("[Info] ====== Mapping Analysis ======")
            val mappingAnalyzer = MappingAnalyzer(sam, refInfo)
            mappingAnalyzer.run()
//          tmpDir.deleteRecursively()
        }
        println("[Success] elapsed time: ${runtime/60000}m${runtime/1000%60}s")
    }
}

fun main(args: Array<String>) = PDMEGA().main(if (args.isEmpty()) arrayOf("--help") else args)

// Processed 240000 reads
// mapped 50000 sequences