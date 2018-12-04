package edu.nus

import me.tongfei.progressbar.ProgressBar
import java.io.File
import java.io.OutputStreamWriter

data class Chromosome(val name:String, val length:Int, var payload:Int, var binCount:Int)

class FastaSplitter(val file: File) {
    private fun chopAndWriteSeqeunce(name:String, seq:StringBuilder, binSize: Int, writer: OutputStreamWriter):Chromosome {
        val chrom = Chromosome(name, seq.length, 0, 0)
        var dropped = 0
        var offset = 0
        while (offset < seq.length){
            while (offset<seq.length && seq[offset] in "Nn") offset++
            if (offset>=seq.length) break
            var nextN = seq.indexOf('N',offset,true)
            if (nextN==-1) nextN = seq.length
            while (nextN - offset > binSize*1.5){
                writer.write(">${name}_${offset}_${offset+binSize}\n${seq.subSequence(offset, offset + binSize)}\n")
                chrom.payload += binSize
                chrom.binCount += 1
                offset += binSize
            }
            if (nextN-offset<binSize*0.5)
                dropped += nextN-offset
            else {
                writer.write(">${name}_${offset}_$nextN\n${seq.subSequence(offset, nextN)}\n")
                chrom.payload += nextN-offset
                chrom.binCount += 1
            }
            offset = nextN
        }
        return chrom
    }
    fun chopAndWriteFasta(binSize: Int, outfile: File):ArrayList<Chromosome> {
        val chromNameSet = mutableSetOf<String>()
        val pb = ProgressBar("Breaking ${file.name}", file.length())
        val writer = outfile.writer()
        var name = ""
        var seq = StringBuilder()
        val list = ArrayList<Chromosome>()
        file.bufferedReader().forEachLine {
            pb.stepBy(it.length.toLong() + 1)
            if (it.isNotEmpty())
            if (it[0] != '>') seq.append(it) else {
                if (name != "") list.add(chopAndWriteSeqeunce(name, seq, binSize, writer))
                name = it.substring(1).takeWhile { c -> !c.isWhitespace() }
                if (chromNameSet.contains(name)) throw Exception("Reference chromosome/contig names after space trimming are indistinguishable")
                chromNameSet.add(name)
                seq = StringBuilder()
            }
        }
        list.add(chopAndWriteSeqeunce(name, seq, binSize, writer))
        writer.close()
        pb.close()
        return list
    }
}
