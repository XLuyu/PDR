package edu.nus

import me.tongfei.progressbar.ProgressBar
import java.io.File
import java.io.OutputStreamWriter

data class ChromInfo(val chrom:String, val length:Int, val binSize:Int, val offset:Int, val binCount:Int)

class FastaSplitter(val file: File) {
    private fun chopAndWriteSeqeunce(name:String, seq:StringBuilder, binSize: Int, writer: OutputStreamWriter):ChromInfo {
        var offset = seq.length % binSize / 2
        while (offset + binSize <= seq.length){
            writer.write(">${name}_$offset\n${seq.subSequence(offset, offset + binSize)}\n")
            offset += binSize
        }
        return ChromInfo(name, seq.length, binSize, seq.length % binSize / 2,seq.length / binSize)
    }
    fun chopAndWriteFasta(binSize: Int, outfile: File):ArrayList<ChromInfo> {
        val chromNameSet = mutableSetOf<String>()
        val pb = ProgressBar("Breaking ${file.name}", file.length())
        val writer = outfile.writer()
        var name = ""
        var seq = StringBuilder()
        val list = ArrayList<ChromInfo>()
        file.bufferedReader().forEachLine {
            pb.stepBy(it.length.toLong() + 1)
            if (it[0] != '>') seq.append(it) else {
                if (name != "") list.add(chopAndWriteSeqeunce(name, seq, binSize, writer))
                name = it.substring(1).takeWhile { !it.isWhitespace() }
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
