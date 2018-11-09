package edu.nus

import java.io.File
import htsjdk.samtools.*
import java.lang.Math.*
import me.tongfei.progressbar.ProgressBar
import kotlinx.coroutines.experimental.*

data class Mapping(val chrom: Int, val start: Int, var end: Int, val contig: String, val alignmentStart: Int, var alignmentEnd: Int){
    fun forward() = alignmentEnd-alignmentStart>=0
    fun getLength() = end-start
    fun pairwiseDistance(other:Mapping):Double {
        if (chrom==other.chrom && contig!=other.contig) return 0.0
        if (chrom!=other.chrom && contig==other.contig) return 0.0
        if (chrom!=other.chrom && contig!=other.contig) return this.getLength().toDouble()*other.getLength()
        return PDIntegeral(this,other).getPD()
    }
}
class MappingSet(): ArrayList<Mapping>() {
    fun chrom() = this[0].chrom
    fun start() = this[0].start
    fun end() = this[0].end
    fun getLength() = this[0].getLength()
    fun sorted(): MappingSet {
        java.util.Collections.sort(this, compareBy({it.contig},{it.alignmentStart}))
        return this
    }
    fun allDownstreams(other:MappingSet): Boolean {
        if (this.size!=other.size) return false
        for ( i in 0 until this.size){
            if (this[i].end!=other[i].start) return false
            if (this[i].forward()!=other[i].forward() || abs(this[i].alignmentEnd-other[i].alignmentStart)>PDMEGA.joinError) return false
        }
        return true
    }
    fun extendDownstreamsBy(other:MappingSet): Boolean {
        if (!allDownstreams(other)) return false
        for ( i in 0 until this.size){
            this[i].alignmentEnd = other[i].alignmentEnd
            this[i].end = other[i].end
        }
        return true
    }
    fun pairwiseDistance(other:MappingSet):Double {
        var best = -1.0
        for (i in this){
            for (j in other)
                best = max(best,i.pairwiseDistance(j))
        }
        return best
    }
}
class BamFileReader(File: File, refInfo: ArrayList<ChromInfo>,
                    private val filter: (SAMRecord) -> Boolean = { record:SAMRecord ->
                        !record.readUnmappedFlag  && record.cigarString.all { it in "0123456789MDI=" } //&& record.mappingQuality>10
                    }){
    private val bamFile = SamReaderFactory.makeDefault().open(File)
    private val iterator = bamFile.iterator()
    private val chromToIdx = refInfo.mapIndexed { index: Int, chromInfo: ChromInfo -> Pair(chromInfo.name,index) }.toMap()
    private var chromIdx = -1
    private var pos = -1
    private var next:Mapping? = getNext()
    private fun getNext(): Mapping? {
        while (true){
            if (!iterator.hasNext()) return null
            val record = iterator.next()
            val tokens = record.readName.split("_")
            val (start,end) = tokens.takeLast(2).map { it.toInt() }
            val idx = chromToIdx[tokens.dropLast(2).joinToString("_")]!!
            if (idx>this.chromIdx || idx==this.chromIdx && start>=this.pos) else throw Exception("Generated Sam file isn't ordered by queryname!")
            this.chromIdx = idx
            this.pos = start
            if (filter(record))
                return Mapping(idx, start, end, record.contig,
                        if (record.readNegativeStrandFlag) record.alignmentEnd else record.alignmentStart,
                        if (record.readNegativeStrandFlag) record.alignmentStart else record.alignmentEnd)
        }
    }
    fun hasNextReadMappings() = next!=null
    fun getNextReadMappingSet(): MappingSet {
        val ms = MappingSet()
        do {
            ms.add(next!!)
            next = getNext()
        } while (next!=null && next!!.chrom==ms.chrom() && next!!.start==ms.start())
        return ms.sorted()
    }
    fun getAllReadMappings(): ArrayList<MappingSet> {
        val all = ArrayList<MappingSet>()
        while (hasNextReadMappings()){
            val ms = getNextReadMappingSet()
            if (all.isNotEmpty() && all.last().extendDownstreamsBy(ms)) else all.add(ms)
        }
        return all
    }
    fun getReadMappingByChrom(): ArrayList<MappingSet> {
        val chrom = next!!.chrom
        val all = ArrayList<MappingSet>()
        while (hasNextReadMappings() && next!!.chrom==chrom){
            val ms = getNextReadMappingSet()
            if (all.isNotEmpty() && all.last().extendDownstreamsBy(ms)) else all.add(ms)
        }
        return all
    }
    fun getAllByChrom(): ArrayList<ArrayList<MappingSet>> {
        val all = ArrayList<ArrayList<MappingSet>>()
        while (hasNextReadMappings())
            all.add(getReadMappingByChrom())
        return all
    }
}
class MappingAnalyzer(bwaSam: File, val refInfo: ArrayList<ChromInfo>) {
    private val binSize = PDMEGA.K
    private val bam = BamFileReader(bwaSam, refInfo)
    fun run() {
        val records = bam.getAllByChrom()
//        println("[Info] load ${records.size} appropriate mappings from SAM file")
        for ( chrom in records){
            println("[Info] chromosome ${refInfo[chrom[0].chrom()]} has ${chrom.size} mapping segments " +
                    "(%.2f%% covered)".format(chrom.sumBy { it.getLength() }.toDouble() * 100 / refInfo[chrom[0].chrom()].length))
        }
        val score = analyzeMappingByChrom(records)
        val totalPos = refInfo.sumByDouble { it.length.toDouble() }
        println("[Success] Absolute Score: $score")
        println("[Success] Ratio Score: ${score/totalPos/totalPos}")
    }

    private fun analyzeMappingByChrom(records: ArrayList<ArrayList<MappingSet>>)= runBlocking<Double> {
        val frecords = records.flatten().toTypedArray()
        var score = 0.0
        val pb = ProgressBar("Pairwise Analysis ", frecords.size.toLong())
        val deferred = frecords.indices.map { i->
            async(CommonPool) {
                var s = (frecords[i].getLength()).toDouble()*(frecords[i].getLength())
                for (j in i + 1 until frecords.size) {
                    try {
                        s += 2 * frecords[i].pairwiseDistance(frecords[j])
                    } catch (e: Exception) {
                        println("$i $j")
                        println(frecords[i])
                        println(frecords[j])
                        throw e
                    }
                }
                s
            }
        }
        deferred.forEach {
            pb.step()
            score += it.await()
        }
        pb.close()
        score
    }
}