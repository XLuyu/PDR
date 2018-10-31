package edu.nus

import java.io.File
import htsjdk.samtools.*
import java.lang.Math.*
import me.tongfei.progressbar.ProgressBar

data class Mapping(val chrom: Int, val pos: Int, val contig: String, var alignmentStart: Int, var alignmentEnd: Int, var binCount:Int = 1){
    fun forward() = alignmentEnd-alignmentStart>=0
    fun pairwiseDistance(other:Mapping):Double {
        if (chrom==other.chrom && contig!=other.contig) return 0.0
        if (chrom!=other.chrom && contig==other.contig) return 0.0
        if (chrom!=other.chrom && contig!=other.contig) return binCount.toDouble()*PDMEGA.binSize*other.binCount*PDMEGA.binSize
        return PDIntegeral(this,other).getPD()
    }
}
class MappingSet(val chrom: Int, val pos: Int): ArrayList<Mapping>() {
    fun getBinCount() = this[0].binCount
    fun sorted(): MappingSet {
        java.util.Collections.sort(this, compareBy({it.chrom},{it.pos}))
        return this
    }
    fun allDownstreams(other:MappingSet): Boolean {
        if (this.size!=other.size) return false
        for ( i in 0 until this.size){
            if (this[i].forward()!=other[i].forward() || abs(this[i].alignmentEnd-other[i].alignmentStart)>50) return false
        }
        return true
    }
    fun extendDownstreamsBy(other:MappingSet): Boolean {
        if (!allDownstreams(other)) return false
        for ( i in 0 until this.size){
            this[i].alignmentEnd = other[i].alignmentEnd
            this[i].binCount += 1
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
    private val chromToIdx = refInfo.mapIndexed { index: Int, chromInfo: ChromInfo -> Pair(chromInfo.chrom,index) }.toMap()
    private var chromIdx = -1
    private var pos = -1
    private var next:Mapping? = getNext()
    private fun getNext(): Mapping? {
        while (true){
            if (!iterator.hasNext()) return null
            val record = iterator.next()
            val splitAt = record.readName.lastIndexOf('_')
            if (splitAt==-1) println(record.toString())
            val chrom = record.readName.substring(0,splitAt)
            val posString = record.readName.substring(splitAt+1)
            val pos = posString.toInt()
            val idx = chromToIdx[chrom]!!
            if (idx>this.chromIdx || idx==this.chromIdx && pos>=this.pos) else throw Exception("Generated Sam file isn't ordered by queryname!")
            this.chromIdx = idx
            this.pos = pos
            if (filter(record))
                return Mapping(idx, pos, record.contig,
                        if (record.readNegativeStrandFlag) record.alignmentEnd else record.alignmentStart,
                        if (record.readNegativeStrandFlag) record.alignmentStart else record.alignmentEnd)
        }
    }
    fun hasNextReadMappings() = next!=null
    fun getNextReadMappingSet(): MappingSet {
        val ms = MappingSet(next!!.chrom, next!!.pos)
        do {
            ms.add(next!!)
            next = getNext()
        } while (next!=null && next!!.chrom==ms.chrom && next!!.pos==ms.pos)
        return ms.sorted()
    }
    fun getAllReadMappings(): ArrayList<MappingSet> {
        val all = ArrayList<MappingSet>()
        while (hasNextReadMappings())
            all.add(getNextReadMappingSet())
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
    private val binSize = PDMEGA.binSize
    private val bam = BamFileReader(bwaSam, refInfo)
    fun run() {
        val records = bam.getAllByChrom()
        println("[Info] load ${records.size} appropriate mappings from SAM file")
        for ( cid in records.indices){
            println("[Info] chromosome ${refInfo[cid]} has ${records[cid].size} mapping segments " +
                    "(%.2f%% covered)".format(records[cid].sumBy { it[0].binCount }.toDouble() * binSize * 100 / refInfo[cid].length))
        }
        val score = analyzeMappingByChrom(records)
        val totalPos = refInfo.sumByDouble { it.binCount*binSize.toDouble() }
        println("Absolute Score: $score")
        println("Ratio Score: ${score/(totalPos*totalPos)}")
    }

    private fun analyzeMappingByChrom(records: ArrayList<ArrayList<MappingSet>>):Double {
        val frecords = records.flatten().toTypedArray()
        var score = 0.0
        val pb = ProgressBar("Pairwise Analysis ", frecords.size.toLong())
        for (i in frecords.indices){
            pb.step()
            score += (frecords[i].getBinCount()*binSize).toDouble()*(frecords[i].getBinCount()*binSize)
            for (j in i+1 until frecords.size) {
                try {
                    score += 2*frecords[i].pairwiseDistance(frecords[j])
                } catch (e: Exception){
                    println("$i $j")
                    println(frecords[i])
                    println(frecords[j])
                    throw e
                }
            }
        }
        pb.close()
        return score
    }
}