package edu.nus

import java.io.File
import htsjdk.samtools.*
import java.lang.Math.*
import me.tongfei.progressbar.ProgressBar
import kotlinx.coroutines.*

data class Mapping(val chrom: Int, val start: Int, var end: Int, val contig: String, val alignmentStart: Int, var alignmentEnd: Int){
    fun forward() = alignmentEnd-alignmentStart>=0
    fun getLength() = end-start
    fun pairwiseDistance(other:Mapping):Double {
        if (chrom==other.chrom && contig!=other.contig) return 0.0
        if (chrom!=other.chrom && contig==other.contig) return 0.0
        if (chrom!=other.chrom && contig!=other.contig) return this.getLength().toDouble()*other.getLength()
        return PDIntegeral(this,other).getPD()
    }
    fun cc50(other: Mapping, threshold: Int): Long {
        if (contig!=other.contig) return 0L
        return CC50Sigma(this, other, threshold).getCC50().toLong()
    }
}
class MappingSet: ArrayList<Mapping>() {
    fun chrom() = this[0].chrom
    fun start() = this[0].start
    fun end() = this[0].end
    fun getLength() = this[0].getLength()
    fun sorted(): MappingSet {
        java.util.Collections.sort(this, compareBy({it.contig},{it.alignmentStart}))
        return this
    }
    private fun allDownstreams(other:MappingSet): Boolean {
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
    fun cc50(other: MappingSet, threshold: Int): Pair<Long, Long> {
        if (chrom()!=other.chrom()) return Pair(0L,0L)
        if (start()>other.start()) return other.cc50(this, threshold)
        val left = max(start()+threshold, other.start())
        val right = min(end()+threshold, other.end())
        val pairnum = max(0L,right-left.toLong())
        var best = 0L
        for (i in this){
            for (j in other)
                best = max(best, i.cc50(j, threshold))
        }
        return Pair(best,pairnum)
    }
}
class BamFileReader(File: File, refInfo: ArrayList<Chromosome>,
                    private val filter: (SAMRecord) -> Boolean = { record:SAMRecord ->
                        !record.readUnmappedFlag  && record.cigarString.all { it in "0123456789MDI=" } //&& record.mappingQuality>10
                    }){
    private val bamFile = SamReaderFactory.makeDefault().open(File)
    private val iterator = bamFile.iterator()
    private val chromToIdx = refInfo.mapIndexed { index: Int, chromInfo: Chromosome -> Pair(chromInfo.name,index) }.toMap()
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
                        if (record.readNegativeStrandFlag) record.alignmentStart-1 else record.alignmentEnd+1)
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
class MappingAnalyzer(bwaSam: File, val refInfo: ArrayList<Chromosome>) {
    private val binSize = PDMEGA.K
    private val bam = BamFileReader(bwaSam, refInfo)
    fun run(): Double {
        val totalPos = refInfo.sumByDouble { it.payload.toDouble() }
        val minLengthToReport = if (PDMEGA.reportLength<0) totalPos*-PDMEGA.reportLength/100 else PDMEGA.reportLength.toDouble()
        val records = bam.getAllByChrom()
        outputAlignment(records)
        for ( chrom in records){
            val ref = refInfo[chrom[0].chrom()]
            if (ref.payload>minLengthToReport)
            println("[Info] $ref " + "(dropped %.2f%%)\t".format((ref.length-ref.payload).toDouble() * 100 / ref.length) +
                    "${chrom.size} Mapping segments " +
                    "(%.2f%% payload covered)".format(chrom.sumBy { it.getLength() }.toDouble() * 100 / ref.payload))
        }
        val score = analyzeMappingByChrom(records)
        println("Genome Payload:     $totalPos\n" +
                "Quality Score Sum:  $score")
        return score/totalPos/totalPos
//        if (!PDMEGA.cc50) return
//        val cc50 = computeCC50(records)
//        println("CC50:           $cc50")
    }
    private fun computeCC50(records: ArrayList<ArrayList<MappingSet>>) = runBlocking {
        operator fun Pair<Long, Long>.plus(other: Pair<Long, Long>) = Pair(first+other.first, second+other.second)
        val frecords = records.flatten().toTypedArray()
        var lowerbound = 0
        var upperbound = refInfo.map { it.length }.max()!!
        val pb = ProgressBar("CC50 binary search", Math.ceil(Math.log(upperbound-lowerbound.toDouble())/Math.log(2.0)).toLong())
        while (lowerbound+2<=upperbound){
            var score = Pair(0L, 0L)
            val mid = (upperbound+lowerbound)/2
            val deferred = frecords.indices.map { i->
                async(Dispatchers.Default) {
                    val self = max(frecords[i].getLength().toLong()-mid, 0L)
                    var ratio = Pair(self, self)
                    for (j in i + 1 until frecords.size) {
                        ratio += frecords[i].cc50(frecords[j], mid)
                    }
                    ratio
                }
            }
            deferred.forEach { score += it.await() }
            println("\n($lowerbound, $mid, $upperbound) Ratio: ${score.first}/${score.second}(${score.first.toDouble()/score.second})")
            pb.step()
            if (score.first.toDouble()/score.second > 0.5) lowerbound = mid else upperbound = mid
        }
        pb.close()
        upperbound
    }
    private fun outputAlignment(records: ArrayList<ArrayList<MappingSet>>) {
        val writer=File(PDMEGA.tmpDir,"MultiAlignment.tab").writer()
        for (chrom in records){
            for (i in chrom){
                writer.write("${i.chrom()}\t${i.start()}\t${i.end()}\t${i.size}\n")
                for (j in i){
                    writer.write("${j.contig}\t${j.alignmentStart}\t${j.alignmentEnd}\n")
                }
            }
        }
        writer.close()
    }
    private fun analyzeMappingByChrom(records: ArrayList<ArrayList<MappingSet>>)= runBlocking {
        val frecords = records.flatten().toTypedArray()
        var score = 0.0
        val deferred = frecords.indices.map { i->
            async(Dispatchers.Default) {
                var s = (frecords[i].getLength()).toDouble()*(frecords[i].getLength())
                for (j in i + 1 until frecords.size) {
                    s += 2 * frecords[i].pairwiseDistance(frecords[j])
                }
                s
            }
        }
        val pb = ProgressBar("PDR Pairwise Analysis ", frecords.size.toLong())
        deferred.forEach {
            pb.step()
            score += it.await()
        }
        pb.close()
        score
    }
}