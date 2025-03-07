#!/bin/bash
while getopts 'fq1:fq2:o:p:d:c:r:t:m:' OPT; do
    case $OPT in
        fq1)
            INPUT_fastq1="$OPTARG";;
        fq2)
            INPUT_fastq2="$OPTARG";;
        o)
            OUTDIR="$OPTARG";;
        p)
            Pipeline="$OPTARG";;
        d) 
            ID="$OPTARG";;
        c)
            COVER="$OPTARG";; 
	    r) 
	        MARKER_RATIO="$OPTARG";;
	    m)
	        MARKER_COVERAGE="$OPTARG";;
	    t)
	        THREAD="$OPTARG";;
        ?)
            echo "Usage: `basename $0` [options] filename"
    esac
done

shift $(($OPTIND - 1))

# sh genome.depth.V6.sh -i ./demo/test.read1_read2.fq -p test -o . -g ./demo/Ground_truth.txt -d  90 -c 80 -r 0.5 -t 4 -m 0.7

PREFIX=$OUTDIR

cat $INPUT_fastq1 $INPUT_fastq2 > $OUTDIR/${PREFIX}.R1R2.fastq

diamond blastx --db $Pipeline/database/PhageMarkerProtein.V6.Diamond.dmnd -q $OUTDIR/${PREFIX}.R1R2.fastq -o $OUTDIR/${PREFIX}.sam --max-target-seqs 1 --outfmt 101 --evalue 1e-6 --unal 0  --id $ID --query-cover $COVER

#sam to bam
samtools view -bSh -T $Pipeline/database/all.324056.cancidate.pc.subdist.n324056.representPR.rmPhageHomo.rmMarkerHomo.addRecovered.rmBacHMM.Shinkage.above3.sameVC.faa $OUTDIR/${PREFIX}.sam > $OUTDIR/${PREFIX}.bam

#sort BAM file
samtools sort -@ $THREAD $OUTDIR/${PREFIX}.bam -o $OUTDIR/${PREFIX}.sorted.bam

##sample index
samtools index -@ $THREAD $OUTDIR/${PREFIX}.sorted.bam $OUTDIR/${PREFIX}.sorted.bai

##reads per marker protein
samtools idxstats  $OUTDIR/${PREFIX}.sorted.bam >$OUTDIR/${PREFIX}.sorted.bam.idxstats

##get marker coverage 
bedtools  genomecov -ibam ${PREFIX}.sorted.bam >${PREFIX}.sorted.bam.coverage

##filter mapped reads
perl $Pipeline/script/filter.coverage.pl ${PREFIX}.sorted.bam.coverage ${PREFIX}.sorted.bam.idxstats  ${PREFIX}.sorted.bam.filter.coverage.idxstats $MARKER_COVERAGE

##calculate relative abudance
python $Pipeline/script/reads_abundance.V6.py $OUTDIR/${PREFIX}.sorted.bam.filter.coverage.idxstats $Pipeline/database/all.324056.cancidate.pc.subdist.n324056.representPR.rmPhageHomo.rmMarkerHomo.addRecovered.rmBacHMM.Shinkage.above3.sameVC $MARKER_RATIO $OUTDIR/${PREFIX}.sorted.bam.idxstats.abundance

##add virus taxonomic information
python $Pipeline/script/add_taxonomy.py $Pipeline/database/20230821_VC_tax_lifestyle_host.txt $OUTDIR/${PREFIX}.sorted.bam.idxstats.abundance $OUTDIR/${PREFIX}.sorted.bam.idxstats.addTaxonomy.final.abundance.txt

echo "The SinProVirP pipline is done!"
