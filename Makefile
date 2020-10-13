SHELL ?= /bin/bash
ART ?= https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier2016.06.05linux64.tgz
DEPTH ?= 2
FRAGMENT_LENGTH ?= 500

# keep intermediate files from pattern rules
.SECONDARY:

READS = $(addprefix reads/sample,$(addsuffix _1.fq,1 2 3))
BAMS = $(addprefix mapping/sample,$(addsuffix .bam,1 2 3))
CONS1 = $(addprefix consensus-1/sample,$(addsuffix .fasta,1 2 3))
CONS2 = $(addprefix consensus-2/sample,$(addsuffix .fasta,1 2 3))

all: $(CONS1) $(CONS2)
reads_all: $(READS)
map_all: $(BAMS)
cons1_all: $(CONS)
cons2_all: $(CONS)

clean:
	rm -rf reads-by-chr/ reads/ mapping/ consensus-1/ consensus-2/
	rm -f reference/reference.fasta.*

### Simulate paire-end reads based on the expected consensus sequences

$(notdir $(ART)):
	wget $(ART)

art_bin_MountRainier:
	tar xzf $(notdir $(ART))

reads-by-chr/%_1.fq reads-by-chr/%_2.fq: samples/%.fasta
	mkdir -p $(dir $@)
	./art_bin_MountRainier/art_illumina \
		-i $^ --out $(subst 1.fq,,$@) --noALN --quiet \
		--rndSeed 5$(shell printf "%d" 0x$$(echo $^ | crc32 /dev/stdin)) --seqSys MSv3 --paired \
		--len 150 --fcov $(DEPTH) --mflen $(FRAGMENT_LENGTH) --sdev 50

reads/%_1.fq: reads-by-chr/%.allele1_1.fq reads-by-chr/%.allele2_1.fq
	mkdir -p $(dir $@)
	cat $^ > $@

reads/%_2.fq: reads-by-chr/%.allele1_2.fq reads-by-chr/%.allele2_2.fq
	mkdir -p $(dir $@)
	cat $^ > $@

### Align simulated reads to known reference

mapping/%.bam.bai: mapping/%.bam
	samtools index $^ $@

mapping/%.bam: mapping/%.sam
	samtools sort $^ > $@

mapping/%.sam: reads/%_1.fq reads/%_2.fq reference/reference.fasta.sa
	mkdir -p $(dir $@)
	bwa mem $(basename $(word 3,$^)) $(word 1,$^) $(word 2,$^) > $@

reference/reference.fasta.sa: reference/reference.fasta
	bwa index $^

### Consensus methods

# Consensus 1: bcftools

consensus-1/%.pileup.vcf.gz: mapping/%.bam mapping/%.bam.bai reference/reference.fasta
	mkdir -p $(dir $@)
	bcftools mpileup --fasta-ref $(word 3,$^) $< -Oz --output $@

consensus-1/%.calls.vcf.gz: consensus-1/%.pileup.vcf.gz
	bcftools call --variants-only --multiallelic-caller $< -Oz --output $@

consensus-1/%.vcf.csi: consensus-1/%.vcf.gz
	bcftools index --force $< 

# What does bcftools norm actually do?  Do I need this? Bypassing for now.
consensus-1/%.calls.norm.vcf.gz: consensus-1/%.calls.vcf.gz reference/reference.fasta
	bcftools norm --fasta-ref $(word 2,$^) $(word 1,$^) -Oz --output $@

consensus-1/%.fasta: consensus-1/%.calls.vcf.gz consensus-1/%.calls.vcf.csi reference/reference.fasta
	bcftools consensus --sample mapping/$*.bam --fasta-ref $(word 3,$^) --iupac-codes $< | seqtk seq -l 0 > $@

# Consensus 2: vcfutils.pl
#
# vcfutils.pl vcf2fq lowercases positions below the minimum DP and/or MQ, and
# sets some to "n" if there isn't a clear allele call.  Zero depth regions
# within the region covered are lowercased, but zero-depth regions before/after
# the region covered are trimmed.

consensus-2/%.fasta: consensus-2/%.fastq
	seqtk seq -a $^ > $@

consensus-2/%.fastq: consensus-2/%.calls.vcf.gz
	zcat $^ | vcfutils.pl vcf2fq -d 3 > $@

consensus-2/%.calls.vcf.gz: consensus-2/%.pileup.vcf.gz
	bcftools call --consensus-caller $< -Oz --output $@

consensus-2/%.pileup.vcf.gz: mapping/%.bam mapping/%.bam.bai reference/reference.fasta
	mkdir -p $(dir $@)
	samtools mpileup -v -f $(word 3,$^) $< > $@
