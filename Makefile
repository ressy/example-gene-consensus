ART = https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier2016.06.05linux64.tgz
DEPTH = 100
SEED = 0
FRAGMENT_LENGTH =  500

# keep intermediate files from pattern rules
.SECONDARY:

READS = $(addprefix reads/sample,$(addsuffix _1.fq,1 2 3))
BAMS = $(addprefix mapping/sample,$(addsuffix .bam,1 2 3))
CONS = $(addprefix consensus/sample,$(addsuffix .fasta,1 2 3))

all: $(CONS)
reads_all: $(READS)
map_all: $(BAMS)
cons_all: $(CONS)

# Simulate paire-end reads based on the expected consensus sequences

$(notdir $(ART)):
	wget $(ART)

art_bin_MountRainier:
	tar xzf $(notdir $(ART))

clean:
	rm -rf reads-by-chr/ reads/ mapping/ consensus/
	rm -f reference/reference.fasta.*

# Simulate paired-end reads for diploid organism

reads-by-chr/%_1.fq reads-by-chr/%_2.fq: samples/%.fasta
	mkdir -p $(dir $@)
	./art_bin_MountRainier/art_illumina \
		-i $^ --out $(subst 1.fq,,$@) --noALN --quiet \
		--rndSeed $(SEED) --seqSys MSv3 --paired \
		--len 150 --fcov $(DEPTH) --mflen $(FRAGMENT_LENGTH) --sdev 50

reads/%_1.fq: reads-by-chr/%.allele1_1.fq reads-by-chr/%.allele2_1.fq
	mkdir -p $(dir $@)
	cat $^ > $@

reads/%_2.fq: reads-by-chr/%.allele1_2.fq reads-by-chr/%.allele2_2.fq
	mkdir -p $(dir $@)
	cat $^ > $@

# Align simulated reads to known reference

mapping/%.bam.bai: mapping/%.bam
	samtools index $^ $@

mapping/%.bam: mapping/%.sam
	samtools sort $^ > $@

mapping/%.sam: reads/%_1.fq reads/%_2.fq reference/reference.fasta.sa
	mkdir -p $(dir $@)
	bwa mem $(basename $(word 3,$^)) $(word 1,$^) $(word 2,$^) > $@

reference/reference.fasta.sa: reference/reference.fasta
	bwa index $^

# Consensus methods

consensus/%.pileup.vcf.gz: mapping/%.bam mapping/%.bam.bai reference/reference.fasta
	mkdir -p $(dir $@)
	bcftools mpileup --fasta-ref $(word 3,$^) $< -Oz --output $@

consensus/%.calls.vcf.gz: consensus/%.pileup.vcf.gz
	bcftools call --variants-only --multiallelic-caller $< -Oz --output $@

consensus/%.vcf.csi: consensus/%.vcf.gz
	bcftools index --force $< 

# What does bcftools norm actually do?  Do I need this? Bypassing for now.
consensus/%.calls.norm.vcf.gz: consensus/%.calls.vcf.gz reference/reference.fasta
	bcftools norm --fasta-ref $(word 2,$^) $(word 1,$^) -Oz --output $@

consensus/%.fasta: consensus/%.calls.vcf.gz consensus/%.calls.vcf.csi reference/reference.fasta
	bcftools consensus --sample mapping/$*.bam --fasta-ref $(word 3,$^) --iupac-codes $< --output $@
