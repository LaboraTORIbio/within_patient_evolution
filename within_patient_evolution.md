# Initial setup

The following directories were created in the current working directory, `Within_patient_evolution`, to store sequencing and analysis data.

```sh
# For section: 'De novo assembling and genomic analysis of E. coli J53 carrying different PVs'
mkdir -p plasmid_variants_J53/reads_trimmed/reads_raw
# For sections: 'Assembly of within-patient evolution strains' and 'Construction of phylogenetic trees'
mkdir -p within_patient_evol_cured/reads_trimmed/reads_raw
mkdir ../Reads_long
mkdir -p ../Reads_WT_pOXA-48/reads_raw
# For storing the GenBank file pOXA-48_K8.gb
mkdir -p ../Closed_sequences/plasmids
```

Download sequencing data from BioProject PRJNA838107 in the specified subdirectories.

| SRA | BioSample | Strain | Library ID | Download in directory |
| --- | --- | --- | --- | --- |
| SRR19213026 | SAMN28405001 | E.coli_J53/PV-free | J53javier_S106 | plasmid_variants_J53/reads_trimmed/reads_raw |
| SRR19213025 | SAMN28405002 | E.coli_J53/PV-A | TCK165_Clon_1_S97 | plasmid_variants_J53/reads_trimmed/reads_raw |
| SRR19213024 | SAMN28405003 | E.coli_J53/PV-B | TCK229_Clon_2_S116 | plasmid_variants_J53/reads_trimmed/reads_raw |
| SRR19213023 | SAMN28405004 | E.coli_J53/PV-C | TCK298_Clon_2_S121 | plasmid_variants_J53/reads_trimmed/reads_raw |
| SRR19213022 | SAMN28405005 | E.coli_J53/PV-D | TCC704_Clon_2_S107 | plasmid_variants_J53/reads_trimmed/reads_raw |
| SRR19213021 | SAMN28405006 | E.coli_J53/PV-F | TCJ57_Clon_1_S100 | plasmid_variants_J53/reads_trimmed/reads_raw |
| SRR19213020 | SAMN28405007 | E.coli_J53/PV-E | TCC288_Clon_3_S103 | plasmid_variants_J53/reads_trimmed/reads_raw |
| SRR19213019 | SAMN28405008 | E.coli_J53/PV-G | TCH53_Clon_2_S110 | plasmid_variants_J53/reads_trimmed/reads_raw |
| SRR19213018 | SAMN28405009 | E.coli_J53/PV-I | TCK163_Clon_11_S113 | plasmid_variants_J53/reads_trimmed/reads_raw |
| SRR19213016 | SAMN28405010 | E.coli_J53/PV-K | J53_pC289clon1_32B5_S206 | plasmid_variants_J53/reads_trimmed/reads_raw |
| SRR19213015 | SAMN28405011 | E.coli_J53/PV-J | TCK26_Clon_2_S119 | plasmid_variants_J53/reads_trimmed/reads_raw |
| SRR19213014 | SAMN28405012 | E.coli_J53/PV-H | TCK57_Clon_1_S122 | plasmid_variants_J53/reads_trimmed/reads_raw |
| SRR19213013 | SAMN28405013 | E.coli_J53/PV-L | TCK273_Clon_1_S120 | plasmid_variants_J53/reads_trimmed/reads_raw |
| SRR19213012 | SAMN28405014 | E.coli_J53/PV-M | TCG14_Clon_1_S135 | plasmid_variants_J53/reads_trimmed/reads_raw |
| SRR19213011 | SAMN28405015 | E.coli_J53/PV-N | TCK153_new_S208 | plasmid_variants_J53/reads_trimmed/reads_raw |
| SRR19213040 | SAMN28404986 | C288 | C288WT_S139 | within_patient_evol_cured/reads_trimmed/reads_raw |
| SRR19213039 | SAMN28404986 | C288 | C288WT_nanopore | ../Reads_long |
| SRR19213028 | SAMN28404987 | C289 | C289WT_S142 | within_patient_evol_cured/reads_trimmed/reads_raw |
| SRR19213017 | SAMN28404987 | C289 | C289WT-31C7_nanopore | ../Reads_long |
| SRR19213010 | SAMN28404988 | C289cured | C289c1_Clon_1_S143 | within_patient_evol_cured/reads_trimmed/reads_raw |
| SRR19213009 | SAMN28404989 | K153 | K153WT_S103 | within_patient_evol_cured/reads_trimmed/reads_raw |
| SRR19213008 | SAMN28404989 | K153 | K153WT-30A7_nanopore | ../Reads_long |
| SRR19213007 | SAMN28404990 | K229 | K229newWT_S102 | within_patient_evol_cured/reads_trimmed/reads_raw |
| SRR19213006 | SAMN28404990 | K229 | K229new-30C6_nanopore | ../Reads_long |
| SRR19213005 | SAMN28404991 | K229cured | K229C5_S89 | within_patient_evol_cured/reads_trimmed/reads_raw |
| SRR19213038 | SAMN28404992 | K163 | K163WT_S99 | within_patient_evol_cured/reads_trimmed/reads_raw |
| SRR19213037 | SAMN28404992 | K163 | K163WT-28F1_nanopore | ../Reads_long |
| SRR19213036 | SAMN28404993 | K165-2 | K165WT_S96 | within_patient_evol_cured/reads_trimmed/reads_raw |
| SRR19213035 | SAMN28404993 | K165-2 | K165WT_nanopore | ../Reads_long |
| SRR19213034 | SAMN28404994 | K165cured | K165c5_S97 | within_patient_evol_cured/reads_trimmed/reads_raw |
| SRR19213033 | SAMN28404995 | K165-3 | 71_S107 | within_patient_evol_cured/reads_trimmed/reads_raw |
| SRR19213032 | SAMN28404996 | K165-4 | 72_S108 | within_patient_evol_cured/reads_trimmed/reads_raw |
| SRR19213031 | SAMN28404997 | K165-5 | 73_S109 | within_patient_evol_cured/reads_trimmed/reads_raw |
| SRR19213030 | SAMN28404998 | K165-1 | 74_S110 | within_patient_evol_cured/reads_trimmed/reads_raw |
| SRR19213029 | SAMN28404999 | K165-6 | 75_S111 | within_patient_evol_cured/reads_trimmed/reads_raw |
| SRR19213027 | SAMN28405000 | K165-7 | 76_S112 | within_patient_evol_cured/reads_trimmed/reads_raw |

Raw Illumina reads from BioProject PRJNA626430 used for constructing phylogenetic trees of *E. coli* and *K. pneumoniae*  are stored in `../Reads_WT_pOXA-48/reads_raw`. The GenBank annotation file for plasmid pOXA-48 (MT441554) is stored in `../Closed_sequences/plasmids` and the `ltrA.fasta` file is placed in the current directory.

# Illumina and Nanopore read preprocessing

Trim Galore v0.6.4 was used to trim Illumina sequences and generate FastQC reports. All read datasets reached good quality.

```sh
cd plasmid_variants_J53/reads_trimmed/reads_raw/

for fq1 in *R1_001.fastq.gz
do
	fq2=${fq1%%R1_001.fastq.gz}"R2_001.fastq.gz"
	name=${fq1/_[A-Z][0-9]*_R1_001.fastq.gz/}
	trim_galore --quality 20 --nextera --length 50 --fastqc --basename MiGS_$name --output_dir .. --paired $fq1 $fq2
done

cd ../../../within_patient_evol_cured/reads_trimmed/reads_raw/

for fq1 in *R1_001.fastq.gz
do
	fq2=${fq1%%R1_001.fastq.gz}"R2_001.fastq.gz"
	name=${fq1/_[A-Z][0-9]*_R1_001.fastq.gz/}
	trim_galore --quality 20 --nextera --length 50 --fastqc --basename MiGS_$name --output_dir .. --paired $fq1 $fq2
done
```

Nanopore reads were filtered with filtlong v0.2.1 with options as suggested by the authors to obtain a subset of high identity reads with a minimum read depth of 85x.

```sh
cd ../.. # From within_patient_evol_cured

filtlong --min_length 1000 --keep_percent 90 --target_bases 500000000 ../../Reads_long/C288WT_nanopore.fastq.gz | gzip > ../../Reads_long/C288WT_nanopore_filt.fastq.gz
filtlong --min_length 1000 --keep_percent 90 --target_bases 500000000 ../../Reads_long/C289WT-31C7_nanopore.fastq.gz | gzip > ../../Reads_long/C289WT-31C7_nanopore_filt.fastq.gz
filtlong --min_length 1000 --keep_percent 90 --target_bases 500000000 ../../Reads_long/K153WT-30A7_nanopore.fastq.gz | gzip > ../../Reads_long/K153WT-30A7_nanopore_filt.fastq.gz
filtlong --min_length 1000 --keep_percent 90 --target_bases 500000000 ../../Reads_long/K229new-30C6_nanopore.fastq.gz | gzip > ../../Reads_long/K229new-30C6_nanopore_filt.fastq.gz
filtlong --min_length 1000 --keep_percent 90 --target_bases 500000000 ../../Reads_long/K163WT-28F1_nanopore.fastq.gz | gzip > ../../Reads_long/K163WT-28F1_nanopore_filt.fastq.gz
filtlong --min_length 1000 --keep_percent 90 --target_bases 500000000 ../../Reads_long/K165WT_nanopore.fastq.gz | gzip > ../../Reads_long/K165WT_nanopore_filt.fastq.gz
```

Reads from BioProject PRJNA626430 stored in `../Reads_WT_pOXA-48/reads_raw` were also trimmed with Trim Galore (--quality 20 --illumina --length 50 --fastqc). Trimmed reads were stored in `../Reads_WT_pOXA-48/trimmed` and were renamed to the corresponding strain names.


# _De novo_ assembling and genomic analysis of *E. coli* J53 carrying different PVs

The *E. coli* J53 strains (carrying different PVs and plasmid-free) used in this study were sequenced to control isogenic conditions. First, the genome of J53 was assembled using SPAdes v3.15.2 and annotated with Prokka v1.14.6.

```sh
cd ../plasmid_variants_J53

spades.py --isolate --cov-cutoff auto -o assemblies_SPAdes/J53javier -1 reads_trimmed/MiGS_J53javier_val_1.fq.gz -2 reads_trimmed/MiGS_J53javier_val_2.fq.gz
quast.py -o assemblies_SPAdes/quast_reports/J53javier --glimmer assemblies_SPAdes/J53javier
prokka --outdir annotation_prokka_J53 --prefix J53 assemblies_SPAdes/J53javier/contigs.fasta
prokka --outdir annotation_prokka_J53_compliant --prefix J53 --compliant assemblies_SPAdes/J53javier/contigs.fasta
```

The *ltrA* gene of pOXA-48_K8 was blasted (BLASTn v2.11.0) against J53 to confirm its absence.

```sh
mkdir blastn_ltrA
makeblastdb -in assemblies_SPAdes/J53javier/contigs.fasta -dbtype nucl
blastn -query ../ltrA.fasta -db assemblies_SPAdes/J53javier/contigs.fasta -outfmt 6 > blastn_ltrA/J53javier.tsv
```

Variant calling was performed to detect mutations in the chromosomes and pOXA-48 PVs. Chromosomal mutations were analyzed following this workflow. First, Snippy v4.6.0 was used to identify variants in the draft genome of J53 by mapping the Illumina reads back to its assembly.

```sh
snippy --report --outdir variants_snippy_ref_J53_map_J53 --ref annotation_prokka_J53_compliant/J53.gbk --R1 reads_trimmed/MiGS_J53javier_val_1.fq.gz --R2 reads_trimmed/MiGS_J53javier_val_2.fq.gz
```

Next, variants in the J53 strains carrying PVs were called with Snippy v4.6.0 and breseq v0.35.6 using as reference the annotated J53 draft genome. Variants that matched the ones found in J53 were discarded since they could constitute assembly errors. From the breseq output, only predicted mutations and unassigned missing coverage (MC) were analyzed, as Illumina data cannot provide good new junction evidence (JC) from draft assemblies. Finally, Snippy was used to call mutations in a reverse approach, mapping the reads of *E. coli* J53 against the assemblies of the PVs carriers. Only mutations that were identified in both directions (direct and reverse) and by both Snippy and breseq were considered. Informatively, these accepted variants generally had higher and more homogeneous read depth than those considered false calls.

To confirm that pOXA-48 did not accumulate mutations during conjugation/transformation, Snippy and breseq were run using the sequence of pOXA-48_K8 (MT441554) as reference. In this case, only mutations called by both programs, as well as MC and JC from breseq, were considered. The sequence of the *ltrA* gene of pOXA-48 was blasted against the assemblies of the PVs carriers to confirm its presence and absence from the respective pOXA-48 variants. The contig containing *ltrA* had similar length in all assemblies and different coverage than the coverage of chromosomal contigs, indicating that the *ltrA* did not move into the chromosome of J53 in any transconjugant/transformant.

Plasmid replicons were detected with ABRicate v1.0.1 using the plasmidfinder database to confirm that no other plasmid was conjugated or mobilized during pOXA-48 conjugation or transformation. Additionally, the Resfinder database of ABRicate v1.0.1 was used to discard the presence of other MGE carrying resistance genes.

```sh
for fq1 in reads_trimmed/*val_1.fq.gz
do
	fq2=${fq1%%1.fq.gz}"2.fq.gz"
	strain=${fq1::-12}
	strain=${strain:14}

	if [[ $strain != *"J53javier"* ]]; then
		# De novo assembly with SPAdes of J53 carrying different PVs
		spades.py --isolate --cov-cutoff auto -o assemblies_SPAdes/$strain -1 $fq1 -2 $fq2

		# Snippy using as reference J53 and mapping short reads of J53 carrying different PVs
		snippy --report --outdir variants_snippy_ref_J53_map_TC/ref_J53_map_$strain --ref annotation_prokka_J53_compliant/J53.gbk --R1 $fq1 --R2 $fq2

		# Snippy using as reference pOXA-48_K8 (accession number MT441554) and mapping short reads of J53 carrying different PVs
		snippy --report --outdir variants_snippy_ref_pOXA-48_K8_map_TC/ref_pOXA-48_K8_map_$strain --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 $fq1 --R2 $fq2

		# Breseq using as reference J53 and mapping short reads of J53 carrying different PVs
		breseq -o variants_breseq_ref_J53_map_TC/ref_J53_map_$strain -c annotation_prokka_J53/J53.gff $fq1 $fq2

		# Breseq using as reference pOXA-48_K8 (accession number MT441554) and mapping short reads of J53 carrying different PVs
		breseq -o variants_breseq_ref_pOXA-48_K8_map_TC/ref_pOXA-48_K8_map_$strain -r ../../Closed_sequences/plasmids/pOXA-48_K8.gb $fq1 $fq2
	fi
done

mkdir plasmids_abricate
mkdir resistance_abricate

for contigs in assemblies_SPAdes/*/contigs.fasta
do
	strain=${contigs:18}
	strain=${strain::-14}

	# Quast quality analysis of assemblies of J53 carrying different PVs
	quast.py -o assemblies_SPAdes/quast_reports/$strain --glimmer $contigs

	# Snippy using as reference the different J53 carrying different PVs and mapping J53
	snippy --report --outdir variants_snippy_ref_TC_map_J53/ref_$strain"_map_J53" --ref $contigs --R1 reads_trimmed/MiGS_J53javier_val_1.fq.gz --R2 reads_trimmed/MiGS_J53javier_val_2.fq.gz

	# Blasting the ltrA gene against J53 carrying different PVs
	makeblastdb -in $contigs -dbtype nucl
	blastn -query ../ltrA.fasta -db $contigs -outfmt 6 > blastn_ltrA/$strain".tsv"

	# Detecting plasmid replicons and resistance genes with ABRicate
	abricate --db plasmidfinder $contigs > plasmids_abricate/$strain".tsv"
	abricate --db resfinder $contigs > resistance_abricate/$strain".tsv"
done

# Summary of ABRicate results
abricate --summary plasmids_abricate/* > plasmids_abricate/summary_plasmids.tsv
abricate --summary resistance_abricate/* > resistance_abricate/summary_resistance.tsv
```


# Analysis of within-patient evolution strains

## Generating closed genomes

Hybrid assemblies of the six strains (C288, C289, K153, K229, K163 and K165) were obtained using Unicycler v0.4.9 with default parameters.

```sh
cd ../within_patient_evol_cured

unicycler -t 18 --spades_path ~/Software/SPAdes-3.13.0-Linux/bin/spades.py -1 reads_trimmed/MiGS_C288WT_val_1.fq.gz -2 reads_trimmed/MiGS_C288WT_val_2.fq.gz -l ../../Reads_long/C288WT_nanopore_filt.fastq.gz -o assemblies_unicycler_hybrid/C288WT
unicycler -t 18 --spades_path ~/Software/SPAdes-3.13.0-Linux/bin/spades.py -1 reads_trimmed/MiGS_C289WT_val_1.fq.gz -2 reads_trimmed/MiGS_C289WT_val_2.fq.gz -l ../../Reads_long/C289WT-31C7_nanopore_filt.fastq.gz -o assemblies_unicycler_hybrid/C289WT-31C7
unicycler -t 18 --spades_path ~/Software/SPAdes-3.13.0-Linux/bin/spades.py -1 reads_trimmed/MiGS_K153_val_1.fq.gz -2 reads_trimmed/MiGS_K153_val_2.fq.gz -l ../../Reads_long/K153WT-30A7_nanopore_filt.fastq.gz -o assemblies_unicycler_hybrid/K153WT-30A7
unicycler -t 18 --spades_path ~/Software/SPAdes-3.13.0-Linux/bin/spades.py -1 reads_trimmed/MiGS_K229_val_1.fq.gz -2 reads_trimmed/MiGS_K229_val_2.fq.gz -l ../../Reads_long/K229new-30C6_nanopore_filt.fastq.gz -o assemblies_unicycler_hybrid/K229new-30C6
unicycler -t 18 --spades_path ~/Software/SPAdes-3.13.0-Linux/bin/spades.py -1 reads_trimmed/MiGS_K163_val_1.fq.gz -2 reads_trimmed/MiGS_K163_val_2.fq.gz -l ../../Reads_long/K163WT-28F1_nanopore_filt.fastq.gz -o assemblies_unicycler_hybrid/K163WT-28F1
unicycler -t 18 --spades_path ~/Software/SPAdes-3.13.0-Linux/bin/spades.py -1 reads_trimmed/MiGS_K165_val_1.fq.gz -2 reads_trimmed/MiGS_K165_val_2.fq.gz -l ../../Reads_long/K165WT_nanopore_filt.fastq.gz -o assemblies_unicycler_hybrid/K165WT
```

Additionally, for the cases where Unicycler was not able to circularize contigs (C288, C289, K153 and K229), long reads were assembled with Flye v2.9 and circularization was confirmed by inspecting the assembly graphs in Bandage.

```sh
flye --nano-raw ../../Reads_long/C288WT_nanopore_filt.fastq.gz --threads 18 --plasmids --out-dir assemblies_flye/C288WT
flye --nano-raw ../../Reads_long/C289WT-31C7_nanopore_filt.fastq.gz --threads 18 --plasmids --out-dir assemblies_flye/C289WT-31C7
flye --nano-raw ../../Reads_long/K153WT-30A7_nanopore_filt.fastq.gz --threads 18 --plasmids --out-dir assemblies_flye/K153WT-30A7
flye --nano-raw ../../Reads_long/K229new-30C6_nanopore_filt.fastq.gz --threads 18 --plasmids --out-dir assemblies_flye/K229new-30C6
```

Consensus sequences of these assemblies were obtained with Medaka v1.4.3.

```sh
medaka_consensus -i ../../Reads_long/C288WT_nanopore_filt.fastq.gz -d assemblies_flye/C288WT/assembly.fasta -o assemblies_flye/C288WT/medaka
medaka_consensus -i ../../Reads_long/C289WT-31C7_nanopore_filt.fastq.gz -d assemblies_flye/C289WT-31C7/assembly.fasta -o assemblies_flye/C289WT-31C7/medaka
medaka_consensus -i ../../Reads_long/K153WT-30A7_nanopore_filt.fastq.gz -d assemblies_flye/K153WT-30A7/assembly.fasta -o assemblies_flye/K153WT-30A7/medaka
medaka_consensus -i ../../Reads_long/K229new-30C6_nanopore_filt.fastq.gz -d assemblies_flye/K229new-30C6/assembly.fasta -o assemblies_flye/K229new-30C6/medaka
```

Several rounds of Pilon v1.24 were performed mapping the trimmed Illumina reads until no further changes in sequence were observed.

```sh
# C288 round 1
bwa index assemblies_flye/C288WT/medaka/consensus.fasta
bwa mem -t 18 assemblies_flye/C288WT/medaka/consensus.fasta reads_trimmed/MiGS_C288WT_val_1.fq.gz reads_trimmed/MiGS_C288WT_val_2.fq.gz | samtools sort -o assemblies_flye/C288WT/medaka/alignments.bam
samtools index assemblies_flye/C288WT/medaka/alignments.bam
pilon --changes --genome assemblies_flye/C288WT/medaka/consensus.fasta --frags assemblies_flye/C288WT/medaka/alignments.bam --output C288WT --outdir assemblies_flye/C288WT/pilon/
# C288 round 2
bwa index assemblies_flye/C288WT/pilon/C288WT.fasta
bwa mem -t 18 assemblies_flye/C288WT/pilon/C288WT.fasta reads_trimmed/MiGS_C288WT_val_1.fq.gz reads_trimmed/MiGS_C288WT_val_2.fq.gz | samtools sort -o assemblies_flye/C288WT/medaka/alignments.bam
samtools index assemblies_flye/C288WT/medaka/alignments.bam
pilon --changes --genome assemblies_flye/C288WT/pilon/C288WT.fasta --frags assemblies_flye/C288WT/medaka/alignments.bam --output C288WT_2 --outdir assemblies_flye/C288WT/pilon/

# C289 round 1
bwa index assemblies_flye/C289WT-31C7/medaka/consensus.fasta
bwa mem -t 18 assemblies_flye/C289WT-31C7/medaka/consensus.fasta reads_trimmed/MiGS_C289WT_val_1.fq.gz reads_trimmed/MiGS_C289WT_val_2.fq.gz | samtools sort -o assemblies_flye/C289WT-31C7/medaka/alignments.bam
samtools index assemblies_flye/C289WT-31C7/medaka/alignments.bam
pilon --changes --genome assemblies_flye/C289WT-31C7/medaka/consensus.fasta --frags assemblies_flye/C289WT-31C7/medaka/alignments.bam --output C289WT-31C7 --outdir assemblies_flye/C289WT-31C7/pilon/
# C289 round 2
bwa index assemblies_flye/C289WT-31C7/pilon/C289WT-31C7.fasta
bwa mem -t 18 assemblies_flye/C289WT-31C7/pilon/C289WT-31C7.fasta reads_trimmed/MiGS_C289WT_val_1.fq.gz reads_trimmed/MiGS_C289WT_val_2.fq.gz | samtools sort -o assemblies_flye/C289WT-31C7/medaka/alignments.bam
samtools index assemblies_flye/C289WT-31C7/medaka/alignments.bam
pilon --changes --genome assemblies_flye/C289WT-31C7/pilon/C289WT-31C7.fasta --frags assemblies_flye/C289WT-31C7/medaka/alignments.bam --output C289WT-31C7_2 --outdir assemblies_flye/C289WT-31C7/pilon/

# K153 round 1
bwa index assemblies_flye/K153WT-30A7/medaka/consensus.fasta
bwa mem -t 18 assemblies_flye/K153WT-30A7/medaka/consensus.fasta reads_trimmed/MiGS_K153_val_1.fq.gz reads_trimmed/MiGS_K153_val_2.fq.gz | samtools sort -o assemblies_flye/K153WT-30A7/medaka/alignments.bam
samtools index assemblies_flye/K153WT-30A7/medaka/alignments.bam
pilon --changes --genome assemblies_flye/K153WT-30A7/medaka/consensus.fasta --frags assemblies_flye/K153WT-30A7/medaka/alignments.bam --output K153WT-30A7 --outdir assemblies_flye/K153WT-30A7/pilon/
# K153 round 2
bwa index assemblies_flye/K153WT-30A7/pilon/K153WT-30A7.fasta
bwa mem -t 18 assemblies_flye/K153WT-30A7/pilon/K153WT-30A7.fasta reads_trimmed/MiGS_K153_val_1.fq.gz reads_trimmed/MiGS_K153_val_2.fq.gz | samtools sort -o assemblies_flye/K153WT-30A7/medaka/alignments.bam
samtools index assemblies_flye/K153WT-30A7/medaka/alignments.bam
pilon --changes --genome assemblies_flye/K153WT-30A7/pilon/K153WT-30A7.fasta --frags assemblies_flye/K153WT-30A7/medaka/alignments.bam --output K153WT-30A7_2 --outdir assemblies_flye/K153WT-30A7/pilon/

# K229 round 1
bwa index assemblies_flye/K229new-30C6/medaka/consensus.fasta
bwa mem -t 18 assemblies_flye/K229new-30C6/medaka/consensus.fasta reads_trimmed/MiGS_K229_val_1.fq.gz reads_trimmed/MiGS_K229_val_2.fq.gz | samtools sort -o assemblies_flye/K229new-30C6/medaka/alignments.bam
samtools index assemblies_flye/K229new-30C6/medaka/alignments.bam
pilon --changes --genome assemblies_flye/K229new-30C6/medaka/consensus.fasta --frags assemblies_flye/K229new-30C6/medaka/alignments.bam --output K229new-30C6 --outdir assemblies_flye/K229new-30C6/pilon/
# K229 round 2
bwa index assemblies_flye/K229new-30C6/pilon/K229new-30C6.fasta
bwa mem -t 18 assemblies_flye/K229new-30C6/pilon/K229new-30C6.fasta reads_trimmed/MiGS_K229_val_1.fq.gz reads_trimmed/MiGS_K229_val_2.fq.gz | samtools sort -o assemblies_flye/K229new-30C6/medaka/alignments.bam
samtools index assemblies_flye/K229new-30C6/medaka/alignments.bam
pilon --changes --genome assemblies_flye/K229new-30C6/pilon/K229new-30C6.fasta --frags assemblies_flye/K229new-30C6/medaka/alignments.bam --output K229new-30C6_2 --outdir assemblies_flye/K229new-30C6/pilon/
# K229 round 3
bwa index assemblies_flye/K229new-30C6/pilon/K229new-30C6_2.fasta
bwa mem -t 18 assemblies_flye/K229new-30C6/pilon/K229new-30C6_2.fasta reads_trimmed/MiGS_K229_val_1.fq.gz reads_trimmed/MiGS_K229_val_2.fq.gz | samtools sort -o assemblies_flye/K229new-30C6/medaka/alignments.bam
samtools index assemblies_flye/K229new-30C6/medaka/alignments.bam
pilon --changes --genome assemblies_flye/K229new-30C6/pilon/K229new-30C6_2.fasta --frags assemblies_flye/K229new-30C6/medaka/alignments.bam --output K229new-30C6_3 --outdir assemblies_flye/K229new-30C6/pilon/
```

Finally, contigs were rotated with Circlator fixstart v1.5.5. The long read assemblies were further quality controlled by inspecting the mappings of Illumina reads in IGV.

```sh
for dir in assemblies_flye/*/; do mkdir -- "$dir/circlator"; done
cd assemblies_flye/C288WT/circlator
circlator fixstart ../pilon/C288WT_2.fasta circlator
cd ../../C289WT-31C7/circlator
circlator fixstart ../pilon/C289WT-31C7_2.fasta circlator
cd ../../K153WT-30A7/circlator
circlator fixstart ../pilon/K153WT-30A7_2.fasta circlator
cd ../../K229new-30C6/circlator
circlator fixstart ../pilon/K229new-30C6_3.fasta circlator
```

Resulting assemblies were copied to `final_closed_assemblies`. Correctness of the assemblies of the pOXA-48 plasmids was confirmed by mapping the short and long reads with BWA-MEM v0.7.17 (14) and minimap2 v2.21 (15) (option -ax map-ont), respectively, to the reference pOXA-48_K8 sequence (MT441554), and by aligning the obtained assemblies to the pOXA-48_K8 reference with minimap2 (option -ax asm5). Alignments were visualized in IGV.

```sh
# Mapping short reads to pOXA-48_K8 (MT441554)

mkdir alignments
bwa index ../../Closed_sequences/plasmids/pOXA-48_K8.fasta

bwa mem -t 18 ../../Closed_sequences/plasmids/pOXA-48_K8.fasta reads_trimmed/MiGS_C288WT_val_1.fq.gz reads_trimmed/MiGS_C288WT_val_2.fq.gz | samtools sort -o alignments/ref_pOXA-48_K8_align_C288_illu.bam
samtools index alignments/ref_pOXA-48_K8_align_C288_illu.bam

bwa mem -t 18 ../../Closed_sequences/plasmids/pOXA-48_K8.fasta reads_trimmed/MiGS_C289WT_val_1.fq.gz reads_trimmed/MiGS_C289WT_val_2.fq.gz | samtools sort -o alignments/ref_pOXA-48_K8_align_C289_illu.bam
samtools index alignments/ref_pOXA-48_K8_align_C289_illu.bam

bwa mem -t 18 ../../Closed_sequences/plasmids/pOXA-48_K8.fasta reads_trimmed/MiGS_K153_val_1.fq.gz reads_trimmed/MiGS_K153_val_2.fq.gz | samtools sort -o alignments/ref_pOXA-48_K8_align_K153_illu.bam
samtools index alignments/ref_pOXA-48_K8_align_K153_illu.bam

bwa mem -t 18 ../../Closed_sequences/plasmids/pOXA-48_K8.fasta reads_trimmed/MiGS_K229_val_1.fq.gz reads_trimmed/MiGS_K229_val_2.fq.gz | samtools sort -o alignments/ref_pOXA-48_K8_align_K229_illu.bam
samtools index alignments/ref_pOXA-48_K8_align_K229_illu.bam

bwa mem -t 18 ../../Closed_sequences/plasmids/pOXA-48_K8.fasta reads_trimmed/MiGS_K163_val_1.fq.gz reads_trimmed/MiGS_K163_val_2.fq.gz | samtools sort -o alignments/ref_pOXA-48_K8_align_K163_illu.bam
samtools index alignments/ref_pOXA-48_K8_align_K163_illu.bam

bwa mem -t 18 ../../Closed_sequences/plasmids/pOXA-48_K8.fasta reads_trimmed/MiGS_K165_val_1.fq.gz reads_trimmed/MiGS_K165_val_2.fq.gz | samtools sort -o alignments/ref_pOXA-48_K8_align_K165_illu.bam
samtools index alignments/ref_pOXA-48_K8_align_K165_illu.bam

# Mapping long reads to pOXA-48_K8 (MT441554)

minimap2 -ax map-ont ../../Closed_sequences/plasmids/pOXA-48_K8.fasta ../../Reads_long/C288WT_nanopore_filt.fastq.gz | samtools sort -o alignments/ref_pOXA-48_K8_align_C288_nano.bam
samtools index alignments/ref_pOXA-48_K8_align_C288_nano.bam

minimap2 -ax map-ont ../../Closed_sequences/plasmids/pOXA-48_K8.fasta ../../Reads_long/C289WT-31C7_nanopore_filt.fastq.gz | samtools sort -o alignments/ref_pOXA-48_K8_align_C289_nano.bam
samtools index alignments/ref_pOXA-48_K8_align_C289_nano.bam

minimap2 -ax map-ont ../../Closed_sequences/plasmids/pOXA-48_K8.fasta ../../Reads_long/K153WT-30A7_nanopore_filt.fastq.gz | samtools sort -o alignments/ref_pOXA-48_K8_align_K153_nano.bam
samtools index alignments/ref_pOXA-48_K8_align_K153_nano.bam

minimap2 -ax map-ont ../../Closed_sequences/plasmids/pOXA-48_K8.fasta ../../Reads_long/K229new-30C6_nanopore_filt.fastq.gz | samtools sort -o alignments/ref_pOXA-48_K8_align_K229_nano.bam
samtools index alignments/ref_pOXA-48_K8_align_K229_nano.bam

minimap2 -ax map-ont ../../Closed_sequences/plasmids/pOXA-48_K8.fasta ../../Reads_long/K163WT-28F1_nanopore_filt.fastq.gz | samtools sort -o alignments/ref_pOXA-48_K8_align_K163_nano.bam
samtools index alignments/ref_pOXA-48_K8_align_K163_nano.bam

minimap2 -ax map-ont ../../Closed_sequences/plasmids/pOXA-48_K8.fasta ../../Reads_long/K165WT_nanopore_filt.fastq.gz | samtools sort -o alignments/ref_pOXA-48_K8_align_K165_nano.bam
samtools index alignments/ref_pOXA-48_K8_align_K165_nano.bam

# Aligning assemblies to pOXA-48_K8 (MT441554)

minimap2 -ax asm5 ../../Closed_sequences/plasmids/pOXA-48_K8.fasta final_closed_assemblies/C288WT.fasta | samtools sort -o alignments/ref_pOXA-48_K8_map_C288flye.bam
samtools index alignments/ref_pOXA-48_K8_map_C288flye.bam

minimap2 -ax asm5 ../../Closed_sequences/plasmids/pOXA-48_K8.fasta final_closed_assemblies/C289WT-31C7.fasta | samtools sort -o alignments/ref_pOXA-48_K8_map_C289flye.bam
samtools index alignments/ref_pOXA-48_K8_map_C289flye.bam

minimap2 -ax asm5 ../../Closed_sequences/plasmids/pOXA-48_K8.fasta final_closed_assemblies/K153WT-30A7.fasta | samtools sort -o alignments/ref_pOXA-48_K8_map_K153flye.bam
samtools index alignments/ref_pOXA-48_K8_map_K153flye.bam

minimap2 -ax asm5 ../../Closed_sequences/plasmids/pOXA-48_K8.fasta final_closed_assemblies/K229new-30C6.fasta | samtools sort -o alignments/ref_pOXA-48_K8_map_K229flye.bam
samtools index alignments/ref_pOXA-48_K8_map_K229flye.bam

minimap2 -ax asm5 ../../Closed_sequences/plasmids/pOXA-48_K8.fasta assemblies_unicycler_hybrid/K229new-30C6/assembly.fasta | samtools sort -o alignments/ref_pOXA-48_K8_map_K229uni.bam
samtools index alignments/ref_pOXA-48_K8_map_K229uni.bam

minimap2 -ax asm5 ../../Closed_sequences/plasmids/pOXA-48_K8.fasta final_closed_assemblies/K163WT-28F1.fasta | samtools sort -o alignments/ref_pOXA-48_K8_map_K163uni.bam
samtools index alignments/ref_pOXA-48_K8_map_K163uni.bam

minimap2 -ax asm5 ../../Closed_sequences/plasmids/pOXA-48_K8.fasta final_closed_assemblies/K165WT.fasta | samtools sort -o alignments/ref_pOXA-48_K8_map_K165uni.bam
samtools index alignments/ref_pOXA-48_K8_map_K165uni.bam
```

Unicycler could not close all contigs in strains C288 and C289. Thus, the Flye assemblies were selected. The assembly of pOXA-48_C288 contained two small deletions and a small insertion not supported by the read mappings to  pOXA-48_K8. Snippy v4.6.0 was used to confirm that these were assembly errors not detected by Pilon, and the .consensus.fa generated was used as the final assembly.

```sh
snippy --report --outdir variants_snippy_ref_C288_map_C288/ --ref final_closed_assemblies/C288WT.fasta --R1 reads_trimmed/MiGS_C288WT_val_1.fq.gz --R2 reads_trimmed/MiGS_C288WT_val_2.fq.gz
cp variants_snippy_ref_C288_map_C288/snps.consensus.fa final_closed_assemblies/C288WT.fasta
```

All contigs in strains K163 and K165 were circularized by Unicycler and thus were selected as the final assemblies. However, the sequence of pOXA-48_K163 had a mutation in the *ltrA* gene (position 15,736 in the reference pOXA-48_K8) that did not appear in the mappings of short and long reads. This is possibly due to other versions of the gene in the chromosome or other plasmids of this strain and could not be identified by Snippy, and therefore was corrected manually to the reference allele. The sequence of pOXA-48_K165 also displayed this mutation and an additional small deletion in the *traC* gene. These mutations are not present in a previous assembly of this plasmid, generated by hybrid assembly with PacBio reads (accession number MT989346). Additionally, the *traC* deletion has been confirmed to be false by Sanger sequencing and is suspected to come from the short read assembly. Therefore, the sequence of this pOXA-48_K165 was replaced by the sequence of MT989346. The hybrid and long read assemblies of K165 lacked a 3.1 Kb contig containing a ColRNAI replicon that was present in the assemblies of K163. However, this contig was assembled with only the short reads of K165, suggesting that Unicycler discarded the contig for not being supported by long reads and that the replicon was lost during the extraction or sequencing protocol for Nanopore. The Illumina reads of K165 were mapped to the assembly of K163 to detect possible mutations in this replicon, but none were detected. Thus, the sequence of the ColRNAI replicon of K163 was added to the final assembly of K165.

```sh
snippy --report --outdir variants_snippy_ref_K163_map_K165/ --ref final_closed_assemblies/K163WT-28F1.fasta --R1 reads_trimmed/MiGS_K165_val_1.fq.gz --R2 reads_trimmed/MiGS_K165_val_2.fq.gz
```

Both assemblies of strain K153 yielded closed contigs, but the Unicycler assembly was discarded since it included the deletion in the *traC* gene (mappings of the sequencing reads against pOXA-48_K8 rejected this deletion). The Flye assembly of pOXA-48_K153 included a one nucleotide insertion in *ltrA* that was not supported by the read alignments and was therefore manually corrected in the fasta file. For strain K229, the Flye assemblies were also selected since Unicycler failed to circularize one contig. However, the Unicycler assembly of pOXA-48_K229 delimited more accurately the positions of the large deletion (in the long read mappings, most reads support the deletion between the positions 33,165-  45,513 bp of the reference pOXA-48_K8), so the pOXA-48_K229 plasmid from the Flye assembly was replaced by the plasmid from the Unicycler assembly. A 777 bp fragment affecting the IS1 upstream *bla*<sub>OXA-48</sub> (positions 9,382-10,160 bp) was missing in the Unicycler and Flye assemblies of both strains. The short and long read mappings supported this deletion in K229, but not in K153. However, the deletion was confirmed for both strains by PCR.

The final closed assemblies used in this study were annotated locally with PGAP v2021-07-01.build5508. The genomes submitted to GenBank were annotated by the NCBI with a newer version of PGAP.


## Variant calling and plasmid and resistance gene content

Breseq was used to identify SNPs and structural variants in the genomes of isolates from the same patient. To discard false positive variant calls due to missassemblies, for each pair of closed strains from the same patient, breseq was run in both directions. This way, variants present in only one of the directions were discarded. Breseq was also run mapping the reads of a strain to its own assembly to discard false positive variants. Since strains K164, K166, K165-1, K165-3, K165-4, K165-5, K165-6 and K165-7 (patient JWC) and K151 and K152 (patient WDV) are not closed, the trimmed reads (BioProjects PRJNA626430 and PRJNA838107) were only mapped to the closed strains from their respective patients. SNPs and structural variants matching previously known false positives in position or equivalent region, respectively, were discarded. Predicted mutations present in only one of the mappings against a reference that should also show against the other reference were also discarded. In some cases, multiple SNPs were identified in the same region as known false positives. These variants were also discarded after observing poor read alignments in these regions.

```sh
# Breseq mapping reads back to their respective closed genomes

breseq -r final_closed_assemblies/C288WT.gbk reads_trimmed/MiGS_C288WT_val_1.fq.gz reads_trimmed/MiGS_C288WT_val_2.fq.gz -o variants_breseq_closed/ref_C288_map_C288
breseq -r final_closed_assemblies/C289WT-31C7.gbk reads_trimmed/MiGS_C289WT_val_1.fq.gz reads_trimmed/MiGS_C289WT_val_2.fq.gz -o variants_breseq_closed/ref_C289_map_C289
breseq -r final_closed_assemblies/K153WT-30A7.gbk reads_trimmed/MiGS_K153_val_1.fq.gz reads_trimmed/MiGS_K153_val_2.fq.gz -o variants_breseq_closed/ref_K153_map_K153
breseq -r final_closed_assemblies/K229new-30C6.gbk reads_trimmed/MiGS_K229_val_1.fq.gz reads_trimmed/MiGS_K229_val_2.fq.gz -o variants_breseq_closed/ref_K229_map_K229
breseq -r final_closed_assemblies/K163WT-28F1.gbk reads_trimmed/MiGS_K163_val_1.fq.gz reads_trimmed/MiGS_K163_val_2.fq.gz -o variants_breseq_closed/ref_K163_map_K163
breseq -r final_closed_assemblies/K165WT.gbk reads_trimmed/MiGS_K165_val_1.fq.gz reads_trimmed/MiGS_K165_val_2.fq.gz -o variants_breseq_closed/ref_K165_map_K165

# Breseq of within-patient cases to the reference closed genomes

breseq -r final_closed_assemblies/C288WT.gbk reads_trimmed/MiGS_C289WT_val_1.fq.gz reads_trimmed/MiGS_C289WT_val_2.fq.gz -o variants_breseq_closed/ref_C288_map_C289
breseq -r final_closed_assemblies/C289WT-31C7.gbk reads_trimmed/MiGS_C288WT_val_1.fq.gz reads_trimmed/MiGS_C288WT_val_2.fq.gz -o variants_breseq_closed/ref_C289_map_C288

breseq -r final_closed_assemblies/K153WT-30A7.gbk reads_trimmed/MiGS_K229_val_1.fq.gz reads_trimmed/MiGS_K229_val_2.fq.gz -o variants_breseq_closed/ref_K153_map_K229
breseq -r final_closed_assemblies/K229new-30C6.gbk reads_trimmed/MiGS_K153_val_1.fq.gz reads_trimmed/MiGS_K153_val_2.fq.gz -o variants_breseq_closed/ref_K229_map_K153
breseq -r final_closed_assemblies/K153WT-30A7.gbk ../../Reads_WT_pOXA-48/trimmed/WTCHG_K151_val_1.fq.gz ../../Reads_WT_pOXA-48/trimmed/WTCHG_K151_val_2.fq.gz -o variants_breseq_closed/ref_K153_map_K151
breseq -r final_closed_assemblies/K153WT-30A7.gbk ../../Reads_WT_pOXA-48/trimmed/WTCHG_K152_val_1.fq.gz ../../Reads_WT_pOXA-48/trimmed/WTCHG_K152_val_2.fq.gz -o variants_breseq_closed/ref_K153_map_K152
breseq -r final_closed_assemblies/K229new-30C6.gbk ../../Reads_WT_pOXA-48/trimmed/WTCHG_K151_val_1.fq.gz ../../Reads_WT_pOXA-48/trimmed/WTCHG_K151_val_2.fq.gz -o variants_breseq_closed/ref_K229_map_K151
breseq -r final_closed_assemblies/K229new-30C6.gbk ../../Reads_WT_pOXA-48/trimmed/WTCHG_K152_val_1.fq.gz ../../Reads_WT_pOXA-48/trimmed/WTCHG_K152_val_2.fq.gz -o variants_breseq_closed/ref_K229_map_K152

breseq -r final_closed_assemblies/K163WT-28F1.gbk reads_trimmed/MiGS_K165_val_1.fq.gz reads_trimmed/MiGS_K165_val_2.fq.gz -o variants_breseq_closed/ref_K163_map_K165
breseq -r final_closed_assemblies/K165WT.gbk reads_trimmed/MiGS_K163_val_1.fq.gz reads_trimmed/MiGS_K163_val_2.fq.gz -o variants_breseq_closed/ref_K165_map_K163
breseq -r final_closed_assemblies/K163WT-28F1.gbk ../../Reads_WT_pOXA-48/trimmed/WTCHG_K164_val_1.fq.gz ../../Reads_WT_pOXA-48/trimmed/WTCHG_K164_val_2.fq.gz -o variants_breseq_closed/ref_K163_map_K164
breseq -r final_closed_assemblies/K163WT-28F1.gbk ../../Reads_WT_pOXA-48/trimmed/WTCHG_K166_val_1.fq.gz ../../Reads_WT_pOXA-48/trimmed/WTCHG_K166_val_2.fq.gz -o variants_breseq_closed/ref_K163_map_K166
breseq -r final_closed_assemblies/K165WT.gbk ../../Reads_WT_pOXA-48/trimmed/WTCHG_K164_val_1.fq.gz ../../Reads_WT_pOXA-48/trimmed/WTCHG_K164_val_2.fq.gz -o variants_breseq_closed/ref_K165_map_K164
breseq -r final_closed_assemblies/K165WT.gbk ../../Reads_WT_pOXA-48/trimmed/WTCHG_K166_val_1.fq.gz ../../Reads_WT_pOXA-48/trimmed/WTCHG_K166_val_2.fq.gz -o variants_breseq_closed/ref_K165_map_K166
breseq -r final_closed_assemblies/K163WT-28F1.gbk reads_trimmed/MiGS_71_val_1.fq.gz reads_trimmed/MiGS_71_val_2.fq.gz -o variants_breseq_closed/ref_K163_map_71
breseq -r final_closed_assemblies/K163WT-28F1.gbk reads_trimmed/MiGS_72_val_1.fq.gz reads_trimmed/MiGS_72_val_2.fq.gz -o variants_breseq_closed/ref_K163_map_72
breseq -r final_closed_assemblies/K163WT-28F1.gbk reads_trimmed/MiGS_73_val_1.fq.gz reads_trimmed/MiGS_73_val_2.fq.gz -o variants_breseq_closed/ref_K163_map_73
breseq -r final_closed_assemblies/K163WT-28F1.gbk reads_trimmed/MiGS_74_val_1.fq.gz reads_trimmed/MiGS_74_val_2.fq.gz -o variants_breseq_closed/ref_K163_map_74
breseq -r final_closed_assemblies/K163WT-28F1.gbk reads_trimmed/MiGS_75_val_1.fq.gz reads_trimmed/MiGS_75_val_2.fq.gz -o variants_breseq_closed/ref_K163_map_75
breseq -r final_closed_assemblies/K163WT-28F1.gbk reads_trimmed/MiGS_76_val_1.fq.gz reads_trimmed/MiGS_76_val_2.fq.gz -o variants_breseq_closed/ref_K163_map_76
breseq -r final_closed_assemblies/K165WT.gbk reads_trimmed/MiGS_71_val_1.fq.gz reads_trimmed/MiGS_71_val_2.fq.gz -o variants_breseq_closed/ref_K165_map_71
breseq -r final_closed_assemblies/K165WT.gbk reads_trimmed/MiGS_72_val_1.fq.gz reads_trimmed/MiGS_72_val_2.fq.gz -o variants_breseq_closed/ref_K165_map_72
breseq -r final_closed_assemblies/K165WT.gbk reads_trimmed/MiGS_73_val_1.fq.gz reads_trimmed/MiGS_73_val_2.fq.gz -o variants_breseq_closed/ref_K165_map_73
breseq -r final_closed_assemblies/K165WT.gbk reads_trimmed/MiGS_74_val_1.fq.gz reads_trimmed/MiGS_74_val_2.fq.gz -o variants_breseq_closed/ref_K165_map_74
breseq -r final_closed_assemblies/K165WT.gbk reads_trimmed/MiGS_75_val_1.fq.gz reads_trimmed/MiGS_75_val_2.fq.gz -o variants_breseq_closed/ref_K165_map_75
breseq -r final_closed_assemblies/K165WT.gbk reads_trimmed/MiGS_76_val_1.fq.gz reads_trimmed/MiGS_76_val_2.fq.gz -o variants_breseq_closed/ref_K165_map_76
```
The strains cured from pOXA-48 were sequenced to control for isogenicity. Variants were called with breseq, mapping the trimmed Illumina reads to their respective reference. The same mutation filtering criteria was applied.

```sh
breseq -r final_closed_assemblies/C288WT.gbk reads_trimmed/MiGS_C288c2_Clon_2_val_1.fq.gz reads_trimmed/MiGS_C288c2_Clon_2_val_2.fq.gz -o variants_breseq_closed/ref_C288_map_C288c2
breseq -r final_closed_assemblies/C289WT-31C7.gbk reads_trimmed/MiGS_C289c1_Clon_1_val_1.fq.gz reads_trimmed/MiGS_C289c1_Clon_1_val_2.fq.gz -o variants_breseq_closed/ref_C289_map_C289c1
breseq -r final_closed_assemblies/K229new-30C6.gbk reads_trimmed/MiGS_K229c5_val_1.fq.gz reads_trimmed/MiGS_K229c5_val_2.fq.gz -o variants_breseq_closed/ref_K229_map_K229c5
breseq -r final_closed_assemblies/K165WT.gbk reads_trimmed/MiGS_K165c5_val_1.fq.gz reads_trimmed/MiGS_K165c5_val_2.fq.gz -o variants_breseq_closed/ref_K165_map_K165c5
```

ABRicate v1.0.1 with the plasmidfinder database was used to confirm no other plasmid was eliminated from the cured strains.

```sh
for fq1 in reads_trimmed/*val_1.fq.gz
do
	fq2=${fq1%%1.fq.gz}"2.fq.gz"
	strain=${fq1::-12}
	strain=${strain:14}
	# De novo assembly with SPAdes
	spades.py --isolate --cov-cutoff auto -o assemblies_SPAdes/$strain -1 $fq1 -2 $fq2
done


mkdir plasmids_abricate
mkdir resistance_abricate

for contigs in assemblies_SPAdes/*/contigs.fasta
do
	strain=${contigs:18}
	strain=${strain::-14}

	abricate --db plasmidfinder $contigs > plasmids_abricate/$strain".tsv"
	abricate --db resfinder $contigs > resistance_abricate/$strain".tsv"
done

for contigs in final_closed_assemblies/*fasta
do
	strain=${contigs:24}
	strain=${strain::-6}

	abricate --db plasmidfinder $contigs > plasmids_abricate/$strain".tsv"
	abricate --db resfinder $contigs > resistance_abricate/$strain".tsv"
done

abricate --summary plasmids_abricate/* > plasmids_abricate/summary_plasmids.tsv
abricate --summary resistance_abricate/* > resistance_abricate/summary_resistance.tsv
```


# Construction of phylogenetic trees

Snippy v4.6.0 was used to find SNPs between all *E. coli* isolates (12 ST10 + 2 ST744 isolates, using C288 as reference), *K. pneumoniae* ST11 (87 isolates, with K153 as reference) and *K. pneumoniae* ST307 (20 isolates, with K163 as reference) from the R-GNOSIS collection. Strain K25 (ST11) was removed from the analysis because the fastq files were truncated.

```sh
# E. coli strains

snippy --ref final_closed_assemblies/C288WT.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_C288_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_C288_val_2.fq.gz --outdir phylogeny_core/snippy_Ec/ref_C288_map_C288
snippy --ref final_closed_assemblies/C288WT.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_C289_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_C289_val_2.fq.gz --outdir phylogeny_core/snippy_Ec/ref_C288_map_C289
snippy --ref final_closed_assemblies/C288WT.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_C165_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_C165_val_2.fq.gz --outdir phylogeny_core/snippy_Ec/ref_C288_map_C165
snippy --ref final_closed_assemblies/C288WT.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_C166_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_C166_val_2.fq.gz --outdir phylogeny_core/snippy_Ec/ref_C288_map_C166
snippy --ref final_closed_assemblies/C288WT.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_C642_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_C642_val_2.fq.gz --outdir phylogeny_core/snippy_Ec/ref_C288_map_C642
snippy --ref final_closed_assemblies/C288WT.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_C643_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_C643_val_2.fq.gz --outdir phylogeny_core/snippy_Ec/ref_C288_map_C643
snippy --ref final_closed_assemblies/C288WT.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_C646_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_C646_val_2.fq.gz --outdir phylogeny_core/snippy_Ec/ref_C288_map_C646
snippy --ref final_closed_assemblies/C288WT.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_C662_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_C662_val_2.fq.gz --outdir phylogeny_core/snippy_Ec/ref_C288_map_C662
snippy --ref final_closed_assemblies/C288WT.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_C667_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_C667_val_2.fq.gz --outdir phylogeny_core/snippy_Ec/ref_C288_map_C667
snippy --ref final_closed_assemblies/C288WT.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_C717_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_C717_val_2.fq.gz --outdir phylogeny_core/snippy_Ec/ref_C288_map_C717
snippy --ref final_closed_assemblies/C288WT.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_C718_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_C718_val_2.fq.gz --outdir phylogeny_core/snippy_Ec/ref_C288_map_C718
snippy --ref final_closed_assemblies/C288WT.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_C728_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_C728_val_2.fq.gz --outdir phylogeny_core/snippy_Ec/ref_C288_map_C728
snippy --ref final_closed_assemblies/C288WT.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_C752_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_C752_val_2.fq.gz --outdir phylogeny_core/snippy_Ec/ref_C288_map_C752
snippy --ref final_closed_assemblies/C288WT.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_N46_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_N46_val_2.fq.gz --outdir phylogeny_core/snippy_Ec/ref_C288_map_N46

# K. pneumoniae ST11

snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K153_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K153_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K153
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K90_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K90_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K90
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K109_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K109_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K109
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K24_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K24_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K24
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K31_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K31_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K31
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K47_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K47_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K47
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K51_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K51_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K51
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K52_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K52_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K52
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K8_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K8_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K8
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K88_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K88_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K88
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K89_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K89_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K89
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K93_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K93_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K93
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K1_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K1_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K1
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K116_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K116_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K116
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K127_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K127_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K127
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K147_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K147_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K147
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K15_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K15_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K15
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K151_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K151_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K151
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K152_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K152_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K152
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K160_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K160_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K160
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K161_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K161_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K161
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K162_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K162_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K162
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K197_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K197_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K197
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K2_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K2_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K2
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K202_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K202_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K202
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K229_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K229_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K229
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K236-1_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K236-1_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K236-1
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K237_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K237_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K237
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K238_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K238_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K238
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K25_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K25_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K25
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K259_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K259_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K259
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K26_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K26_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K26
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K260_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K260_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K260
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K263_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K263_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K263
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K264_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K264_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K264
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K265_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K265_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K265
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K266_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K266_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K266
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K267_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K267_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K267
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K268_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K268_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K268
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K269_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K269_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K269
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K270_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K270_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K270
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K271_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K271_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K271
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K28_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K28_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K28
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K280_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K280_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K280
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K281_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K281_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K281
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K29_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K29_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K29
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K293_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K293_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K293
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K294_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K294_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K294
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K295_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K295_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K295
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K298_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K298_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K298
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K308_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K308_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K308
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K312_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K312_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K312
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K313_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K313_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K313
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K317_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K317_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K317
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K319_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K319_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K319
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K323_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K323_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K323
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K324_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K324_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K324
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K326_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K326_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K326
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K36_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K36_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K36
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K39_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K39_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K39
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K4_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K4_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K4
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K57_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K57_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K57
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_R25_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_R25_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_R25
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_X6_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_X6_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_X6
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_H73_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_H73_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_H73
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K70_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K70_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K70
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K76_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K76_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K76
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_G21_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_G21_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_G21
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_L27_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_L27_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_L27
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K101_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K101_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K101
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_L44_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_L44_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_L44
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_L74_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_L74_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_L74
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_Ñ22_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_Ñ22_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_Ñ22
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_N23_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_N23_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_N23
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_Ñ34_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_Ñ34_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_Ñ34
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_N45_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_N45_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_N45
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_P76_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_P76_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_P76
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_L45_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_L45_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_L45
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K34_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K34_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K34
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_E73_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_E73_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_E73
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_E74_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_E74_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_E74
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_E76_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_E76_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_E76
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_F14_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_F14_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_F14
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K62_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K62_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K62
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_F48_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_F48_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_F48
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K74_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K74_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K74
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_F73_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_F73_val_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_F73

# K. pneumoniae ST307

snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K163_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K163_val_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_K163
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K164_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K164_val_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_K164
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K165_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K165_val_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_K165
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K166_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K166_val_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_K166
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_R50_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_R50_val_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_R50
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K172_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K172_val_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_K172
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K173_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K173_val_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_K173
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K174_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K174_val_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_K174
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K198_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K198_val_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_K198
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K199_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K199_val_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_K199
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K235_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K235_val_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_K235
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K306_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K306_val_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_K306
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K314_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K314_val_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_K314
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K315_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K315_val_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_K315
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K318_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_K318_val_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_K318
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_J61_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_J61_val_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_J61
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_J66_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_J66_val_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_J66
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_N22_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_N22_val_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_N22
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../../Reads_WT_pOXA-48/trimmed/WTCHG_N36_val_1.fq.gz --R2 ../../Reads_WT_pOXA-48/trimmed/WTCHG_N36_val_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_N36
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 reads_trimmed/MiGS_74_val_val_1.fq.gz --R2 reads_trimmed/MiGS_74_val_val_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_74
```

Next, snippy-core was used to find the core genome of each group of alignments. To obtain the core genome for ST11, strain K78 had to be removed for diverging too much from the rest of strains. Recombinant regions were eliminated with Gubbins v3.1.4 and SNPs were extracted from the alignment with snp-sites v2.5.1. Maximum-likelihood trees were constructed with IQ-TREE v1.6.12 from the extracted alignments with the function of automated detection of the best evolutionary model and an ultrafast bootstrap of 1000 optimized by hill-climbing nearest neighbor interchange (NNI) on the corresponding bootstrap alignment.

```sh
cd phylogeny_core
echo ">Reference" > rmid.txt

# E. coli
snippy-core --ref ../final_closed_assemblies/C288WT.gbk --prefix core_Ec snippy_Ec/*
awk '(NR==FNR) { toRemove[$1]; next }/^>/ { p=1; for(h in toRemove) if ( h ~ $0) p=0 }p' rmid.txt core_Ec.full.aln > core_Ec.full_noref.aln
snippy-clean_full_aln core_Ec.full_noref.aln > clean_Ec.full.aln
run_gubbins.py -p gubbins_Ec clean_Ec.full.aln
snp-sites -c gubbins_Ec.filtered_polymorphic_sites.fasta > clean_Ec.core.aln
iqtree -s clean_Ec.core.aln -m MFP -bb 1000 -bnni -pre iqtree_Ec

# K. pneumoniae ST11
snippy-core --ref ../final_closed_assemblies/K153WT-30A7.gbk --prefix core_ST11 snippy_ST11/*
awk '(NR==FNR) { toRemove[$1]; next }/^>/ { p=1; for(h in toRemove) if ( h ~ $0) p=0 }p' rmid.txt core_ST11.full.aln > core_ST11.full_noref.aln
snippy-clean_full_aln core_ST11.full_noref.aln > clean_ST11.full.aln
run_gubbins.py -p gubbins_ST11 clean_ST11.full.aln
snp-sites -c gubbins_ST11.filtered_polymorphic_sites.fasta > clean_ST11.core.aln
iqtree -s clean_ST11.core.aln -m MFP -bb 1000 -bnni -pre iqtree_ST11

# K. pneumoniae ST307
snippy-core --ref ../final_closed_assemblies/K163WT-28F1.gbk --prefix core_ST307 snippy_ST307/*
awk '(NR==FNR) { toRemove[$1]; next }/^>/ { p=1; for(h in toRemove) if ( h ~ $0) p=0 }p' rmid.txt core_ST307.full.aln > core_ST307.full_noref.aln
snippy-clean_full_aln core_ST307.full_noref.aln > clean_ST307.full.aln
run_gubbins.py -p gubbins_ST307 clean_ST307.full.aln
snp-sites -c gubbins_ST307.filtered_polymorphic_sites.fasta > clean_ST307.core.aln
iqtree -s clean_ST307.core.aln -m MFP -bb 1000 -bnni -pre iqtree_ST307
```

Trees were visualized and edited in iTOL (https://itol.embl.de/) and inkscape v0.17.
