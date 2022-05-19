# Initial setup

The following directories were created in the current working directory to store sequencing and analysis data.

```sh
# For section: 'De novo assembling and genomic analysis of E. coli J53 carrying different PVs'
mkdir -p plasmid_variants_J53/reads_trimmed/reads_raw
# For sections: 'Assembly and analysis of pOXA-48 variants in cases of within-patient evolution' and 'Construction of phylogenetic trees'
mkdir -p within_patient_evol_cured/reads_trimmed/reads_raw
mkdir within_patient_evol_cured/reads_nanopore
mkdir -p Reads_pOXA-48_carriers/reads_raw
# For storing the GenBank file pOXA-48_K8.gb ()
mkdir -p closed_sequences/plasmids
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
| SRR19213039 | SAMN28404986 | C288 | C288WT_nanopore | within_patient_evol_cured/reads_nanopore |
| SRR19213028 | SAMN28404987 | C289 | C289WT_S142 | within_patient_evol_cured/reads_trimmed/reads_raw |
| SRR19213017 | SAMN28404987 | C289 | C289WT-31C7_nanopore | within_patient_evol_cured/reads_nanopore |
| SRR19213010 | SAMN28404988 | C289cured | C289c1_Clon_1_S143 | within_patient_evol_cured/reads_trimmed/reads_raw |
| SRR19213009 | SAMN28404989 | K153 | K153WT_S103 | within_patient_evol_cured/reads_trimmed/reads_raw |
| SRR19213008 | SAMN28404989 | K153 | K153WT-30A7_nanopore | within_patient_evol_cured/reads_nanopore |
| SRR19213007 | SAMN28404990 | K229 | K229newWT_S102 | within_patient_evol_cured/reads_trimmed/reads_raw |
| SRR19213006 | SAMN28404990 | K229 | K229new-30C6_nanopore | within_patient_evol_cured/reads_nanopore |
| SRR19213005 | SAMN28404991 | K229cured | K229C5_S89 | within_patient_evol_cured/reads_trimmed/reads_raw |
| SRR19213038 | SAMN28404992 | K163 | K163WT_S99 | within_patient_evol_cured/reads_trimmed/reads_raw |
| SRR19213037 | SAMN28404992 | K163 | K163WT-28F1_nanopore | within_patient_evol_cured/reads_nanopore |
| SRR19213036 | SAMN28404993 | K165-2 | K165WT_S96 | within_patient_evol_cured/reads_trimmed/reads_raw |
| SRR19213035 | SAMN28404993 | K165-2 | K165WT_nanopore | within_patient_evol_cured/reads_nanopore |
| SRR19213034 | SAMN28404994 | K165cured | K165c5_S97 | within_patient_evol_cured/reads_trimmed/reads_raw |
| SRR19213033 | SAMN28404995 | K165-3 | 71_S107 | within_patient_evol_cured/reads_trimmed/reads_raw |
| SRR19213032 | SAMN28404996 | K165-4 | 72_S108 | within_patient_evol_cured/reads_trimmed/reads_raw |
| SRR19213031 | SAMN28404997 | K165-5 | 73_S109 | within_patient_evol_cured/reads_trimmed/reads_raw |
| SRR19213030 | SAMN28404998 | K165-1 | 74_S110 | within_patient_evol_cured/reads_trimmed/reads_raw |
| SRR19213029 | SAMN28404999 | K165-6 | 75_S111 | within_patient_evol_cured/reads_trimmed/reads_raw |
| SRR19213027 | SAMN28405000 | K165-7 | 76_S112 | within_patient_evol_cured/reads_trimmed/reads_raw |

Raw Illumina reads from BioProject PRJNA626430 used for constructing phylogenetic trees of *E. coli* and *K. pneumoniae*  are stored in `Reads_pOXA-48_carriers/reads_raw`. The GenBank annotation file for plasmid pOXA-48 (MT441554) is stored in `closed_sequences/plasmids` and the `ltrA.fasta` file is placed in the current directory.

# Illumina and Nanopore read preprocessing

Trim Galore v0.6.4 was used to trim Illumina sequences and generate FastQC reports. All read datasets reached good quality.

```sh
cd plasmid_variants_J53/reads_trimmed/reads_raw/

for fq1 in *R1_001.fastq.gz
do
	fq2=${fq1%%R1_001.fastq.gz}"R2_001.fastq.gz"
	name=${fq1/_[A-Z][0-9]*_R1_001.fastq.gz/}
	trim_galore --quality 20 --nextera --length 50 --fastqc --basename $name --output_dir .. --paired $fq1 $fq2
done

cd ../../../within_patient_evol_cured/reads_trimmed/reads_raw/

for fq1 in *R1_001.fastq.gz
do
	fq2=${fq1%%R1_001.fastq.gz}"R2_001.fastq.gz"
	name=${fq1/_[A-Z][0-9]*_R1_001.fastq.gz/}
	trim_galore --quality 20 --nextera --length 50 --fastqc --basename $name --output_dir .. --paired $fq1 $fq2
done
```

Nanopore reads were filtered with filtlong v0.2.1 with options as suggested by the authors to obtain a subset of high identity reads with a minimum read depth of 85x.

```sh
cd ../.. # From within_patient_evol_cured

filtlong --min_length 1000 --keep_percent 90 --target_bases 500000000 reads_nanopore/C288WT_nanopore.fastq.gz | gzip > reads_nanopore/C288WT_nanopore_filt.fastq.gz
filtlong --min_length 1000 --keep_percent 90 --target_bases 500000000 reads_nanopore/C289WT-31C7_nanopore.fastq.gz | gzip > reads_nanopore/C289WT-31C7_nanopore_filt.fastq.gz
filtlong --min_length 1000 --keep_percent 90 --target_bases 500000000 reads_nanopore/K153WT-30A7_nanopore.fastq.gz | gzip > reads_nanopore/K153WT-30A7_nanopore_filt.fastq.gz
filtlong --min_length 1000 --keep_percent 90 --target_bases 500000000 reads_nanopore/K229new-30C6_nanopore.fastq.gz | gzip > reads_nanopore/K229new-30C6_nanopore_filt.fastq.gz
filtlong --min_length 1000 --keep_percent 90 --target_bases 500000000 reads_nanopore/K163WT-28F1_nanopore.fastq.gz | gzip > reads_nanopore/K163WT-28F1_nanopore_filt.fastq.gz
filtlong --min_length 1000 --keep_percent 90 --target_bases 500000000 reads_nanopore/K165WT_nanopore.fastq.gz | gzip > reads_nanopore/K165WT_nanopore_filt.fastq.gz
```

Reads from BioProject PRJNA626430 stored in `Reads_pOXA-48_carriers/reads_raw` were also trimmed with Trim Galore (--quality 20 --illumina --length 50 --fastqc). Trimmed reads were stored in `Reads_pOXA-48_carriers` and were renamed to the corresponding strain names.


# _De novo_ assembling and genomic analysis of *E. coli* J53 carrying different PVs

The *E. coli* J53 strains (carrying different PVs and plasmid-free) used in this study were sequenced to control isogenic conditions. First, the genome of J53 was assembled using SPAdes v3.15.2 and annotated with Prokka v1.14.6.

```sh
cd ../plasmid_variants_J53

spades.py --isolate --cov-cutoff auto -o assemblies_SPAdes/J53javier -1 reads_trimmed/J53javier_val_1.fq.gz -2 reads_trimmed/J53javier_val_2.fq.gz
quast.py -o assemblies_SPAdes/quast_reports/J53javier --glimmer assemblies_SPAdes/J53javier
prokka --outdir annotation_prokka_J53 --prefix J53 assemblies_SPAdes/J53javier/contigs.fasta
prokka --outdir annotation_prokka_J53_compliant --prefix J53 --compliant assemblies_SPAdes/J53javier/contigs.fasta
```

The *ltrA* gene of pOXA-48_K8 was blasted (BLASTn v2.11.0) against J53 to confirm its absence.

```sh
mkdir blastn_ltrA
makeblastdb -in assemblies_SPAdes/J53javier/contigs.fasta -dbtype nucl
blastn -query ltrA.fasta -db assemblies_SPAdes/J53javier/contigs.fasta -outfmt 6 > blastn_ltrA/J53javier.tsv
```

Variant calling was performed to detect mutations in the chromosomes and pOXA-48 PVs. Chromosomal mutations were analyzed following this workflow. First, Snippy v4.6.0 was used to identify variants in the draft genome of J53 by mapping the Illumina reads back to its assembly.

```sh
snippy --report --outdir variants_snippy_ref_J53_map_J53 --ref annotation_prokka_J53_compliant/J53.gbk --R1 reads_trimmed/J53javier_val_1.fq.gz --R2 reads_trimmed/J53javier_val_2.fq.gz
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
		snippy --report --outdir variants_snippy_ref_pOXA-48_K8_map_TC/ref_pOXA-48_K8_map_$strain --ref ../closed_sequences/plasmids/pOXA-48_K8.gb --R1 $fq1 --R2 $fq2

		# Breseq using as reference J53 and mapping short reads of J53 carrying different PVs
		breseq -o variants_breseq_ref_J53_map_TC/ref_J53_map_$strain -c annotation_prokka_J53/J53.gff $fq1 $fq2

		# Breseq using as reference pOXA-48_K8 (accession number MT441554) and mapping short reads of J53 carrying different PVs
		breseq -o variants_breseq_ref_pOXA-48_K8_map_TC/ref_pOXA-48_K8_map_$strain -r ../closed_sequences/plasmids/pOXA-48_K8.gb $fq1 $fq2
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
	snippy --report --outdir variants_snippy_ref_TC_map_J53/ref_$strain"_map_J53" --ref $contigs --R1 reads_trimmed/J53javier_val_1.fq.gz --R2 reads_trimmed/J53javier_val_2.fq.gz

	# Blasting the ltrA gene against J53 carrying different PVs
	makeblastdb -in $contigs -dbtype nucl
	blastn -query ltrA.fasta -db $contigs -outfmt 6 > blastn_ltrA/$strain".tsv"

	# Detecting plasmid replicons and resistance genes with ABRicate
	abricate --db plasmidfinder $contigs > plasmids_abricate/$strain".tsv"
	abricate --db resfinder $contigs > resistance_abricate/$strain".tsv"
done

# Summary of ABRicate results
abricate --summary plasmids_abricate/* > plasmids_abricate/summary_plasmids.tsv
abricate --summary resistance_abricate/* > resistance_abricate/summary_resistance.tsv
```


# Assembly and analysis of pOXA-48 variants in cases of within-patient evolution

## Generating closed genomes

Hybrid assemblies of the six strains (C288, C289, K153, K229, K163 and K165) were obtained using Unicycler v0.4.9 with default parameters.

```sh
cd ../within_patient_evol_cured

unicycler -t 18 --spades_path ~/Software/SPAdes-3.13.0-Linux/bin/spades.py -1 reads_trimmed/C288WT_val_1.fq.gz -2 reads_trimmed/C288WT_val_2.fq.gz -l reads_nanopore/C288WT_nanopore_filt.fastq.gz -o assemblies_unicycler_hybrid/C288WT
unicycler -t 18 --spades_path ~/Software/SPAdes-3.13.0-Linux/bin/spades.py -1 reads_trimmed/C289WT_val_1.fq.gz -2 reads_trimmed/C289WT_val_2.fq.gz -l reads_nanopore/C289WT-31C7_nanopore_filt.fastq.gz -o assemblies_unicycler_hybrid/C289WT-31C7
unicycler -t 18 --spades_path ~/Software/SPAdes-3.13.0-Linux/bin/spades.py -1 reads_trimmed/K153_val_1.fq.gz -2 reads_trimmed/K153_val_2.fq.gz -l reads_nanopore/K153WT-30A7_nanopore_filt.fastq.gz -o assemblies_unicycler_hybrid/K153WT-30A7
unicycler -t 18 --spades_path ~/Software/SPAdes-3.13.0-Linux/bin/spades.py -1 reads_trimmed/K229_val_1.fq.gz -2 reads_trimmed/K229_val_2.fq.gz -l reads_nanopore/K229new-30C6_nanopore_filt.fastq.gz -o assemblies_unicycler_hybrid/K229new-30C6
unicycler -t 18 --spades_path ~/Software/SPAdes-3.13.0-Linux/bin/spades.py -1 reads_trimmed/K163_val_1.fq.gz -2 reads_trimmed/K163_val_2.fq.gz -l reads_nanopore/K163WT-28F1_nanopore_filt.fastq.gz -o assemblies_unicycler_hybrid/K163WT-28F1
unicycler -t 18 --spades_path ~/Software/SPAdes-3.13.0-Linux/bin/spades.py -1 reads_trimmed/K165_val_1.fq.gz -2 reads_trimmed/K165_val_2.fq.gz -l reads_nanopore/K165WT_nanopore_filt.fastq.gz -o assemblies_unicycler_hybrid/K165WT
```

Additionally, for the cases where Unicycler was not able to circularize contigs (C288, C289, K153 and K229), long reads were assembled with Flye v2.9 and circularization was confirmed by inspecting the assembly graphs in Bandage.

```sh
flye --nano-raw reads_nanopore/C288WT_nanopore_filt.fastq.gz --threads 18 --plasmids --out-dir assemblies_flye/C288WT
flye --nano-raw reads_nanopore/C289WT-31C7_nanopore_filt.fastq.gz --threads 18 --plasmids --out-dir assemblies_flye/C289WT-31C7
flye --nano-raw reads_nanopore/K153WT-30A7_nanopore_filt.fastq.gz --threads 18 --plasmids --out-dir assemblies_flye/K153WT-30A7
flye --nano-raw reads_nanopore/K229new-30C6_nanopore_filt.fastq.gz --threads 18 --plasmids --out-dir assemblies_flye/K229new-30C6
```

Consensus sequences of these assemblies were obtained with Medaka v1.4.3.

```sh
medaka_consensus -i reads_nanopore/C288WT_nanopore_filt.fastq.gz -d assemblies_flye/C288WT/assembly.fasta -o assemblies_flye/C288WT/medaka
medaka_consensus -i reads_nanopore/C289WT-31C7_nanopore_filt.fastq.gz -d assemblies_flye/C289WT-31C7/assembly.fasta -o assemblies_flye/C289WT-31C7/medaka
medaka_consensus -i reads_nanopore/K153WT-30A7_nanopore_filt.fastq.gz -d assemblies_flye/K153WT-30A7/assembly.fasta -o assemblies_flye/K153WT-30A7/medaka
medaka_consensus -i reads_nanopore/K229new-30C6_nanopore_filt.fastq.gz -d assemblies_flye/K229new-30C6/assembly.fasta -o assemblies_flye/K229new-30C6/medaka
```

Several rounds of Pilon v1.24 were performed mapping the trimmed Illumina reads until no further changes in sequence were observed.

```sh
# C288 round 1
bwa index assemblies_flye/C288WT/medaka/consensus.fasta
bwa mem -t 18 assemblies_flye/C288WT/medaka/consensus.fasta reads_trimmed/C288WT_val_1.fq.gz reads_trimmed/C288WT_val_2.fq.gz | samtools sort -o assemblies_flye/C288WT/medaka/alignments.bam
samtools index assemblies_flye/C288WT/medaka/alignments.bam
pilon --changes --genome assemblies_flye/C288WT/medaka/consensus.fasta --frags assemblies_flye/C288WT/medaka/alignments.bam --output C288WT --outdir assemblies_flye/C288WT/pilon/
# C288 round 2
bwa index assemblies_flye/C288WT/pilon/C288WT.fasta
bwa mem -t 18 assemblies_flye/C288WT/pilon/C288WT.fasta reads_trimmed/C288WT_val_1.fq.gz reads_trimmed/C288WT_val_2.fq.gz | samtools sort -o assemblies_flye/C288WT/medaka/alignments.bam
samtools index assemblies_flye/C288WT/medaka/alignments.bam
pilon --changes --genome assemblies_flye/C288WT/pilon/C288WT.fasta --frags assemblies_flye/C288WT/medaka/alignments.bam --output C288WT_2 --outdir assemblies_flye/C288WT/pilon/

# C289 round 1
bwa index assemblies_flye/C289WT-31C7/medaka/consensus.fasta
bwa mem -t 18 assemblies_flye/C289WT-31C7/medaka/consensus.fasta reads_trimmed/C289WT_val_1.fq.gz reads_trimmed/C289WT_val_2.fq.gz | samtools sort -o assemblies_flye/C289WT-31C7/medaka/alignments.bam
samtools index assemblies_flye/C289WT-31C7/medaka/alignments.bam
pilon --changes --genome assemblies_flye/C289WT-31C7/medaka/consensus.fasta --frags assemblies_flye/C289WT-31C7/medaka/alignments.bam --output C289WT-31C7 --outdir assemblies_flye/C289WT-31C7/pilon/
# C289 round 2
bwa index assemblies_flye/C289WT-31C7/pilon/C289WT-31C7.fasta
bwa mem -t 18 assemblies_flye/C289WT-31C7/pilon/C289WT-31C7.fasta reads_trimmed/C289WT_val_1.fq.gz reads_trimmed/C289WT_val_2.fq.gz | samtools sort -o assemblies_flye/C289WT-31C7/medaka/alignments.bam
samtools index assemblies_flye/C289WT-31C7/medaka/alignments.bam
pilon --changes --genome assemblies_flye/C289WT-31C7/pilon/C289WT-31C7.fasta --frags assemblies_flye/C289WT-31C7/medaka/alignments.bam --output C289WT-31C7_2 --outdir assemblies_flye/C289WT-31C7/pilon/

# K153 round 1
bwa index assemblies_flye/K153WT-30A7/medaka/consensus.fasta
bwa mem -t 18 assemblies_flye/K153WT-30A7/medaka/consensus.fasta reads_trimmed/K153_val_1.fq.gz reads_trimmed/K153_val_2.fq.gz | samtools sort -o assemblies_flye/K153WT-30A7/medaka/alignments.bam
samtools index assemblies_flye/K153WT-30A7/medaka/alignments.bam
pilon --changes --genome assemblies_flye/K153WT-30A7/medaka/consensus.fasta --frags assemblies_flye/K153WT-30A7/medaka/alignments.bam --output K153WT-30A7 --outdir assemblies_flye/K153WT-30A7/pilon/
# K153 round 2
bwa index assemblies_flye/K153WT-30A7/pilon/K153WT-30A7.fasta
bwa mem -t 18 assemblies_flye/K153WT-30A7/pilon/K153WT-30A7.fasta reads_trimmed/K153_val_1.fq.gz reads_trimmed/K153_val_2.fq.gz | samtools sort -o assemblies_flye/K153WT-30A7/medaka/alignments.bam
samtools index assemblies_flye/K153WT-30A7/medaka/alignments.bam
pilon --changes --genome assemblies_flye/K153WT-30A7/pilon/K153WT-30A7.fasta --frags assemblies_flye/K153WT-30A7/medaka/alignments.bam --output K153WT-30A7_2 --outdir assemblies_flye/K153WT-30A7/pilon/

# K229 round 1
bwa index assemblies_flye/K229new-30C6/medaka/consensus.fasta
bwa mem -t 18 assemblies_flye/K229new-30C6/medaka/consensus.fasta reads_trimmed/K229_val_1.fq.gz reads_trimmed/K229_val_2.fq.gz | samtools sort -o assemblies_flye/K229new-30C6/medaka/alignments.bam
samtools index assemblies_flye/K229new-30C6/medaka/alignments.bam
pilon --changes --genome assemblies_flye/K229new-30C6/medaka/consensus.fasta --frags assemblies_flye/K229new-30C6/medaka/alignments.bam --output K229new-30C6 --outdir assemblies_flye/K229new-30C6/pilon/
# K229 round 2
bwa index assemblies_flye/K229new-30C6/pilon/K229new-30C6.fasta
bwa mem -t 18 assemblies_flye/K229new-30C6/pilon/K229new-30C6.fasta reads_trimmed/K229_val_1.fq.gz reads_trimmed/K229_val_2.fq.gz | samtools sort -o assemblies_flye/K229new-30C6/medaka/alignments.bam
samtools index assemblies_flye/K229new-30C6/medaka/alignments.bam
pilon --changes --genome assemblies_flye/K229new-30C6/pilon/K229new-30C6.fasta --frags assemblies_flye/K229new-30C6/medaka/alignments.bam --output K229new-30C6_2 --outdir assemblies_flye/K229new-30C6/pilon/
# K229 round 3
bwa index assemblies_flye/K229new-30C6/pilon/K229new-30C6_2.fasta
bwa mem -t 18 assemblies_flye/K229new-30C6/pilon/K229new-30C6_2.fasta reads_trimmed/K229_val_1.fq.gz reads_trimmed/K229_val_2.fq.gz | samtools sort -o assemblies_flye/K229new-30C6/medaka/alignments.bam
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

Correctness of the assemblies of the pOXA-48 plasmids was confirmed by mapping the short and long reads with BWA-MEM v0.7.17 (14) and minimap2 v2.21 (15) (option -ax map-ont), respectively, to the reference pOXA-48_K8 sequence (MT441554), and by aligning the obtained assemblies to the pOXA-48_K8 reference with minimap2 (option -ax asm5). Alignments were visualized in IGV.

```sh
# Mapping short reads to pOXA-48_K8 (MT441554)

mkdir alignments
bwa index ../closed_sequences/plasmids/pOXA-48_K8.fasta

bwa mem -t 18 ../closed_sequences/plasmids/pOXA-48_K8.fasta reads_trimmed/C288WT_val_1.fq.gz reads_trimmed/C288WT_val_2.fq.gz | samtools sort -o alignments/ref_pOXA-48_K8_align_C288_illu.bam
samtools index alignments/ref_pOXA-48_K8_align_C288_illu.bam

bwa mem -t 18 ../closed_sequences/plasmids/pOXA-48_K8.fasta reads_trimmed/C289WT_val_1.fq.gz reads_trimmed/C289WT_val_2.fq.gz | samtools sort -o alignments/ref_pOXA-48_K8_align_C289_illu.bam
samtools index alignments/ref_pOXA-48_K8_align_C289_illu.bam

bwa mem -t 18 ../closed_sequences/plasmids/pOXA-48_K8.fasta reads_trimmed/K153_val_1.fq.gz reads_trimmed/K153_val_2.fq.gz | samtools sort -o alignments/ref_pOXA-48_K8_align_K153_illu.bam
samtools index alignments/ref_pOXA-48_K8_align_K153_illu.bam

bwa mem -t 18 ../closed_sequences/plasmids/pOXA-48_K8.fasta reads_trimmed/K229_val_1.fq.gz reads_trimmed/K229_val_2.fq.gz | samtools sort -o alignments/ref_pOXA-48_K8_align_K229_illu.bam
samtools index alignments/ref_pOXA-48_K8_align_K229_illu.bam

bwa mem -t 18 ../closed_sequences/plasmids/pOXA-48_K8.fasta reads_trimmed/K163_val_1.fq.gz reads_trimmed/K163_val_2.fq.gz | samtools sort -o alignments/ref_pOXA-48_K8_align_K163_illu.bam
samtools index alignments/ref_pOXA-48_K8_align_K163_illu.bam

bwa mem -t 18 ../closed_sequences/plasmids/pOXA-48_K8.fasta reads_trimmed/K165_val_1.fq.gz reads_trimmed/K165_val_2.fq.gz | samtools sort -o alignments/ref_pOXA-48_K8_align_K165_illu.bam
samtools index alignments/ref_pOXA-48_K8_align_K165_illu.bam

# Mapping long reads to pOXA-48_K8 (MT441554)

minimap2 -ax map-ont ../closed_sequences/plasmids/pOXA-48_K8.fasta reads_nanopore/C288WT_nanopore_filt.fastq.gz | samtools sort -o alignments/ref_pOXA-48_K8_align_C288_nano.bam
samtools index alignments/ref_pOXA-48_K8_align_C288_nano.bam

minimap2 -ax map-ont ../closed_sequences/plasmids/pOXA-48_K8.fasta reads_nanopore/C289WT-31C7_nanopore_filt.fastq.gz | samtools sort -o alignments/ref_pOXA-48_K8_align_C289_nano.bam
samtools index alignments/ref_pOXA-48_K8_align_C289_nano.bam

minimap2 -ax map-ont ../closed_sequences/plasmids/pOXA-48_K8.fasta reads_nanopore/K153WT-30A7_nanopore_filt.fastq.gz | samtools sort -o alignments/ref_pOXA-48_K8_align_K153_nano.bam
samtools index alignments/ref_pOXA-48_K8_align_K153_nano.bam

minimap2 -ax map-ont ../closed_sequences/plasmids/pOXA-48_K8.fasta reads_nanopore/K229new-30C6_nanopore_filt.fastq.gz | samtools sort -o alignments/ref_pOXA-48_K8_align_K229_nano.bam
samtools index alignments/ref_pOXA-48_K8_align_K229_nano.bam

minimap2 -ax map-ont ../closed_sequences/plasmids/pOXA-48_K8.fasta reads_nanopore/K163WT-28F1_nanopore_filt.fastq.gz | samtools sort -o alignments/ref_pOXA-48_K8_align_K163_nano.bam
samtools index alignments/ref_pOXA-48_K8_align_K163_nano.bam

minimap2 -ax map-ont ../closed_sequences/plasmids/pOXA-48_K8.fasta reads_nanopore/K165WT_nanopore_filt.fastq.gz | samtools sort -o alignments/ref_pOXA-48_K8_align_K165_nano.bam
samtools index alignments/ref_pOXA-48_K8_align_K165_nano.bam

# Aligning assemblies to pOXA-48_K8 (MT441554)

minimap2 -ax asm5 ../closed_sequences/plasmids/pOXA-48_K8.fasta final_closed_assemblies/C288WT.fasta | samtools sort -o alignments/ref_pOXA-48_K8_map_C288flye.bam
samtools index alignments/ref_pOXA-48_K8_map_C288flye.bam

minimap2 -ax asm5 ../closed_sequences/plasmids/pOXA-48_K8.fasta final_closed_assemblies/C289WT-31C7.fasta | samtools sort -o alignments/ref_pOXA-48_K8_map_C289flye.bam
samtools index alignments/ref_pOXA-48_K8_map_C289flye.bam

minimap2 -ax asm5 ../closed_sequences/plasmids/pOXA-48_K8.fasta final_closed_assemblies/K153WT-30A7.fasta | samtools sort -o alignments/ref_pOXA-48_K8_map_K153flye.bam
samtools index alignments/ref_pOXA-48_K8_map_K153flye.bam

minimap2 -ax asm5 ../closed_sequences/plasmids/pOXA-48_K8.fasta final_closed_assemblies/K229new-30C6.fasta | samtools sort -o alignments/ref_pOXA-48_K8_map_K229flye.bam
samtools index alignments/ref_pOXA-48_K8_map_K229flye.bam

minimap2 -ax asm5 ../closed_sequences/plasmids/pOXA-48_K8.fasta assemblies_unicycler_hybrid/K229new-30C6/assembly.fasta | samtools sort -o alignments/ref_pOXA-48_K8_map_K229uni.bam
samtools index alignments/ref_pOXA-48_K8_map_K229uni.bam

minimap2 -ax asm5 ../closed_sequences/plasmids/pOXA-48_K8.fasta final_closed_assemblies/K163WT-28F1.fasta | samtools sort -o alignments/ref_pOXA-48_K8_map_K163uni.bam
samtools index alignments/ref_pOXA-48_K8_map_K163uni.bam

minimap2 -ax asm5 ../closed_sequences/plasmids/pOXA-48_K8.fasta final_closed_assemblies/K165WT.fasta | samtools sort -o alignments/ref_pOXA-48_K8_map_K165uni.bam
samtools index alignments/ref_pOXA-48_K8_map_K165uni.bam
```

Unicycler could not close all contigs in strains C288 and C289. Thus, the Flye assemblies were selected. The assembly of pOXA-48_C288 contained two small deletions and a small insertion not supported by the read mappings to  pOXA-48_K8. Snippy v4.6.0 was used to confirm that these were assembly errors not detected by Pilon, and the .consensus.fa generated was used as the final assembly.

```sh
snippy --report --outdir variants_snippy_ref_C288_map_C288/ --ref final_closed_assemblies/C288WT.fasta --R1 reads_trimmed/C288WT_val_1.fq.gz --R2 reads_trimmed/C288WT_val_2.fq.gz
cp variants_snippy_ref_C288_map_C288/snps.consensus.fa final_closed_assemblies/C288WT.fasta
```

All contigs in strains K163 and K165 were circularized by Unicycler and thus were selected as the final assemblies. However, the sequence of pOXA-48_K163 had a mutation in the *ltrA* gene (position 15,736 in the reference pOXA-48_K8) that did not appear in the mappings of short and long reads. This is possibly due to other versions of the gene in the chromosome or other plasmids of this strain and could not be identified by Snippy, and therefore was corrected manually to the reference allele. The sequence of pOXA-48_K165 also displayed this mutation and an additional small deletion in the *traC* gene. These mutations are not present in a previous assembly of this plasmid, generated by hybrid assembly with PacBio reads (accession number MT989346). Additionally, the *traC* deletion has been confirmed to be false by Sanger sequencing and is suspected to come from the short read assembly. Therefore, the sequence of this pOXA-48_K165 was replaced by the sequence of MT989346. The hybrid and long read assemblies of K165 lacked a 3.1 Kb contig containing a ColRNAI replicon that was present in the assemblies of K163. However, this contig was assembled with only the short reads of K165, suggesting that Unicycler discarded the contig for not being supported by long reads and that the replicon was lost during the extraction or sequencing protocol for Nanopore. The Illumina reads of K165 were mapped to the assembly of K163 to detect possible mutations in this replicon, but none were detected. Thus, the sequence of the ColRNAI replicon of K163 was added to the final assembly of K165.

```sh
snippy --report --outdir variants_snippy_ref_K163_map_K165/ --ref final_closed_assemblies/K163WT-28F1.fasta --R1 reads_trimmed/K165_val_1.fq.gz --R2 reads_trimmed/K165_val_2.fq.gz
```

Both assemblies of strain K153 yielded closed contigs, but the Unicycler assembly was discarded since it included the deletion in the *traC* gene (mappings of the sequencing reads against pOXA-48_K8 rejected this deletion). The Flye assembly of pOXA-48_K153 included a one nucleotide insertion in *ltrA* that was not supported by the read alignments and was therefore manually corrected in the fasta file. For strain K229, the Flye assemblies were also selected since Unicycler failed to circularize one contig. However, the Unicycler assembly of pOXA-48_K229 delimited more accurately the positions of the large deletion (in the long read mappings, most reads support the deletion between the positions 33,165-  45,513 bp of the reference pOXA-48_K8), so the pOXA-48_K229 plasmid from the Flye assembly was replaced by the plasmid from the Unicycler assembly. A 777 bp fragment affecting the IS1 upstream *bla*<sub>OXA-48</sub> (positions 9,382-10,160 bp) was missing in the Unicycler and Flye assemblies of both strains. The short and long read mappings supported this deletion in K229, but not in K153. However, the deletion was confirmed for both strains by PCR.

The final closed assemblies used in this study were annotated locally with PGAP v2021-07-01.build5508. The genomes submitted to GenBank were annotated by the NCBI with a newer version of PGAP.


## Variant calling and plasmid and resistance gene content

Breseq was used to identify SNPs and structural variants in the genomes of isolates from the same patient. To discard false positive variant calls due to missassemblies, for each pair of closed strains from the same patient, breseq was run in both directions. This way, variants present in only one of the directions were discarded. Breseq was also run mapping the reads of a strain to its own assembly to discard false positive variants. Since strains K164, K166, K165-1, K165-3, K165-4, K165-5, K165-6 and K165-7 (patient JWC) and K151 and K152 (patient WDV) are not closed, the trimmed reads (BioProjects PRJNA626430 and PRJNA838107) were only mapped to the closed strains from their respective patients. SNPs and structural variants matching previously known false positives in position or equivalent region, respectively, were discarded. Predicted mutations present in only one of the mappings against a reference that should also show against the other reference were also discarded. In some cases, multiple SNPs were identified in the same region as known false positives. These variants were also discarded after observing poor read alignments in these regions.

```sh
# Breseq mapping reads back to their respective closed genomes

breseq -r final_closed_assemblies/C288WT.gbk reads_trimmed/C288WT_val_1.fq.gz reads_trimmed/C288WT_val_2.fq.gz -o variants_breseq_closed/ref_C288_map_C288
breseq -r final_closed_assemblies/C289WT-31C7.gbk reads_trimmed/C289WT_val_1.fq.gz reads_trimmed/C289WT_val_2.fq.gz -o variants_breseq_closed/ref_C289_map_C289
breseq -r final_closed_assemblies/K153WT-30A7.gbk reads_trimmed/K153_val_1.fq.gz reads_trimmed/K153_val_2.fq.gz -o variants_breseq_closed/ref_K153_map_K153
breseq -r final_closed_assemblies/K229new-30C6.gbk reads_trimmed/K229_val_1.fq.gz reads_trimmed/K229_val_2.fq.gz -o variants_breseq_closed/ref_K229_map_K229
breseq -r final_closed_assemblies/K163WT-28F1.gbk reads_trimmed/K163_val_1.fq.gz reads_trimmed/K163_val_2.fq.gz -o variants_breseq_closed/ref_K163_map_K163
breseq -r final_closed_assemblies/K165WT.gbk reads_trimmed/K165_val_1.fq.gz reads_trimmed/K165_val_2.fq.gz -o variants_breseq_closed/ref_K165_map_K165

# Breseq of within-patient cases to the reference closed genomes

breseq -r final_closed_assemblies/C288WT.gbk reads_trimmed/C289WT_val_1.fq.gz reads_trimmed/C289WT_val_2.fq.gz -o variants_breseq_closed/ref_C288_map_C289
breseq -r final_closed_assemblies/C289WT-31C7.gbk reads_trimmed/C288WT_val_1.fq.gz reads_trimmed/C288WT_val_2.fq.gz -o variants_breseq_closed/ref_C289_map_C288

breseq -r final_closed_assemblies/K153WT-30A7.gbk reads_trimmed/K229_val_1.fq.gz reads_trimmed/K229_val_2.fq.gz -o variants_breseq_closed/ref_K153_map_K229
breseq -r final_closed_assemblies/K229new-30C6.gbk reads_trimmed/K153_val_1.fq.gz reads_trimmed/K153_val_2.fq.gz -o variants_breseq_closed/ref_K229_map_K153
breseq -r final_closed_assemblies/K153WT-30A7.gbk ../Reads_pOXA-48_carriers/K151_1.fq.gz ../Reads_pOXA-48_carriers/K151_2.fq.gz -o variants_breseq_closed/ref_K153_map_K151
breseq -r final_closed_assemblies/K153WT-30A7.gbk ../Reads_pOXA-48_carriers/K152_1.fq.gz ../Reads_pOXA-48_carriers/K152_2.fq.gz -o variants_breseq_closed/ref_K153_map_K152
breseq -r final_closed_assemblies/K229new-30C6.gbk ../Reads_pOXA-48_carriers/K151_1.fq.gz ../Reads_pOXA-48_carriers/K151_2.fq.gz -o variants_breseq_closed/ref_K229_map_K151
breseq -r final_closed_assemblies/K229new-30C6.gbk ../Reads_pOXA-48_carriers/K152_1.fq.gz ../Reads_pOXA-48_carriers/K152_2.fq.gz -o variants_breseq_closed/ref_K229_map_K152

breseq -r final_closed_assemblies/K163WT-28F1.gbk reads_trimmed/K165_val_1.fq.gz reads_trimmed/K165_val_2.fq.gz -o variants_breseq_closed/ref_K163_map_K165
breseq -r final_closed_assemblies/K165WT.gbk reads_trimmed/K163_val_1.fq.gz reads_trimmed/K163_val_2.fq.gz -o variants_breseq_closed/ref_K165_map_K163
breseq -r final_closed_assemblies/K163WT-28F1.gbk ../Reads_pOXA-48_carriers/K164_1.fq.gz ../Reads_pOXA-48_carriers/K164_2.fq.gz -o variants_breseq_closed/ref_K163_map_K164
breseq -r final_closed_assemblies/K163WT-28F1.gbk ../Reads_pOXA-48_carriers/K166_1.fq.gz ../Reads_pOXA-48_carriers/K166_2.fq.gz -o variants_breseq_closed/ref_K163_map_K166
breseq -r final_closed_assemblies/K165WT.gbk ../Reads_pOXA-48_carriers/K164_1.fq.gz ../Reads_pOXA-48_carriers/K164_2.fq.gz -o variants_breseq_closed/ref_K165_map_K164
breseq -r final_closed_assemblies/K165WT.gbk ../Reads_pOXA-48_carriers/K166_1.fq.gz ../Reads_pOXA-48_carriers/K166_2.fq.gz -o variants_breseq_closed/ref_K165_map_K166
breseq -r final_closed_assemblies/K163WT-28F1.gbk reads_trimmed/71_val_1.fq.gz reads_trimmed/71_val_2.fq.gz -o variants_breseq_closed/ref_K163_map_71
breseq -r final_closed_assemblies/K163WT-28F1.gbk reads_trimmed/72_val_1.fq.gz reads_trimmed/72_val_2.fq.gz -o variants_breseq_closed/ref_K163_map_72
breseq -r final_closed_assemblies/K163WT-28F1.gbk reads_trimmed/73_val_1.fq.gz reads_trimmed/73_val_2.fq.gz -o variants_breseq_closed/ref_K163_map_73
breseq -r final_closed_assemblies/K163WT-28F1.gbk reads_trimmed/74_val_1.fq.gz reads_trimmed/74_val_2.fq.gz -o variants_breseq_closed/ref_K163_map_74
breseq -r final_closed_assemblies/K163WT-28F1.gbk reads_trimmed/75_val_1.fq.gz reads_trimmed/75_val_2.fq.gz -o variants_breseq_closed/ref_K163_map_75
breseq -r final_closed_assemblies/K163WT-28F1.gbk reads_trimmed/76_val_1.fq.gz reads_trimmed/76_val_2.fq.gz -o variants_breseq_closed/ref_K163_map_76
breseq -r final_closed_assemblies/K165WT.gbk reads_trimmed/71_val_1.fq.gz reads_trimmed/71_val_2.fq.gz -o variants_breseq_closed/ref_K165_map_71
breseq -r final_closed_assemblies/K165WT.gbk reads_trimmed/72_val_1.fq.gz reads_trimmed/72_val_2.fq.gz -o variants_breseq_closed/ref_K165_map_72
breseq -r final_closed_assemblies/K165WT.gbk reads_trimmed/73_val_1.fq.gz reads_trimmed/73_val_2.fq.gz -o variants_breseq_closed/ref_K165_map_73
breseq -r final_closed_assemblies/K165WT.gbk reads_trimmed/74_val_1.fq.gz reads_trimmed/74_val_2.fq.gz -o variants_breseq_closed/ref_K165_map_74
breseq -r final_closed_assemblies/K165WT.gbk reads_trimmed/75_val_1.fq.gz reads_trimmed/75_val_2.fq.gz -o variants_breseq_closed/ref_K165_map_75
breseq -r final_closed_assemblies/K165WT.gbk reads_trimmed/76_val_1.fq.gz reads_trimmed/76_val_2.fq.gz -o variants_breseq_closed/ref_K165_map_76
```
The strains cured from pOXA-48 were sequenced to control for isogenicity. Variants were called with breseq, mapping the trimmed Illumina reads to their respective reference. The same mutation filtering criteria was applied.

```sh
breseq -r final_closed_assemblies/C288WT.gbk reads_trimmed/C288c2_Clon_2_val_1.fq.gz reads_trimmed/C288c2_Clon_2_val_2.fq.gz -o variants_breseq_closed/ref_C288_map_C288c2
breseq -r final_closed_assemblies/C289WT-31C7.gbk reads_trimmed/C289c1_Clon_1_val_1.fq.gz reads_trimmed/C289c1_Clon_1_val_2.fq.gz -o variants_breseq_closed/ref_C289_map_C289c1
breseq -r final_closed_assemblies/K229new-30C6.gbk reads_trimmed/K229c5_val_1.fq.gz reads_trimmed/K229c5_val_2.fq.gz -o variants_breseq_closed/ref_K229_map_K229c5
breseq -r final_closed_assemblies/K165WT.gbk reads_trimmed/K165c5_val_1.fq.gz reads_trimmed/K165c5_val_2.fq.gz -o variants_breseq_closed/ref_K165_map_K165c5
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

snippy --ref final_closed_assemblies/C288WT.gbk --R1 ../Reads_pOXA-48_carriers/C288_1.fq.gz --R2 ../Reads_pOXA-48_carriers/C288_2.fq.gz --outdir phylogeny_core/snippy_Ec/ref_C288_map_C288
snippy --ref final_closed_assemblies/C288WT.gbk --R1 ../Reads_pOXA-48_carriers/C289_1.fq.gz --R2 ../Reads_pOXA-48_carriers/C289_2.fq.gz --outdir phylogeny_core/snippy_Ec/ref_C288_map_C289
snippy --ref final_closed_assemblies/C288WT.gbk --R1 ../Reads_pOXA-48_carriers/C165_1.fq.gz --R2 ../Reads_pOXA-48_carriers/C165_2.fq.gz --outdir phylogeny_core/snippy_Ec/ref_C288_map_C165
snippy --ref final_closed_assemblies/C288WT.gbk --R1 ../Reads_pOXA-48_carriers/C166_1.fq.gz --R2 ../Reads_pOXA-48_carriers/C166_2.fq.gz --outdir phylogeny_core/snippy_Ec/ref_C288_map_C166
snippy --ref final_closed_assemblies/C288WT.gbk --R1 ../Reads_pOXA-48_carriers/C642_1.fq.gz --R2 ../Reads_pOXA-48_carriers/C642_2.fq.gz --outdir phylogeny_core/snippy_Ec/ref_C288_map_C642
snippy --ref final_closed_assemblies/C288WT.gbk --R1 ../Reads_pOXA-48_carriers/C643_1.fq.gz --R2 ../Reads_pOXA-48_carriers/C643_2.fq.gz --outdir phylogeny_core/snippy_Ec/ref_C288_map_C643
snippy --ref final_closed_assemblies/C288WT.gbk --R1 ../Reads_pOXA-48_carriers/C646_1.fq.gz --R2 ../Reads_pOXA-48_carriers/C646_2.fq.gz --outdir phylogeny_core/snippy_Ec/ref_C288_map_C646
snippy --ref final_closed_assemblies/C288WT.gbk --R1 ../Reads_pOXA-48_carriers/C662_1.fq.gz --R2 ../Reads_pOXA-48_carriers/C662_2.fq.gz --outdir phylogeny_core/snippy_Ec/ref_C288_map_C662
snippy --ref final_closed_assemblies/C288WT.gbk --R1 ../Reads_pOXA-48_carriers/C667_1.fq.gz --R2 ../Reads_pOXA-48_carriers/C667_2.fq.gz --outdir phylogeny_core/snippy_Ec/ref_C288_map_C667
snippy --ref final_closed_assemblies/C288WT.gbk --R1 ../Reads_pOXA-48_carriers/C717_1.fq.gz --R2 ../Reads_pOXA-48_carriers/C717_2.fq.gz --outdir phylogeny_core/snippy_Ec/ref_C288_map_C717
snippy --ref final_closed_assemblies/C288WT.gbk --R1 ../Reads_pOXA-48_carriers/C718_1.fq.gz --R2 ../Reads_pOXA-48_carriers/C718_2.fq.gz --outdir phylogeny_core/snippy_Ec/ref_C288_map_C718
snippy --ref final_closed_assemblies/C288WT.gbk --R1 ../Reads_pOXA-48_carriers/C728_1.fq.gz --R2 ../Reads_pOXA-48_carriers/C728_2.fq.gz --outdir phylogeny_core/snippy_Ec/ref_C288_map_C728
snippy --ref final_closed_assemblies/C288WT.gbk --R1 ../Reads_pOXA-48_carriers/C752_1.fq.gz --R2 ../Reads_pOXA-48_carriers/C752_2.fq.gz --outdir phylogeny_core/snippy_Ec/ref_C288_map_C752
snippy --ref final_closed_assemblies/C288WT.gbk --R1 ../Reads_pOXA-48_carriers/N46_1.fq.gz --R2 ../Reads_pOXA-48_carriers/N46_2.fq.gz --outdir phylogeny_core/snippy_Ec/ref_C288_map_N46

# K. pneumoniae ST11

snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K153_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K153_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K153
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K90_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K90_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K90
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K109_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K109_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K109
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K24_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K24_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K24
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K31_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K31_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K31
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K47_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K47_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K47
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K51_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K51_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K51
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K52_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K52_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K52
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K8_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K8_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K8
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K88_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K88_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K88
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K89_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K89_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K89
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K93_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K93_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K93
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K1_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K1_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K1
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K116_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K116_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K116
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K127_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K127_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K127
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K147_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K147_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K147
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K15_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K15_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K15
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K151_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K151_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K151
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K152_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K152_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K152
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K160_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K160_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K160
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K161_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K161_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K161
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K162_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K162_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K162
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K197_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K197_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K197
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K2_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K2_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K2
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K202_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K202_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K202
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K229_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K229_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K229
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K236-1_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K236-1_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K236-1
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K237_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K237_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K237
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K238_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K238_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K238
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K25_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K25_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K25
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K259_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K259_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K259
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K26_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K26_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K26
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K260_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K260_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K260
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K263_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K263_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K263
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K264_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K264_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K264
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K265_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K265_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K265
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K266_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K266_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K266
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K267_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K267_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K267
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K268_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K268_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K268
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K269_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K269_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K269
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K270_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K270_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K270
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K271_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K271_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K271
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K28_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K28_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K28
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K280_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K280_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K280
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K281_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K281_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K281
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K29_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K29_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K29
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K293_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K293_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K293
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K294_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K294_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K294
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K295_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K295_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K295
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K298_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K298_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K298
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K308_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K308_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K308
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K312_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K312_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K312
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K313_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K313_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K313
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K317_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K317_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K317
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K319_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K319_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K319
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K323_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K323_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K323
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K324_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K324_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K324
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K326_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K326_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K326
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K36_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K36_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K36
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K39_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K39_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K39
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K4_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K4_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K4
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K57_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K57_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K57
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/R25_1.fq.gz --R2 ../Reads_pOXA-48_carriers/R25_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_R25
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/X6_1.fq.gz --R2 ../Reads_pOXA-48_carriers/X6_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_X6
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/H73_1.fq.gz --R2 ../Reads_pOXA-48_carriers/H73_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_H73
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K70_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K70_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K70
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K76_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K76_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K76
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/G21_1.fq.gz --R2 ../Reads_pOXA-48_carriers/G21_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_G21
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/L27_1.fq.gz --R2 ../Reads_pOXA-48_carriers/L27_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_L27
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K101_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K101_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K101
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/L44_1.fq.gz --R2 ../Reads_pOXA-48_carriers/L44_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_L44
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/L74_1.fq.gz --R2 ../Reads_pOXA-48_carriers/L74_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_L74
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/Ñ22_1.fq.gz --R2 ../Reads_pOXA-48_carriers/Ñ22_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_Ñ22
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/N23_1.fq.gz --R2 ../Reads_pOXA-48_carriers/N23_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_N23
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/Ñ34_1.fq.gz --R2 ../Reads_pOXA-48_carriers/Ñ34_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_Ñ34
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/N45_1.fq.gz --R2 ../Reads_pOXA-48_carriers/N45_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_N45
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/P76_1.fq.gz --R2 ../Reads_pOXA-48_carriers/P76_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_P76
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/L45_1.fq.gz --R2 ../Reads_pOXA-48_carriers/L45_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_L45
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K34_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K34_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K34
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/E73_1.fq.gz --R2 ../Reads_pOXA-48_carriers/E73_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_E73
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/E74_1.fq.gz --R2 ../Reads_pOXA-48_carriers/E74_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_E74
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/E76_1.fq.gz --R2 ../Reads_pOXA-48_carriers/E76_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_E76
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/F14_1.fq.gz --R2 ../Reads_pOXA-48_carriers/F14_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_F14
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K62_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K62_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K62
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/F48_1.fq.gz --R2 ../Reads_pOXA-48_carriers/F48_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_F48
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/K74_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K74_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_K74
snippy --ref final_closed_assemblies/K153WT-30A7.gbk --R1 ../Reads_pOXA-48_carriers/F73_1.fq.gz --R2 ../Reads_pOXA-48_carriers/F73_2.fq.gz --outdir phylogeny_core/snippy_ST11/ref_K153_map_F73

# K. pneumoniae ST307

snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../Reads_pOXA-48_carriers/K163_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K163_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_K163
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../Reads_pOXA-48_carriers/K164_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K164_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_K164
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../Reads_pOXA-48_carriers/K165_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K165_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_K165
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../Reads_pOXA-48_carriers/K166_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K166_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_K166
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../Reads_pOXA-48_carriers/R50_1.fq.gz --R2 ../Reads_pOXA-48_carriers/R50_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_R50
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../Reads_pOXA-48_carriers/K172_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K172_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_K172
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../Reads_pOXA-48_carriers/K173_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K173_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_K173
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../Reads_pOXA-48_carriers/K174_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K174_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_K174
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../Reads_pOXA-48_carriers/K198_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K198_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_K198
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../Reads_pOXA-48_carriers/K199_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K199_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_K199
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../Reads_pOXA-48_carriers/K235_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K235_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_K235
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../Reads_pOXA-48_carriers/K306_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K306_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_K306
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../Reads_pOXA-48_carriers/K314_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K314_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_K314
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../Reads_pOXA-48_carriers/K315_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K315_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_K315
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../Reads_pOXA-48_carriers/K318_1.fq.gz --R2 ../Reads_pOXA-48_carriers/K318_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_K318
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../Reads_pOXA-48_carriers/J61_1.fq.gz --R2 ../Reads_pOXA-48_carriers/J61_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_J61
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../Reads_pOXA-48_carriers/J66_1.fq.gz --R2 ../Reads_pOXA-48_carriers/J66_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_J66
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../Reads_pOXA-48_carriers/N22_1.fq.gz --R2 ../Reads_pOXA-48_carriers/N22_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_N22
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 ../Reads_pOXA-48_carriers/N36_1.fq.gz --R2 ../Reads_pOXA-48_carriers/N36_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_N36
snippy --ref final_closed_assemblies/K163WT-28F1.gbk --R1 reads_trimmed/74_val_1.fq.gz --R2 reads_trimmed/74_val_2.fq.gz --outdir phylogeny_core/snippy_ST307/ref_K163_map_74
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
