# Cry phylogeny
The repository with working codes and supplementary materials for the study of Cry toxins evolution

## Contents 

This repository contains scripts used for statistical processing, phylogeny reconstruction, inferencing recombination events, and analyzing evolutionary selection in the sequences of 3-domain Cry toxins. Please consult the Methods section in the paper for extra details:

## Figures
Figures from the main text are available in the `pics/` directory. For the description, please consult the Results section of the article.

## Supplementary
All supplementary material is located in the `supplementary/` directory. 
`Supplementary_figures.docx` is a Microsoft Word document with a description of all supplementary figures.
`Supplementary_tables.xlsx` is an Exel table with all supplementary tables.

* The `supplementary/pics/` folder contains individual supplementary figures.
* For convenience, in the `supplementary/tables/` folder individual supporting tables in CSV format are provided. For a detailed description of the tables, please consult the `Supplementary_tables.xlsx` file.

## Data
Analyzed data are included in the `data/` directory. The directory contains sequences and phylogenetic trees. Other essential tables obtained when using data processing scripts are presented in the `data_for_scripts/` directory. 

### Sequences
The `data/fastas` directory contains nucleotide and protein sequences, both full and domain-wise.
The raw set of 3-D Cry protein sequences is provided in the `source_Cry_proteins.fasta` file. The file contains full-length protein sequences inferred from four datasets: BPPRC (Bacterial Pesticidal Protein Resource Center), NCBI IPG database (<i>Bacillus thuringiensis</i>), 467 <i>Bt</i> assemblies from NCBI Assembly, and NCBI Genbank (bacteria non-redundant). Three-domain sequences were obtained with CryProcessor in the “<i>-do</i>” mode. The presented set of sequences is redundant and was further deduplicated.
* The directory `data/fastas/all_sequences` includes non-redundant set from the `source_Cry_proteins.fasta` file. The duplicates with a 100% identity threshold are discarded. Cry toxins absent in the natural strains (presented only in patents) and subjected to artificial mutagenesis are also excluded;</li>
* The directory `data/fastas/ref_sequences` includes reference sequences obtained by applying CD-HIT software with a 95% identity threshold on the non-redundant filtered set of sequences in the `all_sequences` directory.
All the files in two subdirectories follow the same rules for prefixes and postfixes:
* The `nucl_` prefix corresponds to nucleotide sequences;
* The `prot_` prefix marks protein sequences;
* Sequences of individual domains are marked with the `_domain` postfix;
* The `_full` postfix indicates that the sequences contain all regions, including N- and C-terminal ones;
* The `_ processed` postfix indicates that the sequences contain only three main functional domains and are cropped from the beginning of the first domain to the end of the third domain accordingly. 

### Phylogenetic trees
The `data/phylogeny` directory includes reconstructed phylogenies and multiple sequence alignments used for building phylogenetic trees. The directory contains the following subdirectories and files:
1. The `phylogeny/data` subdirectory comprises source files used for phylogeny reconstruction:
* The `phylogeny/data/concat_msa.fasta` file includes concatenated domain-wise alignments;
* The `phylogeny/data/concat_msa_trimmed.fasta` file is the source data for reconstructing reference phylogenetic tree derived from `concat_msa.fasta` using the trimAl utility;
* The `phylogeny/data/concat_msa.fasta` file is a set of partitions for trimmed sequences of the individual domains within the concatenated alignment;
* The `phylogeny/data/domains_msa` subdirectory includes both raw and trimmed alignments of the nucleotide sequences of individual domains;
2. The `phylogeny/init_trees` directory contains both domain-wise trees and the tree inferred from the concatenated partitioned alignment;
3. The `phylogeny/collapsed_trees` comprises domain-wise phylogenetic trees corrected for recombination events with parents and recombinants grouped in case of multifurcation. A detailed description of the filtering procedure is presented in the Methods section of the article.

## Data for scripts
For better reproducibility, the data needed for obtaining all the results described in the article is provided in the `data_for_scripts/` directory. The subdirectories are divided according to the subsections in the Results section of the article. The directory `data_for_scripts/for_viz` contains source data required for producing figures in the main text and supplementary material. The presented data could be generated by scripts in the `scripts/` directory. 

### Data description
In this part of the study, Cry sequences’ properties were analyzed. The following characteristics encompassed the number of putative novel toxins, the lengths, pair-wise identities, and conservation of the sequences as well as the patterns of clusterisations and the distribution of species encoding <i>cry</i> genes.  The respective `data_for_scripts/general_props` directory contains the following files and subdirectories. 
- In the `general_props/domain_sequences` directory the following nucleotide and protein sequences of the individual domains are presented:
  * `domain_sequences/all_seqs_domains` - nucleotide domain sequences form the deduplicated dataset with a 100% identity threshold;
  * `domain_sequences/unique_domain_sequences` - nucleotide domain sequences in the above-described group after domain-wise deduplication;
  * `domain_sequences/ref_nucleotides` - nucleotide domain sequences of the reference clusters obtained by applying CD-HIT software with a 95% identity threshold;
  * `domain_sequences/ref_proteins` - protein sequences of the domains for the above-mentioned reference clusters;
- In `general_props/sequences_with_patents` the processed nucleotide sequences (cropped from the beginning of the first domain to the end of the third domain accordingly) for the sets of sequences with 100%  and 95% identity thresholds are present. Unlike other files, the respective sequences include toxins from patents only, including those subjected to artificial mutagenesis. These sequences were discarded for further analysis;
- In the `general_props/multiclusters` directory the following domain-wie nucleotide and protein sequences of the clusters obtained by applying CD-HIT software with a 95% identity threshold with more than 10 toxins are presented:
  * `nucl` - nucleotide domain sequences for each cluster;
  * `prot` - nucleotide domain sequences for each cluster;
  * `msa_prot` - alignments of the domain-wise protein sequences for each cluster;
  * `nucl ` - codon-wise alignments of the domain-wise nucleotide sequences for each cluster;
- `all_prot_95_full_names.fasta.clstr` - the results of CD-HIT software with clustering patterns applied on full protein Cry sequences with a 95% identity threshold;
- `all_orfs_processed_nucl.fasta`  - nucleotide processed (cropped from the beginning of the first domain to the end of the third domain) sequences form the deduplicated dataset clustered with a 100% identity threshold;
- `ref_nucl_processed.fasta` - the set of processed nucleotide sequences for reference clusters obtained with a 95% identity threshold;
- `all_sequences_nucl.aln.fasta` - multiple sequence alignments of the processed nucleotide sequences from the deduplicated dataset with a 100% identity threshold applied;
- `no_pat_cd_hit_clusters.tsv` - descriptive characteristics of the CD-HIT clusters with a 95% identity threshold (the names of the reference sequence, the number of toxins in the cluster, both raw and unique according to processed nucleotide sequences). Cry toxins absent in the natural strains (presented only in patents) and subjected to artificial mutagenesis are discarded;
- `bt_nomenclature_table.csv` - the names of known Cry proteins with the respective accession numbers (nucleotide and proteins) from BPPRC (Bacterial Pesticidal Protein Resource Center);
- `merged.filtered_events.s70.l3.csv` -  the properties of filtered recombination events. The table here is needed for the annotate_clusters_consistency.py script to assess the presence of sequences from patents in recombination events;
- `clusters_table.csv` - the composition of reference clusters with a 95% identity threshold applied; 
- `all_annotations.tsv` - the data from the IPG database (genome assemblies, species, strain, protein, and nucleotide accessions) for each Cry sequence in the deduplicated dataset with a 100% identity threshold for clusterization;
- `mearged_length_and_id.csv` - the lengths and mean pair-wise identity for individual domains in deduplicated datasets with 100% and 95% identity thresholds, respectively;
- `pairs_identity_no_pat.tsv` - domain-wise pair-wise comparisons of the domain sequences within the reference dataset with a 95% identity threshold applied;
- `full_contat.nwk` - the reference phylogenetic tree inferred from the concatenated partitioned alignment.

### Recombination detection
The following section includes source files used for inferencing and characterizing recombination events, including the filtration procedure, comparing identities between parents and recombinants, reconstructing the recombination graph, etc. 
- `all_domains` - nucleotide domain sequences form the deduplicated dataset with a 100% identity threshold;
- `ref_domains` - nucleotide domain sequences of the reference dataset with a 95% identity threshold for clustering;
- `no_support_trees` - domain-wise phylogenetic trees and the tree based on the concatenated partitioned alignment without supporting values and branch lengths;
- `parsed_trees` - the node-wise structure of phylogenetic trees. In the respective tables, the content of toxins, depth level, and the number of sequences for each node are presented; 
- `collapsed_trees` - domain-wise phylogenetic trees corrected for recombination events with parents and recombinants grouped in case of multifurcation. A detailed description of the filtering procedure is presented in the Methods section of the article;
- `domain_mapings.bed` - coordinates of the domain mappings of the processed sequences (cropped from the beginning of the first domain to the end of the third domain) of the reference dataset with a 95% identity threshold applied; 
- `full_contat.nwk` - the reference phylogenetic tree inferred from the concatenated partitioned alignment;
- `RDP_raw_signals.csv` - the output of the RDP software used for detecting recombination events for the alignment of processed nucleotide sequences;
- `merged.filtered_events.s70.l3.csv` - the characteristics of recombination events inferred from the RDP tool. Presented are the lists of parents and recombinant, coordinates of the breakpoints, p-value levels of detection tests, and pair-wise identity of the domain sequences of parents and children. The list of the events is filtered based on congruence in phylogenetic trees. A detailed scheme for the filtration procedure is present in the Methods section of the article;
- `pairs_identity_no_pat.tsv` - domain-wise pair-wise comparisons of the domain sequences within the reference dataset with a 95% identity threshold applied;
- `unique_multiclusters.tsv` - the domain-wise lists of sequences within the reference clusters obtained using a 95% identity threshold. Unique nucleotide sequences for each domain are selected.

### Analysis of recombination mechanisms
In this section, putative mechanisms of recombination were analyzed. Multiple approaches were applied, namely, assessing the overall homologous recombination rate in genome assemblies with loci coding for Cry toxins, studying the genomic context of <i>cry</i> genes, and comparing sequence identity between parents in regions surrounding recombination breakpoints. The content of the `data_for_scripts/mechanisms` directory goes as follows:
- `MGE_beds` is a directory containing coordinates in the bed format of MGEs (mobile genetic elements) and <i>cry</i> genes within <i>Bacillus thuringiensis</i> genome assemblies. The directory has the following substructure:
  * `MGE_beds/recs` - coordinates of ‘cry’ genes encoding toxins subjected to recombination-driven domain exchanges;
  * `MGE_beds/all_crys` - coordinates of all ‘cry’ genes;
  * `MGE_beds/CRISPRs` - coordinates of CRISPR loci;
  * `MGE_beds/GIs` - coordinates of genomic islands;
  * `MGE_beds/IS` - coordinates of insertion sequences.
Each directory except for ‘recs’ and ‘all_crys’ includes source, crys_inter, and recs_inter subdirectories with the coordinates of genomic features and their intersections with loci encoding Cry toxins and toxins subjected to recombination, respectively. 
- `asmbl_types.csv` - the classification of genome assemblies according to the presence of <i>cry</i> genes and <i>cry</i> loci encoding toxins with recombination signals;
- `bt_assemblies.tsv` - characteristics of genome <i>Bt</i> genome assemblies, including level, serovar (if present in the name of the organism), and links to the FTP repository for downloading;
- `cry_asmbl_stat.tsv` - genomic data for `cry` loci in the <i>Bt</i> assemblies, namely, accession number, protein name, identity with the closest homolog from the BPPRC database, genomic coordinates, strand orientation, and assembly level;
- `Mappings_for_breakpoints_de_novo_search.csv` - coordinates of breakpoints for each parent and recombinant within recombination events as well as the coordinates of the domains. The coordinates are presented both for full and processed (cropped from the beginning of the first domain to the end of the third domain) nucleotide Cry sequences; 
- `merged.filtered_events.s70.l3.csv` - the filtered set of recombination events inferred from the RDP tool based on the congruence in phylogenetic trees. A detailed scheme for the filtration procedure is present in the Methods section of the article;
- `no_pat_cd_hit_clusters.tsv` - descriptive characteristics of the CD-HIT clusters with a 95% identity threshold (the names of the reference sequence, the number of toxins in the cluster, both raw and unique according to processed nucleotide sequences). Cry toxins absent in the natural strains (presented only in patents) and subjected to artificial mutagenesis are discarded;
- `RDP_pre_filtered_fixed.tsv` - properties pf recombination events obtained from the RDP tool with the first step of filtration (events that are marked as dubious are discarded);
- `RDP_raw_signals.csv` - the output of the RDP software used for detecting recombination events for the alignment of processed nucleotide sequences;
- `ref_domains` - nucleotide sequences of the domains for reference representative of the clusters; 
- `Reference_domain_mapppings_full.bed` - coordinates of the domains in full nucleotide sequences in the reference dataset;
- `Reference_domain_mapppings_processed.bed` - coordinates of the domains in processed nucleotide sequences in the reference dataset;
- `ref_all_processed_nucl_aln.fasta` - multiple sequence alignments of the processed nucleotide sequences from the reference dataset;
- `ref_full_nucl_corrected.fasta` - full nucleotide sequences from the reference dataset;
- `ref_nucl_processed.fasta` - processed nucleotide sequences from the reference dataset.

## Scripts
The `scripts/` directory includes all code used for pangenome analysis. The scrips are grouped into two categories, namely, those used for data processing and visualization. These groups are further subdivided following the subsections in the Results section of the article. For convenience, required data are presented in the `data_for_script /` directory. Here, commands to run the scripts are presented. The description of input files is given above. For convenience, the data presented in the `data_for_scripts/` is also divided into subdirectories according to sections in the manuscript. The paths for the input files are given relative to the `scripts/` directory.
