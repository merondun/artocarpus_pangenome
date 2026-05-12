# *Artocarpus altilis* evolutionary genomics (assembly → pangenomics → domestication history)

End-to-end assembly and intraspecific pangenomics for *Artocarpus altilis*, spanning genome QC, assembly, pangenome graph inference, and tracing domestication history with related species.

Core sample metadata live in `samples.info`, with read info and input for [puzzler](https://github.com/merondun/puzzler) assembly in `samples.tsv`.

## Directory map

- [**01_qa_qc_genomescope/**](01_qa_qc_genomescope/) — read/QC and GenomeScope summaries for genome size/heterozygosity context.
- [**02_genome_assembly/**](02_genome_assembly/) — assembly generation and post-processing notes/scripts.
- [**03_whole genome alignments/**](03_whole_genome_alignments/) — initial dotplots against a single reference. 
- [**04_pangenome_inputs/**](04_pangenome_inputs/) — files used for 14-accession primary haplotype pangenome.  


## Overview

This project generates collapsed haploid assemblies for 13 novel breadfruit accessions - in addition to 1 breadfruit genome assembled previously for a [comparative genomics paper](https://github.com/merondun/artocarpus_comparative_genomics), providing 14 breadfruit chromosome-level genomes.

These 14 cultivars are divided into a few 'ploidy groups', based on ploidy and morphological variation. 

| Cultivar Name | sample  | ploidy | chromosomes | flow cytometry size | Ploidy Group         | HiFi Gb inputs | HiFi Length Kb | HiC Gb inputs |
| ------------- | ------- | ------ | ----------- | ------------------- | -------------------- | -------------- | -------------- | ------------- |
| Maafala       | HART001 | 2      | 28          | 891                 | A. altilis 2N        | 59.3           | 11206.3        | 50.5          |
| Huero         | HART030 | 2      | 28          | 917.9               | A. altilis 2N        | 96.6           | 18338.2        | 17.6          |
| Kukumu tasi   | HART050 | 2      | 28          | 883.1               | A. altilis 2N        | 39             | 17753.1        | 18.1          |
| Ulu Fiti      | HART069 | 2      | 28          | 907.1               | A. altilis 2N        | 73             | 12674.4        | 64.5          |
| Hamoa         | HART032 | 3      | 28          | 1367.8              | A. altilis 3N        | 88.3           | 15382.2        | 26.4          |
| Patara        | HART033 | 3      | 28          | 1328.1              | A. altilis 3N        | 139.8          | 14764.9        | 37            |
| Lemai         | HART038 | 3      | 28          | 1331.7              | A. altilis 3N        | 42.6           | 18956.6        | 29.4          |
| Ulu hamoa     | ZZ3     | 2      | 28          | n/a                 | A. altilis hybrid 2N | 29.2           | 18949.4        | 0             |
| Ulu afa elise | ZZ7     | 2      | 28          | n/a                 | A. altilis hybrid 2N | 28.6           | 20738.2        | 0             |
| Ulu afa       | ZZ9     | 2      | 28          | n/a                 | A. altilis hybrid 2N | 30             | 19115.8        | 0             |
| Midolab       | HART046 | 3      | 28          | 1306.8              | A. altilis hybrid 3N | 111.2          | 14132.7        | 37.7          |
| Meinpohnsakar | HART049 | 3      | 28          | 1324                | A. altilis hybrid 3N | 92.9           | 8248.2         | 3.5           |
| Nahnmwal      | HART053 | 3      | 28          | 1320.2              | A. altilis hybrid 3N | 78.2           | 15006          | 32.6          |
| Meikole       | H6      | 4      | 28          | 1790                | A. altilis hybrid 4N | 87.4           | 15742.8        | 35.3          |

Overview of the ploidy groups:

* Seeded diploid (2N) A. altilis (e.g. HART001, HART050)
* Early generation diploid (2N) hybrids (ZZ3, ZZ9)
* Seedless triploid (3N) A. altilis (HART030, HART033, HART053)
* Seedless triploid (3N) hybrids (HART046, HART049)
* Seeded tetraploid (4N) A. altilis (e.g. H6) 

These groups largely follow the morphological classifications from [Jones et al 2013](https://doi.org/10.1007/s10722-012-9824-8) Fig. 6, below:

![Jones et al](/imgs/Jones_etal_2013.png)



Color palette:

| Group                | Color   |
| -------------------- | ------- |
| A. altilis 2N        | #0072B2 |
| A. altilis 3N        | #56B4E9 |
| A. altilis hybrid 2N | #D55E00 |
| A. altilis hybrid 3N | #E69F00 |
| A. altilis hybrid 4N | #009E73 |



## Qs & Cs

Questions or comments reach out to Justin Merondun heritabilities [@] gmail.com or make an issue here. 



