# 04_whole_genome_alignments/

- confirm broad-scale chromosomal synteny and orientation across all accessions using whole‑genome dotplots (HART001 as reference) 

- correct a putative misassembly: invert Chr07 in HART001 based on consistent dotplot inversion + Hi‑C support, ensuring all genomes share a comparable Chr07 direction before pangenome alignment 



Outputs, chromosomal synteny among *Artocarpus altilis* : 

![dotplots](/imgs/20260505_wga.png)

___



## Dotplots against HART001

Using HART001 as the reference, produce dotplots:

```bash
#!/bin/bash

#SBATCH --time=1-00:00:00    
#SBATCH --cpus-per-task=10
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

SAMPLE=$1

# module load minicionda
# source activate puzzler200
module load apptainer

WD=/project/coffea_pangenome/Artocarpus/Pangenome_Paper/dotplots
GENOMES=/project/coffea_pangenome/Artocarpus/Pangenome_Paper/assemblies/unmasked
REF=/project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/primary_asm/init/HART001.Chr07sub.fa

mkdir -p pafs

echo "WGA for ${SAMPLE}"
mashmap -r ${REF} -q ${GENOMES}/${SAMPLE}.chr.fa -t 10 -s 10000 --perc_identity 95 -o pafs/${SAMPLE}.paf 2> pafs/${SAMPLE}.mashmap.log
Rscript ~/apptainer/paf2dotplot.R pafs/${SAMPLE}.paf -r 1e6 -m 1e4 -p 4 -c 1 --sort-by-refid --identity-lower-color 100
```

Submit: 

```bash
cat Samples.list  | xargs -I {} sbatch -J wga_{} 01_WGA_HART001Ref.sh {}
zip wga_dotplots_hart001.zip pafs/*pdf
```

## Invert Chr07 on HART001

Initial dotplots across all 14 accessions against sample HART001 identifies a consistent inversion on Chr07, which after closer inspection of the HiC heat map for HART001, is a putative misassembly. For consistency in this paper, I will reverse complement that contig (confirmed it is a long arm gapped with >50 Ns), so that we have a consistent Chr07 direction. 

![HART001_misassembly](/imgs/HART001_misassembly_inversion.png)

* isolate Chr07 second contig around big N-gap (≥51 Ns) → split into part1 / gap / part2 

* reverse‑complement part2 → rebuild Chr07 with inverted arm (p1 + gap + p2_rc) 

- rename rebuilt contig and verify lengths w/ seqkit 

- replace original Chr07 in HART001.fa with new inverted version → write HART001.Chr07sub.fa 

* Rearrange second contig of Chr07 in HART001. 

```bash
#split 
perl -0777 -pe '
s/^>.*\n//m;        
s/\s+//g;          
if (!/(N{51,})/) { die "ERROR: no N-run >=51 found\n"; }
my $gap = $1;
my ($a,$b) = split(/N{51,}/, $_, 2);
open(F1, ">part1.fa") or die $!;
print F1 ">Chr07_part1\n$a\n";
open(FG, ">gap.fa") or die $!;
print FG ">Chr07_gap\n$gap\n";
open(F2, ">part2.fa") or die $!;
print F2 ">Chr07_part2\n$b\n";
$_="";
' Chr07.fa


seqkit seq -r -p part2.fa > part2.rc.fa
p1=$(perl -ne 'next if /^>/; s/\s+//g; print' part1.fa)
gap=$(perl -ne 'next if /^>/; s/\s+//g; print' gap.fa)
p2rc=$(perl -ne 'next if /^>/; s/\s+//g; print' part2.rc.fa)

{
  echo ">Chr07_invB"
  printf "%s%s%s\n" "$p1" "$gap" "$p2rc"
} > Chr07_invB.fa

#verify
seqkit fx2tab -n -l Chr07.fa Chr07_invB.fa


perl -pe 'if($.==1){s/^>.*/>Chr07/}' Chr07_invB.fa > Chr07.repl.fa

perl -0777 -ne '
BEGIN {
  $repl = do {
    open my $fh, "<", "Chr07.repl.fa" or die $!;
    local $/; <$fh>
  };
}

# split into fasta records
@recs = split(/(?=>)/, $_);
for $r (@recs) {
  next if $r =~ /^\s*$/;

  if ($r =~ /^>(\S+)/) {
    $id = $1;
    if ($id eq "Chr07") {
      print $repl;
      $repl = "";  
    } else {
      print $r;
    }
  }
}
' HART001.fa > HART001.Chr07sub.fa

seqkit fx2tab -n -l HART001.fa | wc -l
seqkit fx2tab -n -l HART001.Chr07sub.fa | wc -l

seqkit fx2tab -n -l HART001.Chr07sub.fa | grep -w '^Chr07'
seqkit fx2tab -n -l Chr07.repl.fa
```

## HART069 "Butterfly"

Sample HART069 looks like it has a truncated Chr09 and an extra large Chr022, which appears as a misassembly in the dotplot against HART001. However, the juicer contact map for this sample does not show a straight story: instead it shows a potential chromosome-arm duplication, as it seems that the portion of chr09 missing ALSO matches chr22. 

This could require a phased genome approach to resolve, but archiving this to remember any issues with chr09 or chr22 in the future. 

![HART069 map](/imgs/HART069_inversion_chr22_chr09.png)
