# small_RNA_analysis

1. Prerequisite: ShortStack (https://github.com/MikeAxtell/ShortStack)
                 Bsmap2 (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-232)
                 python3

2. Purpose: The scripts in this repository are used to analyze argonaute-associated small interfering RNAs (siRNAs) and their roles in RNA-directed DNA methylation in plants.

3. Usage:
   -  Trimmed reads from small RNA seq libraries are mapped to the reference genome by using ShortStack with no more than 2 mismatches

   -  To study the size distribution of immunoprecipitated RNAs, use RNAIP_size_distribution.py:
      
      python RNAIP_size_distribution.py --size --bam YOUR_INPUT_BAM --gff SMALL_RNA_CLUSTER_GFF
      
   -  To simulate small RNA alignments, use siRNA_simulator.py:
      
      python siRNA_simulator.py --genome GENOME_REF_IN_USE --bam YOUR_INPUT_BAM_TO_SHUFFLE --gff3 SMALL_RNA_CLUSTER_GFF --output YOUR_OUTPUT_FILE_NAME --path OUTPUT_DIRECTORY
   
   -  To computationally predict small RNA pairs, use pairs_prediction.py:
      
      python pairs_prediction.py --bam YOUR_INPUT_BAM --gff3 SMALL_RNA_CLUSTER_GFF --standard --output YOUR_OUTPUT_FILE_NAME
      
   -  To calculate mismatch rate per nucleotide position in a given size of small RNAs, use sRNA_mismatches.py
      
      python sRNA_mismatches.py --bam YOUR_INPUT_BAM_TO_SHUFFLE --gff3 SMALL_RNA_CLUSTER_GFF
      
   -  To identify bins of DNA methylation, use bin_methyl.py. You will need to first run bsmap2 first to map reads from whole genome bisulfite sequencing libraries to the reference genome, calculate cytosine methylation across the whole genome and then run 'bin_methyl.py':
      
      python bin_methyl.py -m SUMMARY_OF_CYTOSINE_METHYLATION_(FROM_BSMAP2) -o YOUR_OUTPUT_FILE_NAME -i REFERENCE_GENOME_INDEX -b BIN_SIZE -s MOVING_WINDOW_SIZE
      
      
Contact wangfeng at iu dot edu if you have questions.
