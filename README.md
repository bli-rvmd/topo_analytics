# topo_analytics
Analysis scripts, code snippets, and data templates for TOPO cloning analytics. 
# Analysis workflow:
* QC on .seq files from a TOPO cloning experiments
* Convert QC'ed and cleaned TOPO sequences into .fastq format
* Run CRISPResso2 as AmpliconSeq obtaining analysis report and result files
* Examine %reads getting aligned, %editing efficiency, topmost variants, etc.
* Convert aligned sequences of variants in result file of allelic frequency table into VCF format
* Run VEP for variants functional prediction
* Collate results from VCF and from running VEP to calculate overall %edits and proportion of likely disrupting protein function among edited sequences  
