# SNV_Annotation_Tools
Useful tools for variant annotation

Variant Annotation Tools

[SNP-Nexus] (https://www.snp-nexus.org/v4/)
SNPnexus is a web-based variant annotation tool designed to simplify and assist in the selection and prioritisation of known and novel genomic alterations.
Input: multiple types, including VCF, web-based
Output: annotated files, CSV, TSV, HTML, graphics

[ANNOVAR] (http://annovar.openbioinformatics.org/en/latest/)
command-line tool, supports SNPs, INDELs, CNVs and block substitutions, provides wide variety of annotation techniques, utilizes RefSeq, UCSC Genes, and the Ensembl gene annotation systems; can compare mutations detected in dpSNP or 1000 Genomes Project.
Input: VCF, ANNOVAR input format (simple text-based format); can convert other formats into ANNOVAR input format
Output: VCF (if input VCF), output file with multiple columns, tab-delimited output file

[wANNOVAR] (http://wannovar.usc.edu/)
web-based access to ANNOVAR

[PolyPhen-2] (http://genetics.bwh.harvard.edu/pph2/)
Polymorphism Phenotyping; Web application; predicts impact of amino acid substitution on protein; Calculates Bayes posterior probability
Input: FASTA

[SIFT] (http://sift.jcvi.org/)
predicts how an amino acid substitution will affect protein function; Based on degree of conservation of amino acid residues- collected though PSI-BLAST; Standalone or web app program;
Input: Uniprot ID or Accession, Go term ID, Function name, Species Name or ID, etc

[snpEff] (http://snpeff.sourceforge.net/)
Genetic variant annotation and effect prediction toolbox; integrated with Galaxy, GATK, and GNKO; can annotate SNPs, INDELs, and multiple-nucleotide polymorphisms; categorizes effects into classes by functionality; Standalone or Web app;
Input: VCF, BED
Output: VCF (with new ANN field, also used in ANNOVAR and VEP), HTML summary files

[SnpSIFT] (http://snpeff.sourceforge.net/SnpSift.html)
Filter annotated files; Part of SnpEff main distribution; one variants have been annotated, this can be used to filter your data to find relevant variants

[VAAST 2] (http://www.yandell-lab.org/software/vaast.html)
Variant Annotation, Analysis, and Search Tool; probabilistic search tool for identifying damage genes and the disease causing variants; can score both coding and non-coding variants; Four tools: VAT (Variant annotation tool), VST (Variant Selection Tool), VAAST, pVAAST (for pedigree data); updated April 2015
Input: FASTA, GFF3, GVF
Output: CDR (condenser file), VAAST file (both unique to VAAST)

[VEP] (http://useast.ensembl.org/info/docs/tools/vep/index.html?redirect=no)
Ensembl Variant Effect Predictor; determines effect of variants on genes, transcripts, and protein sequence; uses SIFT and PolyPhen
Input: Coordinates of variants and nucleotide changes; whitespace- separated format, VCF, pileup, HGVS
Output: VCF, JSON, Statistics
