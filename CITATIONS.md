# Citations

The tools and data a release of osmo_refdb is built with, and the references to
cite for them. Cite the ones your work depends on. Each release records the exact
tool versions it was built with in `PROVENANCE.json` and `software_versions.yml`.

UniProt, GMGC and the prevalence-filtered catalogs (Cayman) carry the CC BY 4.0
attribution the release's licence requires (see the README, "Licence and
citations"). The rest are the tools and databases the pipeline runs.

## Workflow and software distribution

- **Nextflow**: the workflow engine.
  Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. *Nat Biotechnol* 35(4):316-319 (2017). [doi:10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820)

- **nf-core**: the pipeline's template, and the institutional profiles from nf-core/configs.
  Ewels PA et al. The nf-core framework for community-curated bioinformatics pipelines. *Nat Biotechnol* 38(3):276-278 (2020). [doi:10.1038/s41587-020-0439-x](https://doi.org/10.1038/s41587-020-0439-x)

- **Bioconda and BioContainers**: the packages and container images the steps run in.
  Grüning B et al. Bioconda: sustainable and comprehensive software distribution for the life sciences. *Nat Methods* 15(7):475-476 (2018). [doi:10.1038/s41592-018-0046-7](https://doi.org/10.1038/s41592-018-0046-7)  
  da Veiga Leprevost F et al. BioContainers: an open-source and community-driven framework for software standardization. *Bioinformatics* 33(16):2580-2582 (2017). [doi:10.1093/bioinformatics/btx192](https://doi.org/10.1093/bioinformatics/btx192)

- **Singularity (Apptainer)**: the container runtime of the HPC profile.
  Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. *PLOS ONE* 12(5):e0177459 (2017). [doi:10.1371/journal.pone.0177459](https://doi.org/10.1371/journal.pone.0177459)

## Sequence search, alignment and annotation

- **DIAMOND**: protein search for the release's database, the cutoff calibration and the translated search of `osmotool profile`.
  Buchfink B, Xie C, Huson DH. Fast and sensitive protein alignment using DIAMOND. *Nat Methods* 12(1):59-60 (2015). [doi:10.1038/nmeth.3176](https://doi.org/10.1038/nmeth.3176)  
  Buchfink B, Reuter K, Drost H-G. Sensitive protein alignments at tree-of-life scale using DIAMOND. *Nat Methods* 18(4):366-368 (2021). [doi:10.1038/s41592-021-01101-x](https://doi.org/10.1038/s41592-021-01101-x)

- **HMMER**: for the family HMMs and their gathering cutoffs (hmmbuild, hmmpress, hmmsearch, hmmscan).
  Eddy SR. Accelerated Profile HMM Searches. *PLoS Comput Biol* 7(10):e1002195 (2011). [doi:10.1371/journal.pcbi.1002195](https://doi.org/10.1371/journal.pcbi.1002195)

- **MAFFT**: the alignments the HMMs are built from.
  Katoh K, Standley DM. MAFFT Multiple Sequence Alignment Software Version 7: Improvements in Performance and Usability. *Mol Biol Evol* 30(4):772-780 (2013). [doi:10.1093/molbev/mst010](https://doi.org/10.1093/molbev/mst010)

- **trimAl**: the trimming of those alignments.
  Capella-Gutiérrez S, Silla-Martínez JM, Gabaldón T. trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. *Bioinformatics* 25(15):1972-1973 (2009). [doi:10.1093/bioinformatics/btp348](https://doi.org/10.1093/bioinformatics/btp348)

- **CD-HIT**: clustering each family's positives and negatives at 90% identity, and the HMM decoys.
  Li W, Godzik A. CD-HIT: a fast program for clustering and comparing large sets of protein or nucleotide sequences. *Bioinformatics* 22(13):1658-1659 (2006). [doi:10.1093/bioinformatics/btl158](https://doi.org/10.1093/bioinformatics/btl158)  
  Fu L, Niu B, Zhu Z, Wu S, Li W. CD-HIT: accelerated for clustering the next-generation sequencing data. *Bioinformatics* 28(23):3150-3152 (2012). [doi:10.1093/bioinformatics/bts565](https://doi.org/10.1093/bioinformatics/bts565)

- **Prodigal**: gene prediction in genomes.
  Hyatt D, Chen G-L, LoCascio PF, Land ML, Larimer FW, Hauser LJ. Prodigal: prokaryotic gene recognition and translation initiation site identification. *BMC Bioinformatics* 11(1):119 (2010). [doi:10.1186/1471-2105-11-119](https://doi.org/10.1186/1471-2105-11-119)

- **OrfM**: open reading frames from reads in `osmotool profile`.
  Woodcroft BJ, Boyd JA, Tyson GW. OrfM: a fast open reading frame predictor for metagenomic data. *Bioinformatics* 32(17):2702-2703 (2016). [doi:10.1093/bioinformatics/btw241](https://doi.org/10.1093/bioinformatics/btw241)

## Read simulation and mapping

- **InSilicoSeq and wgsim**: the simulated reads the benchmarks score. wgsim has no publication: [github.com/lh3/wgsim](https://github.com/lh3/wgsim).
  Gourlé H, Karlsson-Lindsjö O, Hayer J, Bongcam-Rudloff E. Simulating Illumina metagenomic data with InSilicoSeq. *Bioinformatics* 35(3):521-522 (2019). [doi:10.1093/bioinformatics/bty630](https://doi.org/10.1093/bioinformatics/bty630)

- **BWA-MEM and bwa-mem3**: indexing and mapping reads to the gene catalog; bwa-mem3 is a faster implementation of BWA-MEM.
  Li H. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. *arXiv* 1303.3997 (2013). [arXiv:1303.3997](https://arxiv.org/abs/1303.3997)  
  Vasimuddin M, Misra S, Li H, Aluru S. Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems. In *2019 IEEE International Parallel and Distributed Processing Symposium (IPDPS)*, pp. 314-324 (2019). [doi:10.1109/IPDPS.2019.00041](https://doi.org/10.1109/IPDPS.2019.00041)

## Python libraries

- **NumPy, pandas, SciPy and statsmodels**: the analyses and the benchmark metrics.
  Harris CR et al. Array programming with NumPy. *Nature* 585(7825):357-362 (2020). [doi:10.1038/s41586-020-2649-2](https://doi.org/10.1038/s41586-020-2649-2)  
  McKinney W. Data Structures for Statistical Computing in Python. In *Proceedings of the 9th Python in Science Conference*, pp. 56-61 (2010). [doi:10.25080/Majora-92bf1922-00a](https://doi.org/10.25080/Majora-92bf1922-00a)  
  Virtanen P et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nat Methods* 17(3):261-272 (2020). [doi:10.1038/s41592-019-0686-2](https://doi.org/10.1038/s41592-019-0686-2)  
  Seabold S, Perktold J. Statsmodels: Econometric and Statistical Modeling with Python. In *Proceedings of the 9th Python in Science Conference*, pp. 92-96 (2010). [doi:10.25080/Majora-92bf1922-011](https://doi.org/10.25080/Majora-92bf1922-011)

- **scikit-learn and XGBoost**: the benchmark metrics, the BacDive analysis and the NaCl model; XGBoost is the model itself.
  Pedregosa F et al. Scikit-learn: machine learning in Python. *J Mach Learn Res* 12:2825-2830 (2011). [jmlr.org/papers/v12/pedregosa11a.html](https://jmlr.org/papers/v12/pedregosa11a.html)  
  Chen T, Guestrin C. XGBoost: a scalable tree boosting system. In *Proceedings of the 22nd ACM SIGKDD International Conference on Knowledge Discovery and Data Mining*, pp. 785-794 (2016). [doi:10.1145/2939672.2939785](https://doi.org/10.1145/2939672.2939785)

- **Matplotlib**: the benchmark and NaCl model figures.
  Hunter JD. Matplotlib: A 2D Graphics Environment. *Comput Sci Eng* 9(3):90-95 (2007). [doi:10.1109/MCSE.2007.55](https://doi.org/10.1109/MCSE.2007.55)

## Databases and data

- **UniProt**: the reference sequences, and through them the HMMs, the DIAMOND database and the benchmark sets.
  The UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2025. *Nucleic Acids Res* 53(D1):D609-D617 (2025). [doi:10.1093/nar/gkae1010](https://doi.org/10.1093/nar/gkae1010)

- **NCBI and RefSeq**: additional positives from RefSeq for the families that declare RefSeq gene symbols, and the genome assemblies of the BacDive strains from NCBI.
  O'Leary NA et al. Reference sequence (RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation. *Nucleic Acids Res* 44(D1):D733-D745 (2016). [doi:10.1093/nar/gkv1189](https://doi.org/10.1093/nar/gkv1189)  
  Sayers EW et al. Database resources of the National Center for Biotechnology Information in 2025. *Nucleic Acids Res* 53(D1):D20-D29 (2025). [doi:10.1093/nar/gkae979](https://doi.org/10.1093/nar/gkae979)

- **BacDive**: the strain phenotypes (NaCl maximum) behind the NaCl model.
  Schober I et al. BacDive in 2025: the core database for prokaryotic strain data. *Nucleic Acids Res* 53(D1):D748-D756 (2025). [doi:10.1093/nar/gkae959](https://doi.org/10.1093/nar/gkae959)

- **Pfam and InterPro**: curated Pfam models for some families, and the Pfam-A 39.0 families of the random-panel pool (Pfam is distributed through InterPro).
  Mistry J et al. Pfam: The protein families database in 2021. *Nucleic Acids Res* 49(D1):D412-D419 (2021). [doi:10.1093/nar/gkaa913](https://doi.org/10.1093/nar/gkaa913)  
  Blum M et al. InterPro: the protein sequence classification resource in 2025. *Nucleic Acids Res* 53(D1):D444-D456 (2025). [doi:10.1093/nar/gkae1082](https://doi.org/10.1093/nar/gkae1082)

- **GMGC**: the Global Microbial Gene Catalog, the source of the gene catalog.
  Coelho LP et al. Towards the biogeography of prokaryotic genes. *Nature* 601:252-256 (2022). [doi:10.1038/s41586-021-04233-4](https://doi.org/10.1038/s41586-021-04233-4)

- **Cayman and its prevalence-filtered catalogs**: the filtered GMGC gene set the catalog is taken from, and the catalog-based profiling approach and thresholds `osmotool profile --catalog` follows.
  Ducarmon QR et al. Cayman enables large-scale analysis of gut microbiome carbohydrate-active enzyme repertoires. *Nat Microbiol* 11(6):1739-1753 (2026). [doi:10.1038/s41564-026-02318-2](https://doi.org/10.1038/s41564-026-02318-2)  
  Ducarmon Q, Karcher N, Zeller G. Data for "Large-scale computational analyses of gut microbial CAZyme repertoires enabled by Cayman". Zenodo (2024). [zenodo.org/records/10473258](https://zenodo.org/records/10473258) (CC BY 4.0)

