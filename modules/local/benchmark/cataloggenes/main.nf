// The catalog benchmark's two gene sets: the release's proteins back-translated exactly as
// BENCHMARK_SIMULATE builds its read constructs (one fixed codon per amino acid), so a test read
// matches its own gene base for base.
//   A  catalogued      every sequence of the release (train, calib and test positives, both negative
//                      splits, the random background): every test read's own gene is in it
//   B  not catalogued  only what the benchmark does not score (positive.train, positive.calib,
//                      negative.train): CD-HIT at 90% keeps each <= ~90% identical to a test gene,
//                      below profile --catalog's 97% bar
// Checks first that the test constructs' CDS are the back-translated test proteins.
process BENCHMARK_CATALOGGENES {
    // calls: python + pyyaml (simulate_reads.py's protein_to_dna)
    // Seqera container 'simulate' (environment.yml next to this file: its packages)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/insilicoseq_wgsim_python_pyyaml:a4bccdbc78a6f5da' :
        'community.wave.seqera.io/library/insilicoseq_wgsim_python_pyyaml:e5b214e4c49c8d0d' }"
    label 'process_single'
    input:
    path(refs, stageAs: 'refs/*')
    path(families)
    path(reads)          // BENCHMARK_SIMULATE_TEST: the test reads and their constructs

    output:
    tuple val('A'), path('genes_A.fna'), emit: genes_a
    tuple val('B'), path('genes_B.fna'), emit: genes_b
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python

    script:
    """
    python3 - <<'PYEOF'
    import sys
    from pathlib import Path
    sys.path.insert(0, "${moduleDir}/../simulate/templates")
    import simulate_reads as m09

    refs, reads = Path("refs"), Path("${reads}")
    families = m09.load_family_names(Path("${families}"))
    splits = {"A": ["positive.train", "positive.calib", "positive.test", "negative.train", "negative.test"],
              "B": ["positive.train", "positive.calib", "negative.train"]}

    checked = 0
    for fam in families:
        ctx_file = reads / "genomic_context" / f"{fam}.positive.fasta"
        if not ctx_file.exists():
            continue
        ctx = dict(m09.parse_fasta(ctx_file))
        for header, protein in m09.parse_fasta(refs / f"{fam}.positive.test.faa")[:20]:
            cds = ctx[header][m09.FLANK_LEN:len(ctx[header]) - m09.FLANK_LEN]
            if cds != m09.protein_to_dna(protein):
                sys.exit(f"ERROR: the construct of {header} is not protein_to_dna of its protein")
            checked += 1
        if checked >= 100:
            break
    print(f"constructs checked against protein_to_dna: {checked}")

    for scen, kinds in splits.items():
        files = [refs / f"{fam}.{k}.faa" for fam in families for k in kinds]
        if scen == "A":
            files.append(refs / m09.RANDOM_BG_FASTA)
        n = 0
        with open(f"genes_{scen}.fna", "w") as fh:
            for f in files:
                # IDs carry the source file, so a protein in two files stays two genes
                for header, protein in m09.parse_fasta(f):
                    fh.write(f">{f.name[:-4]}|{header.split()[0]}\\n{m09.protein_to_dna(protein)}\\n")
                    n += 1
        print(f"scenario {scen}: {n} catalog genes")
    PYEOF
    """

    stub:
    """
    printf '>a\\nATG\\n' > genes_A.fna
    printf '>b\\nATG\\n' > genes_B.fna
    """
}
