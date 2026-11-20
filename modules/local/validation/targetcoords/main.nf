// Genomic coordinates of the curated target genes of one genome, from Prodigal's headers.
// The ORF ids were confirmed once against independent annotation (assets/real_genome_validation/target_orfs.tsv).
process VALIDATION_TARGETCOORDS {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'biocontainers/python:3.11' }"

    input:
    tuple val(meta), path(proteins)   // osmotool annotate --keep_proteins: <genome>.prodigal.faa
    path(targets)

    output:
    tuple val(meta), path("${meta.id}.gene_regions.tsv"), emit: regions
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python

    script:
    """
    python3 <<'PYEOF'
    # Prodigal's headers carry each ORF's coordinates: >contig_n # start # end # strand # attrs.
    # target_orfs.tsv names the ORF of each target family in this genome.
    import csv

    genome = "${meta.id}"
    with open("${targets}") as fh:
        target_orfs = {r["family"]: r["orf_id"] for r in csv.DictReader(fh, delimiter="\\t")
                       if r["genome"] == genome}
    orf_to_family = {orf: fam for fam, orf in target_orfs.items()}

    found = {}
    with open("${proteins}") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            parts = [p.strip() for p in line[1:].split("#")]
            orf_id = parts[0].split()[0]
            if orf_id in orf_to_family:
                found[orf_to_family[orf_id]] = (orf_id.rsplit("_", 1)[0], int(parts[1]), int(parts[2]))

    missing = set(target_orfs) - set(found)
    if missing:
        print(f"WARNING: no coordinates for {sorted(missing)}: are these the proteins of the annotate run "
              f"that found them?")

    with open(f"{genome}.gene_regions.tsv", "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\\t")
        writer.writerow(["family", "contig", "start", "end"])
        for family, (contig, start, end) in sorted(found.items()):
            writer.writerow([family, contig, start, end])
            print(f"[{family}] {contig}:{start}-{end} ({end - start + 1}bp)")
    PYEOF
    """

    stub:
    """
    printf 'family\\tcontig\\tstart\\tend\\n' > ${meta.id}.gene_regions.tsv
    """
}
