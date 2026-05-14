import re
from collections import defaultdict

def load_slop_map(slop_bed):
    mp = defaultdict(list)
    with open(slop_bed, "r") as f:
        for line in f:
            chrom, start, end, name, score, strand = line.rstrip("\n").split("\t")
            key = f"{start}-{end}-{strand}"
            mp[key].append((chrom, start, end, name, strand))
    return mp

def dreg_windows_to_bed(dreg_file, slop_bed, out_bed):
    seq_re = re.compile(r'^# Sequence:\s*(\d+)-(\d+)\(([+-])\)')
    mp = load_slop_map(slop_bed)

    with open(dreg_file, "r") as infile, open(out_bed, "w") as outfile:
        for line in infile:
            m = seq_re.match(line)
            if not m:
                continue
            start, end, strand = m.groups()
            key = f"{start}-{end}-{strand}"
            for chrom, start, end, name, strand in mp.get(key, []):
                outfile.write(f"{chrom}\t{start}\t{end}\t{name}\t.\t{strand}\n")

if __name__ == "__main__":
    dreg_windows_to_bed("6nrf1_hits.dreg", "4bedtools.slop.output.bed", "7nrf1_hits.bed")