import sys

def convert_to_bed(input_path,output_path):
    with open(input_path,mode='r') as infile,open(output_path,mode='w') as outfile:
        header = infile.readline().strip().split("\t")
        chrom_idx = header.index("chromosome_name")
        start_idx = header.index("transcript_start")
        gene_idx = header.index('external_gene_name')
        strand_idx = header.index('strand')
        valid_chromosomes = {str(i) for i in range(1, 21)}
        for line in infile:
            fields = line.strip().split("\t")
            chrom = fields[chrom_idx]
            start = fields[start_idx]
            gene_name = fields[gene_idx]
            strand = fields[strand_idx]
            end = str(int(start) + 1)
            if chrom in valid_chromosomes:
                bed_line = f"chr{chrom}\t{start}\t{end}\tchr{chrom}@{start}-{end}|{gene_name}\t.\t{"-" if strand == "-1" else "+"}\n"
                outfile.write(bed_line)
        

if __name__ == "__main__":
    convert_to_bed("1human_gene_annotation.tsv", "3output.bed")
