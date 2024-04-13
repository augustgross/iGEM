import sys

MTBREFLENGTH = 4411532
REFGENOMEFILEPATH = "Mtub_H37Rv_NC000962.3.fasta"

def read_vcf_file(file_path):
    try:
        with open(file_path, 'r') as file:
            vcf_array = [line.strip() for line in file if not line.startswith("#")]
        return vcf_array
    except FileNotFoundError as e:
        print(f"File not found: {e}")
        return []

def split_vcf_strings(vcf_array):
    return [line.split("\t") for line in vcf_array]

def process_vcf_data(split_vcf_array):
    position = []
    ref = []
    mut = []
    for data in split_vcf_array:
        position.append(int(data[1]) - 1)  # Adjust for zero-indexing
        ref.append(data[3])
        mut.append(data[4].split(",")[0] if "," in data[4] else data[4])  # Handle multiple mutations
    return position, ref, mut

def calculate_total_offset(ref, mut):
    total_offset = sum(len(m) - len(r) for m, r in zip(mut, ref))
    print(f"totalOffset = {total_offset}")
    return total_offset

def read_reference_genome(file_path):
    try:
        with open(file_path, 'r') as file:
            # Skip the header line when concatenating the genome sequence.
            reference_genome = ''.join([line.strip() for line in file if not line.startswith(">")])
        return reference_genome
    except FileNotFoundError as e:
        print(f"File not found: {e}")
        return ""

def updated_reference_genome(strain, reference_genome, position, ref, mut):
    j = 0
    for k in range(len(position)):
        while j < position[k]:
            strain.append(reference_genome[j])
            j += 1
        for m in mut[k]:
            strain.append(m)
        j += len(ref[k])  # Skip over the reference sequence length in the genome
    strain.extend(reference_genome[j:])  # Append the remaining part of the reference genome
    return strain

def write_char_array_to_file(header, new_fasta, file_path):
    try:
        with open(file_path, 'w') as file:
            file.write(header + '\n')  # Write the header line
            for i, char in enumerate(new_fasta):
                file.write(char)
                if (i + 1) % 70 == 0:  # Formatting lines to be at most 70 characters long
                    file.write('\n')
        print("Data written to file successfully.")
    except IOError as e:
        print(f"Error writing to file: {e}")

def main(vcf_file_path, output_file_path):
    vcf_array = read_vcf_file(vcf_file_path)
    split_vcf_array = split_vcf_strings(vcf_array)
    position, ref, mut = process_vcf_data(split_vcf_array)
    total_offset = calculate_total_offset(ref, mut)
    print(f"Variation Offset Amount = {total_offset}")
    
    # Read the entire reference genome, including the header for the output file.
    try:
        with open(REFGENOMEFILEPATH, 'r') as file:
            header = next(file).strip()  # This assumes the first line is the header.
            reference_genome = ''.join(line.strip() for line in file if not line.startswith(">"))
    except FileNotFoundError as e:
        print(f"File not found: {e}")
        return
    
    strain = updated_reference_genome([], reference_genome, position, ref, mut)
    write_char_array_to_file(header, strain, output_file_path)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_vcf_file> <output_fasta_file>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])