import sys, os

def parse_fasta(filename, chunk_size):
  chunks = []
  with open(filename, 'r') as file:
    header = None
    sequence = ''
    for line in file:
      line = line.strip()
      if line.startswith('>'):
        if header:
          chunks.append((header, sequence))
          sequence = ''
        header = line[1:]
      else:
        sequence += line
  if header and sequence:
    chunks.append((header, sequence))
  print("> Sequences processed.")

  # Split chunks into smaller chunks
  if chunk_size:
    small_chunks = []
    for header, sequence in chunks:
      for i in range(0, len(sequence), chunk_size):
        small_chunks.append((header, sequence[i:i + chunk_size]))
    return small_chunks
  else:
    return chunks

def export_chunks(chunks, file_path):
  curr_header = None
  count = 0

  with open(file_path, "w+") as file:
    for header, sequence in chunks:
      count = count + 1 if curr_header == header else 0
      curr_header = header
      # file_path = f"{folder_path}/{header.replace(' ', '_')}_{count}"
      file.write(f">{header} ({count})\n{sequence}\n")
      # export_fasta(f"{header} ({count}/{len(chunks)})", sequence, file_path)

def export_fasta(header, sequence, file_path):
  with open(file_path, "w+") as file:
    file.write(f">{header}\n{sequence}")

def main():
  if len(sys.argv) != 4:
    print("Usage: python script.py <path_to_FASTA_file> <path_to_chunked_file> <chunk_size>")
    sys.exit(1)

  fasta_file_path = sys.argv[1]
  chunked_file_path = sys.argv[2]
  chunk_size = int(sys.argv[3])

  chunks = parse_fasta(fasta_file_path, chunk_size)
  print("> Writing files...")
  export_chunks(chunks, chunked_file_path)
  print("Done!")

if __name__ == "__main__":
  main()
