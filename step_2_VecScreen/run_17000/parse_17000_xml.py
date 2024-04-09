from Bio.Blast import NCBIXML
import os

HOME = os.environ["HOME"]

folder_path = f"{HOME}/fsl_groups/grp_Bio465_GenomeWars/compute/genomeWars/step_2_VecScreen/run_17000/expected_xml_results"
output_filename = f"{HOME}/fsl_groups/grp_Bio465_GenomeWars/compute/genomeWars/step_2_VecScreen/run_17000/vecscreen_17000_output.tsv"

num_found = 0
count = 0

with open(output_filename, "w+") as output_file:
  output_file.write("query_file\t subject_name\t subject_genome_length\t identities\t gaps\t score\t e_value\t num_alignments\t query_start\t subject_start\t alignment_length\n")
  for file in os.listdir(folder_path):
    count += 1
    print(count, file)
    result_handle = open(f"{folder_path}/{file}", 'r')
    try:
      blast_records = NCBIXML.parse(result_handle)

      for blast_record in blast_records:
        if len(blast_record.alignments) > 0:
          num_found += 1
        for alignment in blast_record.alignments:
          for hsp in alignment.hsps:
            line = f"{file}\t {alignment.title}\t {alignment.length}\t"
            line += f"{hsp.identities}\t {hsp.gaps}\t {hsp.score}\t {hsp.expect}\t {hsp.num_alignments}\t"
            line += f"{hsp.query_start}\t {hsp.sbjct_start}\t {len(hsp.query)}"
            line += "\n"
            output_file.write(line)
    except:
      print("Error parsing:", file)

print(num_found)