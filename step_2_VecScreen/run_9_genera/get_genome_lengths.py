import os

HOME = os.environ["HOME"]

folder_path = f"{HOME}/fsl_groups/grp_Bio465_GenomeWars/compute/genomeWars/step_2_VecScreen/run_9_genera/expected_xml_results"
output_filename = f"{HOME}/fsl_groups/grp_Bio465_GenomeWars/compute/genomeWars/step_2_VecScreen/run_9_genera/vecscreen_9_genera_lengths.tsv"

count = 0

with open(output_filename, "w+") as output_file:
    output_file.write("query_file\t query_genome_length\n")
    for file in os.listdir(folder_path):
        count += 1
        print(f"{count} {file}")
        with open(f"{folder_path}/{file}", 'r') as xml_file:
            for line in xml_file.readlines():
                if line.startswith("  <BlastOutput_query-len>"):
                    items = line.split(">")
                    query_genome_length = items[1].replace("</BlastOutput_query-len", "")
                    new_line = f"{file}\t {query_genome_length}\n"
                    output_file.write(new_line)
