import os 
import csv

# NOTE: Change this to be the correct output dir before running the script
output_dir = "output"

tsv_file = "GUNC.progenomes_2.1.maxCSS_level.tsv"

combined_output = "../../output.tsv"

subdirectories = [d for d in os.listdir(output_dir) if os.path.isdir(os.path.join(output_dir, d))]

for subdirectory in subdirectories:
    print("Current working directory:", os.getcwd())

    subdirectory_path = os.path.join(output_dir, subdirectory)

    os.chdir(subdirectory_path)
    print("Current working directory:", os.getcwd())

    if os.path.exists(tsv_file):
        with open(tsv_file, 'r') as input_file:
            tsv_reader = csv.reader(input_file, delimiter='\t')
            # Skip the first line (header)
            next(tsv_reader)
            # Get the second line
            second_line = next(tsv_reader)

        with open(combined_output, 'a') as output_file:
            tsv_writer = csv.writer(output_file, delimiter='\t')
            tsv_writer.writerow(second_line)

    print("Current working directory:", os.getcwd())

    # Change back to the output_dir
    current_directory = os.getcwd()
    new_directory = os.path.abspath(os.path.join(current_directory, '..', '..'))
    os.chdir(new_directory)




