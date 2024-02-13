#!/usr/bin/python

def split_file_into_chunks(file_path, lines_per_file=5, control_file = 'download_cell_gene_control_file.txt'):
    with open(file_path, 'r') as source_file:
        lines = source_file.readlines()

    total_files = len(lines) // lines_per_file
    if len(lines) % lines_per_file != 0:
        total_files += 1  # Add an extra file if there are remaining lines

    out = []
    for i in range(total_files):
        chunk = lines[i*lines_per_file:(i+1)*lines_per_file]
        output_file_path = f"../data/collection_ids_{i+1}.txt"  # Adjust path as needed
        out.append(output_file_path)
        with open(output_file_path, 'w') as output_file:
            output_file.writelines(chunk)
        print(f"Created {output_file_path}")

    # Now, write the commands to a new file
    commands_file_path = control_file
    with open(commands_file_path, 'w') as commands_file:
        for path in out:
            command_line = f'python iter_download.py {path}\n'
            commands_file.write(command_line)
    
    print(f"Commands file '{commands_file_path}' has been created with all commands.")

# Assuming 'collection_ids.txt' is in the current working directory
source_file_path = '../data/collection_ids.txt'
split_file_into_chunks(source_file_path)