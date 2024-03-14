#!/bin/bash


# List of folder names
names=("leucothrix", "thalassomonas", "velarivirus", "tritimovirus", "dinovernavirus", "bacillarnavirus", "rymovirus", "ignicoccus", "salinimicrobium")

# Iterate over each name
for name in "${names[@]}"; do
    # Convert the name to lowercase
    lowercase_name=$(echo "$name" | tr '[:upper:]' '[:lower:]')

    # Create the folder
    mkdir "$lowercase_name"

    echo "Created folder: $lowercase_name"
done
