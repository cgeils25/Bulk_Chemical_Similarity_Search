# Run the entire pipeline at once

export START_TIME=$(python -c 'import neattime; print(neattime.neattime())') # will be used as a timestamp for the output directories

# download every compound available on PubChem
python -u download_pubchem_compounds.py \
    --output_dir "pubchem_data/full_download_"${START_TIME}

# extract necessary data from PubChem's SDF files and save them
python -u extract_data_from_pubchem_sdf.py \
    --input_dir pubchem_data/full_download_${START_TIME}
    --output_dir extracted_pubchem_data/full_processed_data_${START_TIME}
    --num_processes -1 # use all available cores

# compute the Tanimoto similarity between the extracted data and the comparison dataset
python -u compute_tanimoto_similarity.py \
    --comparison_dataset all_organic_smi.csv \
    --extracted_pubchem_data_dir extracted_pubchem_data/full_processed_data_${START_TIME} \
    --output_dir tanimoto_similarity_results/full_tanimoto_${START_TIME} \
    --num_processes 16
# substitute all_organic_smi.csv with the path to your comparison dataset, and replace 16 with the number of processes you want to use. This
# step is very memory intensive so adjust the number of processes according to your system's capabilities

# filter the compounds based on their tanimoto similarity to the a threshold
python -u filter_tanimoto_results.py \
    --tanimoto_directory tanimoto_similarity_results/full_tanimoto_${START_TIME} \
    --threshold 0.8

echo "Full pipeline done"
