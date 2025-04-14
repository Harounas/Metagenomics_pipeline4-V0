import pandas as pd
import random
from collections import defaultdict
import plotly.express as px
import plotly.io as pio
import os
from .trimmomatic import run_trimmomatic
from .metaspades import run_spades
from .bowtie2 import run_bowtie2
from .kraken2 import run_kraken2
import distinctipy
import numpy as np
import matplotlib.pyplot as plt
import logging

def process_sample(forward, reverse, base_name, bowtie2_index, kraken_db, output_dir, threads, run_bowtie, use_precomputed_reports, use_assembly):
    """
    Processes a metagenomic sample with optional assembly.
    Checks for existing outputs before running each step.
    """
    try:
        if use_precomputed_reports:
            kraken_report = os.path.join(output_dir, f"{base_name}_kraken_report.txt")
            if not os.path.exists(kraken_report):
                raise FileNotFoundError(f"Precomputed Kraken2 report not found: {kraken_report}")
            return kraken_report

        # Define expected output paths
        trimmed_forward = os.path.join(output_dir, f"{base_name}_1_trimmed_paired.fq.gz")
        trimmed_reverse = os.path.join(output_dir, f"{base_name}_2_trimmed_paired.fq.gz")
        bowtie_unmapped_r1 = os.path.join(output_dir, f"{base_name}_1_unmapped.fq.gz")
        bowtie_unmapped_r2 = os.path.join(output_dir, f"{base_name}_2_unmapped.fq.gz")

        # Step 1: Run Trimmomatic
        if not (os.path.exists(trimmed_forward) and os.path.exists(trimmed_reverse)):
            logging.info(f"Running Trimmomatic for sample {base_name}")
            trimmed_forward, trimmed_reverse = run_trimmomatic(
                forward, reverse, base_name, output_dir, threads
            )
        else:
            logging.info(f"Using existing trimmed files for sample {base_name}")

        # Step 2: Optional host depletion
        if run_bowtie:
            if not (os.path.exists(bowtie_unmapped_r1) and os.path.exists(bowtie_unmapped_r2)):
                logging.info(f"Running Bowtie2 host depletion for sample {base_name}")
                unmapped_r1, unmapped_r2 = run_bowtie2(
                    trimmed_forward, trimmed_reverse, base_name, bowtie2_index, output_dir, threads
                )
            else:
                logging.info(f"Using existing Bowtie2 output files for sample {base_name}")
                unmapped_r1, unmapped_r2 = bowtie_unmapped_r1, bowtie_unmapped_r2
        else:
            unmapped_r1, unmapped_r2 = trimmed_forward, trimmed_reverse

        # Step 3: Perform assembly if requested
        if use_assembly:
            contigs_file = os.path.join(output_dir, f"{base_name}_contigs.fasta")
            if not os.path.exists(contigs_file):
                logging.info(f"Performing de novo assembly for sample {base_name}")
                contigs_file = run_spades(unmapped_r1, unmapped_r2, base_name, output_dir, threads)

            kraken_input = contigs_file
        else:
            kraken_input = (unmapped_r1, unmapped_r2)
            logging.info(f"Running Kraken2 directly on reads for sample {base_name}")

        # Step 4: Run Kraken2
        kraken_report = os.path.join(output_dir, f"{base_name}_kraken_report.txt")
        if not os.path.exists(kraken_report):
            kraken_report = run_kraken2(
                kraken_input[0], kraken_input[1] if isinstance(kraken_input, tuple) else None,
                base_name, kraken_db, output_dir, threads
            )
        else:
            logging.info(f"Using existing Kraken2 report for sample {base_name}")

        return kraken_report

    except Exception as e:
        logging.error(f"Error processing sample {base_name}: {e}")
        return None

def generate_sample_ids_csv(kraken_dir):
    """
    Generates a CSV file containing sample IDs extracted from Kraken report filenames.

    Parameters:
    - kraken_dir (str): Path to the directory containing Kraken report files.

    Returns:
    - str: Path to the generated sample_ids.csv file.
    """
    try:
        sample_ids = []
        for fname in os.listdir(kraken_dir):
            if fname.endswith('_kraken_report.txt'):
                sample_id = fname.replace('_kraken_report.txt', '')
                sample_ids.append(sample_id)

        sample_ids_df = pd.DataFrame({'Sample_ID': sample_ids})
        csv_path = os.path.join(kraken_dir, 'sample_ids.csv')
        sample_ids_df.to_csv(csv_path, index=False)
        logging.info(f"Sample IDs written to {csv_path}")
        return csv_path

    except Exception as e:
        logging.error(f"Error generating sample IDs CSV: {e}")
        return None
