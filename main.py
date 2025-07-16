'''
Main file which works as a pipeline for the parsing, processing and analysing of the CRISPR CAS database.

author: U.B.
'''

import os
import subprocess
from crispr_cas_db.db_parser.CrisprDBParser import CrisprDBParser
from crispr_cas_db.processing.CrisprArrays import CrisprArrays
from crispr_cas_db.processing.CrisprRNAs import CrisprRNAs
from crispr_cas_evaluation.analysis.RNAPredictionVisualizer import CRISPRRNAPredictionVisualizer
from crispr_cas_evaluation.analysis.Analyzer import AnalyzerConfig

def main() -> None:
    '''Main function parsing, processing and analysing the CRISPR CAS database'''
    _parse_process_db()
    _remove_old_predictions()
    _prediction_repeats_rnamotifold()
    _prediction_crRNAs_rnamotifold()
    _prediction_repeats_rnamotices()
    _prediction_crRNAs_rnamotices()
    _repeat_analysis()
    _crRNA_analysis()

def _parse_process_db() -> None:
    '''Parses and processes the CRISPR CAS database and stores the sequences to be predicted as fasta-files'''
    parser = CrisprDBParser()
    parser.process_sql_file()
    crispr_arrays = CrisprArrays()
    crispr_arrays.save_repeats_fasta("repeats")
    crispr_RNAs = CrisprRNAs(crispr_arrays)
    crispr_RNAs.save_crRNA_fasta("crRNAs")

def _remove_old_predictions() -> None:
    '''Checks if the "prediction_files" folder exists.
    Deltes the content if it does or creates it if it does not exist.'''
    folder_path = "crispr_cas_evaluation/prediction_files"
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    else:
        for name in os.listdir(folder_path):
            path = os.path.join(folder_path, name)
            if os.path.isfile(path) or os.path.islink(path):
                os.remove(path)

def _prediction_repeats_rnamotifold() -> None:
    '''Executes the prediction via RNAmotiFold for the repeats'''
    script_folder = "./RNAmotiFold"
    cmd = [
        "python3",
        "RNAmotiFold.py",
        "-i", "../crispr_cas_db/fasta_files/repeats.fasta",
        "-o", "../crispr_cas_evaluation/prediction_files/repeats_rnamotifold.csv",
        "-a", "rnamotifold",
        "-s", "--no_update"
    ]
    subprocess.run(cmd, cwd=script_folder, check=True)

def _prediction_crRNAs_rnamotifold() -> None:
    '''Executes the prediction via RNAmotiFold for the crRNAs'''
    script_folder = "./RNAmotiFold"
    cmd = [
        "python3",
        "RNAmotiFold.py",
        "-i", "../crispr_cas_db/fasta_files/crRNAs.fasta",
        "-o", "../crispr_cas_evaluation/prediction_files/crRNAs_rnamotifold.csv",
        "-a", "rnamotifold",
        "-s", "--no_update"
    ]
    subprocess.run(cmd, cwd=script_folder, check=True)

def _prediction_repeats_rnamotices() -> None:
    '''Executes the prediction via RNAmotiCes for the repeats'''
    script_folder = "./RNAmotiFold"
    cmd = [
        "python3",
        "RNAmotiFold.py",
        "-i", "../crispr_cas_db/fasta_files/repeats.fasta",
        "-o", "../crispr_cas_evaluation/prediction_files/repeats_rnamotices.csv",
        "-a", "rnamotices", "--no_update"
    ]
    subprocess.run(cmd, cwd=script_folder, check=True)

def _prediction_crRNAs_rnamotices() -> None:
    '''Executes the prediction via RNAmotiCes for the crRNAs'''
    script_folder = "./RNAmotiFold"
    cmd = [
        "python3",
        "RNAmotiFold.py",
        "-i", "../crispr_cas_db/fasta_files/crRNAs.fasta",
        "-o", "../crispr_cas_evaluation/prediction_files/crRNAs_rnamotices.csv",
        "-a", "rnamotices", "--no_update"
    ]
    subprocess.run(cmd, cwd=script_folder, check=True)

def _repeat_analysis() -> None:
    '''Visualizes the results from the CRISPR Repeat analysis.'''
    config = AnalyzerConfig(
            fasta_path="./crispr_cas_db/fasta_files/repeats.fasta",
            motifold_csv_path="./crispr_cas_evaluation/prediction_files/repeats_rnamotifold.csv",
            motices_csv_path="./crispr_cas_evaluation/prediction_files/repeats_rnamotices.csv"
        )
    visualizer = CRISPRRNAPredictionVisualizer(config, "Repeats")
    visualizer.visualize_all_data()
    visualizer.visualize_subtypes()
    visualizer.visualize_heatmaps()

def _crRNA_analysis() -> None:
    '''Visualizes the results from the mature CRISPR RNA analysis.'''
    config = AnalyzerConfig(
            fasta_path="./crispr_cas_db/fasta_files/crRNAs.fasta",
            motifold_csv_path="./crispr_cas_evaluation/prediction_files/crRNAs_rnamotifold.csv",
            motices_csv_path="./crispr_cas_evaluation/prediction_files/crRNAs_rnamotices.csv"
        )
    visualizer = CRISPRRNAPredictionVisualizer(config, "CRISPR RNA")
    visualizer.visualize_all_data()
    visualizer.visualize_subtypes()
    visualizer.visualize_heatmaps()

if __name__ == "__main__":
    main()