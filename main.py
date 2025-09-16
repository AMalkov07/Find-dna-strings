import argparse
import os
from typing import List, Optional

from utils.data_structures import Config, TelomereSequence, TemplateSwitchData
from data_io.fasta_reader import FastaReader
from analysis.alignment_strategy import AlignmentStrategy
from analysis.template_switching_strategy import TemplateSwitchingStrategy
from data_io.exporters import TemplateSwitchingPrint


def run_analysis(config: Config) -> None:
    # Step 1: Read sequences
    reader = FastaReader(config.fasta_file_path, config.max_ends)
    telomers: List[TelomereSequence]  = reader.parse_fasta()

    # step 2: find pattern
    pattern:str = config.pattern
    if not pattern:
        raise ValueError("no pattern was inputed")

    # Step 3: Analyze sequences
    strategy = config.analysis_strategy
    if strategy == "template_switching":
        analyzer = TemplateSwitchingStrategy(telomers, pattern, config)
    else:
        analyzer = AlignmentStrategy(telomers, pattern, config)
    analysis: List[Optional[TemplateSwitchData]] = analyzer.execute()

    # Step 5: Export results
    exporter = TemplateSwitchingPrint(analysis, telomers, config, pattern)
    exporter.print_analysis()
    #exporter.export_results(pattern, results, graph_path)
            

'''
        
    # Step 2: Find pattern  
    pattern_finder = PatternFinder()
    pattern = pattern_finder.find_consensus_pattern(telomeres)
        
    # Step 4: Create graph
    graph_gen = GraphGenerator()
    graph_path = graph_gen.create_graph(results)
        
    # Step 5: Export results
    exporter = ResultExporter(self.config.output_file)
    exporter.export_results(pattern, results, graph_path)
    '''


def main(args) -> None:
    config = Config(
        fasta_file_path=args.fasta_file_path,
        output_file=args.output,
        analysis_strategy=args.analysis_strategy,
        max_ends=args.maximum_ends,
        pattern=args.pattern,
        maximum_alignment_mutations=args.maximum_alignment_mutations
    )

    # check if file exists
    if not os.path.isfile(config.fasta_file_path):
        raise FileNotFoundError(f"File not found: {config.fasta_file_path}")

    # check if files empty
    if os.path.getsize(config.fasta_file_path) == 0:
        raise ValueError("The FASTA file is empty.")

    # check if valid analysis strategy
    if config.analysis_strategy != "template_switching" and config.analysis_strategy != "alignment":
        raise ValueError("invalid strategy name")

    run_analysis(config)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="program used for finding repeating sequences inside of data, provided in the form of a .fasta file")
    parser.add_argument("fasta_file_path", 
                        help="The path to the file to be processed.")
    parser.add_argument("-o", "--output", 
                        help="Optional output file to save the content.")
    parser.add_argument("-as", "--analysis_strategy", default="template_switching",
                       help="template_switching strategy or alignment strategy for the analysis (default: template_switching)")
    parser.add_argument("--min_length", type=int, default=50,
                        help="The minimum length for a valid repeat sequence (default: 50)")
    parser.add_argument("--graph_dpi", type=int, default=300,
                        help="the dpi of the saved graph (default: 300)")
    parser.add_argument("-go", "--graph_output",
                        help="Optional output file that a graph of the output will be saved too. (if No output is given, no graph will be created)")
    parser.add_argument("-me", "--maximum_ends", type=int, default=32,
                        help="changes the maximum number of chr ends that are expected (default is 32)")
    parser.add_argument("-is", "--ignore_sorting", action="store_true",
                        help="stops the program from attempting to automatically sort the sequences from fasta file based on the sequence headers")
    parser.add_argument("-ia", "--ignore_alignment", action="store_true",
                        help="stops the program from attempting to find imperfecting alighments. This drastically decreases the runtime of the program")
    parser.add_argument("-or", "--original_reference",
                        help="Optional input file with the original reference. New chr ends will be compared against the original reference, and any identical prefixes will not be examined for repeat sequences")
    parser.add_argument("-sp", "--skip_prefix", type=int, default=0,
                        help="skips looking for repeating sequences for the specified number of base pairs in all chr's")
    parser.add_argument("-gt", "--graph_title",
                        help="optional input for graph title (must be used with --graph_output flag)")
    parser.add_argument("-p", "--pattern",
                        help="used to specify the exact pattern to look for instead of the programming trying to find the circle pattern automatically")
    parser.add_argument("-co", "--compare_output",
                        help="used for comparing the output of the program to Ivan's CSV files")
    parser.add_argument("-mam", "--maximum_alignment_mutations", type=int, default=12,
                        help="determines the cutoff for a valid mutation. Default is 1 mutations per 12 bps")
    args = parser.parse_args()

    main(args)
