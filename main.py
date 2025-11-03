import argparse
import os
from typing import List, Optional

from utils.data_structures import Config, TelomereSequence, TemplateSwitchData, AlignmentData
from data_io.fasta_reader import FastaReader
from analysis.alignment_strategy import AlignmentStrategy
from analysis.template_switching_strategy import TemplateSwitchingStrategy
from data_io.graph_alignment import GraphAlignment
from data_io.graph_template_switching import GraphTemplateSwitching
from data_io.template_switching_exporters import TemplateSwitchingPrint
from data_io.alignment_exporters import AlignmentPrint

from analysis.pattern_finder_lcp import pattern_finder_execute

from os.path import splitext

import sys


def run_analysis(config: Config) -> None:
    # Step 1: Read sequences
    reader = FastaReader(config.fasta_file_path, config.max_ends)
    telomers: List[Optional[TelomereSequence]]  = reader.parse_fasta()

    # step 2: find pattern
    pattern: Optional[str] = config.pattern
    if not pattern:
        pattern = pattern_finder_execute(telomers)
    if not pattern:
        raise ValueError("no good pattern was found")

    #sys.exit()

    # Step 3: Analyze sequences
    strategy = config.analysis_strategy
    if strategy == "template_switching":
        print("performing template switching analysis")
        analyzer = TemplateSwitchingStrategy(telomers, pattern, config)
        #template_switch_analysis: List[Optional[TemplateSwitchData]] = analyzer.execute()
        analyzer.execute()

        ## STATS: ##
        #base, ext = splitext(config.output_file)
        #stats_file_name = f"{base}_stats.txt"
        #stats_output_file = open(stats_file_name, 'w')
        #total = 0
        #for telomer in telomers:
            #if telomer and telomer.analysis:
                #total+=len(telomer.analysis.template_switch_event_indexes) - 1
                #print(f"{telomer.chromosome_end_id}\n{len(telomer.analysis.template_switch_event_indexes)-1}\n", file=stats_output_file)
        #print(f"{config.fasta_file_path}: {total}", file=stats_output_file)
                
    else:
        print("performing alignment analysis")
        analyzer = AlignmentStrategy(telomers, pattern, config)
        analyzer.execute()

    # step 4: Graph
    if config.graph_output:
        if strategy == "template_switching":
            grapher = GraphTemplateSwitching(telomers, pattern, config)
            grapher.execute()
        else:
            grapher = GraphAlignment(telomers, pattern, config)
            grapher.execute()
    

    # Step 5: Export results
    if strategy == "template_switching":
        template_switch_exporter = TemplateSwitchingPrint(telomers, config, pattern)
        template_switch_exporter.print_analysis()
    else:
        alignment_exporter = AlignmentPrint(telomers, config, pattern)
        alignment_exporter.print_analysis()
            

def main(args) -> None:
    config = Config(
        fasta_file_path=args.fasta_file_path,
        output_file=args.output,
        analysis_strategy=args.analysis_strategy,
        max_ends=args.maximum_ends,
        pattern=args.pattern,
        maximum_alignment_mutations=args.maximum_alignment_mutations,
        skip_seeding=args.skip_seeding,
        compare_file_path=args.compare_output,
        min_pattern_length=args.min_length,
        graph_output=args.graph_output
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

    # check if compare file path exists
    if config.compare_file_path and not os.path.isfile(config.compare_file_path):
        raise FileNotFoundError(f"File not found: {config.compare_file_path}")

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
    parser.add_argument("-ss", "--skip_seeding", type=bool, default=True,
                        help="determines whether we do a seeding process for alignments in mutagenic zone, or just try to do alignments in the whole mutagenic zone")
    parser.add_argument("--min_length", type=int, default=50,
                        help="The minimum length for a valid pattern (default: 50)")
    parser.add_argument("--max_length", type=int, default=300,
                        help="The maximum length for a valid pattern (default: 300)")
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
    parser.add_argument("-co", "--compare_output", type=str, default=None,
                        help="used for comparing the output of the program to Ivan's CSV files")
    parser.add_argument("-mam", "--maximum_alignment_mutations", type=int, default=12,
                        help="determines the cutoff for a valid mutation. Default is 1 mutations per 12 bps")
    args = parser.parse_args()

    main(args)
