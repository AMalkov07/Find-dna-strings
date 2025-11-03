import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
from typing import List, Optional

from utils.data_structures import TelomereSequence, Config, AlignmentData, ImperfectAlignmentEvent

class GraphAlignment:
    def __init__(self, telomeres: List[Optional[TelomereSequence]], pattern: str, config: Config):
        self.fig, self.ax = plt.subplots(figsize=(16, 4))
        self.telomeres = telomeres
        self.pattern = pattern
        self.config = config
        self.perfect_arrow_distance = len(pattern) - 1  
        self.min_x = 0
        self.max_x = 0
        self._graph_setup()

    def _graph_setup(self):
        self.ax.text(0, 8.5, "chr1", ha='center',
                     va='center', fontsize=10, color='black')
        self.ax.text(0, 8, "chr2", ha='center', va='center',
                     fontsize=10, color='black')
        self.ax.text(0, 7.5, "chr3", ha='center',
                     va='center', fontsize=10, color='black')
        self.ax.text(0, 7, "chr4", ha='center', va='center',
                     fontsize=10, color='black')
        self.ax.text(0, 6.5, "chr5", ha='center',
                     va='center', fontsize=10, color='black')
        self.ax.text(0, 6, "chr6", ha='center', va='center',
                     fontsize=10, color='black')
        self.ax.text(0, 5.5, "chr7", ha='center',
                     va='center', fontsize=10, color='black')
        self.ax.text(0, 5, "chr8", ha='center', va='center',
                     fontsize=10, color='black')
        self.ax.text(0, 4.5, "chr9", ha='center',
                     va='center', fontsize=10, color='black')
        self.ax.text(0, 4, "chr10", ha='center',
                     va='center', fontsize=10, color='black')
        self.ax.text(0, 3.5, "chr11", ha='center',
                     va='center', fontsize=10, color='black')
        self.ax.text(0, 3, "chr12", ha='center',
                     va='center', fontsize=10, color='black')
        self.ax.text(0, 2.5, "chr13", ha='center',
                     va='center', fontsize=10, color='black')
        self.ax.text(0, 2, "chr14", ha='center',
                     va='center', fontsize=10, color='black')
        self.ax.text(0, 1.5, "chr15", ha='center',
                     va='center', fontsize=10, color='black')
        self.ax.text(0, 1, "chr16", ha='center',
                     va='center', fontsize=10, color='black')

        self.ax.set_ylim(.5, 9.5)

        # Remove x/y-axis and add legend
        self.ax.get_yaxis().set_visible(False)
        self.ax.get_xaxis().set_visible(False)

    def graph_perfect_arrow(self, x_axis: int, y_index: float, sign: int, padding: int) -> None:
        x_axis *= sign

        end_point = x_axis+self.perfect_arrow_distance * sign

        tail_length = self.perfect_arrow_distance * .6
        tail_width = .1
        head_length = self.perfect_arrow_distance * .4
        head_width = .2

        if sign == 1:
            bottom_left = x_axis
        else:
            bottom_left = x_axis - tail_length

        tail = Rectangle((bottom_left, y_index - tail_width / 2),
                         tail_length, tail_width, color="teal")
        self.ax.add_patch(tail)

        if sign == -1 and x_axis < self.min_x:
            self.min_x = x_axis
        elif x_axis > self.max_x:
            self.max_x = x_axis

        triangle = Polygon([[x_axis + (tail_length+head_length) * sign, y_index],
                            [x_axis+(tail_length) * sign,
                             y_index-head_width/2],
                            [x_axis+(tail_length) * sign, y_index+head_width/2]],
                           closed=True, color="teal")

        self.ax.add_patch(triangle)

        if self.previous_end_point is not None and abs(x_axis - self.previous_end_point) > padding + 1:
            self.ax.plot([self.previous_end_point+padding/2*sign, x_axis-padding/2*sign], 
                         [y_index, y_index], linewidth=2, color="red", linestyle='-')
        
        self.previous_end_point = end_point

    
    def graph_imperfect_arrow(self, imperfect_alignment: ImperfectAlignmentEvent, x_axis: int, y_index: float, sign: int, padding: int) -> None:
        x_axis *= sign

        insertions = imperfect_alignment.insertion_events
        deletions = imperfect_alignment.deletion_events

        curr_arrow_distance = self.perfect_arrow_distance + len(insertions) - len(deletions)

        end_point = x_axis+curr_arrow_distance * sign

        tail_length = curr_arrow_distance * .6
        tail_width = .1
        head_length = curr_arrow_distance * .4
        head_width = .2

        if sign == 1:
            bottom_left = x_axis
        else:
            bottom_left = x_axis - tail_length

        tail = Rectangle((bottom_left, y_index - tail_width / 2),
                         tail_length, tail_width, color="teal")
        self.ax.add_patch(tail)

        if sign == -1 and x_axis < self.min_x:
            self.min_x = x_axis
        elif x_axis > self.max_x:
            self.max_x = x_axis

        triangle = Polygon([[x_axis + (tail_length+head_length) * sign, y_index],
                            [x_axis+(tail_length) * sign,
                             y_index-head_width/2],
                            [x_axis+(tail_length) * sign, y_index+head_width/2]],
                           closed=True, color="teal")

        self.ax.add_patch(triangle)

        for ins in insertions:
            index = ins[0]
            self.ax.plot([x_axis+index*sign, x_axis+index*sign],
                         [y_index-.1, y_index+.1], color='gold', linestyle='-', lw=1)

        for dele in deletions:
            index = dele[0]
            self.ax.plot([x_axis+index*sign, x_axis+index*sign],
                         [y_index-.1, y_index+.1], color='blue', linestyle='-', lw=1)

        for mis in imperfect_alignment.mismatch_events:
            index = mis[0]
            self.ax.plot([x_axis+index*sign, x_axis+index*sign],
                         [y_index-.1, y_index+.1], color='red', linestyle='-', lw=1)

        if self.previous_end_point is not None and abs(x_axis - self.previous_end_point) > padding + 1:
            self.ax.plot([self.previous_end_point+padding/2*sign, x_axis-padding/2*sign], 
                         [y_index, y_index], linewidth=2, color="red", linestyle='-')
        
        self.previous_end_point = end_point



    def graph_alignment_analysis(self, telomere: TelomereSequence, chr_index: int) -> None:
        if not isinstance(telomere.analysis, AlignmentData):
            raise ValueError("calling graph_alignment_analysis without telomer analysis value being of type AlignmentData")

        sign = chr_index % 2 - 1
        if sign == 0:
            sign += 1
        
        offset = 300
        y_index = (17 - (chr_index//2)/2)
        extra_padding = self.perfect_arrow_distance//2 #mess around with this value

        perfect_alignment_indexes: List[int] = sorted(telomere.analysis.perfect_alignments)
        imperfect_alignments: List[ImperfectAlignmentEvent] = sorted(telomere.analysis.imperfect_alignments, key=lambda x: x.full_telomer_start_index if x.full_telomer_start_index else -1)

        perfect_alignment_n = len(perfect_alignment_indexes)
        imperfect_alignments_n = len(imperfect_alignments)

        perfect_alignment_ptr = 0
        imperfect_alignments_ptr = 0

        counter = 0
        while perfect_alignment_ptr < perfect_alignment_n and imperfect_alignments_ptr < imperfect_alignments_n:
            if perfect_alignment_indexes[perfect_alignment_ptr] < imperfect_alignments[imperfect_alignments_ptr].full_telomer_start_index:
                self.graph_perfect_arrow(perfect_alignment_indexes[perfect_alignment_ptr] + offset + (extra_padding * counter), y_index, sign, extra_padding)
                perfect_alignment_ptr += 1
            else:
                self.graph_imperfect_arrow(imperfect_alignments[imperfect_alignments_ptr], imperfect_alignments[imperfect_alignments_ptr].full_telomer_start_index + offset + (extra_padding * counter), y_index, sign, extra_padding)
                imperfect_alignments_ptr += 1

            counter += 1

        while perfect_alignment_ptr < perfect_alignment_n:
            self.graph_perfect_arrow(perfect_alignment_indexes[perfect_alignment_ptr] + offset + (extra_padding * counter), y_index, sign, extra_padding)
            perfect_alignment_ptr += 1

            counter += 1
        
        while imperfect_alignments_ptr < imperfect_alignments_n:
            self.graph_imperfect_arrow(imperfect_alignments[imperfect_alignments_ptr], imperfect_alignments[imperfect_alignments_ptr].full_telomer_start_index + offset + (extra_padding * counter), y_index, sign, extra_padding)
            imperfect_alignments_ptr += 1

            counter += 1

    def save_graph(self):
        max_x = self.max_x + self.perfect_arrow_distance * 1.5
        min_x = self.min_x - self.perfect_arrow_distance * 1.5
            
        self.ax.set_xlim(min_x, max_x)
        if isinstance(self.config.graph_output, str):
            self.fig.savefig(self.config.graph_output, dpi=300, bbox_inches="tight")
        else:
            raise ValueError("for some reason our config graph_output value is None")

    def execute(self):
        for i, telomere in enumerate(self.telomeres):
            self.previous_end_point = None
            if telomere and telomere.analysis:
                self.graph_alignment_analysis(telomere, i+34)
        self.save_graph()