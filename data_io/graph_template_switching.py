
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
from typing import List, Optional

from utils.data_structures import TelomereSequence, Config, TemplateSwitchData, TemplateSwitchEvent

class GraphTemplateSwitching:
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

    def get_circular_jump_distance(self, i, j) -> int:

        n = len(self.pattern)

        forward = (j - i) % n       # steps forward (wraps around)
        backward = (i - j) % n      # steps backward (wraps around)

        if forward <= backward:
            return forward
        else:
            return (-1 * backward)
        


    def graph_alignment_analysis(self, telomere: TelomereSequence, chr_index: int) -> None:
        if not isinstance(telomere.analysis, TemplateSwitchData):
            raise ValueError("calling graph_alignment_analysis without telomer analysis value being of type TemplateSwitchData")

        offset = 300

        sign = chr_index % 2 - 1
        if sign == 0:
            sign += 1
        
        y_index = (17 - (chr_index//2)/2)

        template_swtiches = telomere.analysis.template_switch_events

        last_end = None
        last_direction = "down"

        for curr_template_switch in template_swtiches:
            #self.ax.plot([(curr_template_switch.telomer_start+offset)*sign, (curr_template_switch.telomer_end+offset)*sign],
                         #[y_index, y_index], linewidth=2, color="teal", linestyle="-")

            line_start = (curr_template_switch.telomer_start+offset) * sign
            line_end = (curr_template_switch.telomer_end+offset) * sign

            if curr_template_switch.is_mutation:
                vert_line_pos = line_start
                self.ax.plot([vert_line_pos, vert_line_pos],
                            [y_index-.1, y_index+.1], color='orange', linestyle='-', lw=1)
                

            else:
                self.ax.plot([line_start, line_end],
                            [y_index, y_index], linewidth=2, color="teal", linestyle="-")

                for i in range(0, len(curr_template_switch.telomer_chunk)//self.perfect_arrow_distance):
                    vert_line_pos = line_start + ((i+1) * self.perfect_arrow_distance) * sign
                    self.ax.plot([vert_line_pos, vert_line_pos],
                                [y_index-.1, y_index+.1], color='blue', linestyle='-', lw=1)

                for ins in curr_template_switch.insertion_events:
                    index = ins[0]
                    self.ax.plot([line_start+index*sign, line_start+index*sign],
                                 [y_index-.1, y_index+.1], color='gold', linestyle='-', lw=1)
                for dele in curr_template_switch.deletion_events:
                    index = dele[0]
                    self.ax.plot([line_start+index*sign, line_start+index*sign],
                                 [y_index-.1, y_index+.1], color='purple', linestyle='-', lw=1)
                for mis in curr_template_switch.mismatch_events:
                    index = mis[0]
                    self.ax.plot([line_start+index*sign, line_start+index*sign],
                                 [y_index-.1, y_index+.1], color='red', linestyle='-', lw=1)


                if sign -1 and line_end < self.min_x:
                    self.min_x = line_end
                elif line_end > self.max_x:
                    self.max_x = line_end

            if last_end is not None and isinstance(last_end, int) and isinstance(curr_template_switch.pattern_start, int):
                jump_distance = self.get_circular_jump_distance(last_end, curr_template_switch.pattern_start)
                if jump_distance < 0:
                    str_jump_distance = str(jump_distance)
                else:
                    str_jump_distance = "+" + str(jump_distance)
                if last_direction == "down":
                    curr_y_index = y_index + .15
                    last_direction = "up"
                else:
                    curr_y_index = y_index - .17
                    last_direction = "down"
                if sign == 1:
                    curr_line_start = line_start - 45
                else:
                    curr_line_start = line_start - 5
                self.ax.text(curr_line_start, curr_y_index, str_jump_distance, va="center", fontsize=3.5)
            
            last_end = curr_template_switch.pattern_end
            offset+=40

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