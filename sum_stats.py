import glob
import re

# Collect all your *_stats.txt files
files = glob.glob("*Output_stats.txt")

# Initialize sums
zone_sum = 0
areas_sum = 0
alignments_sum = 0
ivans_sum = 0
perfect_sum = 0
compared_sum = 0
insertions_sum = 0
deletions_sum = 0
mismatches_sum = 0
given_insertions_sum = 0
given_deletions_sum = 0
given_mismatches_sum = 0

for fname in files:
    with open(fname) as f:
        lines = [line.strip() for line in f][:5]  # now we read first 6 lines

        zone_sum += int(lines[0].split(":")[1].strip())
        areas_sum += int(lines[1].split(":")[1].strip())
        alignments_sum += int(lines[2].split(":")[1].strip())

        # line 4 has two numbers ("X alignments matched perfectly out of Y that were compared")
        parts = lines[3].split()
        perfect_sum += int(parts[0])
        compared_sum += int(parts[6])

# Write the summed output
with open("Sum_stats.txt", "w") as out:
    out.write(f"total chr ends with mutagenic zone found: {zone_sum}\n")
    out.write(f"total chr ends with matching number of mutagenic areas: {areas_sum}\n")
    out.write(f"total chr ends with matching number of alignments: {alignments_sum}\n")
    out.write(f"{perfect_sum} alignments matched perfectly out of {compared_sum} that were compared\n")
