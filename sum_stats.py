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
        lines = [line.strip() for line in f][:7]  # now we read first 6 lines

        zone_sum += int(lines[0].split(":")[1].strip())
        areas_sum += int(lines[1].split(":")[1].strip())
        alignments_sum += int(lines[2].split(":")[1].strip())
        ivans_sum += int(lines[3].split(":")[1].strip())

        # line 5 has two numbers ("X alignments matched perfectly out of Y that were compared")
        parts = lines[4].split()
        perfect_sum += int(parts[0])
        compared_sum += int(parts[6])

        # line 6 has insertions, deletions, mismatches
        nums = re.findall(r"\d+", lines[5])
        insertions_sum += int(nums[0])
        deletions_sum += int(nums[1])
        mismatches_sum += int(nums[2])

        # line 7 has given insertions, deletions, mismatches
        nums = re.findall(r"\d+", lines[6])
        given_insertions_sum += int(nums[0])
        given_deletions_sum += int(nums[1])
        given_mismatches_sum += int(nums[2])

# Write the summed output
with open("SumOutput_stats.txt", "w") as out:
    out.write(f"total chr ends with mutagenic zone found: {zone_sum}\n")
    out.write(f"total chr ends with matching number of mutagenic areas: {areas_sum}\n")
    out.write(f"total chr ends with matching number of alignments: {alignments_sum}\n")
    out.write(f"total chr ends with exact same call as Ivans data: {ivans_sum}\n")
    out.write(f"{perfect_sum} alignments matched perfectly out of {compared_sum} that were compared\n")
    out.write(f"total insertion count: {insertions_sum},  total deletions count: {deletions_sum}, total mismatches count: {mismatches_sum}\n")
    out.write(f"total Ivan insertion count: {given_insertions_sum},  total Ivan deletions count: {given_deletions_sum}, total Ivan mismatches count: {given_mismatches_sum}\n")
