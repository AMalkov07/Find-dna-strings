import matplotlib.pyplot as plt

# Example data for three lines
x1 = [1, 2, 3, 4, 5]
x2 = [1.5, 2.5, 3.5, 4.5, 5.5]      # Offset second line slightly
x3 = [1, 2, 3, 4, 5]                # Third line next to the first
for i in range(len(x3)):
    x3[i] *= -1


# Create a 1D plot with three lines along the x-axis
#plt.hlines(0, min(x1 + x2 + x3) - 1, max(x1 + x2 + x3) + 1, color='gray')  # Base line along the x-axis for reference
plt.plot(x1, [0] * len(x1), 'o-', color='b', label='Line 1')               # First line on x-axis
plt.plot(x2, [-0.02] * len(x2), 'o-', color='g', label='Line 2')           # Second line slightly below the first
plt.plot(x3, [0] * len(x3), 'o-', color='r', label='Line 3')            # Third line slightly above the first

# Add arrows between each point in the first line (pointing right)
for i in range(len(x1) - 1):
    plt.annotate('', xy=(x1[i+1], 0), xytext=(x1[i], 0),
                 arrowprops=dict(arrowstyle='->', color='b', lw=1.5))

# Add arrows between each point in the second line (pointing right)
for i in range(len(x2) - 1):
    plt.annotate('', xy=(x2[i+1], -0.02), xytext=(x2[i], -0.02),
                 arrowprops=dict(arrowstyle='->', color='g', lw=1.5))

# Add arrows between each point in the third line (pointing left)
for i in range(len(x3) - 1):
    plt.annotate('', xy=(x3[i+1], 0), xytext=(x3[i], 0),
                 arrowprops=dict(arrowstyle='->', color='r', lw=1.5))

# Add text "chr1" at coordinates (0, 0)
plt.text(0, 0, "chr1", ha='center', va='center', fontsize=12, color='black')

# Add x-axis label
plt.xlabel("X-axis")

# Set custom y-axis limits to keep lines close
plt.ylim(-0.1, 0.1)

# Remove y-axis and add legend
plt.gca().get_yaxis().set_visible(False)
plt.legend(loc="upper left")

# Save the plot as an image
plt.savefig("1d_graph_with_three_lines.png")  # Adjust the path as needed

# Display the plot
plt.show()