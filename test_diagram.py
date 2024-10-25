import matplotlib.pyplot as plt

# Example data for two lines
x1 = [1, 2, 3, 4, 5]
x2 = [1.5, 2.5, 3.5, 4.5, 5.5]  # Offset second line slightly

# Create a 1D plot (two lines along the x-axis)
plt.hlines(0, min(x1 + x2) - 1, max(x1 + x2) + 1, color='gray')  # Base line along the x-axis for reference
plt.plot(x1, [0] * len(x1), 'o-', color='b', label='Line 1')     # First line on x-axis
plt.plot(x2, [-.2] * len(x2), 'o-', color='g', label='Line 2')  # Second line slightly below the first

# Add arrows between each point in the first line
for i in range(len(x1) - 1):
    plt.annotate('', xy=(x1[i+1], 0), xytext=(x1[i], 0),
                 arrowprops=dict(arrowstyle='->', color='b', lw=1.5))

# Add arrows between each point in the second line
for i in range(len(x2) - 1):
    plt.annotate('', xy=(x2[i+1], -0.2), xytext=(x2[i], -0.2),
                 arrowprops=dict(arrowstyle='->', color='g', lw=1.5))

# Add x-axis label
plt.xlabel("X-axis")

# set custom y-axis limit
plt.ylim(-.5,.5)

# Remove y-axis and add legend
plt.gca().get_yaxis().set_visible(False)
plt.legend(loc="lower left")

# Save the plot as an image
plt.savefig("1d_graph_two_lines.png")  # Adjust the path as needed

# Display the plot
plt.show()