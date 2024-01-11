import matplotlib.pyplot as plt
import numpy as np

# Load data
data = np.loadtxt("plot.dat")
data1 = np.loadtxt("numeric.csv",
                 delimiter=",", skiprows=1)
                 
max_value = data1[:,3].max()


# Increase figure size
fig, ax = plt.subplots(figsize=(3, 3))

# Plot data with customizations
ax.plot(data[:, 0], data[:, 1], linestyle="-", label='Analytical', color='red')
ax.plot(data1[:, 5], data1[:, 3]/max_value, linestyle="-", label='Numeric', color='blue')

plt.legend()

# Add title and labels
plt.xlabel(r"$y$")
plt.ylabel(r"$U^*$")

# Save the plot in high quality as a PNG file
plt.savefig("channel.png", dpi=600, bbox_inches='tight')

