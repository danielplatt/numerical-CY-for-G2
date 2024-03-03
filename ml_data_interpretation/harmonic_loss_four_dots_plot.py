import matplotlib.pyplot as plt
x = [1.27 * 10**(-3), 3.31 * 10**(-4), 4.75 * 10**(-5), 2.42 * 10**(-6)]
y = [3.30 * 10**(-2), 2.88 * 10**(-4), 2.25 * 10**(-2), 8.45 * 10**(-1)]
labels = ['Fermat', 'Random Quintic', 'CICY1', 'CICY2']

fig, ax = plt.subplots()
ax.scatter(x, y)

for i, txt in enumerate(labels):
    ax.annotate(txt, (x[i]+0.00002, y[i]+0.005))

plt.savefig(f'output_images/harmonic_loss_four_dots_plot.png')
plt.show()