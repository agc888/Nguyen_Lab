#### USEFUL PYTHON FUNCTIONS ####


###################### PLOTTING #######################


## Print Colour Palatte Dictionary ##
import matplotlib.pyplot as plt
from math import ceil

def show_palette(palette_dict):
      
    if len(palette_dict) >= 10:
        num_rows = ceil(len(palette) / 5)

        # Create a bar plot to visualize the colors
        fig, axs = plt.subplots(num_rows, 5, figsize=(20, num_rows * 4), constrained_layout=True)

        for i, (label, color) in enumerate(palette.items()):
            row_idx, col_idx = divmod(i, 5)
            axs[row_idx, col_idx].fill_between([0, 1], 0, 1, color=color, label=label)
            axs[row_idx, col_idx].text(0.5, 0.5, color, color='black', ha='center', va='center', fontsize=8)

            axs[row_idx, col_idx].set_xticks([])
            axs[row_idx, col_idx].set_yticks([])
            axs[row_idx, col_idx].spines['top'].set_visible(False)
            axs[row_idx, col_idx].spines['right'].set_visible(False)
            axs[row_idx, col_idx].spines['bottom'].set_visible(False)
            axs[row_idx, col_idx].spines['left'].set_visible(False)

        # Remove empty subplots
        for i in range(len(palette), num_rows * 5):
            row_idx, col_idx = divmod(i, 5)
            fig.delaxes(axs[row_idx, col_idx])

        # Display the plot
        plt.show()
    else:
        fig, ax = plt.subplots(figsize=(10, 1))
        for i, (label, color) in enumerate(palette_dict.items()):
            ax.fill_between([i, i + 1], 0, 1, color=color, label=label)
            ax.text(i + 0.5, 0.5, color, color='black', ha='center', va='center', fontsize=8)

        # Hide the axes
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)

        # Display the plot
        plt.legend(bbox_to_anchor=(1, 1), loc='upper left')
        plt.show()


######################################################