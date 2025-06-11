import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse



def create_hex_grid(rows, cols, radius=1.0):
    """
    Define the shape of the circles for the sphere of clones

    rows: fixed number of row to lot the sphere of clones
    cols: fixed number of columns to lot the sphere of clones
    radius: length of radius for each circle

    """

    positions = []
    for r in range(rows):
        for c in range(cols):
            x = c * np.sqrt(3) * radius + (r % 2) * np.sqrt(3) * radius / 2
            y = r * 1.5 * radius
            positions.append((x, y))
    return positions


def draw_sphere_of_clones(tree_df, palette, output_path):
    """
    Plot VAF Sphere of clones,  one per sample.

    Args:
        tree_df (str): dictionary of samples' dataframes to plot the VAF Heatmap.
        palette (str): Path to clonucopya's color palette (TXT)
        out_path (str): Path to output directory of the project.
        
    """


    # Load Clonucopya's color palette
    with open(palette, 'r') as palette:
        color_palette = palette.readlines()
        colors = [color.strip() for color in color_palette]

    # Load tree inference daframe 
    tree_clones = pd.read_table(tree_df)
    clonal_prev_df = tree_clones[['clone_id', 'sample_id', 'clonal_prev']].drop_duplicates().reset_index(drop=True)
    num_clones = clonal_prev_df['clone_id'].drop_duplicates().tolist()
    sample_ids = clonal_prev_df['sample_id'].drop_duplicates().tolist()

    # Get clonal prevalence for each clone and calculate parent_cell's
    tree_inference = {}
    for sample, group in clonal_prev_df.groupby('sample_id'):
        clone_dict = group.set_index('clone_id')['clonal_prev'].to_dict()
        root_value = 1 - sum(clone_dict.values())
        clone_dict['parent_cell'] = root_value
        clone_dict_rounded = {k: round(v, 2) for k, v in clone_dict.items()}
        tree_inference[sample] = clone_dict_rounded


    for sample in tree_inference.keys():
        sample_id = sample
        proportions = tree_inference[sample]
    
    
        # Create the figure and axis
        fig, ax = plt.subplots(figsize=(6, 6))
        
        # Circle parameters
        circle_radius = 0.45 
        rows, cols = 15, 15
        boundary_radius = 5.0
        boundary_center = (cols * np.sqrt(3) * circle_radius / 2, rows * 1.5 * circle_radius / 2)
        
        # Create grid positions
        positions = create_hex_grid(rows, cols, circle_radius)
        
        # Assign nodes based on proportions
        total_circles = len(positions)
        clone_counts = {clone: int(proportion * total_circles) for clone, proportion in proportions.items()}
    
        
        # Create a list of colors for each node
        
        population_names = sorted(proportions, key=lambda x: (x != 'parent_cell', x))
        clone_color = {clone: colors[i] for i, clone in enumerate(population_names)}
        
        color_list = []
        
        for clone, count in clone_counts.items():
            color_list.extend([clone_color[clone]] * count)
        
        # Fill remaining positions
        remaining = total_circles - len(color_list)
        if remaining > 0:
            color_list.extend(['white'] * remaining)
        
        # Sort the color list to group similar colors together
        sorted_colored_circles = []
        for color in colors:
            sorted_colored_circles.extend([c for c in color_list if c == color])
        
        
        # Draw circles
        for i, (x, y) in enumerate(reversed(positions)):
            # Calculate distance from center
            dist = np.sqrt((x - boundary_center[0])**2 + (y - boundary_center[1])**2)
            # Determine if inside boundary
            inside_boundary = dist <= boundary_radius
            # Only draw circles inside the boundary
            if inside_boundary and i < len(sorted_colored_circles):
                color = sorted_colored_circles[i]
                # Create and add circle
                if color != 'white':
                    circle = plt.Circle((x, y), radius=circle_radius*0.8, color=color, edgecolor='black', linewidth=1)
                    ax.add_patch(circle)
        
        # Add legend
        legend_patches = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color,markersize=10,
                          label=f"{'' if clone == 'parent_cell' else 'clone_'}{clone} ({proportions[clone]*100:.0f}%)")
                          for clone, color in clone_color.items()]
        ax.legend(handles=legend_patches, loc='upper right', bbox_to_anchor=(1.4, 1), borderaxespad=2)
        
        # Set axis properties
        ax.set_aspect('equal')
        ax.set_xlim(0, cols * np.sqrt(3) * circle_radius)
        ax.set_ylim(0, rows * 1.5 * circle_radius)
        ax.axis('off')
        
        plt.title(f"{sample_id}")
        plt.tight_layout()

        outfile = f"{output_path}/{sample_id}_sphere_of_clones.png"
        plt.savefig(outfile)



if __name__ == '__main__':
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--tree_df", action='store', required=True)
    input_parser.add_argument("--palette", action='store', required=True)
    input_parser.add_argument("--output_path", action='store', required=True)

    args = input_parser.parse_args()


    draw_sphere_of_clones(args.tree_df,args.palette, args.output_path)
