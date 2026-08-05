from ete3 import Tree, TreeStyle, NodeStyle, faces
import os
import argparse


def make_branches_thicker(node, new_size):
    node.img_style["hz_line_width"] = new_size
    node.img_style["vt_line_width"] = new_size
    for c in node.children:
        make_branches_thicker(c, new_size)


def draw_subclonal_tree(nwk_file, palette, out_file):
    """
    Plot clonal tree,  one per study.
    Args:
        nwk_file (str): Path to the file of the schema representation of the tree (NWK).
        palette (str): Path to clonucopya's color palette (TXT)
        out_file (str): Path to output file (PNG).
        
    """
    with open(palette, 'r') as palette:
        colors = palette.readlines()
        color_palette = [color.strip() for color in colors]

    t = Tree(nwk_file, format=1)

    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.show_scale = False

    # Set root style
    root_style = NodeStyle()
    root_style["size"] = 15
    root_style["fgcolor"] = color_palette[0]
    root_style["hz_line_type"] = 0
    root_style["hz_line_color"] = "#000"
    t.set_style(root_style)

    t.name = "non-tumor cell"
    
    # Remove root node and keep root title
    if t.name:
        root_label = faces.TextFace(t.name, fsize=12, fgcolor="black")
        t.add_face(root_label, column=0, position="branch-right")
        
    # Set the nodes of the clones
    descendant_palette = color_palette[1:]

    # Clonal nodes, excluding non-tumor cell
    clone_nodes = [node for node in t.traverse() if node != t and node.name]

    # Sort clones
    clone_nodes_sorted = sorted(clone_nodes, key=lambda n: int(n.name))

    # Map clone -> color
    clone_color_map = {
        node.name: descendant_palette[i % len(descendant_palette)]
        for i, node in enumerate(clone_nodes_sorted)
    }

    for node in t.traverse():
        if node == t:
            continue

        custom_style = NodeStyle()
        custom_style["size"] = 15
        custom_style["hz_line_type"] = 0
        custom_style["hz_line_color"] = "#000"

        # Asign a color to each color from the clonucopya's palette
        custom_style["fgcolor"] = clone_color_map[node.name]
        node.set_style(custom_style)

        # Add labels to clone nodes
        if node.name:
            label = faces.TextFace(f"c{node.name}", fsize=12, fgcolor="black")
            node.add_face(label, column=0, position="branch-right")
            
    # Change branch thickness
    make_branches_thicker(t, 1)

    # Set environment variable for offscreen rendering
    os.environ["QT_QPA_PLATFORM"] = "offscreen"

    # Save the tree as a PNG file
    t.render(out_file, w=800, units="px", tree_style=ts)


if __name__ == '__main__':
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--nwk_file", action='store', required=True)
    input_parser.add_argument("--palette", action='store', required=True)
    input_parser.add_argument("--out_file", action='store', required=True)

    args = input_parser.parse_args()


    draw_subclonal_tree(args.nwk_file,args.palette, args.out_file)
