from ete3 import Tree, TreeStyle, NodeStyle, faces
import os
import argparse


def make_branches_thicker(node, new_size):
    node.img_style["hz_line_width"] = new_size
    node.img_style["vt_line_width"] = new_size
    for c in node.children:
        make_branches_thicker(c, new_size)


def draw_subclonal_tree(nwk_file, palette, outfile):
    with open(palette, 'r') as palette:
        colors = palette.readlines()
        color_palette = [color.strip() for color in colors]

    t =  Tree(nwk_file, format=1) 

    # Define basic configuration of the tree
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
    
    # Remove root node and keep root title
    if t.name:
        root_label = faces.TextFace(t.name, fsize=12, fgcolor="black")
        t.add_face(root_label, column=0, position="branch-right")

    # Set the nodes of the clones
    color_index = 1
    
    for node in t.traverse():
        # Pass root node because is has its own style
        if node == t:
            continue
            
        custom_style = NodeStyle()
        custom_style["size"] = 15
        custom_style["hz_line_type"] = 0
        custom_style["hz_line_color"] = "#000"
        
        # Asign a color to each color from the clonucopya's palette
        custom_style["fgcolor"] = color_palette[color_index % len(color_palette)]
        color_index += 1      
        node.set_style(custom_style)
        
        # Add labels to clone nodes
        if node.name:
            label = faces.TextFace(node.name, fsize=12, fgcolor="black")
            node.add_face(label, column=0, position="branch-right")
    
        # Change branch thickness
        make_branches_thicker(t, 1)
        
    # Set environment variable for offscreen rendering
    os.environ["QT_QPA_PLATFORM"] = "offscreen"
    
    # Save the tree as a PNG file
    t.render(outfile, w=800, units="px", tree_style=ts)

if __name__ == '__main__':
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--nwk_file", action='store', required=True)
    input_parser.add_argument("--palette", action='store', required=True)
    input_parser.add_argument("--output_file", action='store', required=True)

    args = input_parser.parse_args()


    draw_subclonal_tree(args.nwk_file,args.palette, args.output_file)
