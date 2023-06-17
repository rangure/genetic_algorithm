import graphviz
import os
import sys

def read_inputs(file_path):
    """
    This function reads two lines from an input file.
    The first line is the circuit vector(int) and the second line is the units vector(float).
    
    Args:
        file_path (str): The path to the input file.
        
    Returns:
        list, list2: Two lists, each containing the values from each line in the input file.
    """
    with open(file_path, 'r') as file:
        data = file.readline().replace(',', ' ')
        data2 = file.readline().replace(',', ' ')
    items = data.split()
    items2 = data2.split()
    list = [int(elem) for elem in items]
    list2 = [round(float(elem), 2) for elem in items2]
    return list, list2


def generate_graph(file_path, graph_name):
    """
    This function generates a directed graph based on the input file and saves it as an SVG image.
    
    Args:
        file_path (str): The path to the input file.
        graph_name (str): The name of the graph (also used as the name of the output image file).
    """
    int_list, int_list2 = read_inputs(file_path)

    # generate a directed graph
    graph = graphviz.Digraph()

    # LR: from left to right
    graph.attr(rankdir='LR')
    #graph.attr(rankdir='LR', size='8,5', nodesep='0.05', ranksep='0.1')

    # Node: shape is 'rectangle'
    graph.attr('node', shape='rectangle')
    index = 1
    graph.edge('Feed',
               'Unit ' + str(int_list[0]),
               color='blue',
               headport='w',
               tailport='e',
               arrowhead='normal')

    # input to feed
    graph.node('100 kg/s of the waste material',
               shape='none',
               width='0',
               height='0')
    graph.edge('100 kg/s of the waste material',
               'Feed', len='0.5',
               color='black',
               headport='w',
               tailport='e',
               arrowhead='normal')
    graph.node('10 kg/s of gerardium', shape='none', width='0', height='0')
    graph.edge('10 kg/s of gerardium',
               'Feed', 
               color='black',
               headport='w',
               tailport='e',
               arrowhead='normal')

    graph.edge('Feed',
               'Unit ' + str(int_list[0]),
               color='blue', len='0.5',
               headport='w',
               tailport='e',
               arrowhead='normal')

    # Get the largest and second largest values in int_list
    copy_list = int_list.copy()
    largest = max(copy_list)
    copy_list.remove(largest)
    second_largest = max(copy_list)
    while (second_largest == largest):
        copy_list.remove(second_largest)
        second_largest = max(copy_list)

    for i in range(int(len(int_list) / 2)):
        node_name = 'Unit ' + str(i)
        if int_list[index] == largest:
            end_node = 'Tailing'
        elif int_list[index] == second_largest:
            end_node = 'Concentrate'
        else:
            end_node = 'Unit ' + str(int_list[index])

        graph.edge(node_name,
                   end_node,
                   label=str(int_list2[index - 1]),
                   color='red',
                   arrowhead='normal')

        if int_list[index + 1] == largest:
            end_node = 'Tailing'
        elif int_list[index + 1] == second_largest:
            end_node = 'Concentrate'
        else:
            end_node = 'Unit ' + str(int_list[index + 1])

        graph.edge(node_name,
                   end_node,
                   label=str(int_list2[index]),
                   color='blue',
                   arrowhead='normal')

        index += 2
         
    str_list = [str(num) for num in int_list]
    vec = ", ".join(str_list)
    vec = "[" + vec + "]"
    graph.node("vector = " + vec, shape='none', width='0', height='0')
    graph.node('', shape='none', width='0', height='0')
    graph.node('Concentrate flow with path',
               shape='none',
               width='0',
               height='0')
    graph.edge('',
               'Concentrate flow with path',
               label='Concentrate flow rate (kg/s)',
               color='red',
               len='0.5')

    graph.node(' ', shape='none', width='0', height='0')
    graph.node('Tailing flow with path', shape='none', width='0', height='0')
    graph.edge(' ',
               'Tailing flow with path',
               label='Tailing flow rate (kg/s)',
               color='blue',
               len='0.5')
    # graph.attr(dpi='300')
    # render: generate graph.
    graph.render(str('output/' + graph_name), cleanup=True, format='svg')


def process_dat_files(directory_path):
    """
    This function processes all .dat files in a specified directory by generating a directed graph for each file.
    
    Args:
        directory_path (str): The path to the directory containing the .dat files.
    """
    dat_files_list = [
        f for f in os.listdir(directory_path) if f.endswith('.dat')
    ]
    dat_files_list.sort()

    # Call the generate_graph function on each .dat file
    for file in dat_files_list:
        file_path = os.path.join(directory_path, file)
        file_name = os.path.splitext(file)[0]
        generate_graph(file_path, file_name)

# Get the directory path from the command line arguments
def pass_command_line_fileNmae(directory_path):
    """
    This function gets the directory path from the command line arguments. If no arguments are provided, it 
    uses the default directory path './out'. It then calls process_dat_files on this directory path.
    
    Args:
        directory_path (str): The default directory path.
    """
    if len(sys.argv) > 1:
        directory_path = sys.argv[1]
    else:
        print("Please provide a directory path as an argument or read data from default path ./out file.")
    process_dat_files(directory_path)
# call function with directory path
directory_path = './out'
pass_command_line_fileNmae(directory_path)

