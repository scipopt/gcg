import pandas as pd
import numpy as np
import os.path as path

def add_node(df, time, node, parent, color):
    """
    add a new node to the dataframe
    :param df: the dataframe
    :param time: timestamp of the operation
    :param node: number of the new node
    :param parent: parent of the new node
    :param color: color of the new node
    :return:
    """
    df.loc[node] = {'time': time,
                    'node_scip': np.NaN,
                    'node_hex': np.NaN,
                    'parent': parent,
                    'color': color,
                    'depth': np.NaN,
                    'var': np.NaN,
                    'primalbound': np.NaN,
                    'dualbound': np.NaN,
                    'nr': np.NaN}

def update_info(df, time, node, info):
    """
    update information of a node, already stored in the dataframe
    :param df: the dataframe
    :param time: timestamp of the operation
    :param node: number of the new node
    :param parent: parent of the new node
    :param color: color of the new node
    :return:
    """
    df.loc[node, 'time'] = time
    # for every information in info set the respective cell in the dataframe
    for key, value in info.iteritems():
        df.loc[node, key] = value
    
def parse_info_string(info_string):
    """
    Takes an info string as given in a vbc-file, parses it and returns a dict containing the information
    :param info_string: info string as given in a vbc-file
    :return info: dict containing the same information
    """
    
    # initialize the dict
    info = {}    
    
    # parse the string
    info_array = info_string.replace('\t',' ').replace('\n',' ').replace('\i',' ').replace('\\t',' ').replace('\\n',' ').replace('\\i',' ').split()
    key = ''
    for text in info_array:
        # everything before a ':' is treated as key
        if text[-1] == ':':
            key = text[:-1]
            info[key] = ''
        # everything after a ':' is treated as value
        else:
            if len(info[key]) > 0:
                info[key] += ' ' + text
            else:
                info[key] += text
                
    # rename some things
    if 'bound' in info:
        info['dualbound'] = float(info.pop('bound'))
    if 'node' in info:
        info['node_scip'] = int(info['node'].split('(')[0])
        info['node_hex'] = info.pop('node').split('(')[1].strip('()')
                
    return info
    
def read(file):
    """
    Reads a vbc-file and stores the information in a pandas dataframe
    :param file: the vbc-file to read
    :return df: the pandas dataframe
    :return info: information about the tree stored in a dict
    """
    filename, ext = path.splitext(path.basename(file))  
    # check if the file exists
    if not path.exists(file):
        print 'no file'
        return None
    
    # initialize variables
    df = pd.DataFrame(columns = ['time', 'node_scip', 'node_hex', 'parent', 'color', 'depth', 'var', 'primalbound', 'dualbound', 'nr'])
    df.index.name = 'node'
    tree_info = {}
    primalbound = np.NaN
    
    # read the lines
    for line in open(file):
        
        # read information about the whole tree from the first lines
        if line.startswith('#'):
            line_array = line[1:].split(':')
            tree_info[line_array[0]] = line_array[1]
            continue

        line_array = line.split()
        
        # first part of the line is the time, formatted as hh:mm:ss.cs
        time = line_array[0].split(':')
        time = int(time[0]) * 3600. + int(time[1]) * 60. + float(time[2])
        
        # second part of the line specifies, which information follows
        ident = line_array[1]
        line_array = line_array[2:]
        if ident == 'A':
            # add information to existing information of a node
            node = int(line_array[0])
            node_info = parse_info_string(' '.join(line_array[1:]))
            node_info['primalbound'] = primalbound
            update_info(df, time, node, node_info)
            continue
        elif ident == 'D':
            # add a new node to the tree (displayed with a delay)
            parent = int(line_array[0])
            node = int(line_array[1])
            color = int(line_array[2])
            add_node(df, time, node, parent, color)
            continue
        elif ident == 'I':
            # set the information of a node
            node = int(line_array[0])
            node_info = parse_info_string(' '.join(line_array[1:]))
            node_info['primalbound'] = primalbound
            update_info(df, time, node, node_info)
            continue
        elif ident == 'L':
            # new global lower bound (primal bound)
            primalbound = float(line_array[0])
            if 'primal_is_upper' in tree_info and tree_info['primal_is_upper']:
                tree_info['primal_is_upper'] = 'Ambiguous'
            elif (not ('primal_is_upper' in tree_info)) or (tree_info['primal_is_upper'] <> 'Ambiguous'):
                tree_info['primal_is_upper'] = False
            continue
        elif ident == 'N':
            parent = int(line_array[0])
            node = int(line_array[1])
            color = int(line_array[2])
            add_node(df, time, node, parent, color)
            # add a new node
            continue
        elif ident == 'P':
            # paint an existing node
            node = int(line_array[0])
            color = int(line_array[1])
            update_info(df, time, node, {'color': color})
            continue
        elif ident == 'U':
            # new global upper bound (primal bound)
            primalbound = float(line_array[0])
            if 'primal_is_upper' in tree_info and not tree_info['primal_is_upper']:
                tree_info['primal_is_upper'] = 'Ambiguous'
            elif (not ('primal_is_upper' in tree_info)) or (tree_info['primal_is_upper'] <> 'Ambiguous'):
                tree_info['primal_is_upper'] = True
            continue
        
    df.apply(pd.to_numeric, errors='ignore')
    df.depth = df.depth.astype('int')
    return df, tree_info