

def rotate_collation(collation):
    """rotates collation 90 degrees anticlockwise, so that each witness is its own list"""
    collation_rotated = [list(reversed(row)) for row in zip(*reversed(collation))]
    return collation_rotated

def hamming(wit1, wit2):
    """calculate the Hamming distance between two transliterations, expressed as the proportion of comparable text that is NOT identical"""
    # internal functions:
    def lacuna_p(a_string):
        """return True if a_string contains the character '[', which indicates it is lacunose"""
        if '[' in a_string:
            return True
        else:
            return False
    def omission_p(a_string):
        """return True if a_string == the character '-', which indicates it is an omission"""
        if '-' == a_string:
            return True
        else:
            return False
    def comparable_parts(wit1, wit2):
        """take two equally long lists of strings; return the parts that can be meaningfully used to calculate their Hamming distance."""
        comparable_parts_list = []
        for position in range(len(wit1)):
            pair = [wit1[position], wit2[position]]
            if not lacuna_p(pair[0]) or not lacuna_p(pair[1]): # disregard lacunae
                if not omission_p(pair[0]) and not omission_p(pair[1]): # disregard shared omissions
                    comparable_parts_list.append(pair)
        return comparable_parts_list
    def identical_parts(list_of_comparable_parts):
        """take list produced by comparable_parts(); return list of identical parts"""
        identical_parts_list = []
        for item in list_of_comparable_parts:
            if item[0] == item[1]:
                identical_parts_list.append(item)
        return identical_parts_list
    # main function body:
    comparable_length = len(comparable_parts(wit1, wit2))
    identical_length = len(identical_parts((comparable_parts(wit1, wit2))))
    hamming_distance = 1 - (identical_length / comparable_length)
    return hamming_distance

def dist_matrix(collation):
    """takes a standard collation; returns a distance matrix"""
    # rotate the collation and strip the reference numbers:
    rot_col = rotate_collation(collation)[1:] # line 0 is just reference numbers, so remove it.
    # initialize empty distance matrix:
    dist_matrix = []
    dist_matrix_top_line = [' '] # top left square should be left blank
    # add all the witness names to the top line of the distance matrix:
    for line in rot_col:
        dist_matrix_top_line.append(line[0])
    # append the constructed top line to the distance matrix:
    dist_matrix.append(dist_matrix_top_line)
    # construct the rest of the distance matrix:
    for this_witness in rot_col:
        h_dists = [this_witness[0]] # first item in table should be witness's name.
        for other_witnesses in rot_col:
            h_dist = hamming(this_witness[1:], other_witnesses[1:]) # position 0 is just the witness name.
            h_dists.append(h_dist)
        dist_matrix.append(h_dists)
    # remove duplicate entries from the matrix:
    deduped_dist_matrix = []
    for line in dist_matrix:
        deduped_line = []
        reached_end = False
        for entry in line:
            if entry == 0.0:
                reached_end = True
            if reached_end == False:
                deduped_line.append(entry)
        deduped_dist_matrix.append(deduped_line)
    return deduped_dist_matrix

def dist_matrix_to_list(d_mat):
    """takes a distance matrix; returns a list where each item has the format (distance, witness1, witness2)"""
    # internal functions:
    def witness_in_column(column_number):
        """takes a column number in the distance matrix; returns the witness name to which it corresponds"""
        column_names_list = d_mat[0]
        return column_names_list[column_number]
    # main body:
    distance_list = [] # initialize empty list of distances
    for row in d_mat[2:]: # ignore top row, which is just witness
                          # names, and second row, which is the first
                          # witness compared to itself.
        wit1 = row[0] # row 0 contains the first witness' name
        for col_num in range(len(row)):
            if col_num != 0: # don't consider the first column number,
                             # as it just contains the witness's name.
                distance = row[col_num]
                wit2 = witness_in_column(col_num)
                distance_list.append([distance, wit1, wit2])
    return sorted(distance_list)

def unaffected_witnesses(dist_list, wit1, wit2):
    """returns all entries in dist_list which do not contain either wit1 or wit2"""
    unaffected = []
    for item in dist_list:
        if not wit1 in item and not wit2 in item:
            unaffected.append(item)
    return unaffected

def affected_witnesses(dist_list, wit1, wit2):
    affected = []
    for item in dist_list:
        if wit1 in item and wit2 in item:
            pass # ignore the item where wit1 and wit2 are both present.
        elif wit1 in item or wit2 in item:
            affected.append(item)
    return affected

def wit_distances(dist_list, wit):
    """returns all distances between wit and all other witnesses mentioned in dist_list"""
    # feed this function a distance list produced by
    # affected_witnesses() to avoid including wit's sibling in the
    # returned list.
    wit_dist_list = []
    for item in dist_list:
        if item[1] == wit:
            distance = [item[0], item[2]]
            wit_dist_list.append(distance)
        if item[2] == wit:
            distance = [item[0], item[1]]
            wit_dist_list.append(distance)
    return wit_dist_list

def average_wit_distance(wit_dist_list):
    """takes output produced by wit_distances; returns the average distance"""
    # for Neighbour-joining, it's important to feed the
    # affected_witnesses list to average_wit_distance, becuase we
    # don't want to include the distance between a witness and its
    # sibling when calculating the average distance:
    num_wits = len(wit_dist_list)
    dists = []
    for item in wit_dist_list:
        dists.append(item[0])
    average_dist = sum(dists) / num_wits
    return average_dist

def dist_between(dist_list, wit1, wit2):
    """returns the distance between wit1 and wit2 as given by dist_list"""
    dist="NOT FOUND" # throw an error if one of the siblings is not in
                     # the list: this should never be returned if
                     # everything is working.
    for item in dist_list:
        if wit1 in item and wit2 in item:
            dist = item[0]
    if wit1 == wit2:
        dist = 0
    return dist

def sib_dist_to_node(dist_list, wit, sibling):
    """takes A FULL DIST_LIST; returns the distance between a witness and its node"""
    # for this calculation, see Saitou 1987: 409, 6a and 6b.
    dist_between_siblings = dist_between(dist_list, wit, sibling)
    affected_witness_dist_list = affected_witnesses(dist_list, wit, sibling)
    wit_dist_to_affecteds = average_wit_distance(wit_distances(affected_witness_dist_list, wit))
    sibling_dist_to_affecteds = average_wit_distance(wit_distances(affected_witness_dist_list, sibling))
    dist_to_node = (dist_between_siblings + wit_dist_to_affecteds - sibling_dist_to_affecteds) / 2
    return dist_to_node

def broadcast_sibling_info(dist_list, wit1, wit2, node_name):
    """receives a dist_list and two witnesses which are known to be siblings; reuturns their distances to their shared node"""
    wit1_dist = sib_dist_to_node(dist_list, wit1, wit2)
    wit2_dist = sib_dist_to_node(dist_list, wit2, wit1)
    sibling_statement = [node_name,[wit1, wit1_dist],[wit2, wit2_dist]]
    return sibling_statement

def list_affected_witnesses(affected_dist_list, sib1, sib2):
    """takes a distance list produced by affected_witnesses; returns a list of the names of the affected witnesses"""
    affected_wits = []
    for row in affected_dist_list:
        if row[1] == sib1 or row[1] == sib2:
            affected_wits.append(row[2])
        elif row[2] == sib1 or row[2] == sib2:
            affected_wits.append(row[1])
    return list(set(affected_wits)) # remove duplicate entries, but return a list rather than a set.

def dists_to_be_averaged(affected_dist_list, witness):
    """takes an affected distance list and a witness; returns the two distances that need averaging"""
    dists_to_be_averaged = []
    for row in affected_dist_list:
        if witness in row:
            dists_to_be_averaged.append(row[0])
    return dists_to_be_averaged

def return_new_node(affected_dist_list, witness, node_name):
    """returns a correctly formatted list for the distance between witness and node_name"""
    dists = dists_to_be_averaged(affected_dist_list, witness)
    average_distance = (dists[0] + dists[1]) / 2
    return [average_distance, node_name, witness]

def return_new_node_dists(affected_dist_list, sib1, sib2, node_name):
    """returns a list of all distances to the new node"""
    wits_to_update = list_affected_witnesses(affected_dist_list, sib1, sib2)
    new_dists = []
    for witness in wits_to_update:
        new_node_dist = return_new_node(affected_dist_list, witness, node_name)
        new_dists.append(new_node_dist)
    return new_dists

def grab_all_witness_names(dist_list):
    """returns a list of all witness names in a distance list"""
    wit_names = []
    for entry in dist_list:
        wit_names.append(entry[1])
        wit_names.append(entry[2])
    return list(set(wit_names))

def neighbour_joiner(dist_list, siblings_list, next_node_num):
    """Takes a SIMPLE, FULL, SORTED distance list, an (initially) empty list, the (starting) node number, and a list of witnesses to place (starting with all witnesses); returns a description of the Neighbour-joined tree"""

    # Identify the next pair of siblings:
    sibling_line = dist_list[0] # the pair of siblings will be the
                                # ones at the top of the sorted
                                # distance list.
    sib1 = sibling_line[1] # sibling_line[0] will be a distance.
    sib2 = sibling_line[2]

    # calculate the sibling info and append it to siblings_list:
    next_siblings = broadcast_sibling_info(dist_list, sib1, sib2, str(next_node_num))
    siblings_list.append(next_siblings)

    # generate the list of unaffected witnesses, to add to the distance list later:
    unaffected_dist_list = (unaffected_witnesses(dist_list, sib1, sib2))

    # generate the distance list for affected witnesses:
    affected_dist_list = affected_witnesses(dist_list, sib1, sib2)

    # generate the new node's distance list:
    new_node_dist_list = return_new_node_dists(affected_dist_list, sib1, sib2, str(next_node_num)) # note stringification of the node number.

    # create the new distance list:
    new_dist_list = []
    new_dist_list.extend(unaffected_dist_list) # note extend, not append.
    new_dist_list.extend(new_node_dist_list) # note extend, not append.

    # SORT THE FINAL LIST, TO PUT THE NEXT PAIR OF SIBLINGS AT THE TOP:
    new_dist_list = sorted(new_dist_list)
    
    # increment the node number for the next time this function is run:
    next_node_num = next_node_num + 1

    # if there are no witnesses left to place, return the description of the siblings:
    if len(new_dist_list) == 1:
        siblings_list.append([new_dist_list[0][1],
[new_dist_list[0][2], new_dist_list[0][0]], [new_dist_list[0][2], new_dist_list[0][0]]])
        # print(siblings_list) # for testing
        # print(new_dist_list)
        return siblings_list
    else:
        # print(siblings_list)
        # print(sib1, sib2)
        # print(new_dist_list)
        return neighbour_joiner(new_dist_list, siblings_list, next_node_num)

def n_j_to_graphviz(n_j_output):
    """receives neighbour-joiner output; returns graphviz code"""
    # internal functions
    def write_line(row):
        """takes a row of n_j_output; returns a line of graphviz code"""
        node = row[0]
        connection_1 = str(row[1][0])
        distance_1 = str(round(row[1][1], 3))
        connection_2 = str(row[2][0])
        distance_2 = str(round(row[2][1], 3))
        graphviz_line_1 = node + ' -- "' + connection_1 + '" [ label="' + distance_1 + '",minlen=30 ];'
        graphviz_line_2 = node + ' -- "' + connection_2 + '" [ label="' + distance_2 + '",minlen=30 ];'
        if graphviz_line_1 == graphviz_line_2: # deal with final entry of the table
            return [graphviz_line_1]
        else:
            return [graphviz_line_1, graphviz_line_2]
    # main body
    graphviz_lines = ['graph G {', 'overlap = false;', 'splines = true']
    for row in n_j_output:
        graphviz_lines.extend(write_line(row))
    graphviz_lines.append('}')
    for row in graphviz_lines:
        print(row)
    return

def nj_topology(collation):
    return neighbour_joiner(dist_matrix_to_list(dist_matrix(collation)), [], 1)

def nj_graphviz(collation):
    return n_j_to_graphviz(nj_topology(collation))
