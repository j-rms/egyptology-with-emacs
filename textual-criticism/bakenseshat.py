# ------------------------------------------------------------------------------
# BAKENSESHAT: tools for Egyptological textual criticism
# ------------------------------------------------------------------------------

import matplotlib
import matplotlib.pyplot
import numpy
import statistics
import math
from itertools import combinations
from itertools import cycle
from itertools import groupby
import re

def rotate_collation(collation):
    """rotates collation 90 degrees anticlockwise, so that each witness is
its own list"""
    collation_rotated = [list(reversed(row)) for row in zip(*reversed(collation))]
    return collation_rotated

# ------------------------------------------------------------------------------
# HAMMING DISTANCE FUNCTIONS
# ------------------------------------------------------------------------------

def hamming(wit1, wit2):
    """calculate the Hamming distance between two transliterations,
expressed as the proportion of comparable text that is NOT
identical"""
    # internal functions:
    def lacuna_p(a_string):
        """return True if a_string contains the character '[', which indicates
it is lacunose"""
        if '[' in a_string:
            return True
        else:
            return False
    def omission_p(a_string):
        """return True if a_string == the character '‑' (non-breaking hyphen, not normal - sign!!!), which indicates it is
an omission"""
        if '‑' == a_string:
            return True
        else:
            return False
    def comparable_parts(wit1, wit2):
        """take two equally long lists of strings; return the parts that can
be meaningfully used to calculate their Hamming distance."""
        comparable_parts_list = []
        for position in range(len(wit1)):
            pair = [wit1[position], wit2[position]]
            if not lacuna_p(pair[0]) and not lacuna_p(pair[1]): # disregard lacunae
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
    if comparable_length == 0:
        comparable_length = 0.0000001 # a comparable length of 0 will produce a divide by zero error when calculating hamming_distance below; have to use an extremely small value instead.
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
            if other_witnesses[0:] == this_witness[0:]:
                h_dists.append("itself")
            else:
                h_dist = hamming(this_witness[1:], other_witnesses[1:]) # position 0 is just the witness name.
                h_dists.append(h_dist)
        dist_matrix.append(h_dists)
    # remove duplicate entries from the matrix:
    deduped_dist_matrix = []
    for line in dist_matrix:
        deduped_line = []
        reached_end = False
        for entry in line:
            if entry == "itself":
                reached_end = True
            if reached_end == False:
                deduped_line.append(entry)
        deduped_dist_matrix.append(deduped_line)
    return deduped_dist_matrix

def full_dist_matrix(collation):
    """takes a standard collation; returns a full distance matrix (i.e. duplicates not removed)"""
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
            if other_witnesses[0:] == this_witness[0:]:
                h_dists.append(0)
            else:
                h_dist = hamming(this_witness[1:], other_witnesses[1:]) # position 0 is just the witness name.
                h_dists.append(h_dist)
        dist_matrix.append(h_dists)
    return dist_matrix

# ------------------------------------------------------------------------------
# NEIGHBOUR-JOINING FUNCTIONS
# ------------------------------------------------------------------------------

def dist_matrix_to_list(d_mat):
    """takes a distance matrix; returns a list where each item has the
format (distance, witness1, witness2)"""
    # internal functions:
    def witness_in_column(column_number):
        """takes a column number in the distance matrix; returns the witness
name to which it corresponds"""
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

def round_dlist(dlist, rounding):
    newlist = []
    for row in dlist:
        newlist.append([round(row[0], rounding), row[1], row[2]])
    return newlist

def rounded_dlist(col, rounding):
    """takes a COLlation and a ROUNDING number; returns an appropriately rounded distance list"""
    return round_dlist(dist_matrix_to_list(dist_matrix(col)), rounding)

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
    """returns all distances between wit and all other witnesses mentioned
in dist_list"""
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
    """takes output produced by wit_distances; returns the average
distance"""
    # for Neighbour-joining, it's important to feed the
    # affected_witnesses list to average_wit_distance, becuase we
    # don't want to include the distance between a witness and its
    # sibling when calculating the average distance:
    # print(wit_dist_list)
    num_wits = len(wit_dist_list)
    dists = []
    for item in wit_dist_list:
        dists.append(item[0])
    # print(dists)
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
    """takes A FULL DIST_LIST; returns the distance between a witness and
its node"""
    # for this calculation, see Saitou 1987: 409, 6a and 6b.
    dist_between_siblings = dist_between(dist_list, wit, sibling)
    affected_witness_dist_list = affected_witnesses(dist_list, wit, sibling)
    wit_dist_to_affecteds = average_wit_distance(wit_distances(affected_witness_dist_list, wit))
    sibling_dist_to_affecteds = average_wit_distance(wit_distances(affected_witness_dist_list, sibling))
    dist_to_node = (dist_between_siblings + wit_dist_to_affecteds - sibling_dist_to_affecteds) / 2
    return dist_to_node

def broadcast_sibling_info(dist_list, wit1, wit2, node_name):
    """receives a dist_list and two witnesses which are known to be
siblings; reuturns their distances to their shared node"""
    wit1_dist = sib_dist_to_node(dist_list, wit1, wit2)
    wit2_dist = sib_dist_to_node(dist_list, wit2, wit1)
    sibling_statement = [node_name,[wit1, wit1_dist],[wit2, wit2_dist]]
    return sibling_statement

def list_affected_witnesses(affected_dist_list, sib1, sib2):
    """takes a distance list produced by affected_witnesses; returns a
list of the names of the affected witnesses"""
    affected_wits = []
    for row in affected_dist_list:
        if row[1] == sib1 or row[1] == sib2:
            affected_wits.append(row[2])
        elif row[2] == sib1 or row[2] == sib2:
            affected_wits.append(row[1])
    return list(set(affected_wits)) # remove duplicate entries, but return a list rather than a set.

def dists_to_be_averaged(affected_dist_list, witness):
    """takes an affected distance list and a witness; returns the two
distances that need averaging"""
    dists_to_be_averaged = []
    for row in affected_dist_list:
        if witness in row:
            dists_to_be_averaged.append(row[0])
    return dists_to_be_averaged

def return_new_node(affected_dist_list, witness, node_name):
    """returns a correctly formatted list for the distance between witness
and node_name"""
    dists = dists_to_be_averaged(affected_dist_list, witness)
    # print(len(dists))
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
    """Takes a SIMPLE, FULL, SORTED distance list, an (initially) empty
list, the (starting) node number, and a list of witnesses to place
(starting with all witnesses); returns a description of the
Neighbour-joined tree"""

    # Identify the next pair of siblings:
    sibling_line = dist_list[0] # the pair of siblings will be the
                                # ones at the top of the sorted
                                # distance list.
    # print("sibling_line: " + str(sibling_line))
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
    # print(new_dist_list)
    # increment the node number for the next time this function is run:
    next_node_num = next_node_num + 1

    # if there are no witnesses left to place, return the description of the siblings:
    if len(new_dist_list) == 1: # or len(new_dist_list) == 0: #not sure if this is needed or not.
        siblings_list.append([new_dist_list[0][1],
[new_dist_list[0][2], new_dist_list[0][0]], [new_dist_list[0][2], new_dist_list[0][0]]])
        # print("SIBLINGS LIST:")
        # print(siblings_list) # for testing
        # print("NEW DIST LIST:")
        # print(new_dist_list)
        return siblings_list
    else:
        # print("SIBLINGS LIST:")
        # print(siblings_list)
        # # print(sib1, sib2)
        # print("NEW DIST LIST:")
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
        graphviz_line_1 = '"' + node + '"' + ' -- "' + connection_1 + '" [ label="' + distance_1 + '",minlen=30,shape=oval ];'
        graphviz_line_2 = '"' + node + '"' + ' -- "' + connection_2 + '" [ label="' + distance_2 + '",minlen=30,shape=oval ];'
        if graphviz_line_1 == graphviz_line_2: # deal with final entry of the table
            return [graphviz_line_1]
        else:
            return [graphviz_line_1, graphviz_line_2]
    # main body
    # graphviz_lines = ['graph G {', 'overlap = false;', 'splines = true']
    graphviz_lines = ['graph G {', 'node [fontname="Charis SIL"]', 'edge [fontname="Charis SIL"]', 'overlap=prism', 'overlap_scaling=2', 'ratio=1', 'splines = true'] # more compact

    for row in n_j_output:
        graphviz_lines.extend(write_line(row))
    graphviz_lines.append('}')
    for row in graphviz_lines:
        print(row)
    return graphviz_lines

def nj_topology(collation):
    return neighbour_joiner(dist_matrix_to_list(dist_matrix(collation)), [], 1)

def nj_graphviz(collation):
    return n_j_to_graphviz(nj_topology(collation))

def nj_pdf(col, file_name, work_name, caption_postscript):
    graph_lines = nj_graphviz(col)
    final_string = ""
    for line in graph_lines:
        final_string = final_string + line + "\n"
    from graphviz import Source
    s = Source(final_string, filename=file_name, format="pdf", engine="neato")
    s.render()
    caption = "#+caption: Neighbour-joiner chain for " + work_name + ". " + caption_postscript
    return caption + '\n' + '#+name: ' + file_name.rsplit('.', maxsplit=1)[0] + '\n#+attr_latex: :placement [t] :width \\textwidth' + '\n' + 'file:' + file_name + ".pdf"


def nj_autoreport(work_name, file_dir, col, lac_threshold):
    autoreport = ""
    file_name = work_name.replace(" ", "_")
    if file_dir[-1] != '/':
        file_dir = file_dir + '/'
    file_path = file_dir + file_name
    # lacunosity graph
    lac_file_name = file_path + "_lac_chart.pdf"
    lac_chartt = lac_chart(col, 'Charis SIL', lac_file_name, '', 4.5, 7.5, work_name)
    autoreport = autoreport + lac_chartt
    # scatter graph
    matplotlib.pyplot.figure(figsize=(6.4,4.8))
    scat_file_name = file_path + "_lact_scatter.pdf"
    stripped_col = strip_lacunose_witnesses(col, lac_threshold)
    scat_strip_amount = 'witness pair: witnesses ≥ ' + str(lac_threshold) + '% lacunose omitted'
    scat_chart = scatter_2_cols(col, stripped_col, 'witness pair: full collation', scat_strip_amount, 'Charis SIL', scat_file_name)
    autoreport = autoreport + "\n \n" + scat_chart
    # neighbour joining full collation
    full_nj_file_name = file_path + "_nj_full"
    full_nj = nj_pdf(col, full_nj_file_name, work_name, "All witnesses.")
    autoreport = autoreport + "\n \n" + full_nj
    # neighbour joining lac-threshold collation
    cut_nj_file_name = file_path + "_nj_cut"
    cut_nj_post_caption = "Witnesses ≥ " + str(lac_threshold) + " lacunose removed."
    cut_nj = nj_pdf(stripped_col, cut_nj_file_name, work_name, cut_nj_post_caption)
    autoreport = autoreport + "\n \n" + cut_nj
    return autoreport

def nj_lac_explore(col, low, high, step, filename):
    """given a collation, an upper and lower threshold, and a step value, generate scatter charts and neighbour-joined chains for each step value"""
    explore_report = ""
    matplotlib.pyplot.figure(figsize=(6.4,4.8))
    for lac_threshold in numpy.arange(low, high, step):
        # generate the stripped collation (omitting all witnesses over lac_threshold):
        stripped_col = strip_lacunose_witnesses(col, lac_threshold)
        # generate the scatter chart comparing the stripped collation with the full collation:
        scat_file_name = filename + "_EXPLORE_" + str(lac_threshold) + "_scat.pdf"
        scat_strip_amount = 'witness pair: witnesses ≥ ' + str(lac_threshold) + '% lacunose omitted'
        scat_chart = scatter_2_cols(col, stripped_col, 'witness pair: full collation', scat_strip_amount, 'Charis SIL', scat_file_name)
        # generate the neighbour-joined chain for the stripped collation:
        cut_nj_file_name = filename + "_EXPLORE_" + str(lac_threshold) + "_nj"
        cut_nj_post_caption = "Witnesses ≥ " + str(lac_threshold) + " lacunose removed."
        cut_nj = nj_pdf(stripped_col, cut_nj_file_name, filename, cut_nj_post_caption)
        # update explore_report:
        explore_report = explore_report + scat_chart + "\n" + cut_nj + "\n \n"
    return explore_report

def rob_one_out(col, filename, workname):
    """given a collation, create nj chains where each chain omits one of the collation's witnesses (robustness check)"""
    namelist = col[0][1:]
    rob_report = str(namelist) + "\n \n"
    # return the full collation
    subset_file_name = filename + "_rob_" + "000"
    subset_post_caption = "No witnesses removed."
    robbed_nj = nj_pdf(col, subset_file_name, workname, subset_post_caption)
    rob_report = rob_report + "\n \n" + robbed_nj
    
    for name in namelist:
        new_namelist = [word for word in namelist if word != name]
        new_col = sub_col_rownums(col, new_namelist) # make a new sub-collation without that witness
        # produce nj chain for that witness subset;
        subset_file_name = filename + "_rob_" + name.replace(" ", "_").replace(".", "_")
        subset_post_caption = name + " removed."
        robbed_nj = nj_pdf(new_col, subset_file_name, workname, subset_post_caption)
        rob_report = rob_report + "\n \n" + robbed_nj
    return rob_report

def col_rep_string(col, str1, str2):
    """in the collation COL, every time the string STR1 appears in a cell, replace it with STR2"""
    new_col = [col[0]]
    for row in col[1:]:
        new_row = []
        new_row.append(row[0])
        for cell in row[1:]:
            contents = str(cell)
            contents = contents.replace(str1, str2)
            new_row.append(contents)
        new_col.append(new_row)
    return new_col


def n_chunk(col, chunks):
    """split a collation into n chunks, and return each chunk as its own sublist, with a sublist containing just the names as the first item in the whole list"""
    chunk_size = round(len(col[1:]) / chunks) # note that this will
                                                # not produce quite
                                                # equal divisions: the
                                                # final chunk could be
                                                # quite a bit larger
                                                # or smaller due to
                                                # rounding.

    chunked_col = [[col[0]]]
    startpoint = 1
    for chunk in range(chunks):
        chunk_start = startpoint
        chunk_end = startpoint + chunk_size
        startpoint = chunk_end
        new_chunk = col[chunk_start:chunk_end]
        chunked_col.append([new_chunk])
    return chunked_col

def remove_chunks(col, chunks):
    """given a collation COL, generate CHUNKS copies with CHUNKS% of the collation removed in each copy, starting with the first CHUNK% of the collation and working through it, and output as a list of lists"""
    col_heading = col[0]
    col_chunks = n_chunk(col, chunks)[1:]
    col_chunks_removed = []
    for ignore_sublist in col_chunks:
        col_a_chunk_removed = [col[0]]
        for sublist in col_chunks:
            if ignore_sublist != sublist:
                for row in sublist:
                    col_a_chunk_removed.extend(row)
        col_chunks_removed.append(col_a_chunk_removed)
    return col_chunks_removed
            
        
def nj_remove_chunks(col, num_chunks, file_name, work_name, caption_postscript):
    """Robustness test: produce neighbour-joined PDFs for COL, incrementally removing NUM_CHUNKS% of the collation"""
    collations = remove_chunks(col, num_chunks)
    counter = 1
    output = ""
    this_file_name = file_name + "_remove_" +str(num_chunks) + "_chunks_NONE_REMOVED"
    full_postscript = "Removing " + str(round(100 / num_chunks)) + "% chunks incrementally: no chunks removed. " + caption_postscript
    output = output + nj_pdf(col, this_file_name, work_name, full_postscript)
    for collation in collations:
        this_file_name = file_name + "_remove_" + str(num_chunks) + "_chunks_num_" + str(counter)
        full_postscript = "Removing " + str(round(100 / num_chunks)) + "% chunks incrementally: chunk " + str(counter) + " of " + str(num_chunks) + " removed. " + caption_postscript
        output = output + "\n \n" + nj_pdf(collation, this_file_name, work_name, full_postscript)
        counter = counter + 1
    return output
        

def return_chunk(col, num_chunks, chunk_to_return):
    """Given a collation, a number of chunks into which to split the collation, and the chunk number to return, return the chunk number headed by a namelist"""
    chunks = n_chunk(col, num_chunks)
    chunk = chunks[chunk_to_return][0] # this works because the first item in the list returned by n_chunk, item 0, is the namelist.
    namelist = chunks[0]
    namelist.extend(chunk)
    return namelist

def type_2_locs_for_chunk(col, num_chunks, chunk_to_return):
    """Return just the type-2 variants for the chunk."""
    namelist = col[0][1:]
    chunk = return_chunk(col, num_chunks, chunk_to_return)
    first_rownum = int(chunk[1][0])
    last_rownum = int(chunk[-1][0])
    type_2s = type_2_locs(sub_col(col, namelist))
    # Now build a list of type 2s that fall between first_rownum and last_rownum, inclusive.
    # don't forget to turn them into integers for comparison!
    returned_type_2s = [chunk[0]]
    for row in type_2s:
        if int(row[0]) >= first_rownum and int(row[0]) <= last_rownum:
            returned_type_2s.append(row)
    return returned_type_2s

# ------------------------------------------------------------------------------
# SEMIAUTOMATED TEXTUAL CRITICISM FUNCTIONS
# ------------------------------------------------------------------------------

def get_wit_text(collation, name):
    """Return a named witness's text from a standard collation table"""
    namerow=collation[0]
    index_number = namerow.index(name)
    col_rotated = rotate_collation(collation)
    return col_rotated[index_number]

def sub_col(collation, name_list):
    """return a sub-collation from COLLATION containing only those witnesses in NAME_LIST"""
    subcol = []
    for name in name_list:
        subcol.append(get_wit_text(collation, name))
    subcol = rotate_collation(subcol) # rotate the collation back to columnar format
    return subcol

def type_2_locs(collation):
    """return a list of all lines with type_2 variations in COLLATION"""
    # internal functions
    def no_lacuna_p(line):
        """returns True if no element of LINE contains '[', i.e. is lacunose"""
        no_lacuna_p = True
        for item in line:
            if '[' in str(item):
                no_lacuna_p = False
        return no_lacuna_p
    def no_insignificant_omission_p(line):
        """returns True if no element of LINE == '‑' (non-breaking hyphen, not normal - sign!!!), i.e. is an insignificant omission"""
        no_insignificant_omission_p = True
        for item in line:
            if "‑" == str(item):
                no_insignificant_omission_p = False
        return no_insignificant_omission_p
    # body
    two_version_lines = []
    linum = 0
    linums = []
    for line in collation:
        number_of_variants = len(set(line))
        if number_of_variants == 2 and no_lacuna_p(line) and no_insignificant_omission_p(line):
            two_version_lines.append(line)
            linums.append(linum)
        linum = linum + 1
    # having established which lines contain two separate readings, we
    # must now only return those lines where each reading is attested
    # at least twice:
    type_2_lines = []
    type_2_linums = []
    entry = 0
    for line in two_version_lines:
        versions = list(set(line))
        if line.count(versions[0]) >= 2 and line.count(versions[1]) >= 2:
            final_line = []
            final_line.append(linums[entry])
            final_line.extend(line)
            type_2_lines.append(final_line)
        entry = entry + 1

    final_table = [0]
    final_table.extend(collation[0])
    final_table = [final_table]
    final_table.extend(type_2_lines)
    return final_table

def witness_combos(witness_list):
    combos = []
    for n in range(1, len(witness_list) + 1): 
        combos += list(combinations(witness_list, n))
    return combos

def pair_combos(combo_list):
    """returns a COMBO_LIST as a list of paired lists"""
    num_wits = len(combo_list[-1]) # the number of witnesses is equal
                                   # to the length of the last item of
                                   # the list.
    wit_set = set(combo_list[-1]) # get the names of all witnesses as a set
    combos = combo_list[:-1] # don't count the last item
    paired_lists = []
    sets_considered = []
    for item in combos:
        this_set = set(item)
        other_set = wit_set.difference(this_set)
        if this_set not in sets_considered:
            pair = [list(this_set), list(other_set)]
            sets_considered.append(other_set)
            paired_lists.append(pair)
    paired_lists.append([list(wit_set),[]])
    return paired_lists

def score_type_2s_by_wit_set(t2_loc_table, divided_wit_set):
    """returns a score indicating how optimally a particular divided witness set divides the readings of a type-2 loc table produced by TYPE_2_LOCS()"""
    loc_table = rotate_collation(t2_loc_table)[1:] # give each witness its own list.
    loc_dict = [] 
    for line in loc_table:
        wit_name = line[0]
        wit_readings = line[1:]
        loc_dict.append((wit_name, wit_readings))
    loc_dict = dict(loc_dict) # transform loc_table into a dictionary,
                              # with witness names as keys.
    wit_set_1 = divided_wit_set[0]
    wit_set_2 = divided_wit_set[1]

    def set_subset(wit_set):
        """returns the correct subset of t2_loc_table for the witness subset"""
        return sub_col(t2_loc_table, wit_set)

    def score_set(wit_set):
        """score a witness set for consistency"""
        if wit_set == []: # if there are no witnesses in the set
            return 0
        set_col = set_subset(wit_set) # get a subset collation table just for this witness set
        num_witnesses = len(set_col[0]) # establish how many witnesses are in the set
        total_score = 0
        for row in set_col[1:]: # go through each row of set_subset apart from the names row
            if len(set(row)) == 1: # if all words in the set are the same
                score = num_witnesses # then the row scores full marks, 1 point per witness.
            else:
                # otherwise, score 1 point per instance of the more frequently represented word:
                word_1 = list(set(row))[0]
                word_2 = list(set(row))[1]
                word_1_freq = row.count(word_1)
                word_2_freq = row.count(word_2)
                if word_1_freq == word_2_freq: # if both words are equally frequent,
                    score = word_1_freq - 1 # let the score be word_1's frequency - 1 (so that a single word by itself scores 0 points)
                elif word_1_freq > word_2_freq:
                    score = word_1_freq - 1
                else:
                    score = word_2_freq - 1
            total_score = total_score + score
        return total_score

    

            # all_words = t2_loc_table[num][1:] # get all instances of the two words for this line
            # word_set = set(all_words) # get the set of these, i.e. just the two words attested across this line
            # word_1 = list(word_set)[0] # get the first word
            # word_2 = list(word_set)[1] # get the second word
            
            
    return score_set(wit_set_1) + score_set(wit_set_2)

def list_all_scores(t2_loc_table):
    witness_set = t2_loc_table[0][1:]
    score_list = []
    for divided_wit_set in pair_combos(witness_combos(witness_set)):
        score = score_type_2s_by_wit_set(t2_loc_table, divided_wit_set)
        # print(divided_wit_set)
        score_list.append([score, divided_wit_set])
    final_list = []
    for row in score_list:
        final_list.append([row[0], row[1][0], row[1][1]])
    return list(reversed(sorted(final_list)))

def list_t2_groupings(witness_list, collation):
    """given a list of witnesses and a collation, return the optimal subdivision for type-2 deviations"""
    type_2_loc_table = type_2_locs(sub_col(collation, witness_list))
    all_scores = list_all_scores(type_2_loc_table)
    return all_scores

def list_groups_to_compare(witness_list):  # 
    """returns a list of every group that must be compared for a given witness list, starting with comparing only 4 witnesses at once, and rising to comparing the total number of witnesses in the group"""
    bifurcation_list = []
    current_group_size = 4
    for group_size in range(current_group_size, len(witness_list) + 1):
        progress_message = "calculating for group size " + str(group_size)
        print(progress_message)
        bifurcation_list += list(combinations(witness_list, current_group_size))
        current_group_size = current_group_size + 1
    return bifurcation_list


def return_winners(witness_list, collation):
    """return the type-2 variation-locations for the witnesses that agree with the highest-scoring bifurcation"""
    witness_texts = rotate_collation(collation)
    def get_witness_text(name):
        for row in witness_texts:
            if row[0] == name:
                return row
        return 'witness name not in collation'
    def group_witness_texts(witness_list):
        """return a list of witness_list's texts, in the order provided"""
        grouped_texts = []
        for witness in witness_list:
            grouped_texts.append(get_witness_text(witness))
        return grouped_texts
    
    def get_variant_list(subgroup, col_row):
        """return a list of all variants for the subgroup at row col_row"""
        variant_list = []
        for name in subgroup:
            variant_list.append(get_witness_text(name)[col_row])
        return variant_list
        
    def count_variants(variant_list):
        """count the number of variant readings in a subgroup for a given row of the collation"""
        return len(set(variant_list))

    # def commoner_variant(variant_list):
    #     """return the more common variant in a variant list containing two variants"""
    #     variants = set(variant_list)
    #     variant1 = variant_list[0]
    #     variant2 = variant_list[1]
    #     variant1_freq = variant_list.count(variant1)
    #     variant2_freq = variant_list.count(variant2)
    #     if variant1_freq == variant2_freq:
    #         return "EQUAL FREQUENCY"
    #     elif variant1_freq > variant2_freq:
    #         return variant1
    #     else:
    #         return variant2
    #     pass

    def process_row_for_subgroup(subgroup, other_subgroup, row_num):
        """do all processes for a subgroup and a row_num"""
        # print(subgroup) # for debug
        group_vars = get_variant_list(subgroup, row_num)
        othergroup_vars = get_variant_list(other_subgroup, row_num)
        how_many_vars = count_variants(group_vars)
        if how_many_vars == 1 and group_vars[0] not in othergroup_vars: # if all readings are identical, and exclusively confined to this subgroup:
            
            row_results = [row_num, subgroup, subgroup] # then all witnesses are winners, since the row agrees with majority opinion.
        else:
            row_results = [row_num, [], subgroup] # then reward no witnesses, as the row does not conform with majority opinion.
        return row_results


    def process_subgroup(subgroup, other_subgroup):
        subgroup_results = []
        for row in variant_table[1:]:
            row_num = row[0]
            subgroup_results.append(process_row_for_subgroup(subgroup, other_subgroup, row_num))
        return subgroup_results


    score_table = list_t2_groupings(witness_list, collation)[:2] # we are only interested in the top two scores
    # print(score_table) # for debug
    if score_table[0][0] == score_table[1][0]:
        return [] # if the two scores are equal, return an empty list.

    variant_table = (type_2_locs(sub_col(collation, witness_list)))
    group1 = score_table[0][1]
    # print(group1) # for debug
    group2 = score_table[0][2]
    # print(group2) # for debug
    scores_to_return = []

    names = witness_list # the order in which the names appear is the same as the order of the witness list.
    set_texts = group_witness_texts(names) # grab the texts in the order used in the variant table.



    group1_results = process_subgroup(group1, group2)
    group2_results = process_subgroup(group2, group1)
    collected_results = group1_results + group2_results

    # # remove all rows for which no witnesses are winners
    # final_results = []
    # if collected_results == []: # if the whole table is empty
    #     final_results = collected_results
    # else:
    #     for row in collected_results:
    #         if row[1] != []:
    #             final_results.append(row)

    return collected_results


def initialize_weighting_table(collation):
    """returns a new weighting table sharing the same structure as the collation table passed to it"""
    weighting_table = [] + [collation[0]] # copy the collation's name row as the first row of the weighting table
    for row in collation[1:]: # and for every row except the first of the collation table:
        new_row = [row[0]]
        for element in row[1:]: # and for every element of the row except for the first (which is the row number)
            new_row.append(0) # append a zero in place of every entry
        weighting_table.append(new_row)
    return weighting_table

def write_winners_to_weighting_table(weighting_table, rounds_table, winners_results):
    """updates WEIGHTING_TABLE and ROUNDS_TABLE with the WINNERS_RESULTS"""
    for row in winners_results:
        rownum = row[0]
        witnesses_who_won = row[1]
        witnesses_who_played = row[2]
        for wit_name in witnesses_who_won:
            wit_index = weighting_table[0].index(wit_name)
            weighting_table[rownum][wit_index] += 1
        for wit_name in witnesses_who_played:
            wit_index = rounds_table[0].index(wit_name)
            rounds_table[rownum][wit_index] += 1

    return # nothing to return, as the weighting_table is permanently edited.

def t2_weighting(collation, witness_set):
    """runs every member of witness_set through RETURN_WINNERS() and returns a weighting table"""
    winners_table = initialize_weighting_table(collation)
    rounds_table = initialize_weighting_table(collation)
    groups_to_compare = list_groups_to_compare(witness_set)
    compare_length_warning = "need to compare " + str(len(groups_to_compare)) + " groups!"
    print(compare_length_warning)
    for group in groups_to_compare:
        winners = return_winners(group, collation)
        write_winners_to_weighting_table(winners_table, rounds_table, winners)

    prediction_table = [winners_table[0]] # make the prediction table's top row the name line.
    for row in zip(winners_table[1:], rounds_table[1:], collation[1:]):
        prediction_row = [row[0][0]] # make first item in row the row number
        for item in zip(row[0][1:], row[1][1:], row[2][1:]): # for all the other items
            predictivity = item
            # if item[1] == 0:
            #     predictivity = 0
            # else:
            #     predictivity = item[0] / item[1] # the number of wins vs. the number of rounds played
            prediction_row.append(predictivity)
        prediction_table.append(prediction_row)
            
            
    return prediction_table

def underline_table_heading(table):
    """puts a line under the first line of a table"""
    output = [table[0]]
    output.append(None)
    output.extend(table[1:])
    return output

def sub_col_rownums(collation, wit_list):
    """return a subcollation of the table with row numbers"""
    subcol = sub_col(collation, wit_list)
    sub_col_rownums = []
    rownum = 0
    for row in subcol:
        new_row = [rownum] + row
        sub_col_rownums.append(new_row)
        rownum = rownum + 1
    return sub_col_rownums

def quartets(wit_list):
    """Returns every possible set of 4 witnesses, for a given witness list"""
    quartets = []
    quartets += list(combinations(wit_list, 4))
    return quartets

def t2_weighting_by_quartets(collation, witness_set):
    """a slimmed-down version of t2_weighting(): exposes every subset's type-2 variation by considering only every possible group of 4 witnesses"""
    winners_table = initialize_weighting_table(collation)
    rounds_table = initialize_weighting_table(collation)
    groups_to_compare = quartets(witness_set)
    compare_length_warning = "need to compare " + str(len(groups_to_compare)) + " groups!"
    print(compare_length_warning)
    for group in groups_to_compare:
        winners = return_winners(group, collation)
        write_winners_to_weighting_table(winners_table, rounds_table, winners)

    prediction_table = [winners_table[0]] # make the prediction table's top row the name line.
    for row in zip(winners_table[1:], rounds_table[1:], collation[1:]):
        prediction_row = [row[0][0]] # make first item in row the row number
        for item in zip(row[0][1:], row[1][1:], row[2][1:]): # for all the other items
            predictivity = item
            prediction_row.append(predictivity)
        prediction_table.append(prediction_row)
    return prediction_table

def sort_collation_by_weight(collation, weighting_table):
    name_row = ['freq'] + collation[0]
    other_rows = []
    for row in zip(weighting_table[1:], collation[1:]):
       total = sum(row[0][1:])
       #updated_row = [total] + row[1]
       updated_row = [total]
       for element in zip(row[0], row[1]):
           updated_row.append(element)
       other_rows.append(updated_row)
    sorted_rows = list(reversed(sorted(other_rows)))
    
    return [name_row] + sorted_rows

# ------------------------------------------------------------------------------
# NEIGHBOUR-JOINING USING QUARTET WEIGHTING
# ------------------------------------------------------------------------------

def quartet_weighted_proportions(col, namelist):
    """return a collation weighted by quartets, but with proportions rather than wins vs. attempts"""
    qweighted_col = t2_weighting_by_quartets(col, namelist)
    prop_col = [qweighted_col[0]]
    for row in qweighted_col[1:]:
        prop_row = [row[0]]
        for cell in row[1:]:
            cell_wins = cell[0]
            cell_tries = cell[1]
            cell_word = cell[2]
            prop_weighted_wins = cell_wins / (cell_tries + 1) # by adding 1 to cell_tries, we avoid a divide by 0 error, and we also reward words more that are extremely predictive.
            new_cell = [prop_weighted_wins, cell_word]
            prop_row.append(new_cell)
        prop_col.append(prop_row)
    return prop_col

# qwp_col = quartet_weighted_proportions(better_col, partial_namelist)

def qw_hamming(wit1, wit2):
    """calculate the Hamming distance between two quartet-weighted
transliterations (extracted from a collation created using
QUARTET_WEIGHTED_PROPORTIONS), expressed as the proportion of
comparable text that is NOT identical

    """
    # internal functions:
    def lacuna_p(a_string):
        """return True if a_string contains the character '[', which indicates
it is lacunose"""
        if '[' in a_string:
            return True
        else:
            return False
    def omission_p(a_string):
        """return True if a_string == the character '‑' (non-breaking hyphen, not normal - sign!!!), which indicates it is
an omission"""
        if '‑' == a_string:
            return True
        else:
            return False

    def comparable_parts(wit1, wit2):
        """take two equally long lists of strings; return the parts that can
be meaningfully used to calculate their Hamming distance."""
        comparable_parts_list = []
        for position in range(len(wit1)):
            pair = [wit1[position][1], wit2[position][1]]
            significance = (wit1[position][0] + wit2[position][0]) / 2.0
            if not lacuna_p(pair[0]) and not lacuna_p(pair[1]): # disregard lacunae
                if not omission_p(pair[0]) and not omission_p(pair[1]): # disregard shared omissions
                    comparable_parts_list.append([pair, significance])
        return comparable_parts_list

    def identical_parts(list_of_comparable_parts):
        """take list produced by comparable_parts(); return list of identical parts"""
        identical_parts_list = []
        for item in list_of_comparable_parts:
            if item[0][0] == item[0][1]:
                identical_parts_list.append(item)
        return identical_parts_list
    
    def calc_length(comparable_parts_list):
        """sums the weightings applied to each row and returns the total"""
        total = 0
        for row in comparable_parts_list:
            total = total + row[1]
        return total
    # main function body:
    comparable_length = calc_length(comparable_parts(wit1, wit2))
    if comparable_length == 0:
        comparable_length = 0.0000001 # a comparable length of 0 will produce a divide by zero error when calculating hamming_distance below; have to use an extremely small value instead.
    identical_length = calc_length(identical_parts((comparable_parts(wit1, wit2))))
    hamming_distance = 1 - (identical_length / comparable_length)
    return hamming_distance

def qw_dist_matrix(collation):
    """takes a collation created with QUARTET_WEIGHTED_PROPORTIONS; returns a distance matrix"""
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
            if other_witnesses[0:] == this_witness[0:]:
                h_dists.append("itself")
            else:
                h_dist = qw_hamming(this_witness[1:], other_witnesses[1:]) # position 0 is just the witness name.
                h_dists.append(h_dist)
        dist_matrix.append(h_dists)
    # remove duplicate entries from the matrix:
    deduped_dist_matrix = []
    for line in dist_matrix:
        deduped_line = []
        reached_end = False
        for entry in line:
            if entry == "itself":
                reached_end = True
            if reached_end == False:
                deduped_line.append(entry)
        deduped_dist_matrix.append(deduped_line)
    return deduped_dist_matrix


def qw_nj_pdf(col, file_name, work_name, caption_postscript):
    namelist = col[0][1:]
    qwp_col = quartet_weighted_proportions(col, namelist)
    dmatrix = qw_dist_matrix(qwp_col)
    dlist = dist_matrix_to_list(dmatrix)
    njed = neighbour_joiner(dlist, [], 1)
    graph_lines = n_j_to_graphviz(njed)
    final_string = ""
    for line in graph_lines:
        final_string = final_string + line + "\n"
    from graphviz import Source
    s = Source(final_string, filename=file_name, format="pdf", engine="neato")
    s.render()
    caption = "#+caption: Neighbour-joiner chain for " + work_name + ". " + caption_postscript
    return caption + '\n' + '#+name: ' + file_name.rsplit('.', maxsplit=1)[0] + '\n#+attr_latex: :placement [t] :width \\textwidth' + '\n' + 'file:' + file_name + ".pdf"

def qw_rob_one_out(col, filename, workname):
    """given a collation, create nj chains where each chain omits one of the collation's witnesses (robustness check), using quartet weighting"""
    namelist = col[0][1:]
    rob_report = str(namelist) + "\n \n"
    # return the full collation
    subset_file_name = filename + "_rob_" + "000"
    subset_post_caption = "No witnesses removed."
    robbed_nj = qw_nj_pdf(col, subset_file_name, workname, subset_post_caption)
    rob_report = rob_report + "\n \n" + robbed_nj
    
    for name in namelist:
        new_namelist = [word for word in namelist if word != name]
        new_col = sub_col_rownums(col, new_namelist) # make a new sub-collation without that witness
        # produce nj chain for that witness subset;
        subset_file_name = filename + "_rob_" + name.replace(" ", "_").replace(".", "_")
        subset_post_caption = name + " removed. Quartet weighting."
        robbed_nj = qw_nj_pdf(new_col, subset_file_name, workname, subset_post_caption)
        rob_report = rob_report + "\n \n" + robbed_nj
    return rob_report


def chain_nolengths(chaincode):
    """Return a CHAINCODE shorn of its length data"""
    chaincode = re.sub('label.*?,', '', chaincode)
    listchain = [line for line in chaincode.splitlines()]
    return listchain



def chain_collapse_nodes(chaincode, node_list, new_node_name):
    """Return a CHAINCODE with NODE_LIST collapsed to a single node called new_node_name"""
    chain_lines = chain_nolengths(chaincode)
    new_chaincode = ""
    for line in chain_lines:
        newline = line
        if "--" in line: # detect whether it's a line you should be replacing
            for node in node_list:
                oldstring = '"' + str(node) + '"'
                newstring = '"' + str(new_node_name) + '"'
                newline = newline.replace(oldstring, newstring)
                # print(newline)
        if newline.count(new_node_name) != 2:
            new_chaincode = new_chaincode + newline + "\n"
            print(newline)
    print(new_chaincode)
    return new_chaincode

def node_collapser(chain_file, new_file, list_of_lists_of_nodes, work_title, caption):
    """Given CHAIN_FILE, return NEW_FILE and a drawer of text for ready insertion into a document, with all the node lists in LIST_OF_LIST_OF_NODES individually collapsed."""
    with open(chain_file, 'r') as file:
        chain = file.read()
        for nodelist in list_of_lists_of_nodes:
            chain = chain_collapse_nodes(chain, nodelist, str(nodelist))
    return gv_string_to_pdf(chain, new_file, work_title, caption, "neato")



# ------------------------------------------------------------------------------
# COLLATION FILLING AND FORMATTING FUNCTIONS
# ------------------------------------------------------------------------------

def fill_collation(collation, refcolnum):
    """takes a collation and fills the empty cells using the values from column <refcolnum>"""
    filled_collation = []
    rownum = 0
    for row in collation:
        refword = row[refcolnum]
        filled_row = []
        cellnum = 0
        for cell in row:
            if cellnum == 0:
                filled_row.append(rownum)
                cellnum = cellnum + 1
                rownum = rownum + 1
            elif cell == '':
                filled_row.append(refword)
            else:
                filled_row.append(cell)
        filled_collation.append(filled_row)
    return filled_collation

def number_collation(collation):
    """takes a collation and renumbers it by row"""
    filled_collation = []
    rownum = 0
    for row in collation:
        filled_row = []
        cellnum = 0
        for cell in row:
            if cellnum == 0:
                filled_row.append(rownum)
                cellnum = cellnum + 1
                rownum = rownum + 1
            else:
                filled_row.append(cell)
        filled_collation.append(filled_row)
    return filled_collation

def ref_collation(raw_collation, filled_collation):
    """takes a raw and a filled collation; returns a reference collation that uses both sets of reference numbers"""
    new_col = []
    rownum = 0
    for row in filled_collation:
        newrow = [raw_collation[rownum][0]]
        newrow.extend(row)
        rownum += 1
        new_col.append(newrow)
    return new_col

# ------------------------------------------------------------------------------
# WRAPPER FUNCTIONS FOR EASE OF USE
# ------------------------------------------------------------------------------

def get_type_2s(collation, witlist):
    """return all type 2 variants within a given collation, for a given list of witnesses"""
    return type_2_locs(sub_col(collation, witlist))

def quartet_weighting(collation, witlist):
    """return the collation weighted by quartets for the given witness list"""
    return sub_col_rownums(t2_weighting_by_quartets(collation, witlist), witlist)

def best_bis(collation, witlist, topnum):
    """return the most parsimonious <topnum> bifurcations for the given witlist, for the given collation"""
    return list_t2_groupings(witlist, collation)[:topnum]

# ------------------------------------------------------------------------------
# DIFF TWO COLLATIONS (witness names must be identical):
# ------------------------------------------------------------------------------

def col_diff(coll1, coll2):
    witnames = coll1[0]
    rownum = 0
    changelist = [['/row/', '/witness/', '/reading/', '', '/normalization/'], None]
    for row in coll1:
        cellnum = 0
        for cell in row:
            coll1_reading = coll1[rownum][cellnum]
            coll2_reading = coll2[rownum][cellnum]
            if coll1_reading != coll2_reading:
                witname = witnames[cellnum]
                # diff = str(rownum) + ": " + witname + ": " + coll1_reading + " → " + coll2_reading
                diff = [rownum, witname, coll1_reading, " → ", coll2_reading]
                changelist.append(diff)
            cellnum = cellnum + 1
        rownum = rownum + 1
    return changelist

# ------------------------------------------------------------------------------
# STRIP A COLLATION OF ALL DETERMINATIVES THAT BEGIN A CELL
# ------------------------------------------------------------------------------
def strip_dets(collation, threshold):
    """takes a collation table; returns the collation table stripped of all lines containing more than <threshold> cells containing determinatives"""
    nodets = [collation[0]] # to contain all rows of the collation containing fewer than <threshold> determinatives
    for row in collation[1:]:
        detcount = 0
        for cell in row:
            # res = any(char.isupper() for char in str(cell)) # with my transliteration scheme, any uppercase character in a cell must be part of a Gardiner sign number.
            # if res == True:
            #safer, just take out transliteration characters on their own row by detecting capital letters at the beginning of the row:
            # print(str(cell)) # for debugging: if you get an index out of range error, probably there's a completely empty cell.
            if str(cell)[0].isupper():
                detcount = detcount + 1
            if "NN" in str(cell):
                detcount = detcount - 1 # remove the penalty for cells containing NN: that just means the name of the manuscript owner.
        if threshold > detcount:
            nodets.append(row)
    return nodets

# ------------------------------------------------------------------------------
# STRIP A COLLATION OF () BRACKETS
# ------------------------------------------------------------------------------
def strip_brackets(collation):
    """strip COLLATION of () brackets (to remove orthographic variants of omission)"""
    stripped = [collation[0]] # start with just the header line.
    for row in collation[1:]:
        newrow = []
        for cell in row:
            cell = re.sub("[\(\)]", "", str(cell))
            newrow.append(cell)
        stripped.append(newrow)
    return stripped

# ------------------------------------------------------------------------------
# RETURN ONLY THOSE CELLS OF A COLLATION THAT CONTAIN DETERMINATIVES
# ------------------------------------------------------------------------------
def only_dets(collation, threshold):
    """takes a collation table; returns the collation table with only those rows containing more than <threshold> cells containing determinatives"""
    nodets = [collation[0]] # to contain all rows of the collation containing fewer than <threshold> determinatives
    for row in collation[1:]:
        detcount = 0
        for cell in row:
            # res = any(char.isupper() for char in str(cell)) # with my transliteration scheme, any uppercase character in a cell must be part of a Gardiner sign number.
            # if res == True:
            #safer, just take out transliteration characters on their own row by detecting capital letters at the beginning of the row:
            # print(str(cell)) # for debugging: if you get an index out of range error, probably there's a completely empty cell.
            if str(cell)[0].isupper():
                detcount = detcount + 1
            if "NN" in str(cell):
                detcount = detcount - 1 # remove the penalty for cells containing NN: that just means the name of the manuscript owner.
        if threshold < detcount:
            nodets.append(row)
    return nodets


# ------------------------------------------------------------------------------
# LACUNOSITY FUNCTIONS
# ------------------------------------------------------------------------------

def lacunosity(text):
    """Takes the text of a collation expressed as a list of words; returns measurements of the text's lacunosity"""
    lacunae_counter = 0
    omission_counter = 0
    for word in text:
        if '[' in word:
            lacunae_counter = lacunae_counter + 1
        if '‑' == word:
            omission_counter = omission_counter + 1
        if  '‑‑' in word: # changed from == to allow tracking of multiple types of significant omissions in a single row (e.g. an isolated omission vs. an omission spanning several lines), by appending a number to ‑‑.
            omission_counter = omission_counter + 1
    percentage = round(lacunae_counter / (len(text) - omission_counter) * 100, 2)
    lacunosity_description = [len(text) - omission_counter, lacunae_counter, percentage]
    return lacunosity_description

def get_lacunosity(collation, witness):
    """Takes a collation and a witness name, and returns the lacunosity of the text as: [<length>, <lines lacunose>, <% lacunose>]"""
    return lacunosity(get_wit_text(collation, witness)[1:]) # because cell 0 is just the witness's name.

def lacunosity_table(collation):
    """returns measurements of lacunosity for a collation's witnesses"""
    wit_list = collation[0][1:]
    lac_table = [['/witness/', '/length/','/lacunose cells/', '/% lacunose/'],None]
    for witness in wit_list:
        line = [witness]
        line.extend(get_lacunosity(collation, witness))
        lac_table.append(line)
    return lac_table

def strip_lacunose_witnesses(collation, threshold):
    """returns a collation stripped of all witnesses whose texts are more than <threshold>% lacunose"""
    lac_table = lacunosity_table(collation)
    lac_list = []
    for row in lac_table[2:]: # because row 0 = titles and row 1 = None (to make horizontal rule)
        lac_list.append([row[0], row[3]]) # append the witness name and the % lacunose to lac_list
    not_too_lacunose = [] # the witness list of witnesses which are not too lacunose
    for entry in lac_list:
        if entry[1] < threshold:
            not_too_lacunose.append(entry[0]) # append the witness name to not_too_lacunose.
    return sub_col_rownums(collation, not_too_lacunose)

def list_lacunose_witnesses(collation, threshold):
    """returns a list of a collation's witnesses which are >= than <threshold>% lacunose"""
    lac_table = lacunosity_table(collation)
    lac_list = []
    for row in lac_table[2:]: # because row 0 = titles and row 1 = None (to make horizontal rule)
        lac_list.append([row[0], row[3]]) # append the witness name and the % lacunose to lac_list
    too_lacunose = [] # the witness list of witnesses which are not too lacunose
    for entry in lac_list:
        if entry[1] >= threshold:
            too_lacunose.append(entry[0]) # append the witness name to not_too_lacunose.
    return too_lacunose

def get_witness_list(collation):
    return collation[0][1:]


def intersection(collation, wit1, wit2):
    "returns a list of the text meaningfully shared by both witnesses"
    subcol = sub_col_rownums(collation, [wit1, wit2])
    intersection = [cell for cell in subcol if cell[1] != '‑' and cell[2] != '‑' and '[' not in cell[1] and '[' not in cell[2]]
    return intersection

def diff_2_wits(collation, wit1, wit2):
    """returns a table of the text that differs between two witnesses"""
    subcol = intersection(collation, wit1, wit2)
    diff = [row for row in subcol if row[1] != row[2]]
    return diff

def diff_table(collation):
    """returns a table showing the amount of text that varies between each witness: firstly the number of actual words, and secondly the proportion of comparable text"""
    witlist = collation[0][1:]
    diff_table = [collation[0][0:-1]]
    for wit in witlist[1:]:
        diff_table_line = [wit]
        for otherwit in witlist[0:-1]:
            if wit == otherwit:
                break
            else:
                num_diffs = len(diff_2_wits(collation, wit, otherwit)[1:])
                wit_text = get_wit_text(collation, wit)[1:]
                otherwit_text = get_wit_text(collation, otherwit)[1:]
                hamming_dist = round(hamming(wit_text, otherwit_text), 3)
                info = [num_diffs, hamming_dist]
                diff_table_line.append(info)
        diff_table.append(diff_table_line)
    return diff_table



def get_mean_units(diff_table):
    """takes a diff_table made with diff_table(); returns the mean number of units (i.e. words) used to produce the Hamming distances in the table"""
    units = []
    for row in diff_table[1:]:
        for cell in row[1:]:
            units.append(cell[0])
    return round(statistics.mean(units), 3)

def get_units(diff_table):
    """takes a diff table and returns a list of the units used to produce Hamming distances in the table, i.e. the number of words that DIFFER from each other for every pair of witnesses"""
    units = []
    for row in diff_table[1:]:
        for cell in row[1:]:
            units.append(cell[0])
    return units

def get_hamming_distances(diff_table):
    """takes a diff table and returns a list of the Hamming distances in the table, i.e. the PROPORTION of words that DIFFER from each other for every pair of witnesses"""
    units = []
    for row in diff_table[1:]:
        for cell in row[1:]:
            units.append(cell[1])
    return units


def intersection_table(collation):
    """returns a table showing the amount of text meaningfully shared by every pair of witnesses in a collation"""
    witlist = collation[0][1:]
    intersection_table = [collation[0][0:-1]]
    for wit in witlist[1:]:
        intersection_table_line = [wit]
        for otherwit in witlist[0:-1]:
            if wit == otherwit:
                break
            else:
                num_intersects = len(intersection(collation, wit, otherwit)[1:])
                intersection_table_line.append(num_intersects)
        intersection_table.append(intersection_table_line)
    return intersection_table

def get_intersection_units(intersection_table):
    """takes an intersection_table made with intersection_table(); returns a list of the numbers of units (i.e. words) that intersect between ever pair of witnesses"""
    units = []
    for row in intersection_table[1:]:
        for cell in row[1:]:
            units.append(cell)
    return units

def get_frequency_of_diff_units(diff_table, intersection_table):
    """returns the frequency with which differing variants occur in every witness pairs' comparable variation places"""
    diff_units = get_units(diff_table)
    inter_units = get_intersection_units(intersection_table)
    freq_list = [round(unit[0] / unit[1], 2) for unit in list(zip(diff_units, inter_units)) if unit[1] != 0]  # have to remove witnesses with 0 comparable material.
    return freq_list

def get_statistics(col, rounding, stripping_threshold, witness_list, population_or_sample):
    """returns a table of statistics about a collation, which could be tagged as 'population' or 'sample': remember to apply lacunosity threshold before calling this!"""
    col_phons = sub_col_rownums(strip_dets(col, stripping_threshold), witness_list)
    col_dets = sub_col_rownums(only_dets(col, stripping_threshold), witness_list)
    # do the above two before the below, so that the stripping threshold for determinatives applies to the whole collation, not just the subcollation.
    col = sub_col_rownums(col, witness_list)

    stat_table = [["variation places","","all data", "phonemic data", "determinative data"], None]
    if population_or_sample == "sample":
        mean_symbol = "x̄"
        stdev_symbol = "s"
        s_div_m_symbol = "s÷x̄"
    else:
        mean_symbol = "μ"
        stdev_symbol = "σ"
        s_div_m_symbol = "σ÷μ"

    comp_units_col = get_intersection_units(intersection_table(col))
    comp_units_phons = get_intersection_units(intersection_table(col_phons))
    comp_units_dets = get_intersection_units(intersection_table(col_dets))

    col_comp_mean = round(statistics.mean(comp_units_col), rounding)
    phons_comp_mean = round(statistics.mean(comp_units_phons), rounding)
    dets_comp_mean = round(statistics.mean(comp_units_dets), rounding)
    stat_table.append(["comparable", mean_symbol, col_comp_mean, phons_comp_mean, dets_comp_mean])

    col_comp_stdev = round(statistics.stdev(comp_units_col), rounding)
    phons_comp_stdev = round(statistics.stdev(comp_units_phons), rounding)
    dets_comp_stdev = round(statistics.stdev(comp_units_dets), rounding)
    stat_table.append(["", stdev_symbol, col_comp_stdev, phons_comp_stdev, dets_comp_stdev])

    stat_table.append(["", s_div_m_symbol, round(col_comp_stdev / col_comp_mean, rounding), round(phons_comp_stdev / phons_comp_mean, rounding), round(dets_comp_stdev / dets_comp_mean, rounding)])

    stat_table.append(None)

    diff_units_col = get_units(diff_table(col))
    diff_units_phons = get_units(diff_table(col_phons))
    diff_units_dets = get_units(diff_table(col_dets))

    col_diff_mean = round(statistics.mean(diff_units_col), rounding)
    phons_diff_mean = round(statistics.mean(diff_units_phons), rounding)
    dets_diff_mean = round(statistics.mean(diff_units_dets), rounding)
    stat_table.append(["differing", mean_symbol, col_diff_mean, phons_diff_mean, dets_diff_mean])

    col_diff_stdev = round(statistics.stdev(diff_units_col), rounding)
    phons_diff_stdev = round(statistics.stdev(diff_units_phons), rounding)
    dets_diff_stdev = round(statistics.stdev(diff_units_dets), rounding)
    stat_table.append(["", stdev_symbol, col_diff_stdev, phons_diff_stdev, dets_diff_stdev])

    stat_table.append(["", s_div_m_symbol, round(col_diff_stdev / col_diff_mean, rounding), round(phons_diff_stdev / phons_diff_mean, rounding), round(dets_diff_stdev / dets_diff_mean, rounding)])

    stat_table.append(None)

    freq_units_col = get_frequency_of_diff_units(diff_table(col), intersection_table(col))
    freq_units_phons = get_frequency_of_diff_units(diff_table(col_phons), intersection_table(col_phons))
    freq_units_dets = get_frequency_of_diff_units(diff_table(col_dets), intersection_table(col_dets))

    col_freq_mean = round(statistics.mean(freq_units_col), rounding)
    phons_freq_mean = round(statistics.mean(freq_units_phons), rounding)
    dets_freq_mean = round(statistics.mean(freq_units_dets), rounding)
    stat_table.append(["freq. differing", mean_symbol, col_freq_mean, phons_freq_mean, dets_freq_mean])

    col_freq_stdev = round(statistics.stdev(freq_units_col), rounding)
    phons_freq_stdev = round(statistics.stdev(freq_units_phons), rounding)
    dets_freq_stdev = round(statistics.stdev(freq_units_dets), rounding)
    stat_table.append(["", stdev_symbol, col_freq_stdev, phons_freq_stdev, dets_freq_stdev])

    stat_table.append(["", s_div_m_symbol, round(col_freq_stdev / col_freq_mean, rounding), round(phons_freq_stdev / phons_freq_mean, rounding), round(dets_freq_stdev / dets_freq_mean, rounding)])
    return stat_table

# ------------------------------------------------------------------------------
# MATPLOTLIB FUNCTIONS
# ------------------------------------------------------------------------------


# Use the histogram functions within a python source block where the
# results are written to a drawer (:results value drawer): this
# enables the filename to be returned cleanly.

def histogram(data, x_label, y_label, font, filename):
    binsize = int(round(2 * (len(data) ** (1. / 3)),0)) # use the Rice Rule to determine bin size: 2 ✕ cube root of the number of data values.
    matplotlib.pyplot.clf()
    matplotlib.rc('axes', axisbelow=False)
    matplotlib.pyplot.style.use('fast')
    matplotlib.rcParams['font.family'] = [font]
    matplotlib.rcParams.update({'font.size': 10})
    matplotlib.pyplot.grid(b=True, which='both', axis='both', linewidth=.5)
    matplotlib.pyplot.hist(data, bins=binsize, linewidth=.5)
    matplotlib.pyplot.xlabel(x_label)
    matplotlib.pyplot.ylabel(y_label)
    matplotlib.pyplot.savefig(filename, format='pdf', bbox_inches="tight")
    return 'file:' + filename

def hamming_distances_histogram(collation, remove_no_comps_p, font, filename, additional_caption_material):
    """Save (to PDF) a histogram of all the Hamming distances in the distance matrix for <collation>"""
    data = get_hamming_distances(diff_table(collation))
    if remove_no_comps_p == True:
        comparable_material = get_intersection_units(intersection_table(collation)) # determine how many units of comparable material there are
        comp_plus_data=list(zip(comparable_material, data))
        stripped_data = [pair[1] for pair in comp_plus_data if pair[0] != 0] # remove distances made on the basis of no comparable material
        data = stripped_data
    histogram(data, 'Hamming distance', 'frequency', font, filename)
    comps_tag = "(witness pairs with no comparable material omitted). " if remove_no_comps_p == True else "(witness pairs with no comparable material included). "
    caption = "Histogram of distance matrix's Hamming distances " + comps_tag + additional_caption_material
    return '#+caption: ' + caption + '\n' + '#+name: ' + filename + '\n#+attr_latex: :placement [t]' + '\n' + 'file:' + filename

def units_for_hamming_distances_histogram(collation, remove_no_comps_p, font, filename, additional_caption_material):
    """Save (to PDF) a histogram of all the units (number of differing readings between each witness pair) used to create the Hamming distances in the distance matrix for <collation>"""
    data = get_units(diff_table(collation))
    if remove_no_comps_p == True:
        comparable_material = get_intersection_units(intersection_table(collation)) # determine how many units of comparable material there are
        comp_plus_data=list(zip(comparable_material, data))
        stripped_data = [pair[1] for pair in comp_plus_data if pair[0] != 0] # remove distances made on the basis of no comparable material
        data = stripped_data
    histogram(data, '№ differing variation places used to calculate Hamming distances', 'frequency', font, filename)
    comps_tag = "(witness pairs with no comparable material omitted). " if remove_no_comps_p == True else "(witness pairs with no comparable material included). "
    caption = "Histogram of the number of differing variation places used to calculate the distance matrix's Hamming distances " + comps_tag + additional_caption_material
    return '#+caption: ' + caption + '\n' + '#+name: ' + filename + '\n#+attr_latex: :placement [t]' + '\n' + 'file:' + filename


def three_histograms(data1, data2, data3, x_label, y_label, font, filename, data1_label, data2_label, data3_label, num_bins):
    binsize1 = int(round(2 * (len(data1) ** (1. / 3)),0)) # use the Rice Rule to determine bin size: 2 ✕ cube root of the number of data values.
    binsize2 = int(round(2 * (len(data2) ** (1. / 3)),0)) # use the Rice Rule to determine bin size: 2 ✕ cube root of the number of data values.
    binsize3 = int(round(2 * (len(data3) ** (1. / 3)),0)) # use the Rice Rule to determine bin size: 2 ✕ cube root of the number of data values.
    matplotlib.pyplot.clf()
    matplotlib.rc('axes', axisbelow=False)
    matplotlib.pyplot.style.use('fast')
    matplotlib.rcParams['font.family'] = [font]
    matplotlib.rcParams.update({'font.size': 10})
    matplotlib.pyplot.grid(b=True, which='both', axis='both', linewidth=.5)
    # matplotlib.pyplot.hist([data1, data2, data3], bins=binsize3, label=[data1_label, data2_label, data3_label])
    #matplotlib.pyplot.hist(data1, bins=binsize1, alpha=1, label=data1_label, edgecolor='blue', linewidth=.5, hatch='OO')
    #matplotlib.pyplot.hist(data2, bins=binsize2, alpha=0.5, label=data2_label, edgecolor='orange', linewidth=.5, hatch='oo')
    #matplotlib.pyplot.hist(data3, bins=binsize3, alpha=0.5, label=data3_label, edgecolor='green', linewidth=.5, hatch='..')
    matplotlib.pyplot.hist([data1, data2, data3], bins=num_bins, alpha=1, label=[data1_label, data2_label, data3_label])
    matplotlib.pyplot.xlabel(x_label)
    matplotlib.pyplot.ylabel(y_label)
    matplotlib.pyplot.legend(loc='best')
    matplotlib.pyplot.savefig(filename, format='pdf', bbox_inches="tight")
    return 'file:' + filename

def hamming_distances_histogram_pda(collation, threshold_dets, threshold_phons, remove_no_comps_p, font, filename, additional_caption_material, num_bins):
    """Save (to PDF) a histogram of all the Hamming distances in the distance matrix for <collation>, for all the data, phonemic, and determinative material"""
    col_all = collation
    col_dets = only_dets(collation, threshold_dets)
    col_phons = strip_dets(collation, threshold_phons)
    data_all = get_hamming_distances(diff_table(col_all))
    data_phon = get_hamming_distances(diff_table(col_phons))
    data_dets = get_hamming_distances(diff_table(col_dets))
    if remove_no_comps_p == True:
        comparable_material_all = get_intersection_units(intersection_table(col_all)) # determine how many units of comparable material there are
        comparable_material_phons = get_intersection_units(intersection_table(col_phons)) # determine how many units of comparable material there are
        comparable_material_dets = get_intersection_units(intersection_table(col_dets)) # determine how many units of comparable material there are
        
        comp_plus_data_all=list(zip(comparable_material_all, data_all))
        comp_plus_data_phons=list(zip(comparable_material_phons, data_phon))
        comp_plus_data_dets=list(zip(comparable_material_dets, data_dets))

        stripped_data_all = [pair[1] for pair in comp_plus_data_all if pair[0] != 0] # remove distances made on the basis of no comparable material
        stripped_data_phons = [pair[1] for pair in comp_plus_data_phons if pair[0] != 0] # remove distances made on the basis of no comparable material
        stripped_data_dets = [pair[1] for pair in comp_plus_data_dets if pair[0] != 0] # remove distances made on the basis of no comparable material

        data_all = stripped_data_all
        data_phon = stripped_data_phons
        data_dets = stripped_data_dets

    three_histograms(data_all, data_phon, data_dets, 'Hamming distance', '№ Hamming distances', font, filename, 'all data', 'phonemic data', 'determinative data', num_bins)
    comps_tag = "(witness pairs with no comparable material omitted). " if remove_no_comps_p == True else "(witness pairs with no comparable material included). "
    caption = "Histogram of distance matrix's Hamming distances for all data, phonemic data, and determinative data " + comps_tag + additional_caption_material
    return '#+caption: ' + caption + '\n' + '#+name: ' + filename + '\n#+attr_latex: :placement [t]' + '\n' + 'file:' + filename


def hamming_distances_histogram_pda_units(collation, threshold_dets, threshold_phons, remove_no_comps_p, font, filename, additional_caption_material, num_bins):
    """Save (to PDF) a histogram of all the Hamming distances in the distance matrix for <collation>, for all the data, phonemic, and determinative material"""
    col_all = collation
    col_dets = only_dets(collation, threshold_dets)
    col_phons = strip_dets(collation, threshold_phons)
    data_all = get_units(diff_table(col_all))
    data_phon = get_units(diff_table(col_phons))
    data_dets = get_units(diff_table(col_dets))
    if remove_no_comps_p == True:
        comparable_material_all = get_intersection_units(intersection_table(col_all)) # determine how many units of comparable material there are
        comparable_material_phons = get_intersection_units(intersection_table(col_phons)) # determine how many units of comparable material there are
        comparable_material_dets = get_intersection_units(intersection_table(col_dets)) # determine how many units of comparable material there are
        
        comp_plus_data_all=list(zip(comparable_material_all, data_all))
        comp_plus_data_phons=list(zip(comparable_material_phons, data_phon))
        comp_plus_data_dets=list(zip(comparable_material_dets, data_dets))

        stripped_data_all = [pair[1] for pair in comp_plus_data_all if pair[0] != 0] # remove distances made on the basis of no comparable material
        stripped_data_phons = [pair[1] for pair in comp_plus_data_phons if pair[0] != 0] # remove distances made on the basis of no comparable material
        stripped_data_dets = [pair[1] for pair in comp_plus_data_dets if pair[0] != 0] # remove distances made on the basis of no comparable material

        data_all = stripped_data_all
        data_phon = stripped_data_phons
        data_dets = stripped_data_dets

    three_histograms(data_all, data_phon, data_dets, '№ differing variation places used to calculate Hamming distances', '№ Hamming distances', font, filename, 'all data', 'phonemic data', 'determinative data', num_bins)
    comps_tag = "(witness pairs with no comparable material omitted). " if remove_no_comps_p == True else "(witness pairs with no comparable material included). "
    caption = "Histogram of the number of differing variation places used to calculate the distance matrix's " + str(len(data_all)) + " Hamming distances for all data, phonemic data, and determinative data " + comps_tag + additional_caption_material
    return '#+caption: ' + caption + '\n' + '#+name: ' + filename.rsplit('.', maxsplit=1)[0] + '\n#+attr_latex: :placement [t]' + '\n' + 'file:' + filename

def scatter_2_cols(col1, col2, label1, label2, font, filename):
    units1 = get_units(diff_table(col1))
    units2 = get_units(diff_table(col2))
    hammings1 = get_hamming_distances(diff_table(col1))
    hammings2 = get_hamming_distances(diff_table(col2))
    matplotlib.pyplot.clf()
    matplotlib.rc('axes', axisbelow=True)
    matplotlib.pyplot.style.use('fast')
    matplotlib.rcParams['font.family'] = [font]
    matplotlib.rcParams.update({'font.size': 10})
    matplotlib.pyplot.grid(b=True, which='both', axis='both', linewidth=.5)
    matplotlib.pyplot.scatter(units1, hammings1, s=2, label=label1)
    matplotlib.pyplot.scatter(units2, hammings2, s=.5, label=label2)
    matplotlib.pyplot.xlabel('№ differing variation places used to calculate Hamming distances')
    matplotlib.pyplot.ylabel('Hamming distance')
    matplotlib.pyplot.legend(loc='best')
    matplotlib.pyplot.savefig(filename, format='pdf', bbox_inches='tight')
    caption = "Scatter graph:  Hamming distances for " + label1 + " vs. " + label2 + "."
    return '#+caption: ' + caption + '\n' + '#+name: ' + filename.rsplit('.', maxsplit=1)[0] + '\n#+attr_latex: :placement [t]' + '\n' + 'file:' + filename




def scatter_3_cols(col1, col2, col3, label1, label2, label3, font, filename, caption_postscript, plot_lines_p):
    units1 = get_units(diff_table(col1))
    units2 = get_units(diff_table(col2))
    units3 = get_units(diff_table(col3))    
    hammings1 = get_hamming_distances(diff_table(col1))
    hammings2 = get_hamming_distances(diff_table(col2))
    hammings3 = get_hamming_distances(diff_table(col3))
    matplotlib.pyplot.clf()
    matplotlib.rc('axes', axisbelow=True)
    matplotlib.pyplot.style.use('fast')
    matplotlib.rcParams['font.family'] = [font]
    matplotlib.rcParams.update({'font.size': 10})
    matplotlib.pyplot.grid(b=True, which='both', axis='both', linewidth=.5)
    if plot_lines_p == True:
        max_units = max(units1 + units2 + units3)
        lines_x = list(zip(units1, units2, units3))
        lines_y = list(zip(hammings1, hammings2, hammings3))
        lines_xy_coords = list(zip(lines_x, lines_y))
        cycol = cycle('bgrcmk')
        alpha = cycle([0.2, 0.4, 0.6, 0.8])
        for coords in lines_xy_coords:
            # distance = math.dist(coords)
            x_coords = list(coords[0])
            y_coords = list(coords[1])
            x_coords.append(x_coords[0])
            y_coords.append(y_coords[0])
            matplotlib.pyplot.plot(x_coords, y_coords, color=next(cycol), alpha=0.5, zorder=0, linewidth=.5)
    matplotlib.pyplot.scatter(units1, hammings1, s=12, label=label1)
    matplotlib.pyplot.scatter(units2, hammings2, s=6, label=label2)
    matplotlib.pyplot.scatter(units3, hammings3, s=4, label=label3)
    matplotlib.pyplot.xlabel('№ differing variation places used to calculate Hamming distances')
    matplotlib.pyplot.ylabel('Hamming distance')
    matplotlib.pyplot.legend(loc='best')
    matplotlib.pyplot.savefig(filename, format='pdf', bbox_inches='tight')
    caption = "Scatter graph:  Hamming distances for " + label1 + " vs. " + label2 + " vs. " + label3 + ". Coloured triangles connect Hamming distances referencing the same witness pairs. " + caption_postscript
    return '#+caption: ' + caption + '\n' + '#+name: ' + filename.rsplit('.', maxsplit=1)[0] + '\n#+attr_latex: :placement [t]' + '\n' + 'file:' + filename

def scatter_3_cols_mean_vp(col1, col2, col3, label1, label2, label3, font, filename, caption_postscript, plot_lines_p):
    units1 = get_units(diff_table(col1))
    units2 = get_units(diff_table(col2))
    units3 = get_units(diff_table(col3))
    mean_units = [statistics.mean([float(unitlist[0]), float(unitlist[1]), float(unitlist[2])]) for unitlist in list(zip(units1, units2, units3))]
    hammings1 = get_hamming_distances(diff_table(col1))
    hammings2 = get_hamming_distances(diff_table(col2))
    hammings3 = get_hamming_distances(diff_table(col3))
    matplotlib.pyplot.clf()
    matplotlib.rc('axes', axisbelow=True)
    matplotlib.pyplot.style.use('fast')
    matplotlib.rcParams['font.family'] = [font]
    matplotlib.rcParams.update({'font.size': 10})
    matplotlib.pyplot.grid(b=True, which='both', axis='both', linewidth=.5)
    if plot_lines_p == True:
        max_units = max(mean_units)
        lines_x = mean_units
        lines_y = list(zip(hammings1, hammings2, hammings3))
        lines_xy_coords = sorted(list(zip(lines_x, lines_y)),  key=lambda item: item[1])
        cycol = cycle('bgrcmk')
        alpha = cycle([0.2, 0.4, 0.6, 0.8])
        for coords in lines_xy_coords:
            x_coords = [coords[0], coords[0], coords[0]]
            y_coords = list(coords[1])
            distance = max(y_coords) - min(y_coords)
            x_coords.append(x_coords[0])
            y_coords.append(y_coords[0])
            matplotlib.pyplot.plot(x_coords, y_coords, color=next(cycol), alpha=max(distance, .1), zorder=0, linewidth=3)
    matplotlib.pyplot.scatter(mean_units, hammings1, s=3, label=label1)
    matplotlib.pyplot.scatter(mean_units, hammings2, s=2, label=label2)
    matplotlib.pyplot.scatter(mean_units, hammings3, s=1, label=label3)
    matplotlib.pyplot.xlabel('Mean № differing variation places used to calculate Hamming distances')
    matplotlib.pyplot.ylabel('Hamming distance')
    matplotlib.pyplot.legend(loc='best')
    matplotlib.pyplot.savefig(filename, format='pdf', bbox_inches='tight')
    caption = "Scatter graph:  Hamming distances for " + label1 + " vs. " + label2 + " vs." + label3 + ". " + caption_postscript
    return '#+caption: ' + caption + '\n' + '#+name: ' + filename.rsplit('.', maxsplit=1)[0] + '\n#+attr_latex: :placement [t]' + '\n' + 'file:' + filename


def scatter_3_cols_scaled(col1, col2, col3, label1, label2, label3, font, filename, caption_postscript, plot_lines_p):
    units1 = get_units(diff_table(col1))
    units2 = get_units(diff_table(col2))
    units3 = get_units(diff_table(col3))
    units1 = [round(float(unit) / max(units1), 4) for unit in units1]
    units2 = [round(float(unit) / max(units2), 4) for unit in units2]
    units3 = [round(float(unit) / max(units3), 4) for unit in units3]
    hammings1 = get_hamming_distances(diff_table(col1))
    hammings2 = get_hamming_distances(diff_table(col2))
    hammings3 = get_hamming_distances(diff_table(col3))
    matplotlib.pyplot.clf()
    matplotlib.rc('axes', axisbelow=True)
    matplotlib.pyplot.style.use('fast')
    matplotlib.rcParams['font.family'] = [font]
    matplotlib.rcParams.update({'font.size': 10})
    matplotlib.pyplot.grid(b=True, which='both', axis='both', linewidth=.5)
    if plot_lines_p == True:
        max_units = max(units1 + units2 + units3)
        lines_x = list(zip(units1, units2, units3))
        lines_y = list(zip(hammings1, hammings2, hammings3))
        lines_xy_coords = list(zip(lines_x, lines_y))
        cycol = cycle('bgrcmk')
        alpha = cycle([0.2, 0.4, 0.6, 0.8])
        for coords in lines_xy_coords:
            # distance = math.dist(coords)
            x_coords = list(coords[0])
            y_coords = list(coords[1])
            x_coords.append(x_coords[0])
            y_coords.append(y_coords[0])
            matplotlib.pyplot.plot(x_coords, y_coords, color='gray', alpha=0.2, zorder=0, linewidth=.3)
    matplotlib.pyplot.scatter(units1, hammings1, s=3, label=label1)
    matplotlib.pyplot.scatter(units2, hammings2, s=2, label=label2)
    matplotlib.pyplot.scatter(units3, hammings3, s=1, label=label3)
    matplotlib.pyplot.xlabel('Proportion of differing variation places used to calculate Hamming distances')
    matplotlib.pyplot.ylabel('Hamming distance')
    matplotlib.pyplot.legend(loc='best')
    matplotlib.pyplot.savefig(filename, format='pdf', bbox_inches='tight')
    caption = "Scatter graph: Hamming distances for " + label1 + " vs. " + label2 + " vs. " + label3 + ". " + caption_postscript
    return '#+caption: ' + caption + '\n' + '#+name: ' + filename.rsplit('.', maxsplit=1)[0] + '\n#+attr_latex: :placement [t]' + '\n' + 'file:' + filename



def lac_chart(collation, font, filename, caption_postscript, width, height, work_name):
    """Returns a bar chart giving lacunosity values for <col>'s witnesses"""
    lac_table = lacunosity_table(collation)
    witnesses = [row[0] for row in lac_table[2:]]
    lengths = [row[1] for row in lac_table[2:]]
    lac_cells = [row[2] for row in lac_table[2:]]
    pc_lacs = [row[3] for row in lac_table[2:]]
    unlac_lengths = [pair[0] - pair[1] for pair in list(zip(lengths, lac_cells))]
    matplotlib.pyplot.clf()
    matplotlib.pyplot.figure(figsize=(width,height))
    matplotlib.rc('axes', axisbelow=True)
    matplotlib.pyplot.style.use('fast')
    matplotlib.rcParams['font.family'] = [font]
    matplotlib.rcParams.update({'font.size': 10})
    matplotlib.pyplot.grid(b=True, which='both', axis='both', linewidth=.5)
    matplotlib.pyplot.barh(range(len(unlac_lengths)), unlac_lengths, 0.8, 0, label='extant variation places')
    matplotlib.pyplot.barh(range(len(lac_cells)), lac_cells, 0.8, unlac_lengths, label='lacunose variation places')
    matplotlib.pyplot.yticks(range(len(witnesses)), witnesses, rotation='horizontal')
    for index,data in enumerate(lengths):
        if pc_lacs[index] != 0.0:
            matplotlib.pyplot.text(x = data, y = index, s = f"{pc_lacs[index]}%", color='black', ha='right', va='center_baseline', fontsize='small')
    #matplotlib.pyplot.legend(loc='best')
    matplotlib.pyplot.legend(loc='upper center', bbox_to_anchor=(0.5, 1.07), ncol=2)
    matplotlib.pyplot.xlabel('№ variation places')
    matplotlib.pyplot.savefig(filename, format='pdf', bbox_inches='tight')
    caption = 'Witness lacunosity for ' + work_name + '.'
    return '#+caption: ' + caption + '\n' + '#+name: ' + filename.rsplit('.', maxsplit=1)[0] + '\n#+attr_latex: :placement [t]' + '\n' + 'file:' + filename
    return witnesses

def two_histograms(data1, data2, x_label, y_label, font, filename, data1_label, data2_label, num_bins):
#    binsize1 = int(round(2 * (len(data1) ** (1. / 3)),0)) # use the Rice Rule to determine bin size: 2 ✕ cube root of the number of data values.
#    binsize2 = int(round(2 * (len(data2) ** (1. / 3)),0)) # use the Rice Rule to determine bin size: 2 ✕ cube root of the number of data values.
    matplotlib.pyplot.clf()
    matplotlib.rc('axes', axisbelow=False)
    matplotlib.pyplot.style.use('fast')
    matplotlib.rcParams['font.family'] = [font]
    matplotlib.rcParams.update({'font.size': 10})
    matplotlib.pyplot.grid(b=True, which='both', axis='both', linewidth=.5)
    matplotlib.pyplot.hist([data1, data2], bins=num_bins, alpha=1, label=[data1_label, data2_label])
    matplotlib.pyplot.xlabel(x_label)
    matplotlib.pyplot.ylabel(y_label)
    matplotlib.pyplot.legend(loc='best')
    matplotlib.pyplot.savefig(filename, format='pdf', bbox_inches="tight")
    return 'file:' + filename

def compare_two_histograms(col1, col2, remove_no_comps_p, font, filename, additional_caption_material, col1_label, col2_label, bins):
    """Save (to PDF) a histogram of all the Hamming distances in the distance matrix for <collation>, for all the data, phonemic, and determinative material"""
    data_col1 = get_hamming_distances(diff_table(col1))
    data_col2 = get_hamming_distances(diff_table(col2))
    if remove_no_comps_p == True:
        comparable_material_col1 = get_intersection_units(intersection_table(col1)) # determine how many units of comparable material there are
        comparable_material_col2 = get_intersection_units(intersection_table(col2)) # determine how many units of comparable material there are
        
        comp_plus_data_col1=list(zip(comparable_material_col1, data_col1))
        comp_plus_data_col2=list(zip(comparable_material_col2, data_col2))

        stripped_data_col1 = [pair[1] for pair in comp_plus_data_col1 if pair[0] != 0] # remove distances made on the basis of no comparable material
        stripped_data_col2 = [pair[1] for pair in comp_plus_data_col2 if pair[0] != 0] # remove distances made on the basis of no comparable material

        data_col1 = stripped_data_col1
        data_col2 = stripped_data_col2

    two_histograms(data_col1, data_col2, 'Hamming distance', '№ Hamming distances', font, filename, col1_label, col2_label, bins)
    comps_tag = "(witness pairs with no comparable material omitted). " if remove_no_comps_p == True else "(witness pairs with no comparable material included). "
    caption = "Histogram of distance matrix's Hamming distances for " + col1_label + " vs. " + col2_label + " " + comps_tag + additional_caption_material
    return '#+caption: ' + caption + '\n' + '#+name: ' + filename.rsplit('.', maxsplit=1)[0] + '\n#+attr_latex: :placement [t]' + '\n' + 'file:' + filename

# ------------------------------------------------------------------------------
# MORE PARSIMONY FUNCTIONS
# ------------------------------------------------------------------------------

def find_witdetails(dates_table, witname, stringstoremove_list):
    """return the details for the correct manuscript, provided with a given witness name: all strings in <stringstoremove_list> will be removed before searching"""
    stripped_witname = witname
    for string in stringstoremove_list:
        stripped_witname = stripped_witname.replace(string, '')
    return [line for line in dates_table if stripped_witname in line[0]][0]

def make_witness_order_table(dates_table, lacunosity_table, max_lac, witlist, stringstoremove_dates_table, stringstoremove_lac_table):
    """return the optimal order for analysing witnesses, given a dates_table (produced using function in witness database org file), a lacunosity_table (produced using bakenseshat.lacunosity_table(col), the <max_lac> (maximum tolerable lacunosity), the <witlist> a list of witnesses, and <stringstoremove_list_dates_table> and <stringstoremove_lac_table>, two lists of strings to remove from <witlist> so that the names are found in <dates_table> and <lacunosity_table> respectively.)"""
    witorder = []
    for name in witlist:
        witline = [name]
        witdates = find_witdetails(dates_table, name, stringstoremove_dates_table)
        print(witdates[0:2])
        witline.extend(witdates[0:2])
        witline.extend(find_witdetails(lacunosity_table, name, stringstoremove_lac_table)[3:])
        witorder.append(witline)
    sorted_witorder = sorted(witorder, key=lambda i: i[2], reverse=True)
    sorted_witorder = sorted(sorted_witorder, key=lambda i: i[3] if (i[3] >= max_lac) else 0.0 , reverse=False)
    table_headers = [['witness', 'manuscript', 'mean date (ʙ.ᴄ.)', 'lacunosity'], None]
    table_headers.extend(sorted_witorder)
    return table_headers

def strip_unweighted(weighted_collation):
    """remove all rows of a weighted collation (created using t2_weighting_by_quartets(sub_col_rownums(col, witlist))) where no's weight has been calculated"""
    new_collation = [weighted_collation[0]]
    for row in weighted_collation[1:]:
        no_calcs = True
        for cell in row[1:]:
            if cell[1] != 0:
                no_calcs = False
        if no_calcs == False:
            new_collation.append(row)
    return new_collation

def t2_weighting_by_quartets_stripped_unweighted(collation, witlist):
    """return a quartet-weighted collation for witness list <witlist> with all rows removed where no weighting has been applied"""
    weighted_col = t2_weighting_by_quartets(sub_col_rownums(collation, witlist), witlist)
    stripped_col = strip_unweighted(weighted_col)
    return stripped_col

def added_type_2_locs(collation, old_witlist, new_witlist):
    """returns the type 2 variation places added in old_witlist which are not present in new_witlist"""
    new_table = t2_weighting_by_quartets_stripped_unweighted(collation, new_witlist)
    old_table = t2_weighting_by_quartets_stripped_unweighted(collation, old_witlist)
    old_type2_locs = set([row[0] for row in old_table])
    all_new_type2_locs = set([row[0] for row in new_table])
    added_type_2_locs = sorted(all_new_type2_locs - old_type2_locs)
    added_type_2_loc_list = [row for row in new_table if row[0] in added_type_2_locs]
    return [new_table[0]] + (added_type_2_loc_list)


def list_variants(wbq_str_unw, index):
    """returns a list of variants for a given row number of the collation, in a t2_weighting_by_quartets_stripped_unweighted table"""
    line_of_interest = [row for row in wbq_str_unw if row[0] == index][0]
    # line_of_interest_stripped = [item[2] for item in line_of_interest[1:] if item[2] not in ["[...]", "‑", "‑‑"]] # do not include fully lacunose witnesses here, or witnesses with significant or insignificant omissions.
    line_of_interest_stripped = [item[2] for item in line_of_interest[1:]]
    witline = wbq_str_unw[0][1:]
    reads = set(line_of_interest_stripped) # the set of available readings
    ret_list = sorted([[read,[]] for read in reads])
    wit_reads = list(zip(witline, line_of_interest_stripped)) # list of each individual witness's reading, in the format (witness, reading)
    for item in wit_reads:
        witname = item[0]
        witread = item[1]
        for item in ret_list:
            if witread == item[0]:
                item[1].append(witname)
    str_ret_list = ""
    for cell in ret_list:
        if cell[0] not in ["[...]", "‑", "‑‑"]:
            print(cell)
            witness_string = ""
            for witness in cell[1]:
                witness_string = witness_string + witness + " "
            string_cell = "▌" + " *" + str(cell[0]) + "* " + ": " + witness_string + "▐"
            str_ret_list = str_ret_list + string_cell
    return str_ret_list

def next_wit_report(collation, old_witlist, new_witlist):
    """returns a blank next witness report for the collation (i.e. for the last witness in the collation)"""
    stripped_quartets = t2_weighting_by_quartets_stripped_unweighted(collation, new_witlist)
    wit_report = [[row[0], row[-1]] for row in stripped_quartets]
    wit_report.insert(1, None)
    wit_report[0].insert(2, 'add?')
    wit_report[0].insert(3, 'topology')
    wit_report[0].insert(4, 'goes like')
    wit_report[0].insert(5, '[...]')
    wit_report[0].insert(6, 'variants')
    wit_report[0].insert(7, '‑‑')
    wit_report[0].insert(8, '‑')


    new_table = stripped_quartets
    old_table = t2_weighting_by_quartets_stripped_unweighted(collation, old_witlist)
    old_type2_locs = set([row[0] for row in old_table])
    all_new_type2_locs = set([row[0] for row in new_table])
    added_type_2_locs = sorted(all_new_type2_locs - old_type2_locs)


    for row in wit_report[2:]:
        if row[0] in added_type_2_locs:
            # print(row)
            row.append('ADD')
        else:
            row.append('')

    total_goes_likes = []
    total_lacunoses = []
    for index, row in enumerate(wit_report[2:]):
        goes_likes = []
        lacunoses = []
        sig_oms = []
        insig_oms = []
        word_to_match = row[1][2] # find the word to match
        equiv_collation_row = new_table[index + 1]
        for cell_index, cell in enumerate(equiv_collation_row[1:-1]):
            if cell[2] == word_to_match:
                matched_witness_name = new_table[0][cell_index + 1]
                goes_likes.append(matched_witness_name)
                total_goes_likes.extend([str(matched_witness_name)])
            if cell[2] == '[...]':
                matched_witness_name = new_table[0][cell_index + 1]
                lacunoses.append(matched_witness_name)
                total_lacunoses.extend([str(matched_witness_name)])
            if cell[2] == '‑':
                matched_witness_name = new_table[0][cell_index + 1]
                insig_oms.append(matched_witness_name)
            if cell[2] == '‑‑':
                matched_witness_name = new_table[0][cell_index + 1]
                sig_oms.append(matched_witness_name)
        row.append('')
        row.append(" ".join(goes_likes))
        row.append(" ".join(lacunoses))
        row.append(list_variants(stripped_quartets, row[0]))
        row.append(" ".join(sig_oms))
        row.append(" ".join(insig_oms))

    goes_likes_counts = []
    for item in set(total_goes_likes):
        item_count = total_goes_likes.count(item) - total_lacunoses.count(item)
        goes_likes_counts.append([item, item_count])

    goes_like_counts = sorted(goes_likes_counts, key=lambda x: x[1], reverse=True)
    wit_report.append(None)
    wit_report.append(['', '', '', '', 'counts'])
    for count in goes_like_counts:
        wit_report.append(['', '', '', '', count])
                    
    return wit_report

def reorder_wit_report(wit_report):
    """Reorders a witness report according to the topology indications"""
    new_tab = [wit_report[0], None]
    pre_ordered_tab = []
    for row in wit_report[1:]:
        if row != None:
            pre_ordered_tab.append(row)
    sorted_tab = sorted(pre_ordered_tab, key=lambda x: x[3]) # reorder according to topology row
    new_tab.extend(sorted_tab)
    return new_tab

# ------------------------------------------------------------------------------
# COLLATION SEARCH FUNCTIONS
# ------------------------------------------------------------------------------

def str_in_row_p(col, string_list):
    """return all rows containing at least one instance of a string as part of the contents of at least one cell"""
    subcol = [col[0]]
    for row in col[1:]:
        include_row_p = False
        for cell in row:
            for string in string_list:
                if string in str(row):
                    include_row_p = True
        if include_row_p == True:
            subcol.append(row)
    return subcol

def detect_empties(col):
    """return a list of the lines in which a collation has empty cells"""
    rownums = []
    for row in col[1:]:
        for cell in row:
            if str(cell) == "":
                rownums.append(row[0])
    return rownums

# ------------------------------------------------------------------------------
# MISC FUNCTIONS
# ------------------------------------------------------------------------------

def gv_string_to_pdf(gv_string, file_name, work_name, caption_postscript, render_engine):
    """produce a pdf from a string containing graphviz code"""
    from graphviz import Source
    s = Source(gv_string, filename=file_name, format="pdf", engine=render_engine)
    s.render()
    caption = "#+caption: Neighbour-joiner chain for " + work_name + ". " + caption_postscript
    return caption + '\n' + '#+name: ' + file_name.rsplit('.', maxsplit=1)[0] + '\n#+attr_latex: :placement [t] :width \\textwidth' + '\n' + 'file:' + file_name + ".pdf"

def wits_agree_others_dont(col, witlist):
    """return a table showing where the witnesses in WITLIST all agree
    with each other, but the other witnesses have at least
    one competing reading"""
    agreements = [col[0]]
    coldict = {witness: col[0].index(witness) for witness in col[0][1:]}
    for row in col[1:]:
        compare_list = []
        for witness in witlist:
            compare_list.append(row[coldict[witness]])
        if len(set(compare_list)) == 1 and len(set(row[1:])) > 1:
            agreements.append(row)
    return agreements

def compare_diffs(col, main_witness, others_witlist):
    """Compare the variation places of MAIN_WITNESS in COL(lation)
    with those in OTHERS_WITLIST"""
    witlist = [main_witness] + others_witlist
    others_witlist.append(main_witness)
    print(witlist)
    subcol_rownums = sub_col_rownums(col, witlist)
    toprow = subcol_rownums[0].append(main_witness + " goes like:")
    combos = list(combinations(witlist, 2))
    variation_places = []
    for combo in combos:
        combovars = wits_agree_others_dont(subcol_rownums, combo)
        variation_places = variation_places + combovars
    sorted_variation_places = sorted(variation_places, key=lambda x: x[0])
    set_varplaces = list(row for row,_ in groupby(sorted_variation_places))

    vps_compared = []

    def is_comparable(variant):
        """returns True if variant is comparable, False, if not"""
        if "[" not in variant and variant != "‑":
            return True
        else:
            return False

    mainwit_index = set_varplaces[0].index(main_witness)
    for row in set_varplaces:
        goeslikes = ""
        mainwit_reading = row[mainwit_index]
        reading_index = 0
        for reading in row:
            if reading == mainwit_reading and is_comparable(reading):
                witname = set_varplaces[0][reading_index]
                if witname != main_witness:
                    goeslikes = goeslikes + " " + witname
            reading_index = reading_index + 1
        row.append(goeslikes)
    return set_varplaces


# ------------------------------------------------------------------------------
# MDS FUNCTIONS
# ------------------------------------------------------------------------------
import sklearn
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import numpy as np


def sk_dmatrix(distance_matrix):
    dm = distance_matrix
    skdm = []
    for row in dm[1:]: # bypass top row of names
        newrow = []
        for cell in row[1:]: # bypass name in first cell
            if cell == "itself":
                newrow.append(0)
            else:
                newrow.append(cell)
        skdm.append(newrow)
    return sklearn.utils.check_array(skdm) # return a sklearn-suitable array.

def MDS_figgen(x, y, xlabel, ylabel, filename, workname, caption_postscript, col, stress, coords, comps, metric_p, normstress_p):
    """helper function: generates MDS figures for MDS_2d_plot and MDS_3d_plot"""
    namelist = col[0][1:]

    # generate the figure:
    plt.clf()
    matplotlib.rcParams['font.family'] = 'Charis SIL'
    matplotlib.rcParams.update({'font.size' : 10})
    matplotlib.rcParams["figure.figsize"] = (8, 8)
    matplotlib.pyplot.grid(b=True, which='both', axis='both', linewidth=.5)
    matplotlib.pyplot.xlabel(xlabel)
    matplotlib.pyplot.ylabel(ylabel)
    plt.plot(x, y, 'o', color='black')
    plt.gca().set_aspect('equal', adjustable='box')

    for num, name in enumerate(namelist):
        plt.annotate(name, (x[num]+0.004, y[num]+0.004))
        # plt.show()
    # save the figure:
    plt.savefig(filename, format='pdf')

    if metric_p == True:
        metric_or_nonmetric = "Metric "
    else:
        metric_or_nonmetric = "Nonmetric "

    if normstress_p == True:
        normstressornot = "Normalized"
    else:
        normstressornot = "Raw"
    caption = metric_or_nonmetric + "MDS plot for " + workname + ". " + str(comps) + " components. " + normstressornot + " stress = " + str(stress) + ". " + caption_postscript
    return '#+caption: ' + caption + '\n' + '#+name: ' + filename.rsplit('.', maxsplit=1)[0] + '\n#+attr_latex: :placement [t] :width \\textwidth' + '\n' + 'file:' + filename


def MDS_2d_plot(col, metric_p, normstress_p, workname, caption_postscript, filename):
    """return a Multidimensional Scaling two-dimensional graph of the Hamming distance matrix of COL, where METRIC is either True or False"""

    # set up MDS for taking a precompiled distance matrix:
    comps = 2
    from sklearn.manifold import MDS
    # mds = MDS(dissimilarity='precomputed', random_state=0, metric=metric_p, n_components=comps)
    mds = MDS(dissimilarity='precomputed', random_state=0, metric=metric_p, n_components=comps, normalized_stress=normstress_p)
    dm = full_dist_matrix(col) # get a full distance matrix of COL
    X = np.array(sk_dmatrix(dm))

    # generate the coords:
    coords = mds.fit_transform(X)
    # nb: if you increase n_components above, this will increase the number of dimensions here!
    
    stress = round(mds.stress_, 3) # calculate the stress value, i.e. how much the coordinates fail to represent reality.

    # rotate the data:
    # this part taken from: https://scikit-learn.org/stable/auto_examples/manifold/plot_mds.html#example-manifold-plot-mds-py
    from sklearn.decomposition import PCA
    clf = PCA(n_components=comps)
    coords = clf.fit_transform(coords)

    # assign coordinates to x and y lists:
    x = [cell[0] for cell in coords]
    y = [cell[1] for cell in coords]

    returned = MDS_figgen(x, y, 'Component 1', 'Component 2', filename + ".pdf", workname, caption_postscript, col, stress, coords, comps, metric_p, normstress_p)
    return returned


def MDS_3d_plot(col, metric_p, normstress_p, workname, caption_postscript, filename):
    """return a Multidimensional Scaling three-dimensional graph of the Hamming distance matrix of COL, where METRIC is either True or False"""

    # set up MDS for taking a precompiled distance matrix:
    comps = 3
    from sklearn.manifold import MDS
    # mds = MDS(dissimilarity='precomputed', random_state=0, metric=metric_p, n_components=comps)
    mds = MDS(dissimilarity='precomputed', random_state=0, metric=metric_p, n_components=comps, normalized_stress=normstress_p)

    dm = full_dist_matrix(col) # get a full distance matrix of COL
    X = np.array(sk_dmatrix(dm))

    # generate the coords:
    coords = mds.fit_transform(X)
    # nb: if you increase n_components above, this will increase the number of dimensions here!
    
    stress = round(mds.stress_, 3)

    # rotate the data:
    # this part taken from: https://scikit-learn.org/stable/auto_examples/manifold/plot_mds.html#example-manifold-plot-mds-py
    from sklearn.decomposition import PCA
    clf = PCA(n_components=comps)
    coords = clf.fit_transform(coords)


    # assign coordinates to x and y lists:
    x = [cell[0] for cell in coords]
    y = [cell[1] for cell in coords]
    z = [cell[2] for cell in coords]

    returned_xy = MDS_figgen(x, y, 'Component 1', 'Component 2', filename + "_12.pdf", workname, "Components 1 and 2. " + caption_postscript, col, stress, coords, comps, metric_p, normstress_p)
    returned_zy = MDS_figgen(y, z, 'Component 2', 'Component 3', filename + "_23.pdf", workname, "Components 2 and 3. " + caption_postscript, col, stress, coords, comps, metric_p, normstress_p)
    returned_xz = MDS_figgen(x, z, 'Component 1', 'Component 3', filename + "_13.pdf", workname, "Components 1 and 3. " + caption_postscript, col, stress, coords, comps, metric_p, normstress_p)
    return returned_xy + "\n\n" + returned_zy + "\n\n" + returned_xz

def MDS_3d_plot_just_return_coords(col, metric_p, normstress_p):
    """return a Multidimensional Scaling three-dimensional graph of the Hamming distance matrix of COL, where METRIC is either True or False"""

    # set up MDS for taking a precompiled distance matrix:
    comps = 3
    from sklearn.manifold import MDS
    # mds = MDS(dissimilarity='precomputed', random_state=0, metric=metric_p, n_components=comps)
    mds = MDS(dissimilarity='precomputed', random_state=0, metric=metric_p, n_components=comps, normalized_stress=normstress_p)

    dm = full_dist_matrix(col) # get a full distance matrix of COL
    X = np.array(sk_dmatrix(dm))

    # generate the coords:
    coords = mds.fit_transform(X)
    # nb: if you increase n_components above, this will increase the number of dimensions here!
    
    stress = round(mds.stress_, 3)

    # rotate the data:
    # this part taken from: https://scikit-learn.org/stable/auto_examples/manifold/plot_mds.html#example-manifold-plot-mds-py
    from sklearn.decomposition import PCA
    clf = PCA(n_components=comps)
    coords = clf.fit_transform(coords)


    # assign coordinates to x and y lists:
    x = [cell[0] for cell in coords]
    y = [cell[1] for cell in coords]
    z = [cell[2] for cell in coords]
    return [x, y, z, stress]

def MDS_interactive_3d_plot(collation, metric_p, scaled_stress_p):
    # adapted from https://matplotlib.org/stable/gallery/mplot3d/surface3d.html
    coords = MDS_3d_plot_just_return_coords(collation, metric_p, scaled_stress_p)
    # rotate the data:
    # this part taken from: https://scikit-learn.org/stable/auto_examples/manifold/plot_mds.html#example-manifold-plot-mds-py
    from sklearn.decomposition import PCA
    clf = PCA(n_components=3)
    zippedcoords=list(zip(coords[0], coords[1], coords[2]))
    zippedcoords = clf.fit_transform(zippedcoords)
    # print(zippedcoords)
    
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import axes3d
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # X = coords[0]
    # Y = coords[1]
    # Z = coords[2]
    X = [coord[0] for coord in zippedcoords]
    Y = [coord[1] for coord in zippedcoords]
    Z = [coord[2] for coord in zippedcoords]
    namelist = collation[0][1:]
    for item in list(zip(X, Y, Z, namelist)):
        ax.plot(item[0], item[1], item[2], marker='o', label='foo')
        ax.text(item[0], item[1], item[2], item[3], color='black')
    plt.show(block=False)
    return 

def two_nearest(names, dists):
    """given NAMES, an ordered list of names, and DISTS, a similarly ordered list of distances, return the names and distances corresponding to the smallest and next-smallest distances"""
    nd = list(zip(names, dists))
    new_nd = [pair for pair in nd if pair[1] != 0]
    new_nd.sort(key=lambda x: x[1]) # sort from smallest to largest distance
    smallest = new_nd[0] # nd[0] will be the witness itself, with a distance of 0.
    next_smallest = new_nd[1]
    return [smallest, next_smallest]


def MDS_interactive_3d_plot_lined(collation, metric_p, scaled_stress_p, workname):
    # parts adapted from https://matplotlib.org/stable/gallery/mplot3d/surface3d.html
    # another part adapted from: https://scikit-learn.org/stable/auto_examples/manifold/plot_mds.html#example-manifold-plot-mds-py

    coords = MDS_3d_plot_just_return_coords(collation, metric_p, scaled_stress_p)
    # rotate the data:
    from sklearn.decomposition import PCA
    clf = PCA(n_components=3)
    zippedcoords=list(zip(coords[0], coords[1], coords[2]))
    zippedcoords = clf.fit_transform(zippedcoords)
    # print(zippedcoords)
    
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import axes3d
    matplotlib.rcParams['font.family'] = 'Charis SIL'

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # X = coords[0]
    # Y = coords[1]
    # Z = coords[2]
    X = [coord[0] for coord in zippedcoords]
    Y = [coord[1] for coord in zippedcoords]
    Z = [coord[2] for coord in zippedcoords]
    namelist = collation[0][1:]
    distmatrix = full_dist_matrix(collation)[1:]
    for item in list(zip(X, Y, Z, namelist, distmatrix)):

        itemname=item[3]

        twonearest=two_nearest(namelist, item[4][1:])
        # print(twonearest)
        
        nearestname=twonearest[0][0]
        nextnearestname=twonearest[1][0]
        # print(itemname, ": ", nearestname, nextnearestname)

        # stress = round(mds.stress_, 3)
        ax.set_title(workname + ": stress = " + str(coords[3]))
        # ax.set_title(workname)
        for nitem in list(zip(X, Y, Z, namelist)):
            if nitem[3] == nearestname:
                nearestx = nitem[0]
                nearesty = nitem[1]
                nearestz = nitem[2]
            if nitem[3] == nextnearestname:
                nextnearestx = nitem[0]
                nextnearesty = nitem[1]
                nextnearestz = nitem[2]
        # print(nearestx, nearesty, nearestz, nextnearestx, nextnearesty, nextnearestz)
        
        ax.plot([item[0],nearestx], [item[1], nearesty], [item[2], nearestz], color='red', linestyle='solid', alpha=0.5)
        ax.plot([item[0],nextnearestx], [item[1], nextnearesty], [item[2], nextnearestz], color='orange', linestyle='solid', alpha=0.3)

        ax.plot(item[0], item[1], item[2], marker='o')
        ax.text(item[0], item[1], item[2], item[3], color='black')
    plt.show(block=False)
    return 

def remove_witnesses(col, witlist):
    """Return collation COL shorn of all witnesses in WITLIST"""
    namelist = col[0][1:]
    for wit in witlist:
        namelist.remove(wit)
    subcol = sub_col_rownums(col, namelist)
    return subcol

# ------------------------------------------------------------------------------
# INTERACTIVE 3D NEIGHBOUR-JOINER
# ------------------------------------------------------------------------------

def nj_post_dist_matrix(col):
    """return a distance matrix including all nodes created by neighbour-joiner, with the distances in the initital topology added to it
    """
    # get the initial distance matrix:
    initial_dist_matrix = full_dist_matrix(col)
    # get the neighbour_joiner topology:
    dlist = dist_matrix_to_list(dist_matrix(col))
    topology = neighbour_joiner(dlist, [], 1)
    # figure out how many nodes are needed:
    num_nodes = len(topology) - 1 # because the last item on the topology will also refer to the last node.
    node_names = [row[0] for row in topology][:-1] # again, because the last name should be identical to the penultimate one.
    # create the new distance matrix, with empty values for the nodes:
    new_dm_headline = initial_dist_matrix[0]
    new_dm = [new_dm_headline + node_names]
    witness_filler_list = ['fill_me' for item in node_names]
    for row in initial_dist_matrix[1:]:
        new_dm.append(row + witness_filler_list)
    node_filler_list = ['fill_me' for item in initial_dist_matrix[0][1:]]
    for node in node_names:
        new_dm.append([node] + node_filler_list + witness_filler_list)
    for row in topology:
        newnode = row[0]
        child1 = row[1][0]
        child1_distance = row[1][1]
        child2 = row[2][0]
        child2_distance = row[2][1]
        # insert distances between the children and the node into the distance matrix:
        newnode_index = new_dm[0].index(newnode)
        child1_index = new_dm[0].index(child1)
        child2_index = new_dm[0].index(child2)
        # insert distances between new_node and its children:
        new_dm[newnode_index][child1_index] = child1_distance
        new_dm[child1_index][newnode_index] = child1_distance
        new_dm[newnode_index][child2_index] = child2_distance
        new_dm[child2_index][newnode_index] = child2_distance
    return new_dm



def find_distance(me, otherwit, dmatrix, topology):
    """for ME, find the distance to OTHERWIT using the DMATRIX (output of NJ_POST_DIST_MATRIX()) and TOPOLOGY (output of NEIGHBOUR_JOINER())"""
    import statistics
    me_index = dmatrix[0].index(me)
    otherwit_index = dmatrix[0].index(otherwit)
    if dmatrix[me_index][otherwit_index] != 'fill_me':
        distance = dmatrix[me_index][otherwit_index]
        # print("  distance for " + otherwit + " is " + str(distance))
        return distance
    else:
        wits_to_mean = []
        for row in topology:
            if row[0] == otherwit:
                wits_to_mean.append(row[1][0])
                wits_to_mean.append(row[2][0])
        wits_to_mean = set(wits_to_mean)
        # print("no distance for " + otherwit + ": need to calculate using " + str(wits_to_mean))
        distances = []
        for wit in wits_to_mean:
            distances.append(find_distance(me, wit, dmatrix, topology))
            mean_distance = statistics.mean(distances)
            # print("returning means for " + str(wits_to_mean) + ": " + str(distances) + " = " + str(mean_distance))
        return mean_distance

# now, for each witness in dmatrix, I want to find the distance to each node, if not already known.
def nj_complete_dist_matrix(col):
    """given COL, a collation, return a complete distance matrix including neighbour-joiner nodes"""
    dmatrix = nj_post_dist_matrix(col)
    dlist = dist_matrix_to_list(dist_matrix(col))
    topology = neighbour_joiner(dlist, [], 1)
    witlist = col[0][1:]
    nodelist = dmatrix[0][len(witlist) + 1:]

    # I need to add 0's for where it's a distance to itself.
    for row in dmatrix:
        name1 = row[0]
        for item, name2 in zip(row, dmatrix[0]):
            if name1 == name2:
                dmatrix[dmatrix[0].index(name1)][dmatrix[0].index(name2)] = 0

    for wit in witlist:
        wit_index = dmatrix[0].index(wit)
        for node in nodelist:
            node_index = dmatrix[0].index(node)
            dist = find_distance(wit, node, dmatrix, topology)
            dmatrix[wit_index][node_index] = dist
            dmatrix[node_index][wit_index] = dist

    for wit in nodelist:
        wit_index = dmatrix[0].index(wit)
        for node in nodelist:
            node_index = dmatrix[0].index(node)
            dist = find_distance(wit, node, dmatrix, topology)
            dmatrix[wit_index][node_index] = dist
            dmatrix[node_index][wit_index] = dist

    return dmatrix


def nj_3d(col, metric_p, scaled_stress_p, workname):
    from sklearn.manifold import MDS
    import numpy as np
    mds = MDS(dissimilarity='precomputed', random_state=0, metric=metric_p, n_components=3, normalized_stress=scaled_stress_p)

    dm = nj_complete_dist_matrix(col)
    final_matrix = np.array(sk_dmatrix(dm))
    coords = mds.fit_transform(final_matrix)
    coords

    stress = round(mds.stress_, 3)

    # rotate the data:
    # this part taken from: https://scikit-learn.org/stable/auto_examples/manifold/plot_mds.html#example-manifold-plot-mds-py
    from sklearn.decomposition import PCA
    clf = PCA(n_components=3)
    coords = clf.fit_transform(coords)


    # assign coordinates to x and y lists:
    x = [cell[0] for cell in coords]
    y = [cell[1] for cell in coords]
    z = [cell[2] for cell in coords]
    # [x, y, z, stress]

    namelist = dm[0][1:]

    import matplotlib
    matplotlib.rcParams['font.family'] = 'Charis SIL'
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import axes3d
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title(workname + ": stress = " + str(stress))

    dlist = dist_matrix_to_list(dist_matrix(col))
    topology = neighbour_joiner(dlist, [], 1)
    for row in topology:
        parent = row[0]
        child1 = row[1][0]
        child2 = row[2][0]
        parent_coords = coords[namelist.index(parent)]
        child1_coords = coords[namelist.index(child1)]
        child2_coords = coords[namelist.index(child2)]
        ax.plot([parent_coords[0], child1_coords[0]], [parent_coords[1], child1_coords[1]], [parent_coords[2], child1_coords[2]], color='black', linestyle='solid', alpha=0.5)
        ax.plot([parent_coords[0], child2_coords[0]], [parent_coords[1], child2_coords[1]], [parent_coords[2], child2_coords[2]], color='black', linestyle='solid', alpha=0.5)

    witnesses = col[0][1:]
    firstnode = len(witnesses) + 1
    nodes = dm[0][firstnode:]
    # print(witnesses, nodes)
    
    for item in list(zip(x, y, z, namelist)):
        if item[3] in witnesses:
            ax.plot(item[0], item[1], item[2], marker='o')
        else:
            ax.plot(item[0], item[1], item[2], color='gray')
        ax.text(item[0], item[1], item[2], item[3], color='black')
    
    return plt.show(block=False)

def median_witness(col):
    """return a sorted list from median witness to most outlying witness, with distances"""
    all_dists = nj_complete_dist_matrix(col)
    sum_dists = [[row[0], sum(row[1:])] for row in all_dists[1:]]
    return sorted(sum_dists, key=lambda x: x[1])
