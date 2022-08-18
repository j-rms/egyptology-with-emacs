# ------------------------------------------------------------------------------
# BAKENSESHAT: tools for Egyptological textual criticism
# ------------------------------------------------------------------------------

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
                if entry == 0.0:
                    deduped_line.append(0.0001) # you can't have 0.0 without getting a divide by zero error when calculating neighbour joining; have to use an extremely small value instead.
                else:
                    deduped_line.append(entry)
        deduped_dist_matrix.append(deduped_line)
    return deduped_dist_matrix

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
        print("SIBLINGS LIST:")
        print(siblings_list) # for testing
        print("NEW DIST LIST:")
        print(new_dist_list)
        return siblings_list
    else:
        print("SIBLINGS LIST:")
        print(siblings_list)
        # print(sib1, sib2)
        print("NEW DIST LIST:")
        print(new_dist_list)
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
            if '[' in item:
                no_lacuna_p = False
        return no_lacuna_p
    def no_insignificant_omission_p(line):
        """returns True if no element of LINE == '‑' (non-breaking hyphen, not normal - sign!!!), i.e. is an insignificant omission"""
        no_insignificant_omission_p = True
        for item in line:
            if "‑" == item:
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

from itertools import combinations
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
    changelist = []
    for row in coll1:
        cellnum = 0
        for cell in row:
            coll1_reading = coll1[rownum][cellnum]
            coll2_reading = coll2[rownum][cellnum]
            if coll1_reading != coll2_reading:
                witname = witnames[cellnum]
                diff = str(rownum) + ": " + witname + ": " + coll1_reading + " → " + coll2_reading
                changelist.append(diff)
            cellnum = cellnum + 1
        rownum = rownum + 1
    return changelist
