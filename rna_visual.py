# Pavithra
# April 2019
# rna_visual.py
##############
# Use: Modify dot bracket notation according to bootstrap confidence, unfavorable short base-pair runs, or reactivity differences.
###############


"""
Important: dot_bracket_modify() with corresponding bpp and reactivity data matrices must be run in MATLAB first!
           Example (in MATLAB, not in terminal): dot_bracket_modify(bpp_1M7, d_1M7)
Once the bpp.txt and reaccs_only.txt file are generated, rna_visual.py can be run multiple times.
rna_visual.py depends on the .txt files formatted from the MATLAB script. Without them, will not work.



IMPORTANT: Make sure to delete Dot_Bracket_Sequences.txt before running, if you want a fresh run/no previous dot bracket sequences in that file. This script APPENDS to that file.


Sample run on terminal:
python rna_visual.py '(((....((((((..)...........((((((...((((((.....(((.(((((((((.(((((((((....))))))))).(((((....))))).)))).....))))))))....))))))..)).)))))))..((((...((((((((.....))))))))..))))..................((((((.....))))))......................' GGCCAAAGGCGUCGAGUAGACGCCAACAACGGAAUUGCGGGAAAGGGGUCAACAGCCGUUCAGUACCAAGUCUCAGGGGAAACUUUGAGAUGGCCUUGCAAAGGGUAUGGUAAUAAGCUGACGGACAUGGUCCUAACCACGCAGCCAAGUCCUAAGUCAACAGAUCUUCUGUUGAUAUGGAUGCAGUUCAAAACCAAACCGUCAGCGAGUAGCUGACAAAAAGAAACAACAACAACAAC --check_hairpins=True --check_confidence=True --check_reactivity=True --ones=False --twos=True --confidence_cutoff=0.7 --reactivity_cutoff=1.0 --bpp_textfile=bpp.txt --react_textfile=reaccs_only.txt

Arguments needed in command line:
python rna_visual.py '[Dot bracket sequence]' [RNA sequence] --check_hairpins=[True or False] --check_reactivity=[True or False] --ones=False --twos=True --confidence_cutoff[value between 0 and 1] --reactivity_cutoff=[some value] --bpp_textfile=[name of bpp txt file] --react_textfile=[name of reactivity txt file]


What each argument means...
--check_hairpins=
    If value is True, the given dot bracket sequence will be modified according to whether a single base-pair is found i.e. ....(.).... is found or a double base-pair i.e. ....((.)).....
    ASSUMPTION: either --ones= or --twos= must be set to True, or else no modifications/checks will be made.
--check_confidence=
    If value is True, the given dot bracket sequence will be modified according to the provided cutoff value
    ASSUMPTION: --confidence_cutoff=must be set to some threshold value in [0,1.0] Example value could be 0.50. S
    #Example: Suppose confidence value at a given base-pair at locations 12 and 101 have a confidence value of 0.40. In this case, dot bracket sequence will be modified to .. at those respective positions since is below the cutoff.
--check_reactivity=
    If value is True, the given dot bracket sequence will be modified according to the provided cutoff value for the maximum difference between two base-paired nucleotides.
    #Example: A and C in sequence at positions 145 and 25 are base paired. A has a reactivity value of 2.45 but C has one of 0.5. This difference in values exceeds the threshold of 1.0, so the base-pair will be modified to .. at the respective positions in dot bracket sequence
    ASSUMPTION: --reactivity_cutoff=must be set to some threshold value of choice. Example value could be 1.0.
--ones=
    If value is True, dot bracket sequence will be modified according to whether a single base-pair is found i.e. ....(.)....
--twos=
    If value is True, dot bracket sequence will be modified according to whether a single base-pair is found i.e. ....((.))....
--confidence_cutoff=
    must be set to some threshold value in [0,1.0] Example value could be 0.50.
    #Example: Suppose confidence value at a given base-pair at locations 12 and 101 have a confidence value of 0.40. In this case, dot bracket sequence will be modified to .. at those respective positions since is below the cutoff.
--reactivity_cutoff=
    --reactivity_cutoff=must be set to some threshold value of choice. Example value could be 1.0.
    #Example: A and C in sequence at positions 145 and 25 are base paired. A has a reactivity value of 2.45 but C has one of 0.5. This difference in values exceeds the threshold of 1.0, so the base-pair will be modified to .. at the respective positions in dot bracket sequence
--bpp_textfile=
    must be set to the name of the bpp.txt file created using dot_bracket_modify.m
    #ASSUMPTION: must use dot_bracket_modify.m to create text file so that parsing data from this .txt file is possible
--react_textfile=
    must be set to the name of reactivity.txt file created using dot_bracket_modify.m
    #ASSUMPTION: must use dot_bracket_modify.m to create text file so that parsing data from this .txt file is possible

"""

import itertools


class CommandLine():
    '''
    Handles user input arguments.
    '''

    def __init__(self, inOpts=None):

        import argparse
        self.parser = argparse.ArgumentParser(description='user specifications for rna_construct.py',
                                              epilog='This tool will help you generate a (hopefully optimal) RNA construct with buffer and hairpin regions.',
                                              add_help=True,  # default is True
                                              prefix_chars='-',
                                              usage='%(prog)s [options] -option1[default] <input >output'
                                              )
        self.parser.add_argument('dot_bracket_sequence', action='store', help='dot bracket sequence')
        self.parser.add_argument('rna_sequence', action='store', help='rna')
        self.parser.add_argument('-check_hairpins', '--check_hairpins', action='store', help='Do you want to modify dot bracket sequence if hairpins or loops with stems of length 1 or 2 are found?')
        self.parser.add_argument('-check_confidence', '--check_confidence', action='store', help='Do you want to modify dot bracket sequence if the bootstrapped probability value for base pairing is below a certain value?')
        self.parser.add_argument('-check_reactivity', '--check_reactivity', action='store', help='Do you want to modify dot bracket sequence if reactivity values between base-paired nucleotides differ by a certain level of reactivity?')
        self.parser.add_argument('-ones', '--ones', action='store', help='Example: Changes ...((....))... to ..............')
        self.parser.add_argument('-twos', '--twos', action='store', help='Example: Changes ...(....)... to .............. ')
        self.parser.add_argument('-confidence_cutoff', '--confidence_cutoff', action='store', default=False, help='Example: 0.60 in probability')
        self.parser.add_argument('-reactivity_cutoff', '--reactivity_cutoff', action='store', default=False, help='Example: A reactivity difference of 0.9')
        self.parser.add_argument('-bpp_textfile', '--bpp_textfile', action='store', help='Example: 0.60 in probability')
        self.parser.add_argument('-react_textfile', '--react_textfile', action='store', help='Example: A reactivity difference of 0.9')
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)


def map_notation(dot_bracket_sequence):
    """Goal: Find positions of bases pairing with one another. To do this, use stack to map each ( to corresponding ).
       Input: Dot bracket sequence (ASSUMPTION: is a string)
       Output: Returns a dictionary with corresponding base pair positions. Example: {0:101, 25:212}
    """
    stack = []
    position_map = dict()
    for index in range(len(dot_bracket_sequence)):
        if dot_bracket_sequence[index] == '(':
            stack.append(index)
        else:
            if dot_bracket_sequence[index] == ')':
                position_map[stack.pop()] = index
    return position_map
    """
    for (a, b) in list(position_map.items()):
        print("{0}:{1}".format(rna[a], rna[b]))
    """


def bootstrap_check(file_name, length, threshold, bp_positions, dot_bracket_sequence):
    """Goal: Modify dot bracket notation so that base paired nucleotides all have confidence values above the threshold.
       Inputs: bpp.txt file, length of RNA, threshold value for confidence, dictionary of base pair positions (ASSUMPTION: from map_notation() function), dot bracket sequence(ASSUMPTION: is a string)
       Output: Returns the modified/unchanged dot bracket sequence, and prints to file named 'Dot_Bracket_Sequences.txt'
    """
    f = open('Dot_Bracket_Sequences.txt', 'a')
    row_vectors = length * [0]
    with open(file_name, 'r') as file:
        index = 0
        for key, group in itertools.groupby(file, lambda x: x.endswith('\n')):
            for entry in group:
                row_vectors[index] = entry
                index += 1

    i = -1
    changes = []

    for position in range(len(row_vectors)):
        tracker = -1
        for i in range(0, len(row_vectors[position]) - 4, 5):
            tracker += 1
            if float(row_vectors[position][i:i + 4]) < threshold and float(row_vectors[position][i:i + 4]) != 0.00 and (tracker in bp_positions or tracker in bp_positions.values()):
                if tracker in bp_positions and bp_positions[tracker] == position:
                    changes.append((float(row_vectors[position][i:i + 4]), tracker, position))
                else:
                    for key, value in bp_positions.items():
                        if key == position and bp_positions[position] == tracker:
                            changes.append((float(row_vectors[position][i:i + 4]), position, tracker))
                            break
    changes.sort(key=lambda x: x[1])
    dots_list = list(dot_bracket_sequence)
    for tupl in changes:
        dots_list[tupl[1]] = '.'
        dots_list[tupl[2]] = '.'
    if len(changes) > 0:
        f.write("\nModified Dot Bracket Sequence from considering bootstrap probabilites:\n {0}\n".format(''.join(dots_list)))
        f.write("Previous Dot Bracket Sequence:\n {0}\n".format(dot_bracket_sequence))
        print("\nModified Dot Bracket Sequence from considering bootstrap probabilites:\n {0}\n".format(''.join(dots_list)))
        print("Previous Dot Bracket Sequence:\n {0}\n".format(dot_bracket_sequence))
        return(''.join(dots_list))

    else:
        print("No modifications. Bootstrap probabilities at base-paired sites above specified threshold.")
        f.write("No modifications. Bootstrap probabilities at base-paired sites above specified threshold.")
        f.write("Dot Bracket Sequence:\n {0}\n".format(dot_bracket_sequence))
        print("Dot Bracket Sequence:\n {0}\n".format(dot_bracket_sequence))
        return(dot_bracket_sequence)

    """                
    print(changes)
    for (a, b, c) in changes:
        print("{2} {0}:{1}".format(a, rna[b], rna[c]))  # 0:4 5:9 10:14
    """


def reactivity_check(file_name, length, threshold, bp_positions, dot_bracket_sequence):
    """Goal: Modify dot bracket notation so that base paired nucleotides all have reactivity differences lower than the cutoff.
        #Example: if base at position 15 base-paired with base at position 75 has value of 1.5 and position 75 has value of 3.0, reactivity difference is 1.5.
        Inputs: reactivity.txt file, length of RNA, cutoff value for reactivity difference, dictionary of base pair positions (ASSUMPTION: from map_notation() function), dot bracket sequence(ASSUMPTION: is a string)
        Output: Returns the modified/unchanged dot bracket sequence, and prints to file named 'Dot_Bracket_Sequences.txt'
    """
    f = open('Dot_Bracket_Sequences.txt', 'a')
    check = False
    reactivities = length * [0]
    dots_list = list(dot_bracket_sequence)
    file = open(file_name, "r")
    index = 0
    for line in file.readlines():
        l = line.split(' ')
        for entry in l:
            if entry and entry[0] == '+':
                reactivities[index] = float(entry[1:len(entry)])
                index += 1
            if entry and entry[0] == '-':
                reactivities[index] = float(entry)
                index += 1

    for i in bp_positions:
        if abs(reactivities[bp_positions[i]] - reactivities[i]) > threshold:
            check = True
            dots_list[i] = '.'
            dots_list[bp_positions[i]] = '.'

    if check == True:
        f.write("\nModified Dot Bracket Sequence from considering reactivity values:\n {0}\n".format(''.join(dots_list)))
        f.write("Previous Dot Bracket Sequence:\n {0}\n".format(dot_bracket_sequence))
        print("\nModified Dot Bracket Sequence from considering reactivity values:\n {0}\n".format(''.join(dots_list)))
        print("Previous Dot Bracket Sequence:\n {0}\n".format(dot_bracket_sequence))
    else:

        print("No modifications. No asymmetric reactivity.")
        print("Dot Bracket Sequence:\n {0}\n".format(dot_bracket_sequence))
        f.write("No modifications. No asymmetric reactivity.")
        f.write("Dot Bracket Sequence:\n {0}\n".format(dot_bracket_sequence))

        return(dot_bracket_sequence)
    return(''.join(dots_list))


def hairpin_check(positions, position_map, twos_check, ones_check, dot_bracket_sequence):
    """Goal: Modify dot bracket notation so that short runs of base pairing are removed.
        #Example: modified according to whether a single base-pair is found i.e. ....(.).... is found or a double base-pair i.e. ....((.)).....
        Inputs: list of base pair positions (ASSUMPTION: from map_notation() and main()),dictionary of base pair positions (ASSUMPTION: from map_notation() function), ones_check from command line, twos_check from command line, dot bracket sequence(ASSUMPTION: is a string)
        Output: Returns the modified/unchanged dot bracket sequence, and prints to file named 'Dot_Bracket_Sequences.txt'
    """
    f = open('Dot_Bracket_Sequences.txt', 'a')
    dots_list = list(dot_bracket_sequence)
    i = -1
    check = False
    while (i + 1) < len(positions):
        count = 0
        i += 1
        current = i
        while((i + 1) < len(positions) and positions[i] + 1 == positions[i + 1]):
            count += 1
            i += 1
        if ones_check and count < 1:
            check = True
            dots_list[positions[current]:position_map[positions[current]] + 1] = '.' * (position_map[positions[current]] + 1 - positions[current])
        if twos_check and count < 2 and count != 0:
            check = True
            dots_list[positions[current]:position_map[positions[current]] + 1] = '.' * (position_map[positions[current]] + 1 - positions[current])
        count = 0
    if check:
        f.write("\nModified Dot Bracket Notation from considering short loop or hairpin base pairing sites:\n {0}\n".format(''.join(dots_list)))
        f.write("Previous Dot Bracket Notation:\n {0}\n".format(dot_bracket_sequence))
        print("\nModified Dot Bracket Notation from considering short loop or hairpin base pairing sites:\n {0}\n".format(''.join(dots_list)))
        print("Previous Dot Bracket Notation:\n {0}\n".format(dot_bracket_sequence))
        return(''.join(dots_list))
    else:
        f.write("\nNo modifications made to dot bracket notation. No such hairpin loops or bulges detected.")
        f.write("Dot Bracket Sequence:\n {0}\n".format(dot_bracket_sequence))
        print("\nNo modifications made to dot bracket notation. No such hairpin loops or bulges detected.")
        print("Dot Bracket Sequence:\n {0}\n".format(dot_bracket_sequence))
        return(dot_bracket_sequence)


def main():

    # Maps base pair positions, and stores in both dictionary and list of one of the bases in each pair.
    # Example: position_map = {0:15, 1:21}
    # Example: positions= [0,1]
    myCommandLine = CommandLine()
    position_map = map_notation(myCommandLine.args.dot_bracket_sequence)
    dots_list = list(myCommandLine.args.dot_bracket_sequence)
    positions = list(position_map.keys())
    positions.sort()

    # According to command line arguments, modifies dot bracket sequence according to user specifications.
    length = len(myCommandLine.args.rna_sequence)
    dot_bracket_sequence = myCommandLine.args.dot_bracket_sequence
    if myCommandLine.args.check_hairpins:
        dot_bracket_sequence = hairpin_check(positions, position_map, myCommandLine.args.twos, myCommandLine.args.ones, dot_bracket_sequence)
    if myCommandLine.args.check_confidence:
        dot_bracket_sequence = bootstrap_check('bpp.txt', length, float(myCommandLine.args.confidence_cutoff), position_map, dot_bracket_sequence)
    if myCommandLine.args.check_reactivity:
        dot_bracket_sequence = reactivity_check('reaccs_only.txt', length, float(myCommandLine.args.reactivity_cutoff), position_map, dot_bracket_sequence)

    """
    position_map = map_notation('....(.).....((..)).....')
    dots_list = list('....(.).....((..)).....')
    positions = list(position_map.keys())
    positions.sort()

    hairpin_check(positions, position_map, False, True, '....(.).....((..)).....')
    """


main()
