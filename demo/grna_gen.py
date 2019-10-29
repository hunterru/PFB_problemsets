#! /usr/bin/env python3


import argparse, re, datetime, time

def cas2dict(cas_prot_tsv,cas):

    # import a cas protein tsv to get cas specific information
    with open(cas_prot_tsv, 'r') as cas_file:
        for line in cas_file:
            # prep the lines to work with by removing irrelevant white space
            if line[0] != '#':
                if line[-1] == '\n':
                    line = line.rstrip()

                # populate the dictionary of cas protein information
                cas_param = line.split('\t')
                if cas_param[0] == cas:
                    cas_prot = cas_param[0]
                    pam = cas_param[1]
                    pam_orient = cas_param[2]
                    grna_len = cas_param[3]
                    proximal_nt = cas_param[4].split(',')
                    middle_nt = cas_param[5].split(',')
                    distal_nt = cas_param[6].split(',')
                    mm_toler = int(cas_param[7])
                    prox_toler = int(cas_param[8])
                    dint_toler = cas_param[9]
                    prox_weight = float(cas_param[10])
                    mid_weight = float(cas_param[11])
                    dis_weight = float(cas_param[12])

                    # split the boundaries in the cas protein tsv into ints
                    for i in range(2):
                        proximal_nt[i] = int(proximal_nt[i])
                        middle_nt[i] = int(middle_nt[i])
                        distal_nt[i] = int(distal_nt[i])

                    cas_list = [pam, pam_orient, grna_len, proximal_nt, middle_nt, \
                            distal_nt, mm_toler, prox_toler, dint_toler, prox_weight, \
                            mid_weight, dis_weight]

    return cas_list


def fasta2dict(fasta_file):
    fasta_dict = {}

    with open(fasta_file, 'r') as fasta:
        str_fasta = ''
        for line in fasta:
            # rstrips if it's not a sequence line and isn't a line with only \n
            if line[0] != '>' and line[0] != '\n':
                line = line.rstrip()
            # concatenates a new line character if it is the first sequence line
            if line[0] == '>' and str_fasta != '':
                str_fasta = str_fasta + '\n'
            # concatenates the prep string with the prepped line from the fasta file
            str_fasta = str_fasta + line

        # extracts 0) seq ID 1) description (or the \n if there is non) 2) sequence
        extracted = re.findall(r'(>\S+)([^\n]\s\w+.+\n|\n)(\S+[^>])', str_fasta)

        # adds a new dictionary for each gene
        for index in range(len(extracted)):
            gene = extracted[index][0]
            gene = gene[1:]

            # remove description if it is not actually present
            descrip = extracted[index][1]
            descrip = descrip[1:]
            seq = extracted[index][2]
            if seq[-1] is '\n' or '\r':
                seq = seq.rstrip()

            # prepares dictionaries for each gene with description, seq, rvcmpl_seq, \
            # and codons
            fasta_dict[gene] = {}
            if descrip == '\n':
                fasta_dict[gene]['description'] = None
            else:
                fasta_dict[gene]['description'] = descrip
            fasta_dict[gene]['sequence'] = seq.upper()

    return fasta_dict


# populates a list of dictionaries of all possible gRNAs for the sequences submitted
# remove any that don't fit our criteria
def grna_finder(fasta_dict,cas,gc_upper_bound,gc_lower_bound):

    # dict for translating ambiguous nucleotides
    ambiguous_nt = {
            'N': 'ATCGYRBDHV',
            'R': 'AG',
            'Y': 'CT',
            'B': 'CTGY',
            'D': 'ATGR',
            'H': 'ATCY',
            'V': 'ACGR'
            }

    # dict for translating sequences to their reverse complement, including ambiguous nucleotides
    reverse_nt = {
            'N': 'ATCGYRBDHV',
            'R': 'TC',
            'Y': 'GA',
            'B': 'GACVR',
            'D': 'TACY',
            'H': 'GATR',
            'V': 'CGTY',
            'G': 'C',
            'A': 'T',
            'C': 'G',
            'T': 'A'
            }

    # creating an empty list to populate with dictionaries of gRNAs
    grna_dict_list = []
    index = -1
    # for each cas protein in the input tsv, gather the info in each column
    pam = cas[0]
    pam_orientation = cas[1]
    grna_length = int(cas[2])

    # populating PAM strings for concatenation
    pam_re = ''
    rev_pam_re = ''

    # for each character in the inputted PAMs
    for char in pam:
        # if it is a standard nucleotide, add it to 'pam_re'
        if char in ['A','T','C','G']:
            pam_re += char
            rev_pam_re = reverse_nt[char] + rev_pam_re
        # if it isn't a standard nucleotide, refer to the ambiguous dictionary to obtain the \
        # nucleotides the nonstandard nucleotide is composed
        else:
            pam_re += '[' + ambiguous_nt[char] + ']'
            rev_pam_re = '[' + reverse_nt[char] + ']' + rev_pam_re

    # create a dictionary for counting duplicate gRNAs
    grna_counter = {}

    # for each sequence in the fasta dictionary
    for seq_name in fasta_dict:
        sequence = fasta_dict[seq_name]['sequence']
        # if the PAM sequence is on the 3' end of the gRNA
        if pam_orientation == "3'":
            # match a regular expression containing word characters equivalent to the gRNA len \
            # and obtain the PAM nucleotides from the regular expression ready converted sequence
            grna_matches = re.finditer(r'(?=(\w{' + str(grna_length) + '}' + pam_re + '))', sequence)
            # match a regular expression containing the prepared reverse PAM regular expression \
            # and obtain word characters afterward up to the gRNA len
            rev_grna_matches = re.finditer(r'(?=('+ rev_pam_re +'\w{' + str(grna_length) + '}))', sequence)
            # populate the match information from the iterable regular expression into a tuple \
            # index 0 is the start index of the gRNA PAM combo, 1 is the end index, 2 is the gRNA \
            # 3 is the PAM sequence - remember this is with respect to the orientation of the PAM \
            # relative to the gRNA
            match_info = [(match.start(1), match.end(1), match.group(1)[:grna_length], match.group(1)[grna_length:]) for match in grna_matches]
            # same for the reverse sequence
            rev_match_info = [(rev_match.start(1), rev_match.end(1), rev_match.group(1)[len(pam):], rev_match.group(1)[:len(pam)]) for rev_match in rev_grna_matches]
        # if the PAM sequence is on the 5' end of the gRNA
        elif pam_orientation == "5'":
            # this is essentially the reciprocal of the above 'if' command
            grna_matches = re.finditer(r'(?=(' + pam_re + '\w{' + str(grna_length) + '}))', sequence)
            rev_grna_matches = re.finditer(r'(?=(\w{' + str(grna_length) + '}'+rev_pam_re+'))', sequence)
            match_info = [(match.start(1), match.end(1), match.group(1)[len(pam):], match.group(1)[:len(pam)]) for match in grna_matches]
            rev_match_info = [(rev_match.start(1), rev_match.end(1), rev_match.group(1)[:grna_length], rev_match.group(1)[grna_length:]) for rev_match in rev_grna_matches]

        print('\n'+str(len(match_info)+len(rev_match_info)),'potential gRNAs.')
        # create a new index to populate index values for deletion from the compiled dict
        # we will be deleting any gRNA that is a complete match cuz who wants a 100% off \
        # target ( if it's you, then whatever dude, go figure your life out )
        gc_removal_count = 0
        at_removal_count = 0
        for grna in match_info:
            g_count = grna[2].count('G')
            c_count = grna[2].count('C')
            gc_count = g_count + c_count
            gc_content = gc_count/len(grna[2])
            if gc_content > gc_upper_bound:
                gc_removal_count += 1
                continue
            if gc_content < gc_lower_bound:
                at_removal_count += 1
                continue
            # add 1 to the index for each iteration
            index += 1
            # if the gRNA is a key in the counter dictionary
            # append the index of the current seq to the gRNA counter
            if grna_counter.get(grna[2]):
                grna_counter[grna[2]].append([index,'+'])
            else:
            # otherwise, add a new key for the gRNA and populate it with the current index
                grna_counter[grna[2]] = [[index,'+']]
            # append a dictionary to the list of gRNAs and populate it with .bed specific \
            # parameters
            grna_dict_list.append({})
            grna_dict = grna_dict_list[-1]
            grna_dict['chrom'] = seq_name
            grna_dict['chromStart'] = grna[0]
            grna_dict['chromEnd'] = grna[1]
            grna_dict['grna_seq'] = grna[2]
            grna_dict['strand'] = '+'
            grna_dict['pam'] = grna[3]
           # grna_dict['name'] = cas

        # now continue with the reverse sequences
        for grna in rev_match_info:
            # removing values based on GC bounds
            g_count = grna[2].count('G')
            c_count = grna[2].count('C')
            gc_count = g_count + c_count
            gc_content = gc_count/len(grna[2])
            if gc_content > gc_upper_bound:
                gc_removal_count += 1
                continue
            if gc_content < gc_lower_bound:
                at_removal_count += 1
                continue

            index += 1
            grna_dict_list.append({})
            grna_dict = grna_dict_list[-1]
            grna_dict['chrom'] = seq_name
            grna_dict['chromStart'] = grna[0]
            grna_dict['chromEnd'] = grna[1]

            # time to create the reverse complements of each sequence
            # these will still be reported as the start and stop indices of the sense \
            # sequence in the .bed file, though these are the antisense sequences
            before_com_seq = list(grna[2])
            for nt_index in range(len(before_com_seq)):
                before_com_seq[nt_index] = reverse_nt[before_com_seq[nt_index]]
            com_grna_seq = ''.join(before_com_seq)
            rev_grna_seq = com_grna_seq[::-1]

            # now we've created the reverse complements, we can check the dictionary \
            # for same sequence gRNAs and add their indices if there is overlap
            if grna_counter.get(rev_grna_seq):
                grna_counter[rev_grna_seq].append([index,'-'])
            else:
                grna_counter[rev_grna_seq] = [[index,'-']]

            grna_dict['grna_seq'] = rev_grna_seq
            grna_dict['strand'] = '-'
            before_com_pam = list(grna[3])
            # get the antisense for the PAM as well

            for nt_index in range(len(before_com_pam)):
                before_com_pam[nt_index] = reverse_nt[before_com_pam[nt_index]]
            com_pam_seq = ''.join(before_com_pam)
            rev_pam_seq = com_pam_seq[::-1]
            grna_dict['pam'] = rev_pam_seq
            #grna_dict['name'] = ca

        # delete 100% gRNA matches by first populating a list of deletion indices
        print('\nPopulating list of gRNAs that 100% overlap... \ncuz who needs them, amiright?\n')
        deletion_index = []
        for grna_seq in grna_counter:
            seq_count = grna_counter[grna_seq]
            if len(seq_count) > 1:
                for del_index in seq_count:
                    deletion_index.append(del_index[0])
        # so as to not break the index-list relationship, sort the indices from \
        # highest to lowest
        deletion_index.sort(reverse=True)

        print('Removed',gc_removal_count,'gRNAs that exceed GC upper boundary:',gc_upper_bound)
        print('Removed',at_removal_count,'gRNAs that exceed GC lower boundary:',gc_lower_bound)

            # delete the largest index first and work through the rest
    for index in deletion_index:
        del grna_dict_list[index]

    print('Removed',len(deletion_index),'gRNAs that have 100% similarity to another.')
    print('\n'+str(len(grna_dict_list)),'valid gRNAs to score against.')

    return grna_dict_list


# this function scores gRNAs by comparing each to all others for mismatches within the constraints
# in cas protein .tsv
def mismatch_scoring(grna_dict_list,cas,min_score,min_index,max_index):
    print('Scoring gRNAs...',str(len(grna_dict_list)**2),'iterations - ignore if index boundary specified.')

    # simplifying boundary values from cas protein information list
    prox_start = cas[3][0]
    prox_end = cas[3][1]
    mid_start = cas[4][0]
    mid_end = cas[4][1]
    dis_start = cas[5][0]
    dis_end = cas[5][1]

    # for each guide RNA, give it a starting score of 100
    for index1 in range(len(grna_dict_list)):
        grna1 = grna_dict_list[index1]
        if grna1['chromStart'] < min_index:
            continue
        if grna1['chromEnd'] > max_index and str(max_index) != '-1':
            continue
        score = 100

        # if the orientation of the PAM sequence is 5', sort the proximal, middle \
        # and distal sequences accordingly, same with 3'
        if cas[1] == "5'":
            prox1 = grna1['grna_seq'][prox_start:prox_end]
            mid1 = grna1['grna_seq'][mid_start:mid_end]
            dis1 = grna1['grna_seq'][dis_start:dis_end]
        elif cas[1] == "3'":
            temp_rev = grna1['grna_seq'][::-1]
            prox1_rev = grna1['grna_seq'][prox_start:prox_end]
            mid1_rev = grna1['grna_seq'][mid_start:mid_end]
            dis1_rev = grna1['grna_seq'][dis_start:dis_end]
            prox1 = prox1_rev[::-1]
            mid1 = mid1_rev[::-1]
            dis1 = dis1_rev[::-1]

        # for all other gRNAs, if the score is 0 at the start, just stop it
        for index2 in range(len(grna_dict_list)):
            if score <= 0:
                break
            grna2 = grna_dict_list[index2]

            # this 'if' checks to make sure its not comparing a sequence to itself
            if grna1 != grna2:
                # getting the proximal, middle, and distal sequences for the second gRNA
                if cas[1] == "5'":
                    prox2 = grna2['grna_seq'][prox_start:prox_end]
                    mid2 = grna2['grna_seq'][mid_start:mid_end]
                    dis2 = grna2['grna_seq'][dis_start:dis_end]
                elif cas[1] == "3'":
                    temp_rev = grna2['grna_seq'][::-1]
                    prox2_rev = grna2['grna_seq'][prox_start:prox_end]
                    mid2_rev = grna2['grna_seq'][mid_start:mid_end]
                    dis2_rev = grna2['grna_seq'][dis_start:dis_end]
                    prox2 = prox2_rev[::-1]
                    mid2 = mid2_rev[::-1]
                    dis2 = dis2_rev[::-1]

                # consecutive nucleotide mismatches can completely nullify an off-target \
                # so let's setup a counter to check for that
                consec_nt_count = 0

                # also, let's make a counter for proximal mismatches
                prox_mm_count = 0
                # check each nucleotide in each sequence and compare to see if they mismatch
                # if they do add one to the proximal mismatch counter, add 1 to the \
                # consecutive nucleotide mismatch counter
                for index in range(len(prox1)):
                    if prox1[index] != prox2[index]:
                        prox_mm_count += 1
                        consec_nt_count += 1
                        # if that second nucelotide mismatch counter is > 1 and the cas \
                        # protein specifications say that is not an offtarget, then break
                        if consec_nt_count > 1 and cas[8] == 'n':
                            break
                    # if a mismatch isn't detected, reset the consecutive nt counter
                    else:
                        consec_nt_count = 0

                # now that loop is finished/if we broke out of it due to consecutive nt \
                # move on to the next sequence to check for mismatches
                if consec_nt_count > 1 and cas[8] == 'n':
                    continue

                # if the proximal mismatch counter is above the threshold for this portion \
                # continue on to the next sequence as well because this isn't a mismatch
                if prox_mm_count > int(cas[7]):
                    continue

                # create a mid sequence mismatch counter and proceed the same as the proximal
                mid_mm_count = 0
                for index in range(len(mid1)):
                    if mid1[index] != mid2[index]:
                        mid_mm_count += 1
                        consec_nt_count += 1
                        if consec_nt_count > 1 and cas[8] == 'n':
                            break
                    else:
                        consec_nt_count = 0
                if consec_nt_count > 1 and cas[8] == 'n':
                    continue

                # now the same for the distal
                dis_mm_count = 0
                for index in range(len(dis1)):
                    if dis1[index] != dis2[index]:
                        dis_mm_count += 1
                        consec_nt_count += 1
                        if consec_nt_count > 1 and cas[8] == 'n':
                            break
                    else:
                        consec_nt_count = 0
                if consec_nt_count > 1 and cas[8] == 'n':
                    continue

                # if the total mismatches is above the threshold in the cas protein tsv \
                # this sequence isn't a valid off target for scoring, so continue
                if prox_mm_count + mid_mm_count + dis_mm_count > int(cas[6]):
                    continue

                # the score to subtract is (100/(zd+ym+xp)) where 'p' = proximal mismatches
                # 'm' = mid mismatches, 'd' = distal mismatches, and their coefficients are
                # their respective mismatch weights derived from the cas protein tsv
                subtract_score = 100/(cas[9]*prox_mm_count+cas[10]*mid_mm_count+cas[11]*dis_mm_count)
                score -= subtract_score

        # if the score is less than 0, set the value for that gRNA to 0, otherwise set
        # it to the score
        if score < 0:
            grna_dict_list[index1]['mismatch_score'] = 0
        else:
            grna_dict_list[index1]['mismatch_score'] = score

    # if the minimal score is greater than 0, check if each value is less than it
    # if it is, chuck that gRNA out
    if min_score > 0:
        print('\nPopulating list of gRNAs below minimum score ('+str(min_score)+')')
        del_list = []
        for index in grna_dict_list:
            if grna_dict_list[index]['mismatch_score'] < min_score:
                del_list.append(index)
        print('Removing gRNAs below minimum score.')
        for del_index in del_list:
            del grna_dict_list[del_index]

    return grna_dict_list


def grna_dict2bed(fasta,hits):

    # prepare the header for the output bed string
    output_string = '#chrom\tchromStart\tchromEnd\tguideSeq\tstrand\tpam\tmismatch_score'
    # name the output file by concatenating the following to the end of the fasta
    with open(fasta+'_grna_hits.bed', 'w') as output_file:
        for grna in hits:
            if grna.get('mismatch_score'):
                output_string += '\n'
                output_string += grna['chrom']
                output_string += '\t'+str(grna['chromStart'])
                output_string += '\t'+str(grna['chromEnd'])
                output_string += '\t'+grna['grna_seq']
                output_string += '\t'+grna['strand']
                output_string += '\t'+grna['pam']
               # output_string += '\t'+grna['name'] 
                output_string += '\t'+str(grna['mismatch_score'])
        output_file.write(output_string)


def exit(logo):

    print("\n\n\nYou're data is",end='')
    time.sleep(1)
    print('.',end='')
    time.sleep(1)
    print('.',end='')
    time.sleep(1)
    print('.\n\n\n')

    with open(logo, 'r') as cat:
        for line in cat:
            print(line,end='')


def denovoGuideRnaAnno(cas_file,cas_prot,fasta,gc_lower,gc_upper,min_index,max_index,score):

    start = time.time()
    start_time = datetime.datetime.now()
    print('Execution begun:',str(start_time))

    cas_list = cas2dict(cas_file,cas_prot)
    fasta_dict = fasta2dict(fasta)
    hits = grna_finder(fasta_dict,cas_list,float(gc_upper),float(gc_lower))
    scored = mismatch_scoring(hits,cas_list,score,int(min_index),int(max_index))
    grna_dict2bed(fasta,scored)

    end_time = datetime.datetime.now()
    end = time.time()

    print('\n',str(len(scored)),'gRNA analyzed.')
    print('Execution finished:',str(end_time))
    print('Execution time:',str(end-start))

    while True:
        try:
            with open('cat.txt','r') as cat_file:
                print()
                exit(cat_file)
            break
        except FileNotFoundError:
            break

def main(args):

	parser = argparse.ArgumentParser(description='Imports fasta, extracts possible gRNAs')
	parser.add_argument('-f','--fasta',required=True,help='Fasta file to extract gRNAs from')
	parser.add_argument('-s','--score',default=0,help='Lower boundary for off target scores, default = 0')
	parser.add_argument('-c','--cas',required=True,help='Cas protein to extract from tab delimitted list of Cas protein parameters - must exactly match name')
	parser.add_argument('--cas_list',required=True,help='Tab delimitted list of Cas protein parameters')
	parser.add_argument('--gc_upper',default=.8,help='GC content upper bound (in decimal) of gRNA sequences to not consider, default = .8')
	parser.add_argument('--gc_lower',default=.1,help='GC content lower bound (in decimal) of gRNA sequences to not consider, defulat = .1')
	parser.add_argument('--min_index',default=0,help='Only calculates scores for sequences starting at this index - can save a lot of time.')
	parser.add_argument('--max_index',default=-1,help='Only calculates scores for sequences ending at this index - can save a lot of time.')
	args = parser.parse_args()


	start = time.time()
	start_time = datetime.datetime.now()
	print('Execution begun:',str(start_time))

	cas_list = cas2dict(args.cas_list,args.cas)
	fasta_dict = fasta2dict(args.fasta)
	hits = grna_finder(fasta_dict,cas_list,float(args.gc_upper),float(args.gc_lower))
	scored = mismatch_scoring(hits,cas_list,args.score,int(args.min_index),int(args.max_index))
	grna_dict2bed(args.fasta,scored)

	end_time = datetime.datetime.now()
	end = time.time()

	print('\n',str(len(scored)),'gRNA analyzed.')
	print('Execution finished:',str(end_time))
	print('Execution time:',str(end-start))

	exit('cat.txt')


if __name__ == '__main__':
    main(args)
