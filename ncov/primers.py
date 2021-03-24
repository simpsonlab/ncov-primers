"""
Functions to handle and process primers
"""

import os
import csv
import re
import sys
import pybedtools
import tempfile
import shutil

def import_bed_file(bed):
    """
    Import a BED file containing regions and return a matrix containing those
    regions.

    Arguments:
        * bed: full path to a BED file to process
    
    Return Value:
        Returns a list of primers from the BED file
    """
    primers = list()
    with open(bed, 'r') as fh:
        reader = csv.reader(fh, delimiter='\t')
        for line in reader:
            primers.append(line)
    return primers


def get_amplicon_id(primer, index=3, delimiter='_', idloc=1):
    """
    Extract the amplicon id from the name of the format:
        "virus_idloc__location"
    """
    return str(primer[index]).split(delimiter)[idloc]


def get_primer_direction(primer, delimiter='_', left='left', right='right'):
    """
    Check if the primer is on the left or the right side of the amplicon.
    """
    primer = delimiter.join([str(item) for item in primer])
    try:
        if left in primer.lower():
            return 'LEFT'
        elif right in primer.lower():
            return 'RIGHT'
    except:
        print('Cannot determine LEFT or RIGHT primer side')


def create_primer_info(id):
    """
    Create the data structure for the primer info
    """
    primer = dict()
    primer[id] = dict()
    primer[id]['chr'] = 'MN908947.3'
    primer[id]['qual'] = 60
    primer[id].update(create_forward_primer_info())
    primer[id].update(create_reverse_primer_info())
    return primer


def create_reverse_primer_info():
    """
    Create the data structure for the primer dictionary
    """
    right = {'right' : {'start': set(), 'end': set(), 'strand': '-'}}
    return right


def create_forward_primer_info():
    """
    Create the data structure for the primer dictionary.
    """
    left = {'left' : {'start': set(), 'end': set(), 'strand': '+'}}
    return left


def add_primer(primers, primer):
    """
    Add a primer to the primer dictionary. 
    """
    id = get_amplicon_id(primer=primer)
    if id not in primers.keys():
        primers[id] = dict()
        primers[id]['chr'] = 'MN908947.3'
        primers[id]['qual'] = 60
        primers[id].update(create_reverse_primer_info())
        primers[id].update(create_forward_primer_info())
    primer_direction = get_primer_direction(primer=primer)
    if primer_direction == 'LEFT':
        primers[id]['left']['start'].add(primer[1])
        primers[id]['left']['end'].add(primer[2])
        primers[id]['left']['strand'] = primer[5]
    elif primer_direction == 'RIGHT':
        primers[id]['right']['start'].add(primer[1])
        primers[id]['right']['end'].add(primer[2])


def create_amplicon(primers, type, path='.'):
    """
    Create an amplicon data structure.
    """
    amplicons = list()
    primer_ids = list()
    if type == 'unique_amplicons':
        amplicons = create_unique_amplicon(primers=primers,
                                           path=path,
                                           first=first_primer_id,
                                           last=last_primer_id,
                                           is_unique=True)
    elif type == 'full' or type == 'no_primers':
        for id in primers:
            try:
                if type == 'full':
                    amplicons.append(create_full_amplicon(id=id,
                                                          primer=primers[id],
                                                          first=first_primer_id,
                                                          last=last_primer_id))
                elif type == 'no_primers':
                    amplicons.append(create_no_primer_amplicon(id=id, primer=primers[id]))
            except:
                print(f'Skipping primer {id}...')
    else:
        sys.exit('Invalid amplicon type...')
    return amplicons


def create_full_amplicon(id, primer, first=None, last=None, is_unique=False):
    """
    Create the full amplicon including primers.  For those with
    alternative primers, the amplicon will include the largest
    region possible.

    Arguments:
        * id: the id number of the primer
        * primer: a dictionary containg the primer information
        * first: the id for the first amplicon
        * last: the id for the last amplicon
        * is_unique: a Boolean to apply unique amplicons

    Return Value:
        Returns a list containing primer details including:
            * chr
            * start
            * end
            * name
            * score
            * strand (+)
    """
    amplicon = list()
    amplicon.append(primer['chr'])
    if is_unique:
        if str(id) == str(first):
            amplicon.append(str(max(primer['left']['end'])))
            amplicon.append(str(min(primer['right']['end'])))
        elif str(id) == str(last):
            amplicon.append(str(min(primer['left']['start'])))
            amplicon.append(str(max(primer['right']['start'])))
        else:
            amplicon.append(str(min(primer['left']['start'])))
            amplicon.append(str(max(primer['right']['end'])))
    else:
        amplicon.append(str(min(primer['left']['start'])))
        amplicon.append(str(max(primer['right']['end'])))
    amplicon.append(f'nCoV-2019_{id}')
    amplicon.append(str(primer['qual']))
    amplicon.append('+')
    return amplicon


def create_no_primer_amplicon(id, primer, first, last):
    """
    Create an amplicon without primers.

    Arguments:
        * id: the id/index of the amplicons
        * primer: a primer dictionary
        * first: id of the first primer
        * last: id of the last primer
    
    Return Value:
        Returns an amplicon as a list
    """
    amplicon = list()
    amplicon.append(primer['chr'])
    amplicon.append(str(max(primer['left']['end'])))
    amplicon.append(str(min(primer['right']['start'])))
    amplicon.append(f'nCoV-2019_{id}')
    amplicon.append(str(primer['qual']))
    amplicon.append('+')
    return amplicon


def write_amplicon_to_bed(amplicons, outfile):
    """
    Write the amplicon list as a BED file.

    Arguments:
        * amplicons: a list of amplicons
        * outfile: full path to a file to write the amplicons
    """
    path = os.path.dirname(outfile)
    if path == '':
        path = '.'
    if not os.path.exists(path):
        os.makedirs(path)
    with open(outfile, 'w') as ofh:
        for amplicon in amplicons:
            ofh.write('\t'.join(amplicon))
            ofh.write('\n')
    ofh.close()
    return 0


def create_unique_amplicon(primers, path='.'):
    """
    Create unique amplicon regions by using the full amplicons
    and removing any overlapping regions from other amplicons.
    """
    amplicons = list()
    unique_amplicons = list()
    primer_ids = list()
    for primer_id in primers.keys():
        primer_ids.append(primer_id)
    last_primer_id = max(primer_ids)
    first_primer_id = min(primer_ids)
    #amplicons = create_amplicon(primers=primers, type='full')
    for id in primers:
        amplicons.append(create_full_amplicon(id=id,
                                              primer=primers[id],
                                              first=first_primer_id,
                                              last=last_primer_id,
                                              is_unique=True))
    for index, amplicon in enumerate(amplicons):
        # create a list of amplicons to remove
        tmp_amplicon = list()
        tmp_amplicon.append(amplicon)
        tmp_amplicons = amplicons.copy()
        del(tmp_amplicons[index])
        id_file = f'{index}.bed'
        outfile = '/'.join([path, str(index), id_file])
        write_amplicon_to_bed(amplicons=tmp_amplicon, outfile=outfile)
        amplicon_bed = pybedtools.BedTool(outfile)
        id_remove_file = f'{index}.remove.bed'
        remove_outfile = '/'.join([path, str(index), id_remove_file])
        write_amplicon_to_bed(amplicons=tmp_amplicons, outfile=remove_outfile)

        # use bedtool subtract to remove all other overlapping amplicons from
        # a given amplicon
        amplicon_remove_bed = pybedtools.BedTool(remove_outfile)
        amplicon_unique_file = f'{index}.unique.bed'
        amplicon_unique_outfile = '/'.join([path, str(index), amplicon_unique_file])
        amplicon_bed.subtract(amplicon_remove_bed).saveas(amplicon_unique_outfile)

        # read the unique amplicon files and store the amplicons in a list
        with open(amplicon_unique_outfile, 'r') as amp_file:
            reader = csv.reader(amp_file, delimiter='\t')
            for line in reader:
                unique_amplicons.append(line)
    return unique_amplicons


def main(file, path='./qc_primers'):
    primers = import_bed_file(bed=file)
    primer_dict = dict()
    for primer in primers:
        add_primer(primers=primer_dict, primer=primer)
    return primer_dict

