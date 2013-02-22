#!/Users/kalindbl/bin/virtualenv/biopython/bin/python

''' 

Take an xml file of MDS-IES annotations and produce a new XML file with these added annotations:

- whether the MDS annotation is "valid," i.e. covers the MAC chromosome completely with only pointer overlap
- which MAC contigs come from multiple MIC loci
- which MAC contigs come from scrambled loci
- if scrambled, the sequence of MDS at the MIC locus

Voted least efficient script of the year by everyone ever.

'''

from lxml import etree
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser(description='annotate MDS/IES descriptors with scrambling information')
parser.add_argument('input', help='input MDS/IES database, xml format')
parser.add_argument('contig_lengths', help='text file containing the lengths of MAC contigs, sans telomeres')
parser.add_argument('--coverage_cutoff', help='number of BP on the MAC chromosome that can be unaccounted for before it\'s called invalid.', type=int, default=50)
parser.add_argument('--overlap_cutoff', help='percent of an MDS that must be covered by another before it\'s called invalid.', type=float, default=0.8)
parser.add_argument('output', help='file to which to write (xml) output')
args = parser.parse_args()

def validate_contig(mds_map, contig_length):
    '''
        Takes the MDS structure of a contig and determines whether it is adequately covered, without overlap. Returns true if contig is valid.
    '''

    # First, we want the MDS annotations to extend all along the chromosome
    if mds_map[0][0] > args.coverage_cutoff or contig_length - mds_map[-1][1] > args.coverage_cutoff:
        # This does run into the problem that because the
        # MDS were sorted with respect to final endpoint, the first
        # on the list may not have the lowest coordinate.
        # However, if the endpoint of the actual first MDS is
        # lower than the endpoint of the second MDS, the two overlap
        # and should be discarded anyway.
        return False

    # Are there gaps greater than the coverage cutoff
    # within the sequence? Are there overlaps?
    # Move pairwise over the MDS...
    for pair in zip(mds_map, mds_map[1:]):
        # Determine gap
        gap = pair[1][0] - pair[0][1]
        if gap > args.coverage_cutoff:
            return False
        if gap < 0 and abs(gap) > (0.8 * pair[0][1] - pair[0][0]):
            return False

    # If all that passed, the contig is valid
    return True

def get_scrambling(mac_mds_map, mic_mds_map):
    '''
        Take MDS structure for a MAC contig at a particular MIC locus and determine whether it's out of order. Returns the scrambling pattern if yes, or false otherwise.
    '''

    # Because neither the MAC nor the MIC is necessarily
    # in order in the XML, this is going to be a clumsy
    # sorting operation

    # Label both sequences with indices
    # (We only need the end coordinate for each MDS,
    # since that's what we're going to sort on anyway.)
    mac = zip([tup[1] for tup in mac_mds_map], xrange(len(mac_mds_map)))
    mic = zip([tup[1] for tup in mic_mds_map], xrange(len(mic_mds_map)))

    # Now sort both 5'-3'
    mac.sort(key=lambda x: x[0])
    mic.sort(key=lambda x: x[0])

    # If their indices no longer match, the contig is scrambled...
    new_mac_idx = [tup[1] for tup in mac]
    new_mic_idx = [tup[1] for tup in mic]

    if new_mac_idx != new_mic_idx:
        # If they're mirrors, that just means the MAC contig is backwards...
        if sorted(new_mac_idx) != new_mic_idx:
            # ...and the order of scrambling is the order of
            # the sorted MIC indices sorted by the sorted MAC indices! Whew!
            return ':'.join([str(new_mic_idx[idx]) for idx in new_mac_idx])

    return False

mac_contigs = defaultdict(list)
mic_loci = defaultdict(list)
multi_locus_contigs = []
valid_maps = []
scrambled_maps = {}

# Make the length dict
lengths = {line.split('\t')[0]:int(line.split('\t')[1].strip()) for line in open(args.contig_lengths)}

# Pull in the XML annotations
parser = etree.XMLParser(remove_blank_text=True)
tree = etree.parse(args.input, parser)

maps = tree.getiterator('map')

for map in maps:
    map_id = map.get('id')
    mic = []
    mac = []
    contig = map[1].text
    mds = map.getiterator('mds')
    for indv_mds in mds:
        mic.append((int(indv_mds[0].get('start')), int(indv_mds[0].get('end'))))
        mac.append((int(indv_mds[1].get('start')), int(indv_mds[1].get('end'))))

    # Mark this as a scrambled map if appropriate
    scrambling = get_scrambling(mac, mic)
    if scrambling:
       scrambled_maps[map_id] = scrambling

    if len(mic_loci[contig]):
        # flag this as a multi-locus contig
        if contig not in multi_locus_contigs:
            multi_locus_contigs.append(contig)

    mac_contigs[contig] += mac
    mic_loci[contig] += mic

# Now go through and validate all contigs
for contig in mac_contigs:
    # Since MDS may be coming from several separate maps, they need to be sorted by coordinate
    mac_contigs[contig].sort(cmp=lambda x, y: cmp(x[1], y[1]))
    if validate_contig(mac_contigs[contig], lengths[contig]):
        valid_maps.append(contig)

# Annnnd go through the tree *again* to add appropriate annotations for output
for map in tree.getiterator('map'):
    if map.get('id') in scrambled_maps:
        map.insert(2, etree.Element('scrambled', order=scrambled_maps[map.get('id')]))
    if map[1].text in valid_maps:
        map.insert(2, etree.Element('valid'))

# Finally, write this monstrosity back to a file
outfile = open(args.output, 'w')
outfile.write(etree.tostring(tree.getroot(), xml_declaration=True, pretty_print=True))
outfile.close()
