#!/Users/kalindbl/bin/virtualenv/biopython/bin/python

''' 

Take an XML database of rearrangement maps and create a text file containing the scrambling order for those that are a) valid and b)scrambled.

'''

from lxml import etree
import argparse

parser = argparse.ArgumentParser(description='extract valid scrambling patterns')
parser.add_argument('input', help='input MDS/IES database, xml format')
parser.add_argument('output', help='file to which to write output')
args = parser.parse_args()

# Output file
outfile = open(args.output, 'w')

# Open XML annotations
tree = etree.parse(args.input).getroot()

# Select all valid, scrambled maps
valid_scrambled = tree.xpath('/mapSet/map[valid and scrambled]')

# A valid contig may have multiple scrambled maps,
# so each one needs to be identified by *both* MIC and MAC contig
for valid_map in valid_scrambled:
    mic = valid_map.xpath('./mic/text()')[0]
    mac = valid_map.xpath('./mac/text()')[0]

    scrambling = valid_map.xpath('./scrambled/@order')[0]

    outfile.write('{0}:{1}\t{2}\n'.format(mic, mac, scrambling))

outfile.close()
