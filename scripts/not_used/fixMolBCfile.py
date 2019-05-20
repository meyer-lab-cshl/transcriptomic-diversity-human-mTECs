import sys, optparse, re

parser = optparse.OptionParser()
parser.add_option('--infile', dest = 'infile', type="string", nargs = 1)
parser.add_option('--outfile', dest = 'outfile', type="string", nargs = 1)

options, args = parser.parse_args()

infile = open(options.infile,'r')
outfile=open(options.outfile,'w')
outfile.write('read_name\tmolecular_bc\n')

line = infile.readline()
while line!='':
	name=re.split('\t', re.sub('[\r\n]','',line))[0]
	name2=re.split(' ', name)[0]
	bc=re.split('\t', re.sub('[\r\n]','',line))[1]
	outfile.write( name2 + '\t' + bc + '\n' )
	line = infile.readline()

outfile.close()
infile.close()
