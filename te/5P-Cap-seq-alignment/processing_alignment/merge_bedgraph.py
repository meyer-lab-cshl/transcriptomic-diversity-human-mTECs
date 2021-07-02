import pandas as pd
import optparse
import pdb

parser = optparse.OptionParser()
parser.add_option('--plus', dest = 'plus', type="string", nargs = 1)
parser.add_option('--minus', dest = 'minus', type="string", nargs = 1)
parser.add_option('--outfile', dest = 'outfile', type="string", nargs = 1)
parser.add_option('--spikeins', dest = 'spikeins', type="string", nargs = 1,
        default="M24537,X04603,X17013")

options, args = parser.parse_args()

required="outfile plus minus".split()
for r in required:
    if options.__dict__[r] is None:
        parser.error("Parameter %s required but not provided"%r)

plus = pd.read_csv(options.plus, sep='\t', names = ["Chr", "Start", "Coverage"])
minus = pd.read_csv(options.minus, sep='\t',
        names = ["Chr", "Start", "Coverage"])
minus['Coverage'] *= -1

both = pd.concat([plus, minus], ignore_index=True)
both = both.sort_values(by=['Chr', 'Start'])

both['End'] = both['Start']+1

# Remove spike in "chromosomes"
if options.spikeins != "":
    spikeins = options.spikeins.split(",")
    for sp in spikeins:
        both = both[~both['Chr'].str.startswith(sp)]

# Keep only standard chromosomes
both = both[both['Chr'].str.startswith('chr')]

both = both[['Chr', 'Start', 'End', 'Coverage' ]]

#save dataframe
both.to_csv(options.outfile, sep='\t', index=False, header=False)
