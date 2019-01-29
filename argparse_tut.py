import argparse
import random as rand



parser=argparse.ArgumentParser(description='Compute (a-b)^q')
parser.add_argument('a',  default=1,
	help='1st number', type=float)
parser.add_argument('b',
	help='2nd number', type=float)
parser.add_argument('-q', '--q', default=1,  help='exponent (q=1 by default)', type=float)
parser.add_argument("-v", "--verbose", help="increase output verbosity",
                    action="store_true")

##Parsing command line arguments.
args=parser.parse_args()
if args.verbose:
	print "a: {0} ".format(args.a)+"b: {0} ".format(args.b)+"q: {0} ".format(args.q)+"(a-b)^q: {0}".format((args.a-args.b)**args.q)
else:
	print (args.a-args.b)**args.q
