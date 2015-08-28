EXAMPLE = '''
Example run:

    python generate-abc-file.py ii-nh3

where:
    ii-nh3  --- is the prefix for which the output will be produced
'''

import sys
import smallFunctions

try:
    prefix = sys.argv[1]
except IndexError:
    print EXAMPLE
    sys.exit(1)

for angle in smallFunctions.drange(-179.00,-59.5,0.5):
    print "python forceAnalysis.py %.2f %s &> wyput_%s_%.2f" % (angle,prefix,prefix,angle)



for angle in smallFunctions.drange(-179.00,-59.5,0.5):
    print "python bootstrappingError.py %.2f %s &> wyput_%s_%.2f" % (angle,prefix,prefix,angle)
