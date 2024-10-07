import optparse
import sys
import argparse


def menu(args):
    parser2 = argparse.ArgumentParser(prog='Utils')
    parser2.add_argument("--xx", nargs='+', action='append')
    parser2.add_argument("--zz", nargs='+', action='append')
    argsxx = parser2.parse_args()
    print(argsxx.xx)

'''
    # create OptionParser object
    parser = optparse.OptionParser()
    # add options
    parser.add_option('--sample', dest='samples',
                      type='string',
                      action='append',
                      default=[],
                      help='sample ids')
    parser.add_option('--case', dest='cases',
                      type='string',
                      action='append',
                      help='case ids')
    parser.add_option('-i', '--input_file', dest='i_file',
                      type='string',
                      action="append",
                      help='specify the n''th input files')
    parser.add_option('-m', '--mutations_file', dest='mutations',
                      type='string',
                      help='specify the mutations file')
    parser.add_option("--cnv_file", dest='cnv_file',
                      type='string',
                      action="append",
                      help='specify the n''th cnv files')

    (options, args) = parser.parse_args()
    if options.samples is not None:
        print(options.samples)
    if options.cases is not None:
        print(options.cases)
    if options.i_file is not None:
        print(options.i_file)
    if options.mutations is not None:
        print(options.mutations)
    if options.cnv_file is not None:
        print(options.cnv_file)'''


if __name__ == '__main__':
    menu(sys.argv[1:])
