import argparse
import shutil

if __name__ == '__main__':
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('-i', help='Input files', nargs="+")
    parser.add_argument('-o', help='Output files', nargs="+")

    args = parser.parse_args()
    for i, o in zip(args.i, args.o):
        shutil.move(i, o)
