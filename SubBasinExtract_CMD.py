import argparse
import os
import shutil

import SubBasinExtract


parser = argparse.ArgumentParser('SubbasinExtraction_CMD.py')
parser.add_argument('-d', "--dem", required=True, help='dem path')
parser.add_argument('-w', "--wsp", required=False, default="./", type=str, help='path of workspace.default ./')
parser.add_argument('-p', '--pshp', required=True, type=str,
                    help='shpfile path of the position of hydrological station.')
parser.add_argument('-r', '--rad', required=False, default=1000, type=int,
                    help='search radius in snapping pour point.default 1000')
parser.add_argument('-a', '--acc', required=False, default=100, type=int,
                    help='Threshold of waster accumulation above which river comes into being.default 100 ')
parser.add_argument('-l', '--luse', required=False, type=str, help='Landuse data file of USGS classfication standard.')
parser.add_argument('-s', '--sub', required=False, type=str, help='Parameter file of river in subbasin. ')
args = parser.parse_args()
if os.path.exists(args.wsp):
    shutil.rmtree(args.wsp)
    os.mkdir(args.wsp)
else:
    os.mkdir(args.wsp)

print("   ____ ____  _____ _   _ __  __")
print(" / ___| __ )| ____| | | |  \/  |")
print("| |  _|  _ \|  _| | |_| | |\/| |")
print("| |_| | |_) | |___|  _  | |  | |")
print(" \____|____/|_____|_| |_|_|  |_|")
print("                                ")
print("   SUBBASIN EXTRACTION PROGRAM   ")
print("---------------------------------")

SubBasinExtract.Main(args.dem, args.pshp, args.wsp, args.rad, args.acc, landusePath=args.luse, subdat=args.sub)
