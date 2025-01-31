# !/usr/bin/env python3
import sys, argparse

# Retrieve the arguments
args = argparse.ArgumentParser(description='Offline preprocessing of the data')
startTime = args.add_argument('--startTime', type=str, help='Start time of the data to be processed')
endTime = args.add_argument('--endTime', type=str, help='End time of the data to be processed')
paths = args.add_argument('--paths', type=str, help='Paths of the data to be processed', default='.')
args = args.parse_args()

# Check if the arguments are correct
if args.startTime is None:
    print('Please provide a start time')
    sys.exit(1)
if args.endTime is None:
    print('Please provide an end time')
    sys.exit(1)

# Test the arguments
print('Start time:', args.startTime)
print('End time:', args.endTime)
print('Paths:', args.paths)
