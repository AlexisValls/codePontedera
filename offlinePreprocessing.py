#!/usr/bin/env python3
import sys, argparse, os

def linkFolder(source, dest):
    # Recursively link the source folder to the destination folder, keeping the same structure
    if not os.path.exists(dest):
        os.makedirs(dest)
    for item in os.listdir(source):
        if os.path.isdir(source+'/'+item):
            linkFolder(source+'/'+item, dest+'/'+item)
        elif os.path.isfile(source+'/'+item):
            os.system('ln ' + source+'/'+item + ' ' + dest+'/'+item)

# Define the arguments
args = argparse.ArgumentParser(description='Offline preprocessing of the data')
startTime = args.add_argument('--startTime', type=float, help='Start time of the data to be processed')
endTime = args.add_argument('--endTime', type=float, help='End time of the data to be processed')
deltaTime = args.add_argument('--deltaTime', type=float, help='Time step between two snapshots')
paths = args.add_argument('--paths', type=str, help='Paths of the data to be processed separated by commas : path1,path2,path3', default='.')

# Retrieve the arguments
args = args.parse_args()
startTime = args.startTime
endTime = args.endTime
deltaTime = args.deltaTime
paths = args.paths.split(',')

# Check if the arguments are correct
if startTime is None:
    print('Please provide a start time')
    sys.exit(1)
if endTime is None:
    print('Please provide an end time')
    sys.exit(1)
if deltaTime is None:
    print('Please provide a time step')
    sys.exit(1)

# Cleaning the ITHACAoutput/Offline folder
if not os.path.exists('ITHACAoutput/Offline'):
    os.makedirs('ITHACAoutput/Offline')
elif os.path.exists('ITHACAoutput/Offline'):
    os.system('rm -rf ITHACAoutput/Offline/*')

# Wanted times
tmpTimeDirs = []
current = startTime
nDecimal = len(str(deltaTime).split('.')[1])
for i in range(int((endTime-startTime)/deltaTime)+1): # First list with .0 for integers
    tmp = str(current).split('.')
    tmpTimeDirs.append(tmp[0]+'.'+tmp[1][:nDecimal])
    current += deltaTime
timeDirs = []
for timeDir in tmpTimeDirs: # Final list without .0 for integers
    if float(timeDir) == int(float(timeDir)):
        timeDirs.append(str(int(float(timeDir))))
    else:
        timeDirs.append(timeDir)

# Preprocessing the data
for i in range(len(paths)):
    # If it is the first path, also link 0, system and constant folders
    if i == 0:
        linkFolder(paths[i]+'/0', 'ITHACAoutput/Offline/0')
        linkFolder(paths[i]+'/system', 'ITHACAoutput/Offline/system')
        linkFolder(paths[i]+'/constant', 'ITHACAoutput/Offline/constant')
    
    # Link the wanted times
    for j in range(len(timeDirs)):
        linkFolder(paths[i]+'/'+timeDirs[j], 'ITHACAoutput/Offline/' + str(j+1+len(timeDirs)*i))
