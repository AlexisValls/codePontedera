#!/usr/bin/env python3
import sys, argparse, os

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
times = [(startTime + i*deltaTime) for i in range(int((endTime-startTime)/deltaTime)+1)]
nSnap = len(times)

# Preprocessing the data
for i in range(len(paths)):
    # If it is the first path, also link 0, system and constant folders
    if i == 0:
        os.system('ln -s ' + paths[i] + '/0 ITHACAoutput/Offline/0')
        os.system('ln -s ' + paths[i] + '/constant ITHACAoutput/Offline/constant')
        os.system('ln -s ' + paths[i] + '/system ITHACAoutput/Offline/system')
    
    # Retrieve the time directories
    timeDirs = os.listdir(paths[i])
    timeDirs = [timeDirs[j] for j in range(len(timeDirs)) if timeDirs[j].replace('.', '', 1).isdigit()] # Remove the non-digit directories


    # Link the wanted times
    tempTimes = times.copy()
    for timeDir in timeDirs:
        if float(timeDir) in tempTimes:
            index=times.index(float(timeDir))
            os.system('ln -s ' + paths[i] + '/' + timeDir + ' ITHACAoutput/Offline/' + str(index+1+(nSnap*i)))
            tempTimes.remove(float(timeDir))

    # If some times are missing, raise an error
    if len(tempTimes) != 0:
        print('Missing times in ' + paths[i] + ' : ' + str(tempTimes))
        sys.exit(1)
