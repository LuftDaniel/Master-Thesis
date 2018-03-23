import numpy as np
import datetime
import os

__this_files_path = os.path.realpath(__file__)
__this_files_dir  = os.path.dirname(__this_files_path)
outputfolder = os.path.join(__this_files_dir,
                            "testordner")
os.mkdir(outputfolder)

testfile = open(os.path.join(outputfolder, "testcsv"), 'w' )
testfile.write("was geht ab \n hi scooter")
print(outputfolder)