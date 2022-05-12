
from distutils.command.config import config
from genericpath import isdir
from posixpath import dirname
import subprocess
import os
from sys import argv


def execute_with(simulator, trace, ssd_config, temperature, peCycle, others):
    if others != "" :
        dir_name = './log' + '/temp' + temperature + '/PE' +  peCycle + '/' + others + '/' + trace
    else :
        dir_name = './log' + '/temp' + temperature + '/PE' +  peCycle + '/' + trace

    if(not os.path.isdir(dir_name)):
        os.makedirs(dir_name)

    subprocess.run([simulator,'./config/' + trace + '.cfg', './simplessd/config/' + ssd_config, dir_name])

if len(argv) < 2:
    exit()

trace = argv[1]
temperature = "25"
peCycle = "3000"
ssd_config = "temp"+ temperature + "_pe" + peCycle + ".cfg"
others = ""
print("########################################################################################")
print("Trace : " + trace + " temperature : " + temperature + ", PE : " + peCycle)
print("########################################################################################")

execute_with("./simplessd-standalone", trace, ssd_config, temperature, peCycle, others)

print("########################################################################################")
print()

trace = argv[1]
temperature = "25"
peCycle = "5000"
ssd_config = "temp"+ temperature + "_pe" + peCycle + ".cfg"
others = ""
print("########################################################################################")
print("Trace : " + trace + " temperature : " + temperature + ", PE : " + peCycle)
print("########################################################################################")

execute_with("./simplessd-standalone", trace, ssd_config, temperature, peCycle, others)

print("########################################################################################")
print()

if len(argv) < 3:
    exit()

trace = argv[2]
temperature = "25"
peCycle = "3000"
ssd_config = "temp"+ temperature + "_pe" + peCycle + ".cfg"
others = ""
print("########################################################################################")
print("Trace : " + trace + " temperature : " + temperature + ", PE : " + peCycle)
print("########################################################################################")

execute_with("./simplessd-standalone", trace, ssd_config, temperature, peCycle, others)

print("########################################################################################")
print()

trace = argv[2]
temperature = "25"
peCycle = "5000"
ssd_config = "temp"+ temperature + "_pe" + peCycle + ".cfg"
others = ""
print("########################################################################################")
print("Trace : " + trace + " temperature : " + temperature + ", PE : " + peCycle)
print("########################################################################################")

execute_with("./simplessd-standalone", trace, ssd_config, temperature, peCycle, others)

print("########################################################################################")
print()
