import subprocess
from sys import argv
#out=subprocess.check_output('ls')
def execute(trace):
    simulator_config = [(trace + "_3000_50.cfg"), (trace + "_5000_50.cfg")]
    ssd_config = [("temp50_pe3000.cfg"), ("temp50_pe5000.cfg")]

    for i in range(len(simulator_config)):
        print('./config/' + simulator_config[i], './simplessd/config/' + ssd_config[i])
        subprocess.run(['./simplessd-standalone','./config/' + simulator_config[i], './simplessd/config/' + ssd_config[i], './log'])
        #output_file = open('test_output.txt', "at")
        #output_file.write("####################### " + simulator_config[i] + " #################################\n")
        #output_file.write(out.decode('ascii'))
        #output_file.write("####################### " + "simulation " + str(i+1) + "/" + str(len(simulator_config)) + " ended" + " #################################\n\n")
        #output_file.close()

def execute_with(simulator, simulator_config, ssd_config):
    subprocess.run([simulator,'./config/' + simulator_config, './simplessd/config/' + ssd_config, './log'])
    #out = subprocess.check_output([simulator,'./config/' + simulator_config, './simplessd/config/' + ssd_config, './log'])
    '''
    output_file = open(simulator_config + '_output.txt', "at")
    output_file.write("######################### " + simulator_config + " ####################################\n")
    output_file.write(out.decode('ascii'))
    output_file.write("####################### " + simulator_config  + " ended" + " #################################\n\n")
    output_file.close()
    '''


#execute (argv[1])

#execute_with("./simplessd-standalone", argv[1] + "_3000_50.cfg", "3000PE_50.cfg")
execute_with("./simplessd-standalone", argv[1] + "_5000_50.cfg", "5000PE_50.cfg")
#print(argv[1] + "_" + argv[2] + ".cfg", "temp25_pe" + argv[2] + ".cfg")