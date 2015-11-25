import sys
import os
import numpy as np
import datetime
import json

ITERATIONS = 2 
MIN_NUM_NEURONS = 0
MAX_NUM_NEURONS = 2001
STEP_SIZE = 1000

MIN_TIMESTEP = 0.000
MAX_TIMESTEP = 0.011
TIMESTEP_STEP_SIZE = 0.01


sim_time = 10 # in sec
sim_time_build = 0.001 # in sec

data_filename = sys.argv[1]

if os.path.exists(data_filename):
    os.remove(data_filename)


data ={"num_neurons": [], "time": [], "type": [], "iteration": [], "timestep": []}


for timestep in np.arange(MIN_TIMESTEP, MAX_TIMESTEP, TIMESTEP_STEP_SIZE):
    for num_neurons in np.arange(MIN_NUM_NEURONS, MAX_NUM_NEURONS, STEP_SIZE):

        if timestep == 0:
            timestep = 0.001
        if num_neurons == 0:
            num_neurons = 1

        print "\n\n\n\n\ RUNNING", num_neurons, "NEURONS WITH TIMESTEP", timestep, "\n\n\n\n"
    
        music_base_config = \
                    "[sensor]\n\
                      binary=../ros_sensor_adapter\n\
                      args=\n\
                      np=1\n\
                      music_timestep="+ str(timestep) +"\n\
                      ros_topic=/jubot/laserscan\n\
                      message_type=Laserscan\n\
                      sensor_update_rate=30\n\
                    [connect]\n\
                      binary=../connect_adapter\n\
                      args=\n\
                      np=1\n\
                      music_timestep="+ str(timestep) +"\n\
                    [encoder]\n\
                      binary=../rate_encoder\n\
                      args=\n\
                      np=1\n\
                      music_timestep="+ str(timestep) +"\n\
                      rate_min=1\n\
                      rate_max=2\n\
                    [decoder]\n\
                      binary=../linear_readout_decoder\n\
                      args=\n\
                      np=1\n\
                      music_timestep="+ str(timestep) +"\n\
                      music_acceptable_latency="+ str(timestep) +"\n\
                      tau=0.03\n\
                    [command]\n\
                      binary=../ros_command_adapter\n\
                      args=\n\
                      np=1\n\
                      music_timestep="+ str(timestep) +"\n\
                      ros_topic=/jubot/cmd_vel\n\
                      message_mapping_filename=float_mapping.dat\n\
                      command_rate=20\n\
                    sensor.out->connect.in[640]\n\
                    connect.out->encoder.in[" + str(num_neurons) +"]\n\
                    encoder.out->decoder.in[" + str(num_neurons) +"]\n\
                    decoder.out->command.in[2]"
    
        music_config_build = "stoptime=" + str(sim_time_build) + "\n"\
                             + music_base_config
    
        music_config_build_file = open("config_build.music", 'w+')
        music_config_build_file.writelines(music_config_build)
        music_config_build_file.close()
    
        music_config_run = "stoptime=" + str(sim_time) + "\n"\
                           + music_base_config
    
        music_config_run_file = open("config_run.music", 'w+')
        music_config_run_file.writelines(music_config_run)
        music_config_run_file.close()
    
    
        for it in range(ITERATIONS):
    
            start = datetime.datetime.now()
            os.system("mpirun \-np 5 music config_build.music ")
            end = datetime.datetime.now()
    
            dt_build = end - start
            build_time = dt_build.seconds + dt_build.microseconds / 1000000.
    
            data["num_neurons"].append(num_neurons)
            data['type'].append("build-time")
            data['iteration'].append(it)
            data["time"].append(build_time)
            data["timestep"].append(timestep)
    
            start = datetime.datetime.now()
            os.system("mpirun \-np 5 music config_run.music ")
            end = datetime.datetime.now()
            
            dt_run = end - start
            run_time = dt_run.seconds + dt_run.microseconds / 1000000.
    
            data["num_neurons"].append(num_neurons)
            data['type'].append("total-time")
            data['iteration'].append(it)
            data["time"].append(run_time)
            data["timestep"].append(timestep)
    
            rtf = sim_time / (run_time - build_time)
    
            data["num_neurons"].append(num_neurons)
            data['type'].append("real-time factor")
            data['iteration'].append(it)
            data["time"].append(rtf)
            data["timestep"].append(timestep)
        
        if os.path.exists("config_run.music"):
            os.remove("config_run.music")
        if os.path.exists("config_build.music"):
            os.remove("config_build.music")



data_file = open(data_filename, "w+")
json.dump(data, data_file)
data_file.close()


    
        
    

