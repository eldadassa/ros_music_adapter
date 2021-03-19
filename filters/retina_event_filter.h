#include <music.hh>
#include <mpi.h>

#include <queue>
//#include <cmath>
//#include <unistd.h>
//#include <iostream>
//#include <fstream>
#include "sys/time.h"
//#include <json/json.h>

#define DEBUG_OUTPUT false 

const double DEFAULT_TIMESTEP = 1e-3;

class RetinaFilter : MUSIC::EventHandlerGlobalIndex{
    public:
        void init(int argc, char** argv);
        void runMUSIC();
        void finalize();

    private:
	class Event {
	public:
	  double t;
	  int id;
	  Event (double t_, int id_) : t (t_), id (id_) { }
	  bool operator< (const Event& other) const { return t > other.t; }
	};

        MPI::Intracomm comm;
        MUSIC::Runtime* runtime;
        double stoptime;
        double timestep;
        double acceptable_latency;
        int input_sensor_xdim;
        int input_sensor_ydim;
        int fovea_dim;
        int periphery_downsample_factor;

        uint periphery_xdim;
        uint periphery_ydim;
        int size_spike_data_in;
        int size_spike_data_out;
        uint num_spikes0;

        MUSIC::EventInputPort* port_in;
        MUSIC::EventOutputPort* port_out;
        std::priority_queue<Event> spikes;

        void initMUSIC(int argc, char** argv);
        void operator() (double t, MUSIC::GlobalIndex id );
};
