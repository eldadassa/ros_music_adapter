#include "retina_event_filter.h"

int main(int argc, char** argv)
{

    RetinaFilter retina_filter;
    retina_filter.init(argc, argv);
    retina_filter.runMUSIC();
    retina_filter.finalize();

}

void RetinaFilter::init(int argc, char** argv)
{
    std::cout << "initializing retina filter" << std::endl;

    timestep = DEFAULT_TIMESTEP;
    input_sensor_xdim = 16;
    input_sensor_ydim = 16;
    periphery_downsample_factor = 1;
    fovea_dim = 4;

    num_spikes0 = 0;

    // init MUSIC to read config
    initMUSIC(argc, argv);
}

void RetinaFilter::initMUSIC(int argc, char** argv)
{
    MUSIC::Setup* setup = new MUSIC::Setup (argc, argv);

    setup->config("stoptime", &stoptime);
    setup->config("music_timestep", &timestep);
    setup->config("input_sensor_xdim", &input_sensor_xdim);
    setup->config("input_sensor_ydim", &input_sensor_ydim);
    setup->config("periphery_downsample_factor", &periphery_downsample_factor);
    setup->config("fovea_dim", &fovea_dim);

    acceptable_latency = timestep;
    periphery_xdim = input_sensor_xdim >> periphery_downsample_factor;
    periphery_ydim = input_sensor_ydim >> periphery_downsample_factor;

    port_in = setup->publishEventInput("in");
    port_out = setup->publishEventOutput("out");

    comm = setup->communicator ();
    int rank = comm.Get_rank ();
    int nProcesses = comm.Get_size ();
    if (nProcesses > 1)
    {
        std::cout << "ERROR: num processes (np) not equal 1" << std::endl;
        comm.Abort(1);
    }

    // get dimensions of sensory data and spike data
    if (port_in->hasWidth() && port_out->hasWidth())
    {
        size_spike_data_in = port_in->width();
        size_spike_data_out = port_out->width();

        if (input_sensor_xdim*input_sensor_ydim != size_spike_data_in)
        {
            std::cout << "ERROR: node arguments (sensor dimensions) don't match input port width" << std::endl;
            comm.Abort(1);
        }

        if (periphery_xdim*periphery_ydim + fovea_dim*fovea_dim != size_spike_data_out)
        {
            std::cout << "ERROR: node arguments (fovea/periphery dimensions) don't match output port width" << std::endl;
            comm.Abort(1);
        }
    }
    else
    {
        std::cout << "ERROR: Port-width not defined" << std::endl;
        comm.Abort(1);
    }

    // map linear index to event in port
    MUSIC::LinearIndex l_index_in(0, size_spike_data_in);
    port_in->map(&l_index_in, this, acceptable_latency, 1);

    // map linear index to event out port
    MUSIC::LinearIndex l_index_out(0, size_spike_data_out);
    port_out->map(&l_index_out, MUSIC::Index::GLOBAL, 1);

    MPI::COMM_WORLD.Barrier();
    runtime = new MUSIC::Runtime (setup, timestep);
}

void RetinaFilter::runMUSIC()
{
    std::cout << "running retina filter" << std::endl;

    int num_spikes_transfered = 0;

    struct timeval start;
    struct timeval end;
    gettimeofday(&start, NULL);

    int fovea_min_x = input_sensor_xdim/2 - fovea_dim/2;
    int fovea_max_x = input_sensor_xdim/2 + fovea_dim/2;
    int fovea_min_y = input_sensor_ydim/2 - fovea_dim/2;
    int fovea_max_y = input_sensor_ydim/2 + fovea_dim/2;
    uint periphery_n = periphery_xdim*periphery_ydim;

    double t_spike;
    int id;
    uint receptorx;
    uint receptory;
    uint periphry_x;
    uint periphry_y;
    uint periphry_index;
    double t = runtime->time();
    for (t = runtime->time (); t < stoptime; t = runtime->time ())
    {
        double next_t = t + timestep;
        while (!spikes.empty () && spikes.top ().t < next_t)
        {
            t_spike = spikes.top ().t;
            id = spikes.top ().id;
            receptorx = id % input_sensor_xdim;
            receptory = id / input_sensor_xdim;

            periphry_x = receptorx >> periphery_downsample_factor;
            periphry_y = receptory >> periphery_downsample_factor;
            periphry_index = periphry_y*periphery_xdim + periphry_x;

            port_out->insertEvent(t_spike+0.002, MUSIC::GlobalIndex(periphry_index));

            if (receptorx >= fovea_min_x && receptorx < fovea_max_x &&
                receptory >= fovea_min_y && receptory < fovea_max_y)
            {
               uint fovea_indx = periphery_n + (receptory-fovea_min_y)*fovea_dim + (receptorx-fovea_min_x);
               port_out->insertEvent(t_spike+0.002, MUSIC::GlobalIndex(fovea_indx));
            }

            spikes.pop (); // remove spike from queue
            ++num_spikes_transfered;
        }


//#if DEBUG_OUTPUT

//#endif
        runtime->tick();
    }

    gettimeofday(&end, NULL);
    unsigned int dt = (end.tv_sec*1000000 + end.tv_usec)  - (start.tv_sec*1000000 + start.tv_usec);
    std::cout << "retina filter: total simtime: " << dt/1000000.0 << " received spikes " << num_spikes0 << " spikes decoded " << num_spikes_transfered << std::endl;
}

void RetinaFilter::finalize(){
    runtime->finalize();
    delete runtime;
}

void RetinaFilter::operator () (double t, MUSIC::GlobalIndex id){
    // Decoder: add incoming spikes to map
    num_spikes0++;
    //std::cout <<"decoder: got event"<<std::endl;

    spikes.push (Event (t + acceptable_latency, id));
}
