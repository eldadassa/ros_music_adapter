#include "ros_sensor_adapter.h"

#include "rtclock.h"

#include <math.h>
#include <regex>

#define PI 3.14159265

static void*
ros_thread(void* arg)
{
    RosSensorAdapter* ros_adapter = static_cast<RosSensorAdapter*>(arg);
    ros_adapter->runROS();
}

int
main(int argc, char** argv)
{

    RosSensorAdapter ros_adapter;
    ros_adapter.init(argc, argv);

    MPI::COMM_WORLD.Barrier();
    // If sensor_update_rate and timestep match to a relative
    // precision of 0.1%, lump the ROS and MUSIC event loops
    // together.
    if (ros_adapter.ratesMatch (0.001))
    {
      ros_adapter.runROSMUSIC();
    }
    else
    {
      pthread_t t;
	    pthread_create (&t, NULL, ros_thread, &ros_adapter);

    	ros_adapter.runMUSIC();
    	pthread_join(t, NULL);
    }

    ros_adapter.finalize();

}


bool
RosSensorAdapter::ratesMatch (double precision)
{
    return std::abs (sensor_update_rate * timestep - 1.) < precision;
}


void
RosSensorAdapter::init(int argc, char** argv)
{
    std::cout << "initializing ROS sensor adapter" << std::endl;

    timestep = DEFAULT_TIMESTEP;
    sensor_update_rate = DEFAULT_SENSOR_UPDATE_RATE;
    ros_node_name = DEFAULT_ROS_NODE_NAME;
    rtf = DEFAULT_RTF;

    pthread_mutex_init(&data_mutex, NULL);

    // MUSIC before ROS to read the config first!
    initMUSIC(argc, argv);
    initROS(argc, argv);
}


void
RosSensorAdapter::initROS(int argc, char** argv)
{
    //std::cout << "initROS" << std::endl;

    ros::init(argc, argv, ros_node_name);
    ros::start();

    ros::NodeHandle n;
    switch (msg_type)
    {
        case Laserscan:
            subscriber = n.subscribe(ros_topic, 1000, &RosSensorAdapter::laserscanCallback, this);
            break;
        case Twist:
            subscriber = n.subscribe(ros_topic, 1000, &RosSensorAdapter::twistCallback, this);
            break;
        case Float64MultiArray:
            subscriber = n.subscribe(ros_topic, 1000, &RosSensorAdapter::float64MultiArrayCallback, this);
            break;
        case LinkStates:
            subscriber = n.subscribe(ros_topic, 1000, &RosSensorAdapter::gazeboLinkStatesCallback, this);
            break;
    }
}

void
RosSensorAdapter::initMUSIC(int argc, char** argv)
{
  //std::cout << "initMUSIC" << std::endl;

  setup = new MUSIC::Setup (argc, argv);

  setup->config("ros_topic", &ros_topic);
  setup->config("stoptime", &stoptime);
  setup->config("music_timestep", &timestep);
  setup->config("sensor_update_rate", &sensor_update_rate);
  setup->config("ros_node_name", &ros_node_name);
  setup->config("rtf", &rtf);

    std::string _msg_type;
    setup->config("message_type", &_msg_type);

    //std::cout << "_msg_type: "<<_msg_type<< std::endl;//debug

    if (_msg_type.compare("Laserscan") == 0){
        msg_type = Laserscan;
    }
    else if (_msg_type.compare("Twist") == 0){
        msg_type = Twist;
    }
    else if (_msg_type.compare("FloatArray") == 0){
        msg_type = Float64MultiArray;
    }
    else if (_msg_type.compare("LinkStates") == 0){
        msg_type = LinkStates;
        std::string link_names, link_infs;
        setup->config("link_names", &link_names);
        setup->config("link_information", &link_infs);

        link_name_vec = split(link_names, ",");
        std::vector<std::string> link_inf_strs = split(link_infs, ",");

        //std::cout << "link_names: "<<link_name_vec[0]<<" "<<link_name_vec[1]<<std::endl; //debug
        //std::cout << "link_information: "<<link_inf_strs[0]<<" "<<link_inf_strs[1]<< std::endl; //debug

        for (int i = 0; i < link_inf_strs.size(); ++i)
        {
           if (link_inf_strs[i].compare("position") == 0){
             //std::cout << "position1"<<std::endl;//debug
              link_inf_vec.push_back(Position);
              //std::cout << "position2"<<std::endl;//debug
           }
           else if (link_inf_strs[i].compare("azz") == 0){
             //std::cout << "azz1"<<std::endl;//debug
              link_inf_vec.push_back(Azz);
              //std::cout << "azz2"<<std::endl;//debug
           }
           else {
             std::cout << "ERROR: link information unknown" << std::endl;
             finalize();
           }
        }
    }
    else
    {
        std::cout << "ERROR: msg type unknown" << std::endl;
        finalize();
    }

    //std::cout << "publishContOutput" << std::endl;//debug
    MUSIC::ContOutputPort* port_out = setup->publishContOutput ("out");

    comm = setup->communicator ();
    int rank = comm.Get_rank ();
    int nProcesses = comm.Get_size ();
    if (nProcesses > 1)
    {
        std::cout << "ERROR: num processes (np) not equal 1" << std::endl;
        comm.Abort(1);
    }


    if (port_out->hasWidth ())
    {
        datasize = port_out->width ();
    }
    else
    {
        std::cout << "ERROR: Port-width not defined" << std::endl;
        comm.Abort (1);
    }

    data = new double[datasize];
    for (unsigned int i = 0; i < datasize; ++i)
    {
        data[i] = 0.;
    }

    // Declare where in memory to put data
    MUSIC::ArrayData dmap (data,
      		 MPI::DOUBLE,
      		 0,
      		 datasize);
    port_out->map (&dmap, 1);
}

void
RosSensorAdapter::runROSMUSIC()
{
    std::cout << "running sensor adapter with update rate of " << sensor_update_rate << std::endl;
    RTClock clock( 1. / (sensor_update_rate * rtf) );

    ros::spinOnce();
    runtime = new MUSIC::Runtime (setup, timestep);

    for (int t = 0; runtime->time() < stoptime; t++)
    {

#if DEBUG_OUTPUT
        std::cout << "ROS Sensor Adapter: ";
        for (int i = 0; i < datasize; ++i)
        {
            std::cout << data[i] << " ";
        }
        std::cout << std::endl;
#endif

        clock.sleepNext();
        ros::spinOnce();
        runtime->tick();
    }

    std::cout << "sensor: total simtime: " << clock.time () << " s" << std::endl;
}

void
RosSensorAdapter::runROS()
{
    RTClock clock( 1. / (sensor_update_rate * rtf) );

    // wait until first sensor update arrives
    while (ros::Time::now().toSec() == 0.)
    {
        clock.sleepNext();
    }

    ros::Time stop_time = ros::Time::now() + ros::Duration(stoptime/rtf);

    ros::spinOnce();
    for (ros::Time t = ros::Time::now(); t < stop_time; t = ros::Time::now())
    {
#if DEBUG_OUTPUT
        std::cout << "ROS Sensor Adapter: ";
        for (int i = 0; i < datasize; ++i)
        {
            std::cout << data[i] << " ";
        }
        std::cout << std::endl;
#endif

        clock.sleepNext();
        ros::spinOnce();
   }
}

void
RosSensorAdapter::runMUSIC()
{
    std::cout << "running sensor adapter with update rate of " << sensor_update_rate << std::endl;
    RTClock clock(timestep / rtf);

    runtime = new MUSIC::Runtime (setup, timestep);

    for (int t = 0; runtime->time() < stoptime; t++)
    {
        clock.sleepNext();
    	pthread_mutex_lock(&data_mutex);
        runtime->tick();
	    pthread_mutex_unlock(&data_mutex);
    }

    std::cout << "sensor: total simtime: " << clock.time () << " s" << std::endl;
}

void
RosSensorAdapter::laserscanCallback(const sensor_msgs::LaserScanConstPtr& msg)
{
    pthread_mutex_lock(&data_mutex);
    for (unsigned int i = 0; i < msg->ranges.size(); ++i)
    {
      // scale data between -1 and 1
      // TODO: catch exception if ranges.size not width of port
      if (isinf(msg->ranges.at(i))){
	data[i] = 1.;
      }
      else{
	data[i] = ((msg->ranges.at(i) - msg->range_min) / (msg->range_max - msg->range_min) ) * 2 - 1;
      }
    }
    pthread_mutex_unlock(&data_mutex);
}

void
RosSensorAdapter::twistCallback(const geometry_msgs::Twist msg)
{
    pthread_mutex_lock(&data_mutex);

    data[0] = msg.linear.x;
    data[1] = msg.angular.z;
    for (unsigned int i = 0; i < 2; ++i) // Twist msg has 2 dimensions
    {
        // limit data between -1 and 1
        if (data[i] > 1)
            data[i] = 1;
        else if (data[i] < -1)
            data[i] = -1;

    }

    pthread_mutex_unlock(&data_mutex);
}

void
RosSensorAdapter::float64MultiArrayCallback(const std_msgs::Float64MultiArray msg)
{
    pthread_mutex_lock(&data_mutex);

    for (unsigned int i = 0; i < datasize; ++i)
    {
        data[i] = msg.data[i];
    }

    pthread_mutex_unlock(&data_mutex);
}

void
RosSensorAdapter::gazeboLinkStatesCallback(const gazebo_msgs::LinkStates &msg)
{
    pthread_mutex_lock(&data_mutex);

    int data_indx = 0;
    for (int j = 0; j < link_name_vec.size(); ++j)
    {
       int link_indx = -1;
       for (int i = 0; i < msg.name.size(); ++i)
       {
          if (msg.name[i] == link_name_vec[j])
              //ROS_INFO("msg name: [%s]\n", msg.name[i].c_str());
              link_indx = i;
       }

       if (link_indx == -1)
       {
         //  ROS_ERROR_STREAM("Failed to find link "<<link_name<<".");
         std::cout <<"Failed to find link "<<link_name_vec[j]<<"."<<std::endl;
         comm.Abort (1);
          //return;
       }

       if (link_inf_vec.size() > j) {

         switch (link_inf_vec[j])
         {
           case Position:
              if (data_indx+2 < datasize) {
                 data[data_indx++] = msg.pose[link_indx].position.x;
                 data[data_indx++] = msg.pose[link_indx].position.y;
                 data[data_indx++] = msg.pose[link_indx].position.z;
               }
              break;
           case Azz:
              //ROS_INFO("indx: %d", link_indx);
              //ROS_INFO("orientation- x: %f, y: %f, z: %f, w: %f", msg.pose[link_indx].orientation.x, msg.pose[link_indx].orientation.y, msg.pose[link_indx].orientation.z, msg.pose[link_indx].orientation.w);

              //Calculating azz - caculating rotation angle assuming the rotation is around the z axis
              double theta = 2.0 * atan2( sqrt(pow(msg.pose[link_indx].orientation.x,2) +
              pow(msg.pose[link_indx].orientation.y,2) + pow(msg.pose[link_indx].orientation.z,2)),
              msg.pose[link_indx].orientation.w) - PI/2;

              if (data_indx < datasize)
                 data[data_indx++] = theta;

              break;
         }
       }
  }

  pthread_mutex_unlock(&data_mutex);
}

void RosSensorAdapter::finalize(){
    runtime->finalize();
    delete runtime;
}

std::vector<std::string>RosSensorAdapter::split(const std::string str, const std::string regex_str)
{
    std::regex regexz(regex_str);
    std::vector<std::string> list(std::sregex_token_iterator(str.begin(), str.end(), regexz, -1),
                                  std::sregex_token_iterator());
    return list;
}
