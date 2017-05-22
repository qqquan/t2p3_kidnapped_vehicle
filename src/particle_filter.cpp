/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 		Updated by QH
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <map>
#include <vector>
#include <unordered_set>
#include <algorithm> 


#include "particle_filter.h"
#include "helper_functions.h"

void ParticleFilter::init(double x, double y, double theta, double std_sigma[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	const int kInitalWeight = 1;	
	const int kSearchScale = 1; //enlarge particle distribution coverage, to counter other noises from the real-world
  	num_particles = 128;

	std::default_random_engine gen;
	std::normal_distribution<double> N_x_init(x, std_sigma[0]*kSearchScale);
	std::normal_distribution<double> N_y_init(y, std_sigma[1]*kSearchScale);
	std::normal_distribution<double> N_theta_init(theta, std_sigma[2]);

	for(int i =0; i<num_particles; i++)
	{
		const double x_prtcl = N_x_init(gen);
		const double y_prtcle = N_y_init(gen);
		const double theta_prtcle = N_theta_init(gen);

		Particle prtcle = {i, x_prtcl, y_prtcle, theta_prtcle, kInitalWeight};

		particles.push_back(prtcle);

	}

	is_initialized=true;



}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: since particle data is already with gaussian noise, why should we have additional random Gaussian noise at prediction step
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	std::default_random_engine gen;
	std::normal_distribution<double> noise_dist_x(0, std_pos[0]);
	std::normal_distribution<double> noise_dist_y(0, std_pos[1]);
	std::normal_distribution<double> noise_dist_theta(0, std_pos[2]);	

	for (auto& a_prtcle : particles)
	{
		double px       = a_prtcle.x;
		double py       = a_prtcle.y;
		double v        = velocity;
		double yaw      = a_prtcle.theta;
		// double yaw_rate = yaw_rate;


		double px_new       = 0.0;
		double py_new       = 0.0;
		double yaw_new      = 0.0;


		//if (fabs(yaw_rate )>0.001)
		if ((yaw_rate )!=0)
		{
         px_new       = px + v/yaw_rate*( sin(yaw + yaw_rate*delta_t) - sin(yaw));
         py_new       = py + v/yaw_rate*(-cos(yaw + yaw_rate*delta_t) + cos(yaw));
		}
		else
		{
		 px_new       = px + v*cos(yaw)*delta_t;
		 py_new       = py + v*sin(yaw)*delta_t; 
		}
      

		yaw_new      = yaw + yaw_rate*delta_t;


		a_prtcle.x = px_new + noise_dist_x(gen);
		a_prtcle.y = py_new + noise_dist_y(gen);
		a_prtcle.theta = yaw_new + noise_dist_theta(gen);



	}

}


//predicted means the map data in the range
void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	for(int i=0; i<observations.size(); i++)
	{
	  auto& obs = observations[i];

      if (predicted.size() >0)
      {
        obs.id = predicted[0].id; //assume the first map landmark is the closest
        double min_dist = dist(predicted[0].x, predicted[0].y, obs.x, obs.y);
        for (unsigned int i = 1; i < predicted.size(); ++i)
        {
          double new_dist = dist(predicted[i].x, predicted[i].y, obs.x, obs.y);
          if (new_dist < min_dist)
          {
            obs.id = predicted[i].id;
            min_dist = new_dist;
          }
        }

        //remove obs that is too far from any landmark in sight
        if (min_dist>10)//TODO: use macro or create new init state 
        {
        	// observations.erase( observations.begin()+i); 
        	obs.id = -1;
        }

      }
      else
      {
      	std::cout <<"!!!!! Error: no predicted observations."<<std::endl;
      }



	}

}

 

static double calc_multivar_gaussian(double x, double y, double mean_x, double mean_y, double sigma_x, double sigma_y)
{

	double dx =  x - mean_x;
	double dy =  y - mean_y;
    // double x_diff = pow((x - mean_x), 2) / (2.0*sigma_x*sigma_x);
    // double y_diff = pow((y - mean_y), 2) / (2.0*sigma_y*sigma_y);
    double nomalizer = 1.0 / (2.0*M_PI*sigma_x*sigma_y);
    double nominator = exp(-( (dx*dx)/(2.0*sigma_x*sigma_x)  + (dy*dy)/(2.0*sigma_y*sigma_y) ));
	return nomalizer*nominator;

}

//TODO: clarify if the std_landmark[] is about position and bearing, or about positions at x and y axises
void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html



	//map data 
	std::map<int, Map::single_landmark_s> landmark_idx_map;		
	for( const auto& l : map_landmarks.landmark_list)
	{
		 landmark_idx_map.insert({l.id_i, l}); 
	}

	//calculate weight for each particle
	double total_weight = 0.0;
	for(auto& p : particles)
	{
	
		// convert sensor data into map coordinates with respect the particle position
		std::vector<LandmarkObs> mapped_obs_list;
		for(const auto& obs : observations)
		{
			LandmarkObs mapped_obs = {0};

			mapped_obs.id = obs.id;
			mapped_obs.x = obs.x*cos(p.theta)   - obs.y*sin(p.theta)  + p.x;
			mapped_obs.y = obs.x*sin(p.theta)   + obs.y*cos(p.theta)  + p.y;

			mapped_obs_list.push_back(mapped_obs);
		}		

		//find the map data in the range, and covert the map data into LandmarkObs class
		std::vector<LandmarkObs> predicted_observations_list;
		for(const auto& l : map_landmarks.landmark_list)
		{
			if (dist(l.x_f, l.y_f, p.x, p.y) < sensor_range)
			{
				predicted_observations_list.push_back(LandmarkObs{ l.id_i,l.x_f,l.y_f });
			}
		}

		std::vector<Map::single_landmark_s> DEBUG_detected_landmarks_list;

		if(predicted_observations_list.size() > 0)
		{
			//find the landmarks
	    	dataAssociation(predicted_observations_list, mapped_obs_list);

			//calculate weight based on the detected/found landmarks
	    	double new_weight = 1;
	    	for(const auto& mapped_obs : mapped_obs_list)
	    	{
	    		if (-1 != mapped_obs.id)
	    		{
		    		Map::single_landmark_s map_landm = landmark_idx_map[mapped_obs.id];

		    		// std::cout<<"Debug:: "<<"Landmark="<<map_landm.id_i << ", "<<map_landm.x_f<<", "<<map_landm.y_f<<";  ";
		    		// std::cout<< "Mapped Obs="<<mapped_obs.x<<", "<<mapped_obs.y<<std::endl;
		    		DEBUG_detected_landmarks_list.push_back(map_landm);

					//calculate gaussian
					new_weight *= calc_multivar_gaussian(	mapped_obs.x,
														mapped_obs.y,
														map_landm.x_f,
														map_landm.y_f,
														std_landmark[0],
														std_landmark[1]
														);
	    		}
	    		else
	    		{
	    			//no landmark found for this observation point
	    		}

	    	}


			total_weight += new_weight; 		

			p.weight = new_weight;
		}


	}

	//normalize the weight
	if(0 != total_weight)
	{
		for (auto& p : particles)
		{	
			p.weight  /= total_weight;
		}	
	}

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	std::vector<double> weight_list;

	for(const auto& p:particles)
	{
		weight_list.push_back(p.weight);
	}

	double weight_sum = std::accumulate(weight_list.begin(), weight_list.end(), 0.0);
	if (0 != weight_sum)
	{


	  std::default_random_engine gen;
	  std::discrete_distribution<int> sampled_idx_distribution({weight_list.begin(), weight_list.end()});

	  std::vector<Particle> particles_old = particles; 

	  for (int i = 0; i < particles.size(); ++i)
	  {
		  int new_idx = sampled_idx_distribution(gen);

		  particles[i] = particles_old[new_idx];

	  }
	}
	else
	{
		std::cout<<"!!!!!!!!!!!!!!! Err: All weights are zero."<<std::endl;
	}





}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
