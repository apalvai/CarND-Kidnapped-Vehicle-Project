/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <random>
#include <map>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    
    cout << "Initializing particle filter..." << endl;
    cout << "x: " << x << " stdev: " << std[0] << endl;
    cout << "y: " << y << " stdev: " << std[1] << endl;
    cout << "theta: " << theta << " stdev: " << std[2] << endl;
    
    num_particles = 100;
    
    std::default_random_engine x_generator;
    std::default_random_engine y_generator;
    std::default_random_engine theta_generator;
    
    std::normal_distribution<double> x_dist(x, std[0]);
    std::normal_distribution<double> y_dist(y, std[1]);
    std::normal_distribution<double> theta_dist(theta, std[2]);
    
    for (int i=0; i<num_particles; i++) {
        Particle particle;
        particle.id = i;
        particle.x = x_dist(x_generator);
        particle.y = y_dist(y_generator);
        particle.theta = theta_dist(theta_generator);
        particle.weight = 1.0;
        
        particles.push_back(particle);
    }
    
    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
    
    cout << "prediction..." << endl;
    
    std::default_random_engine x_generator;
    std::default_random_engine y_generator;
    std::default_random_engine theta_generator;
    
    double c1 = 0.0;
    if (yaw_rate == 0.0) {
        cout << "***** yaw_rate is zero *****" << endl;
    } else {
        c1 = velocity / yaw_rate;
    }
    
    for (int i=0; i<particles.size(); i++) {
        Particle particle = particles[i];
        
        double x = particle.x;
        double y = particle.y;
        double theta = particle.theta;
        
        if (c1 == 0.0) {
            x += velocity * cos(theta);
            y += velocity * sin(theta);
        } else {
            double theta_updated = theta + yaw_rate * delta_t;
            
            x += c1 * (sin(theta_updated) - sin(theta));
            y += c1 * (cos(theta) - cos(theta_updated));
            theta = theta_updated;
        }
        
        std::normal_distribution<double> x_dist(x, std_pos[0]);
        std::normal_distribution<double> y_dist(y, std_pos[1]);
        std::normal_distribution<double> theta_dist(theta, std_pos[2]);
        
        particle.x = x_dist(x_generator);
        particle.y = y_dist(y_generator);
        particle.theta = theta_dist(theta_generator);
        
        particles[i] = particle;
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

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
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
    
    cout << "updating weights..." << endl;
    
    double sigma_dist = std_landmark[0];
    double sigma_theta = std_landmark[1];
    
    // cout << "sensor_range: " << sensor_range << endl;
    // cout << "sigma_dist: " << sigma_dist << endl;
    // cout << "sigma_theta: " << sigma_theta << endl;
    
    double var_dist = pow(sigma_dist, 2);
    
    // cout << "map_landmarks: " << map_landmarks.landmark_list.size() << endl;
    // cout << "observations: " << observations.size() << endl;
    
    for (int i=0; i<particles.size(); i++) {
        
        // std::vector<int> associations;
        // std::vector<double> sense_x;
        // std::vector<double> sense_y;
        
        Particle particle = particles[i];
        
        double theta = particle.theta;
        double weight = particle.weight;
        
        for (int j=0; j<observations.size(); j++) {
            LandmarkObs observation = observations[j];
            
            // transform observations from car's local coordinate system to map's coordinate system
            double obs_x = particle.x*cos(theta) - particle.y*sin(theta) + observation.x;
            double obs_y = particle.x*sin(theta) + particle.y*cos(theta) + observation.y;
            
            // Find nearest neighbor
            Map::single_landmark_s closestLandmark;
            double minDistance;
            
            for (int k=0; k<map_landmarks.landmark_list.size(); k++) {
                Map::single_landmark_s landmark = map_landmarks.landmark_list[k];
                
                double distance = sqrt(pow((landmark.x_f-obs_x),2) + pow((landmark.y_f-obs_y),2));
                
                // begin by assigning the first distance as minDistance or
                // the computed distance is less than the minDistance
                if (k == 0 || distance < minDistance) {
                    minDistance = distance;
                    closestLandmark = landmark;
                }
            }
            
            // sense_x.push_back(obs_x);
            // sense_y.push_back(obs_y);
            // associations.push_back(closestLandmark.id_i);
            
            // compute weight from observation
            double c1 = 1/(2*var_dist);
            double c2 = c1/M_PI;
            double diff_x2 = pow((obs_x-closestLandmark.x_f),2);
            double diff_y2 = pow((obs_y-closestLandmark.y_f),2);
            double c3 = diff_x2+diff_y2;
            double c4 = c1*c3;
            double w_obs = c2 * exp(-c4);
            // cout << "c1: " << c1 << " c3: " << c3 << " c4: " << c4 << " w_obs: " << w_obs << endl;
            weight = weight * w_obs;
        }
        
        // Update associations, sense_x, sense_y
        // particle = SetAssociations(particle, associations, sense_x, sense_y);
        
        // update particle weight
        // cout << "particle weight: " << weight << endl;
        particle.weight = weight;
        
        // update partciles vector with the updated particle
        particles[i] = particle;
    }
    
    double cumm_weight = 0.0;
    for (int i=0; i<particles.size(); i++) {
        Particle particle = particles[i];
        cumm_weight += particle.weight;
    }
    cout << "cummulative weight: " << cumm_weight << endl;
    
    // normalize weights
    for (int i=0; i<particles.size(); i++) {
        Particle particle = particles[i];
        particle.weight = particle.weight/cumm_weight;
        particles[i] = particle;
    }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    
    cout << "resampling..." << endl;
    
    // weights
    std::vector<double> weights;
    for (int i=0; i<particles.size(); i++) {
        Particle particle = particles[i];
        weights.push_back(particle.weight);
    }
    
    // resampling
    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<> d(weights.begin(), weights.end());
    std::map<int, int> m;
    
    for(int n=0; n<particles.size(); ++n) {
        ++m[d(gen)];
    }
    std::vector<Particle> resampled;
    int index = 0;
    for(auto p : m) {
        std::cout << p.first << " generated " << p.second << " times\n";
        
        for (int i=0; i<p.second; i++) {
            Particle oldParticle = particles[p.first];
            
            Particle newParticle;
            newParticle.id = index;
            newParticle.x = oldParticle.x;
            newParticle.y = oldParticle.y;
            newParticle.theta = oldParticle.theta;
            newParticle.weight = 1.0;
            resampled.push_back(newParticle);
            index ++;
        }
    }
    particles = resampled;
    
    for (int i=0; i<particles.size(); i++) {
        Particle particle = particles[i];
        cout << "particle " << particle.id << ": " << particle.x << ", " << particle.y << endl;
    }
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
