---
title:  "Particle Filter"
categories: post
mathjax: true
---
![particle_filter](https://github.com/SeokLeeUS/seokleeus.github.io/raw/master/_images/_Particle_filter/particle_filter_gif.gif)
## Summary:
 - This project is to implement 2D particle filter in C++. The particle filter will have an initial condition (GPS data) along with map and observation/controls data..

## Particle filter simulation result (please click the below thumbnail):
[![video result](https://img.youtube.com/vi/phv0iByFB4E/hqdefault.jpg)](https://youtu.be/phv0iByFB4E) 

## Why particle filter?

|filter variance  |State space  |belief        |efficiency |in robotics|
|:---             |:-           |:-            |:-         |:-         |
|Histogram filter |discrete     |multi modal   |exponential|approximate|
|Kalman filter    |continuous   |uni modal     |quadratic  |approximate|
|Particle filter  |continuous   |multimodal    |?          |approximate|

## A flow chart describes particle filter impmentation:

![Particlefiiter_flowchart 1](https://github.com/SeokLeeUS/seokleeus.github.io/raw/master/_images/_Particle_filter/particle_filter_algorithm_flowchart.png)

## Steps for implementation:

|sequence|algorithm flow chart|
|:--     |:--                 |
|1       |initialization      |
|2       |prediction          |
|3       |weight updates      |
|4       |resampling          |

### 1. initialization
- Estimate position from GPS input
- This is the most practical way to initialize the particles and generate real time output

### 2. prediction
- Add control input(yaw rate& velocity) for all particles

### 3. update
- Update particle weights using map landmark positions and feature measurements

### 4. resampling
- resample M times drawing a prticle i proportional to its weight. 

### 5. return new particle set to prediciton
- new set of particles representing the Bayes filter posterior probability. 

## Project implementation tips

### 1. initialization 

- sample Gaussian ditribution to determine particles
- consider Gaussian sensor noise around initial GPS poisition and heading estimates
- C++ standard library normal distribution and C++ standard library random engine function 
- C++ random library supports normal distribtion class. The class receives mean (mean()) and standard deviation (stddev()) as inputs. 
- Example:

```
std::default_random_engine eng;
std::normal_distribution<double> dist(5.0,2.0);
// mean 5, standard deviation 2
array<int,10> tmp = {0, }; 
//initialize with zero and set the array size of 10
for (int i = 0; i < 10000; ++i) {
  double number = dist(eng);
  if((number>= 0.0) &&(number <10.0)) ++tmp[int(number)];
}

```
[Refererence](https://blog.naver.com/drvoss/220340760979)

- implementation: 

```
/*
function taking GPS position and initial heading as input.
output: print out to the terminal 3 samples from a normal distribution with mean equal to the GPS position
and initial heading measurements and standard deviation of 2 m for the x and y position and 0.05 radians for the heading of the car.

*/

#include <iostream>
#include <random>

using namespace std;
using std::normal_distribution;

void printSamples(double gps_x, double gps_y, double theta);

int main() {

	// set GPS provided state of the car
	double gps_x = 4983;
	double gps_y = 5029;
	double theta = 1.201;

	printSamples(gps_x, gps_y, theta);
	return 0;

}

void printSamples(double gps_x, double gps_y, double theta) {

	std::default_random_engine gen;
	double std_x, std_y, std_theta;

	/*
	stadard deviation of x: 2m
	standard deviation of y: 2m
	standard deviation of theta: 0.05 radians
	*/

	std_x = 2;
	std_y = 2;
	std_theta = 0.05;

	normal_distribution<double> dist_x(gps_x, std_x);
	normal_distribution<double> dist_y(gps_y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);

	for (int i = 0; i < 3; ++i) {
		double sample_x, sample_y, sample_theta;

		// TODO: sample from these normal distributions like this:
		sample_x = dist_x(gen);
		sample_y = dist_y(gen);
		sample_theta = dist_theta(gen);
		// where "gen" is the random



		cout << "sample" << i + 1 << "" << sample_x << "" << sample_y << sample_theta << endl;

	}


	return;
}
```

### 2. Prediction

- utilize motion model to predict where the vehicle would be at the next step. 

- Here is the excerpt of motion model equation for bicycle model:

![motion_model1](https://github.com/SeokLeeUS/seokleeus.github.io/raw/master/_images/_Particle_filter/motion_model.png)

### 3. Update

![update_step1](https://github.com/SeokLeeUS/seokleeus.github.io/raw/master/_images/_Particle_filter/update_step.png)

### 4.Homogeneous transformation 
![Homogeneous_transform](https://github.com/SeokLeeUS/seokleeus.github.io/raw/master/_images/_Particle_filter/homogeneous_transform.png)

```
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <math.h>
using namespace std;
int main() {
	double x_part, y_part, x_obs, y_obs, theta;
	x_part = 4;
	y_part = 5;
	x_obs = 2;
	y_obs = 2;
	theta = -M_PI / 2;

	double x_map;
	x_map = x_part + (cos(theta) * x_obs) - (sin(theta) * y_obs); //x_obs,y_obs are car coordinates. 

	double y_map;
	y_map = y_part + (sin(theta) * x_obs) + (cos(theta) * y_obs);

	cout << int(round(x_map)) << ',' << int(round(y_map)) << endl;

	return 0;
}
```

### 5. Association 
After transforming observation into map coordinate, this step associates the transformed observation with the nearest landmark. 

![Association](https://github.com/SeokLeeUS/seokleeus.github.io/raw/master/_images/_Particle_filter/association.png)

[Reference](http://ais.informatik.uni-freiburg.de/teaching/ws09/robotics2/pdfs/rob2-11-dataassociation.pdf)

### 6. Particle weight

![Particle_weight](https://github.com/SeokLeeUS/seokleeus.github.io/raw/master/_images/_Particle_filter/particle_weights_solution.png)

```
#define _USE_MATH_DEFINES
#include "multivariate_gaussian_particle_weight.h"
#include <cmath>

double multiv_prob(double sig_x, double sig_y, double x_obs, double y_obs, double mu_x, double mu_y) {

	double gauss_norm;
	gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);

	double exponent;
	exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2))) + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));

	double weight;
	weight = gauss_norm * exp(-exponent);

	return weight;
}
```
### 7. Resampling:

1. Resampling wheel. what is it?
[Video reference #1](https://www.youtube.com/watch?v=wNQVo6uOgYA)
2. How to implement in python?
[Video reference #2](https://www.youtube.com/watch?v=aHLslaWO-AQ)
3. How to implement in c++?

```
void ParticleFilter::resample() {

	/**   TODO: Resample particles with replacement with probability proportional
		* to their weight.
		* NOTE : You may find std::discrete_distribution helpful here.
		* http ://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
		* https://robotics.stackexchange.com/questions/479/particle-filters-how-to-do-resampling
		* http://cecas.clemson.edu/~ahoover/ece854/lecture-notes/lecture-pf.pdf
		* https://www.youtube.com/watch?v=aHLslaWO-AQ
		*

	*/
	
	for (int i = 0; i < particles.size(); ++i)
	{
		cout << "weight:" << particles[i].weight << endl;
	}
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count(); // creating seed
	std::default_random_engine gen_resample(seed);
	std::uniform_real_distribution<> ranN(0, 1);

	int index;
	int N = particles.size();
	index = int(ranN(gen_resample) * N);
	double beta = 0.0;
	double mw = 0.0;
	std::vector<Particle> p3;

	for (int j = 0; j < N; j++)
	{
		mw = max(particles[j].weight, mw);

	}

	double max_weight = mw; //obtain maximum weight

	for (int i = 0; i < N; i++)
	{
		beta = beta + ranN(gen_resample) * 2 * max_weight;
		while (particles[index].weight < beta) {

			beta = beta - particles[index].weight;
			index = (index + 1) % N;
		}
		p3.push_back(particles[index]);
	}
	particles = p3;
}
```

### 8. Result:

[![video result](https://img.youtube.com/vi/phv0iByFB4E/hqdefault.jpg)](https://youtu.be/phv0iByFB4E) 

### 9. Reference:

- All figures shown here is taken from Udacity's particle filter course material 
