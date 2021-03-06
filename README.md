# Extended Kalman Filters Project

[![Udacity - Self-Driving Car NanoDegree](https://s3.amazonaws.com/udacity-sdc/github/shield-carnd.svg)](http://www.udacity.com/drive)

**Shameless plug:** While you're here, you may want to have a look at my
[libfixkalman](https://github.com/sunsided/libfixkalman) library for
fixed-point Kalman filters in embedded systems.

This project utilizes both a regular and extended Kalman filter to estimate the state 
of a moving object of interest with noisy [LIDAR](https://en.wikipedia.org/wiki/Lidar) and 
[RADAR](https://en.wikipedia.org/wiki/Radar) measurements. 
Passing the project requires obtaining [RMSE](https://en.wikipedia.org/wiki/Root-mean-square_deviation) values that are lower than the tolerance 
outlined in the project rubric. 

This project involves the Term 2 Simulator which can be downloaded 
[here](https://github.com/udacity/self-driving-car-sim/releases).

### WebSocket communication

This repository includes two files that can be used to set up and install 
[uWebSocketIO](https://github.com/uWebSockets/uWebSockets) for either Linux or Mac 
systems. For windows you can use either Docker, VMware, or even 
[Windows 10 Bash on Ubuntu](https://www.howtogeek.com/249966/how-to-install-and-use-the-linux-bash-shell-on-windows-10/) 
to install uWebSocketIO. Please see [this concept in the classroom](https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/0949fca6-b379-42af-a919-ee50aa304e6a/lessons/f758c44c-5e40-4e01-93b5-1a82aa4e044f/concepts/16cf4a78-4fc7-49e1-8621-3450ca938b77) 
for the required version and installation scripts.

Once the install for uWebSocketIO is complete, the main program can be built and run 
by doing the following from the project top directory.

1. `mkdir build`
2. `cd build`
3. `cmake ..`
4. `make`
5. `./ExtendedKF`

Tips for setting up your environment can be found [here](https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/0949fca6-b379-42af-a919-ee50aa304e6a/lessons/f758c44c-5e40-4e01-93b5-1a82aa4e044f/concepts/23d376c7-0195-4276-bdf0-e02f1f3c665d).
Note that the files that need to be completed to accomplish the project are 
`src/FusionEKF.cpp`, `src/FusionEKF.h`, `kalman_filter.cpp`, `kalman_filter.h`, `tools.cpp`, and `tools.h`.
The `main.cpp` has already been filled out, but feel free to modify it.

Here is the main protocol that `main.cpp` uses for uWebSocketIO in communicating with the simulator:

#### Input values provided by the simulator to the c++ program

- `["sensor_measurement"]` => the measurement that the simulator observed (either lidar or radar)

#### Output values provided by the c++ program to the simulator

- `["estimate_x"]` <= kalman filter estimated position x
- `["estimate_y"]` <= kalman filter estimated position y
- `["rmse_x"]`
- `["rmse_y"]`
- `["rmse_vx"]`
- `["rmse_vy"]`

---

## Other Important Dependencies

* cmake >= 3.5
  * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1 (Linux, Mac), 3.81 (Windows)
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools](https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make` 
   * On windows, you may need to run: `cmake .. -G "Unix Makefiles" && make`
4. Run it: `./ExtendedKF `

## Generating Additional Data

This is optional!

If you'd like to generate your own radar and lidar data, see the
[utilities repo](https://github.com/udacity/CarND-Mercedes-SF-Utilities) for
Matlab scripts that can generate additional data.
