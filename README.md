# 3-DOF Rocket Simulation
MATLAB simulation of a rocket flight for the PSP flight dynamics team 
## Overview
This simulation models the flight of a rocket in 2D space, accounting 
for thrust, drag, and gravitational forces, but not wind. The program calculates the 
complete trajectory from launch to impact and provides visualizations of the 
flight path, velocities, and forces over time.
## Features
- Physics: Implements thrust, drag, and gravitational forces
- Numerical Integration: Uses Euler method with second-order position
correction (goes to acceleration term)
- Multiple plots showing trajectory, velocities, and forces
- Animated Flight Path: Real-time animation with rotating rocket model
- Calculates apogee, range, flight time, and maximum speed
## Simulation Parameters
### Default Configuration
- Thrust: 2000 N (burn time is 5 seconds)
- Mass: 40 kg
- Rocket Diameter: 0.1 m
- Drag Coefficient: 0.75
- Launch Angle: 85° from horizontal (I chose this angle because no angle
was specified in the onboarding assignment, and 85 degrees is almost
vertical, so the flight time can be longer, but the rocket is still angled 
enough to show some interesting physics)
- Time Step: 0.01 seconds (just to get lots of data points)
### Physical Constants
- Gravity: 9.81 m/s² (this is standard)
- Air Density: 1.225 kg/m³ (at sea level)
## How to Run
1. Install MATLAB 2025a
2. Save the code as a `.m` file (ex. `rocket_simulation.m`)
3. Run the script in MATLAB:
```matlab rocket_simulation```
## Output
### Console Output
The simulation displays:
- Initial conditions
- Key flight parameters (apogee, range, flight time, max speed)
- Animation status updates
### Generated Plots
1. Figure 1 - Trajectory Plot: plots x by y -> shows 2D flight path with marked apogee and launch point
2. Figure 2 - Velocity Components: Three subplots showing horizontal velocity, vertical velocity, and total speed vs time
3. Figure 3 - Force Analysis: Thrust, drag, and weight forces throughout the flight
4. Figure 4 - Animated Flight: Real-time animation with rotating rocket and trajectory trail
## Physics Implementation
### Forces Modeled
- Thrust Force: F_thrust = 2000 N (during the 5 second burn time only)
- Drag Force: F_drag = 0.5 × ρ × v² × A × C_d (this equation was given,
and drag always opposes motion, so I made sure to include that in the code)
- Gravitational Force: F_gravity = m × g (constant downward)
### Coordinate System
- X-axis: Horizontal distance (positive = moving to the right)
- Y-axis: Vertical altitude (positive = upward)  
- Origin: Launch point (0, 0)
### Integration Method
The simulation uses the Euler method with second-order correction:
x(i) = x(i-1) + v(i-1)*dt + 0.5*a*dt² (this is a kinematics equation x = xo + vot +1/2 at^2
v(i) = v(i-1) + a*dt (this is a kinematics equation vf = vo+at) 
i-1 means the previous data point
## Customization
### Modifying Parameters
To change simulation parameters, edit the values in the "Given Parameters" section:
thrust = 2000;        % Thrust force in Newtons
mass = 40;           % Rocket mass in kg
burnTime = 5;        % Burn duration in seconds
diameter = 0.1;      % Rocket diameter in meters
dragCoeff = 0.75;    % Drag coefficient (dimensionless)
launchAngle = 85;    % Launch angle in degrees (it would be cool to play
around with the angle)
I want to add that to the rotating code that was provided, I changed the
dimensions of the rocket because given the scales of my plots, it didn't
make sense for the rocket to be so small because I couldn't see it, so I
made it a lot bigger)
## Added Key Features 
### Automatic Ground Detection
The simulation automatically terminates when the rocket hits the ground (y ≤ 0), ensuring realistic flight times.
### Dynamic Drag Calculation
Drag force is calculated based on instantaneous velocity and opposes the direction of motion.
### Rocket Orientation
The animated rocket rotates to align with its velocity vector, providing realistic visual representation.
### Performance Optimization
The animation skips frames (every 5th point) for a smoother visual.
## Results Interpretation
### Typical Output Values
- Apogee: Maximum altitude reached
- Range: Horizontal distance traveled
- Flight Time: Total time from launch to impact
- Max Speed: Peak velocity during flight
### Understanding the Plots
- Velocity plots show the transition from powered to unpowered flight
- Force plots show the impact of each force throughout the flight
- Trajectory plot shows the complete flight path with key points marked
## Limitations
- 2D simulation rather than 3D
- Constant air density (not realistic)
- Rigid body assumption (no structural dynamics)
- Point mass approximation (not realistic)
- Constant drag coefficient
## Future Enhancements
Potential improvements could include:
- Variable air density with altitude
- Wind effects
- Multi-stage rockets
- 3D trajectory capability
## Dependencies
- MATLAB (2025a version)
- No additional toolboxes required
- Built-in functions: `polyshape`, `rotate`, `translate`
## Troubleshooting
### Common Issues
- Earlier, the animation was not visible, I fixed this by increasing the 
rocket size parameters `w` and `h` in the given code.
- Simulation runs too long: Increase `dt`%
- Memory issues: Reduce simulation resolution or duration
### Performance Tips
- Use larger time steps for faster simulation (less accurate)
- Reduce animation frame rate by increasing the skip value in the animation loop
- Close unnecessary figure windows to save memory
- Use larger time steps for faster simulation (less accurate)
- Reduce animation frame rate by increasing the skip value in the animation loop
- Close unnecessary figure windows to save memory
