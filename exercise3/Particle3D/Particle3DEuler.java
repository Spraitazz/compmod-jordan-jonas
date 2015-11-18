import java.io.*;

/**
 * Computer modelling exercise 3 - the euler integration method file output
 * program.
 * 
 * @author Louis Headley
 * @author Jonas Paulavicius
 *
 */
public class Particle3DEuler {

	public static void main(String[] args) throws IOException {
		
		//the output file for the t,x,y,z data is the 1st commandline argument
		String outFile = args[0];
		
		//the integer number of steps is the 2nd argument
		int numstep = Integer.valueOf(args[1]);
		
		//the printwriter for data output
		PrintWriter output = new PrintWriter(new FileWriter(outFile));
		
		//a separate printwriter for the total energy file
		PrintWriter outputEnergy = new PrintWriter(new FileWriter("EulerEnergy.out"));
		
		//initialize the first particle, the orbiting particle
		Vector3D particlePos = new Vector3D(1,0,0);
		Vector3D particleVel = new Vector3D(0,1,0);		
		Particle3D particle = new Particle3D(particlePos, particleVel, 0.0001, "orbiting particle");
		
		//initialize the second particle, the fixed mass
		Vector3D massPos = new Vector3D(0,0,0);
		Vector3D massVel = new Vector3D(0,0,0);
		Particle3D mass = new Particle3D(massPos, massVel, 1, "fixed mass");
		
		//initial force on the particle due to the mass
		Vector3D force = Particle3D.graviForce(particle, mass);
		
		//algorithm parameters			
		double dt = 0.1;		
		double t = 0;
	    double energy;					   
		
	    //the data for the particle at t = 0
		output.printf("%s %10.2f\n", particle, t);
		
		//begin euler integration
		for (int i = 0; i < numstep; i++) {
			
			/*
			 * this follows the euler algorithm - increment position due to
			 * current velocity, then increment the force at the new position,
			 * increase the velocity due to the new force and continue
			 */
			
			particle.jumpPosition(dt);
			
			force = Particle3D.graviForce(particle, mass);
			
			particle.jumpVelocity(dt, force);
			
			t = t + dt;
			
			//output the data at t = t + dt
			output.printf("%s %10.2f\n", particle, t);
			
			/*
			 * the total energy of the particle - kinetic energy + gravitational
			 * potential energy due to the mass
			*/ 
			energy = particle.kineticEnergy() + Particle3D.graviPotential(particle, mass);
			
			//output the energy into the relevant file
			outputEnergy.printf("%10.2f %10.7f \n", t, energy);
			
		}
		
		//close the output files
		output.close();
		outputEnergy.close();
		
	}

}
