import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;


/**
 * Computer modelling exercise 3 - the implementation of the symplectic
 * verlet integration algorithm for a particle orbiting a fixed mass
 * 
 * @author Louis Headley
 * @author Jonas Paulavicius
 *
 */
public class Particle3DVerlet {

	public static void main(String[] args) throws IOException {
		
		//the output file is the 1st commandline argument
		String outFile = args[0];
		PrintWriter output = new PrintWriter(new FileWriter(outFile));
		
		//the integer number of steps is the 2nd argument
		int numstep = Integer.valueOf(args[1]);
		
	    //the output file for total energy
		PrintWriter outputEnergy = new PrintWriter(new FileWriter("VerletEnergy.out"));
		
		//initialize the small particle and the fixed mass
		Vector3D particlePos = new Vector3D(1,0,0);
		Vector3D particleVel = new Vector3D(0,1,0);		
		Particle3D particle = new Particle3D(particlePos, particleVel, 0.0001, "orbiting particle");
		
		Vector3D massPos = new Vector3D(0,0,0);
		Vector3D massVel = new Vector3D(0,0,0);
		Particle3D mass = new Particle3D(massPos, massVel, 1, "fixed mass");
		
		//calculate the initial force on the particle
		Vector3D force = Particle3D.graviForce(particle, mass);
		
		
		//algorithm parameters		
		double dt = 0.1;		
		double t = 0;
		double energy;
		
		//output the initial data at t = 0
		output.printf("%s %10.2f\n", particle, t);	
		
		for (int i = 0; i < numstep; i++) {
			
			//jump the position according to the force and timestep
			particle.jumpPosition(dt, force);
			
			//update the force for the new position
			Vector3D forceNew = Particle3D.graviForce(particle, mass);
			
			//get average of the two forces	(old and new)		
			Vector3D averageForce = Vector3D.vectorAdd(force, forceNew);
			averageForce.scalarDivide(2);
			
			//jump the velocity by the average force
			particle.jumpVelocity(dt, averageForce);
			
			/*
			 * set the force to the force at the new position, copy to avoid
			 * referencing problems
			 */
			
			//set the force in the next loop to the new force
			force = new Vector3D(forceNew);
			
			//increment time by the timestep
			t = t + dt;
			
			//print the new data to the output file
			output.printf("%s %10.2f\n", particle, t);
			
			//calculate and output the total energy
			energy = particle.kineticEnergy() + Particle3D.graviPotential(particle, mass);
			outputEnergy.printf("%10.2f %10.7f \n", t, energy);
			
		}
		
		//close the output files
		output.close();
		outputEnergy.close();

	}

}
