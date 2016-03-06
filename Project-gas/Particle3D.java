import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;

/**
 * 
 * Computer modelling project B "Liquid & gas simulations" 
 * This is a class representing a particle with a mass,
 * position and velocity in 3 dimensions. It also includes the relevant methods
 * for calculation of forces, energies and jumps in integration methods, as well
 * as the appropriate LJ methods to support minimum image convention, periodic
 * boundary conditions and operations on an array of particles
 * 
 * @author Jonas Paulavicius
 * @author Louis Headley
 *
 */

public class Particle3D {
	
	private Vector3D position;
	private Vector3D velocity;
	private double mass;
	private String label;
	
	/**
	 * @param position Vector3D, the position of the particle
	 */
	public void setPosition(Vector3D position){
		this.position = position;
	}
	
	/**
	 * @param velocity Vector3D, the velocity of the particle
	 */
	public void setVelocity(Vector3D velocity){
		this.velocity = velocity;
	}
	
	/**
	 * @param mass double, the mass of the particle
	 */
	public void setMass(double mass) {
		this.mass = mass;
	}
	
	/**
	 * @param label String, the label for the particle
	 */
	public void setLabel(String label) {
		this.label = label;
	}
	
	/**
	 * @return the Vector3D representing the position of the particle
	 */
	public Vector3D getPosition() {
		return position;
	}
	
	/**
	 * @return the Vector3D representing the velocity of the particle
	 */
	public Vector3D getVelocity() {
		return velocity;
	}
	
	/**
	 * @return the mass of the particle
	 */
	public double getMass() {
		return mass;
	}
	
	/**
	 * @return the label of the particle
	 */
	public String getLabel() {
		return label;
	}
	
	/**
	 * Creates a new particle with the given position, velocity, mass and label
	 * 
	 * @param position Vector3D
	 * @param velocity Vector3D
	 * @param mass double
	 * @param label String
	 */
	public Particle3D(Vector3D position, Vector3D velocity, double mass, String label) {
		this.position = position;
		this.velocity = velocity;
		this.mass = mass;
		this.label = label;
	}
	
	/**
	 * Creates a new particle with undefined parameters
	 * 
	 */
	public Particle3D() {
		this.position = null;
		this.velocity = null;
		this.mass = Double.NaN;
		this.label = null;
	}
	
	public String toString() {
		//using 3 decimal point format
		
		String formattedFloats = String.format("%10.3f %10.3f %10.3f", position.getX(), position.getY(), position.getZ());
	    return label + " " + formattedFloats;
		
	}
	
	/**
	 * Outputs the label and the x, y, z coordinates in the appropriate vmd
	 * format in a single line using the given FileWriter
	 * 
	 * @param outFile FileWriter, the initialised file writer for the vmd output file
	 * @throws IOException
	 */
	public void toVMD(PrintWriter outFile) throws IOException {
		outFile.println(this.toString());
		
	}
	
	
	/**
	 * Creates a new particle from a scanner which reads 7 doubles and a string
	 * in that order and creates a particle where the first 3 doubles are the
	 * x,y,z components of position, the next 3 are the components of velocity,
	 * the last double is mass and the following string is the label.
	 * All components are expected to be input
	 * 
	 * @param scanner Scanner, the scanner from which to take data for the particle
	 * @throws Exception 
	 */
	public Particle3D(Scanner scanner) throws Exception{
		//expecting exactly 7 doubles - check if true
		double[] vecComp = new double[7];
		for (int i = 0; i < 7; i++) {
			if (scanner.hasNextDouble()) {
				vecComp[i] = scanner.nextDouble();
			} else {
				//not a double here - throw error
				throw new Exception("The input format is incorrect");
			}			
		}
		
		position = new Vector3D(vecComp[0],vecComp[1],vecComp[2]);
		velocity = new Vector3D(vecComp[3],vecComp[4],vecComp[5]);
		mass = vecComp[6];
		label = scanner.next();
	}
	
	/**
	 * @return the kinetic energy of the particle - T = 0.5*m*v^2
	 */
	public double kineticEnergy() {
		return 0.5 * mass * velocity.magnitudeSquared();
	}
	
	/**
	 * Increments the velocity by timestep dt due to some given force,
	 * used in both Verlet and Euler integration 
	 * dv = (F / m) * dt
	 * v = v + dv
	 * 
	 * @param dt double, the timestep of integration
	 * @param force Vector3D, the force acting on the particle
	 */
	public void jumpVelocity(double dt, Vector3D force) {		
		Vector3D forceEffect = new Vector3D(force);
		forceEffect.scalarMultiply(dt/mass);
		velocity = Vector3D.vectorAdd(velocity, forceEffect);		
	}
	
	
	/**
	 * Increments the velocity of all particles using the jumpVelocity() method
	 * 
	 * @param dt double, the timestep of integration
	 * @param particles Particle3D[], the array holding all the particles
	 * @param forces Vector3D[], the array of forces acting on each particle
	 */
	public static void jumpVelocityArray(double dt, Particle3D[] particles, Vector3D[] forces) {
		for (int i = 0; i < particles.length; i++) {			
			particles[i].jumpVelocity(dt, forces[i]);
		}
	}
	
	
	/**
	 * Calculates the LJ force on each particle by making a copy list of the 
	 * particle array without that particular particle and looping through this list to
	 * calculate the force contribution from each particle except itself
	 * 
	 * @param particles Particle3D[], the array holding all the particles
	 * @param boxSize Vector3D, a vector holding the dimensions of the simulation box
	 * @param cutoffDistance double, the distance at which the LJ potential is cut off
	 * @return an array holding the force on each particle
	 */
	public static Vector3D[] forceArrayLJ(Particle3D[] particles, Vector3D boxSize, double cutoffDistance) {
		
		Vector3D[] forces = new Vector3D[particles.length];
		
		for (int i = 0; i < particles.length; i++) {
			Particle3D p = particles[i];
			Vector3D force = new Vector3D(0,0,0);
			ArrayList <Particle3D> temp = new ArrayList<Particle3D>(Arrays.asList(particles));
			temp.remove(p);
			for (Particle3D p2: temp) {
				force = Vector3D.vectorAdd(force, Particle3D.LJForce(p, p2, boxSize, cutoffDistance));
			}			
			forces[i] = force;
		}		
		return forces;		
	}
	
	/**
	 * Increment the position of the particle by timestep dt due to a force
	 * r = r + v*dt + 0.5 * (F / m) * dt^2
	 * 
	 * @param dt double, the timstep
	 * @param force Vector3D, the force applied
	 */
	public void jumpPosition(double dt, Vector3D force) {
		Vector3D copyVelocity = new Vector3D(velocity);
		copyVelocity.scalarMultiply(dt);
		Vector3D forceEffect = new Vector3D(force);
		forceEffect.scalarMultiply(dt*dt*0.5/mass);
		position = Vector3D.vectorAdd(position, copyVelocity);
		position = Vector3D.vectorAdd(position, forceEffect);
	}

		
	/**
	 * Increments the position of all particles through the jumpPosition() method 
	 * and takes into account Periodic Boundary Conditions (PBC) to return each 
	 * particle to the correct position in the box in case it steps out. Each direction
	 * is independent, and is treated separately, where the mod (%) operator here
	 * can give negative values.	 * 
	 * 
	 * @param dt double, the timestep of integration
	 * @param particles Particle3D[], the array holding all the particles
	 * @param forces Vector3d[], the array of forces acting on each particle
	 * @param boxSize Vector3D, a vector holding the dimensions of the simulation box
	 */
	public static void jumpPositionArrayVerlet(double dt, Particle3D[] particles, Vector3D[] forces, Vector3D boxSize) {
		
		for (int i = 0; i < particles.length; i++) {	
			
			particles[i].jumpPosition(dt, forces[i]);			
					
			double newX = particles[i].getPosition().getX() % boxSize.getX();
			if (newX < 0) {
				newX = newX + boxSize.getX();
			}
			
			double newY = particles[i].getPosition().getY() % boxSize.getY();
			if (newY < 0) {
				newY = newY + boxSize.getY();
			}
			
			double newZ = particles[i].getPosition().getZ() % boxSize.getZ();
			if (newZ < 0) {
				newZ = newZ + boxSize.getZ();
			}
			
			Vector3D new_position = new Vector3D();
			new_position.setX(newX);
			new_position.setY(newY);
			new_position.setZ(newZ);			
			particles[i].setPosition(new_position);		
			
		}		
	}
	
	
	/**
	 * Calculates the total kinetic, potential and the total energy of the LJ fluid,
	 * where the potential is calculated by going from i=1 to N and from j=i+1 to N
	 * so that each interaction is taken into account in the least number of steps.
	 * 
	 * @param particles Particle3D[], the array holding all the particles
	 * @param boxSize Vector3D, a vector holding the dimensions of the simulation box
	 * @param cutoffDistance double, the distance at which the LJ potential is cut off
	 * @return a double array containing the total kinetic and potential energy, respectively
	 */
	public static double[] energyArray(Particle3D[] particles, Vector3D boxSize, double cutoffDistance) {
		
		double totPotential = 0;
		double totKinetic = 0;		
		
		for (int i = 0; i < particles.length-1; i++) {						
			for (int j = i+1; j < particles.length; j++) {
				totPotential = totPotential + Particle3D.LJPotential(particles[i], particles[j], boxSize, cutoffDistance);
			}						
			totKinetic = totKinetic + particles[i].kineticEnergy();
		}
		
		//missing last one as i loop stops at length-2
		totKinetic = totKinetic + particles[particles.length-1].kineticEnergy();
		
		double[] energies = new double[2];
		energies[0] = totKinetic;
		energies[1] = totPotential;		
		return energies;
	}
	
		
	/**
	 * Returns the relative separation of the two particles r = r1 - r2
	 * 
	 * @param p1 Particle3D, the first particle
	 * @param p2 Particle3D, the second particle
	 * @return the relative separation between the first and the second particle.
	 */
	public static Vector3D relativeSeparation(Particle3D p1, Particle3D p2) {
		Vector3D separation = Vector3D.vectorSubtract(p1.getPosition(), p2.getPosition());
		return separation;
	}
	
	
	/**
	 * The new relative separation method used in Molecular Dynamics simulation that
	 * also includes the minimum image convention (MIC) - the particles only 
	 * interact with the closest image of their neighbour in the PBC space.
	 * 
	 * @param p1 Particle3D, the first particle
	 * @param p2 Particle3D, the second particle
	 * @param boxSize Vector3D, a vector holding the dimensions of the simulation box
	 * @return
	 */
	public static Vector3D relativeSeparationMIC(Particle3D p1, Particle3D p2, Vector3D boxSize) {
		
		Vector3D separation = relativeSeparation(p1, p2);			
		
		if (Math.abs(separation.getX())>boxSize.getX()/2) {			
			double new_x = boxSize.getX() - Math.abs(separation.getX());					
			//a sign flip occurs as the closer image will be in the opposite direction
			separation.setX(-1D * new_x);			
		}
		
		if (Math.abs(separation.getY())>boxSize.getY()/2) {			
			double new_y = boxSize.getY() - Math.abs(separation.getY());			
			separation.setY(-1D * new_y);		
		}
		
		if (Math.abs(separation.getZ())>boxSize.getZ()/2) {			
			double new_z = boxSize.getZ() - Math.abs(separation.getZ());		
			separation.setZ(-1D * new_z);
		}
		
		return separation;
	}
	
	//TESTING ONLY
	/*
	public static double minSeparation(Particle3D[] particles, Vector3D boxSize) {
		double minSeparation = 1000;
		for (int i = 0; i < particles.length-1; i++) {
			for (int j = i+1; j < particles.length; j++) {
				if (relativeSeparationMIC(particles[i], particles[j], boxSize).magnitude() < minSeparation) {
					minSeparation = relativeSeparationMIC(particles[i], particles[j], boxSize).magnitude();
				}
			}
			
		}
		return minSeparation;
	}
	*/
	
	/**
	 * Calculates the force on particle 1 due to particle 2 due to a
	 * Leonard-Jones force for idealised particles in reduced units
	 * F = 48(r^-13 - (1/2)r^-7)
	 * MIC is used together with a cutoff distance to reduce calculation
	 * 
	 * @param p1 Particle3D, the first particle
	 * @param p2 Particle3D, the second particle
	 * @param boxSize Vector3D, a vector holding the dimensions of the simulation box
	 * @param cutoffDistance double, the distance at which the LJ potential is cut off
	 * @return the force on particle p1 due to particle p2
	 */	
	public static Vector3D LJForce(Particle3D p1, Particle3D p2, Vector3D boxSize, double cutoffDistance){
		Vector3D force = Particle3D.relativeSeparationMIC(p2, p1, boxSize);
		double r = force.magnitude();
		
		if (r <= cutoffDistance) {				
			force.scalarDivide(r);				
			double magnitude = 48*(Math.pow(r, -13) - Math.pow(r, -7)/2);
			force.scalarMultiply(magnitude);			
		} else {
			force = new Vector3D(0,0,0);
		}
		
		return force;
	}
	
	/**
	 * Calculates the potential energy of particle 1 due to particle 2 given a
	 * Leonard-Jones potential for idealised particles in reduced units
	 * U = 4(r^-12 - r^-6)
	 * It uses MIC when calculating the relative separation of the two particles
	 * and includes a discontinuity correction to allow for a cutoff distance to be
	 * used in the simulation
	 * 
	 * @param p1 Particle3D, the first particle
	 * @param p2 Particle3D, the second particle
	 * @param boxSize Vector3D, a vector holding the dimensions of the simulation box
	 * @param cutoffDistance double, the distance at which the LJ potential is cut off
	 * @return the potential at particle p1 due to particle p2
	 */	
	public static double LJPotential(Particle3D p1, Particle3D p2, Vector3D boxSize, double cutoffDistance) {
		double r = Particle3D.relativeSeparationMIC(p2, p1, boxSize).magnitude();
		double potential = 0;		
		
		if (r <= cutoffDistance) {
			double discCorr = 4*(Math.pow(cutoffDistance, -12) - Math.pow(cutoffDistance, -6));
			potential = 4*(Math.pow(r, -12) - Math.pow(r, -6)) - discCorr;
		} 
		
		return potential;
	}
	
	
	/**
	 * Calculates the Radial Distribution Function (RDF) values for the gas by 
	 * sampling at discrete relative separations and stores these values in an array
	 * that is updated throughout the simulation to yield a time-average.
	 * 
	 * g(R) = (V/N)< (sum for i not equal to j) delta(Rij - R) >
	 * where g(R) is the value of the RDF at R, the <> denote a time average (an average
	 * over snapshots from the system over time), and delta is the dirac delta function,
	 * with Rij the distance between the i-th and j-th particles. The R only needs to run
	 * from 0 to cutoffdistance in discrete steps that depend on the ratio of the cutoff
	 * distance to the length of the RDF array
	 * 
	 * This method has no return, but updates the RDF array that is held in the 
	 * main method of the simulation
	 * 
	 * @param particles Particle3D[], the array holding all the particles
	 * @param boxSize Vector3D, a vector holding the dimensions of the simulation box
	 * @param deltaPeaks double[], an array holding all the RDF values
	 * @param cutoffDistance double, the distance at which the LJ potential is cut off
	 */
	public static void calculateRDF(Particle3D[] particles, Vector3D boxSize, double[] deltaPeaks, double cutoffDistance) {
				
		double lengthUnit = deltaPeaks.length / cutoffDistance;
		
		//i=0 to N-1, j=i+1 to N looping ensures the j not equal to i condition is satisfied
		for (int i = 0; i < particles.length-1; i++) {						
			for (int j = i+1; j < particles.length; j++) {
				//only put a peak at Rij = R, taking the dirac delta function
				//as 1 at Rij = R and 0 elsewhere
				double deltaPeak = Particle3D.relativeSeparationMIC(particles[i], particles[j], boxSize).magnitude() * lengthUnit;
				int arrayPos = (int) Math.round(deltaPeak);
				if (arrayPos < deltaPeaks.length) {
					//averaging is done automatically by updating the same RDF array 
					//2 because the same applies for both particles
					deltaPeaks[arrayPos] = deltaPeaks[arrayPos] + 2;
				}				
			}			
		}		
	}
	
	/**
	 * Prints the RDF values to the file rdf.out that is in the same folder as all
	 * of the simulation files. 
	 * 
	 * @param deltaPeaks double[], an array holding all the RDF values
	 * @param boxSize Vector3D, a vector holding the dimensions of the simulation box
	 * @param particleNo int, the number of particles in the simulation
	 * @param cutoffDistance double, the distance at which the LJ potential is cut off
	 * @throws IOException
	 */
	public static void printRDF(double[] deltaPeaks, Vector3D boxSize, int particleNo, double cutoffDistance) throws IOException {
		
		double volume = boxSize.getX()*boxSize.getY()*boxSize.getZ();
		double prefactor = volume/particleNo;
		
		double lengthUnit = cutoffDistance / deltaPeaks.length;
		PrintWriter rdfOutput = new PrintWriter(new FileWriter("rdf.out"));
		for (int i = 0; i < deltaPeaks.length; i++) {
			rdfOutput.printf("%10.5f %10.5f\n", i*lengthUnit, deltaPeaks[i]*prefactor);
		}
		rdfOutput.close();		
	}
	

	public static void main(String[] args) {			
	}
}