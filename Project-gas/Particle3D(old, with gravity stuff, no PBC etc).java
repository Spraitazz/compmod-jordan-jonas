import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;


/**
 * 
 * Computer modelling exercise 3 - a class representing a particle with a mass,
 * position and velocity in 3 dimensions. It also includes the relevant methods
 * for calculation of forces, energies and jumps in integration methods
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
		//copy necessary as Vector3D force is passed by reference
		Vector3D forceEffect = new Vector3D(force);
		forceEffect.scalarMultiply(dt/mass);
		velocity = Vector3D.vectorAdd(velocity, forceEffect);		
	}
	
	
	public static void jumpVelocityArray(double dt, Particle3D[] particles, Vector3D[] forces) {
		for (int i = 0; i < particles.length; i++) {			
			particles[i].jumpVelocity(dt, forces[i]);
		}
	}
	
	//separately so as to differentiate between old and new forces (no star wars puns)
	//could use alternative algorithm here to achieve ~ 2x speedup, but
    //both still O(n^2)
	public static Vector3D[] forceArrayLJ(Particle3D[] particles, Vector3D boxSize) {
		Vector3D[] forces = new Vector3D[particles.length];
		for (int i = 0; i < particles.length; i++) {
			Particle3D p = particles[i];
			//start with a 0 vector
			Vector3D force = new Vector3D(0,0,0);
			//temporary arraylist with all particles except p
			ArrayList <Particle3D> temp = new ArrayList<Particle3D>(Arrays.asList(particles));
			temp.remove(p);
			//add up the forces on p due to all particles except p
			for (Particle3D p2: temp) {
				force = Vector3D.vectorAdd(force, Particle3D.LJForce(p, p2, boxSize));
			}
			forces[i] = force;
		}
		
		return forces;		
	}

	
	/**
	 * Increment the position by timestep dt given the current velocity
	 * r = r + v*dt
	 * 
	 * @param dt double, the timestep
	 */
	public void jumpPosition(double dt) {
		Vector3D copyVelocity = new Vector3D(velocity);
		copyVelocity.scalarMultiply(dt);
		position = Vector3D.vectorAdd(position, copyVelocity);
	}
	
	//take into account PBC
	public static void jumpPositionArrayVerlet(double dt, Particle3D[] particles, Vector3D[] forces, Vector3D boxSize) {
		
		for (int i = 0; i < particles.length; i++) {			
			particles[i].jumpPosition(dt, forces[i]);
			
			//if stepped out of bounds, updates to correct position, else nothing changes
			//each direction independent
			//mod (%) well-defined here - returns negative if first number negative, else positive
			
			Vector3D new_position = new Vector3D();
			
			//(x % n) - (x < 0 ? n : 0); from http://stackoverflow.com/questions/4403542/how-does-java-do-modulus-calculations-with-negative-numbers
			//should be a plus though?
			
			
			double newX = (particles[i].getPosition().getX() % boxSize.getX()) + (particles[i].getPosition().getX() % boxSize.getX() < 0 ? boxSize.getX() : 0);
			double newY = (particles[i].getPosition().getY() % boxSize.getY()) + (particles[i].getPosition().getY() % boxSize.getY() < 0 ? boxSize.getY() : 0);
			double newZ = (particles[i].getPosition().getZ() % boxSize.getZ()) + (particles[i].getPosition().getZ() % boxSize.getZ() < 0 ? boxSize.getZ() : 0);

			new_position.setX(newX);
			new_position.setY(newY);
			new_position.setZ(newZ);
			
			particles[i].setPosition(new_position);
		
		}
		
	}
	
	//[0] - kinetic [1] - potential [2] - total
	public static double[] energyArray(Particle3D[] particles, Vector3D boxSize) {
		
		double totalPotential = 0;
		double totalKinetic = 0;
		
		for (Particle3D p: particles) {			
			double thisPotential = 0;
			//temporary arraylist with all particles except p
			ArrayList <Particle3D> temp = new ArrayList<Particle3D>(Arrays.asList(particles));
			temp.remove(p);
			//add up the potentials at p due to all particles except p
			for (Particle3D p2: temp) {
				thisPotential = thisPotential + Particle3D.LJPotential(p, p2, boxSize);
			}
			totalPotential = totalPotential + thisPotential;
			totalKinetic = totalKinetic + p.kineticEnergy();
		}
		
		double[] energies = new double[3];
		energies[0] = totalKinetic;
		energies[1] = totalPotential;
		energies[2] = totalKinetic + totalPotential;
		
		return energies;
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
	
	//this includes the minimum image convention MIC
	public static Vector3D relativeSeparationMIC(Particle3D p1, Particle3D p2, Vector3D boxSize) {
		//same, but scaled down if necessary
		Vector3D separation = relativeSeparation(p1, p2);
		
		//absolute value to check, sign flips if moved (i think)
		//each direction again independent
		if (Math.abs(separation.getX())>boxSize.getX()/2) {
			
			double new_x = 0;
			
			
			//this should work in finding the smaller distance to image in the particular direction
			if (p1.getPosition().getX() > p2.getPosition().getX()) {
				new_x = boxSize.getX() - (p1.getPosition().getX() - p2.getPosition().getX());
			} else {
				new_x = boxSize.getX() - (p2.getPosition().getX() - p1.getPosition().getX());
			}
			
			//sign flip, make a sketch to see why
			separation.setX(-1 * Math.signum(separation.getX())*new_x);
			
		}
		
		if (Math.abs(separation.getY())>boxSize.getY()/2) {
			
			double new_y = 0;		
			
			if (p1.getPosition().getY() > p2.getPosition().getY()) {
				new_y = boxSize.getY() - (p1.getPosition().getY() - p2.getPosition().getY());
			} else {
				new_y = boxSize.getY() - (p2.getPosition().getY() - p1.getPosition().getY());
			}
			
			separation.setY(-1 * Math.signum(separation.getY())*new_y);
			
		}
		
		if (Math.abs(separation.getZ())>boxSize.getZ()/2) {
			
			double new_z = 0;
			
			if (p1.getPosition().getZ() > p2.getPosition().getZ()) {
				new_z = boxSize.getZ() - (p1.getPosition().getZ() - p2.getPosition().getZ());
			} else {
				new_z = boxSize.getZ() - (p2.getPosition().getZ() - p1.getPosition().getZ());
			}
			
			separation.setZ(-1 * Math.signum(separation.getZ())*new_z);
		}
		
		
		
		
		return separation;
	}
	
	
	
	/**
	 * Returns the gravitational force on particle 1 due to particle 2 with
	 * the constant G taken to be 1
	 * F = (m1 * m2) / r^2 in the unit r direction
	 * 
	 * @param p1 Particle3D, the first particle
	 * @param p2 Particle3D, the second particle
	 * @return the force due to gravity on particle p1 due to particle p2 due 
	 * to a Newtonian gravitational force
	 */
	public static Vector3D graviForce(Particle3D p1, Particle3D p2) {
		Vector3D force = Particle3D.relativeSeparation(p2, p1);
		force.scalarDivide(Math.pow(force.magnitude(), 3));
		force.scalarMultiply(p1.getMass()*p2.getMass());
		return force;
	}
	
	
	/**
	 * Calculates the gravitational potential energy of p1 due to p2,
	 * V = -(m1 * m2) / r where G is again 1
	 * 
	 * @param p1 Particle3D, the first particle
	 * @param p2 Particle3D, the second particle
	 * @return the potential at particle p1 due to particle p2 due to a
	 * Newtonian gravitational force
	 */
	public static double graviPotential(Particle3D p1, Particle3D p2) {
		return -1 * p1.getMass() * p2.getMass() / Particle3D.relativeSeparation(p2, p1).magnitude();
	}
	
	
	/**
	 * Calculates the force on particle 1 due to particle 2 due to a
	 * Leonard-Jones potential for idealised particles in reduced units
	 * F = 48(r^-13 - (1/2)r^-7)
	 * 
	 * @param p1 Particle3D, the first particle
	 * @param p2 Particle3D, the second particle
	 * @return the force on particle p1 due to particle p2
	 */
	//CHANGED TO INCLUDE MINIMUM IMAGE CONVENTION
	//cut off distance used
	public static Vector3D LJForce(Particle3D p1, Particle3D p2, Vector3D boxSize, double cutoffDistance){
		Vector3D force = Particle3D.relativeSeparationMIC(p2, p1, boxSize);
		double r = force.magnitude();
		force.scalarDivide(r);
		double magnitude = 48*(Math.pow(r, -13)-Math.pow(r, -7)/2);
		force.scalarMultiply(magnitude);
		return force;
	}
	
	/**
	 * Calculates the potential energy of particle 1 due to particle 2 given a
	 * Leonard-Jones potential for idealised particles in reduced units
	 * U = 4(r^-12 - r^-6)
	 * 
	 * @param p1 Particle3D, the first particle
	 * @param p2 Particle3D, the second particle
	 * @return the potential at particle p1 due to particle p2
	 */
	//changed to include MIC
	//use cutoff distance
	public static double LJPotential(Particle3D p1, Particle3D p2, Vector3D boxSize, double cutoffDistance) {
		double r = Particle3D.relativeSeparationMIC(p2, p1, boxSize).magnitude();
		return 4*(Math.pow(r, -12) - Math.pow(r, -6));
	}
	
	//sample by units of 1/x of lengthscale
	//r from 0 to 3 length scales
	public static void calculateRDF(Particle3D[] particles, Vector3D boxSize, double[] deltaPeaks) {
		
		double volume = boxSize.getX()*boxSize.getY()*boxSize.getZ();
		//may not work as a good length scale
		double lengthScale = Math.pow(volume/particles.length, 1/3);
		//150 steps of this unit - make checks to determine where the delta function peaks
		double lengthUnit = lengthScale*3 / deltaPeaks.length;
		//double prefactor = volume / particles.length;
		for (int i = 0; i < particles.length; i++) {
			Particle3D p = particles[i];
			
			//temporary arraylist with all particles except p
			ArrayList <Particle3D> temp = new ArrayList<Particle3D>(Arrays.asList(particles));
			temp.remove(p);
			
			for (Particle3D p2: temp) {
				double deltaPeak = Particle3D.relativeSeparationMIC(p, p2, boxSize).magnitude()/lengthUnit;
				int arrayPos = (int) Math.round(deltaPeak);
				if (arrayPos < deltaPeaks.length) {
					//averaging by adding
					deltaPeaks[arrayPos] = deltaPeaks[arrayPos] + 1;
				}
				
				
			}
			
		}
		
	}
	
	public static void printRDF(double[] deltaPeaks, Vector3D boxSize, int particleNo) throws IOException {
		
		double volume = boxSize.getX()*boxSize.getY()*boxSize.getZ();
		double prefactor = volume/particleNo;
		double lengthScale = Math.pow(volume/particleNo, 1/3);
		//150 steps of this unit - make checks to determine where the delta function peaks
		double lengthUnit = lengthScale*3 / deltaPeaks.length;
		PrintWriter rdfOutput = new PrintWriter(new FileWriter("rdf.out"));
		for (int i = 0; i < deltaPeaks.length; i++) {
			rdfOutput.printf("%10.5f %10.5f\n", i*lengthUnit, deltaPeaks[i]*prefactor);
		}
		rdfOutput.close();
		
		
	}
	

	public static void main(String[] args) {	
		
	}
}