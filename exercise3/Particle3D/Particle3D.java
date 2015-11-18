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
	
	public String toString(){
		//using 3 decimal point format
		
		String formattedFloats = String.format("%10.3f %10.3f %10.3f", position.getX(), position.getY(), position.getZ());
	    return label + " " + formattedFloats;
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
	
	
	
	/**
	 * Returns the gravitational force on particle 1 due to particle 2 with
	 * the constant G taken to be 1
	 * F = (m1 * m2) / r^2 in the unit r direction
	 * 
	 * @param p1 Particle3D, the first particle
	 * @param p2 Particle3D, the second particle
	 * @return the force due to gravity on particle 1 due to particle 2
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
	 * @return the potential at particle p1 due to particle p2
	 */
	public static double graviPotential(Particle3D p1, Particle3D p2) {
		return -1 * p1.getMass() * p2.getMass() / Particle3D.relativeSeparation(p2, p1).magnitude();
	}
	

	public static void main(String[] args) {	
		
	}
}