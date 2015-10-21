/**
 * Computer Modelling, Exercise 3: Particle1D class example.
 *
 * @author A. Hermann
 * @author A. R. Turner
 * @version "10/2015"
 *
 */
public class Particle1D {

    /* ******************************************
     * Properties
     ********************************************/
    private double mass;
    private double position;
    private double velocity;

    // Setters and Getters
    
    /** Get the position of a particle.
     *
     * @return a double representing the position.
     */
    public double getPosition() { return position; }

    /** Get the velocity of a particle.
     *
     * @return a double representing the velocity.
     */
    public double getVelocity() { return velocity; }

    /** Get the mass of a particle.
     *
     * @return a double representing the mass.
     */
    public double getMass()     { return mass; }

    /** Set the position of a particle.
     *
     * @param p a double representing the position.
     */
    public void setPosition(double p) { this.position = p; }
    
    /** Set the velocity of a particle.
     *
     * @param v a double representing the velocity.
     */
    public void setVelocity(double v) { this.velocity = v; }
    
    /** Set the mass of a particle.
     *
     * @param m a double representing the mass.
     */
    public void setMass(double m)     { this.mass = m; }

    /* ******************************************
     * Constructors
     ********************************************/
    
    /** Default constructor. Sets all properties to "not a number" 
     * to indicate that they are uninitialised.
     */
    public Particle1D() {
        this.setMass(Double.NaN);
        this.setPosition(Double.NaN);
        this.setVelocity(Double.NaN);
    }

    /** Explicit constructor. Constructs a new Particle1D with
     * explicitly given position, velocity, and mass.
     *
     * @param m a double that defines the mass.
     * @param p a double that defines the position.
     * @param v a double that defines the velocity.
     */
    public Particle1D(double m, double p, double v) {
        this.setMass(m);
        this.setPosition(p);
        this.setVelocity(v);
    }

    /* ******************************************
     * toString Method
     ********************************************/
    
    /** Returns a String representation of Particle1D.
     * Used to print a Particle1D instance using the "%s"
     * format identifier.
     */
    public String toString() {
        return "{x=" + this.getPosition() + ",v=" + this.getVelocity() + 
	       "m=" + this.getMass() + "}";
    }

    /* ******************************************
     * Instance Methods
     ********************************************/
    
    /** Returns the kinetic energy of a Particle1D,
     * calculated as 1/2*m*v^2.
     *
     * @return a double that is the kinetic energy.
     */
    public double kineticEnergy() { return 0.5*mass*velocity*velocity; }

    /** Time integration support: evolve the velocity
     * according to dv = f/m * dt.
     *
     * @param dt a double that is the timestep.
     * @param force a double that is the current force on the particle.
     */
    public void leapVelocity(double dt, double force) {
        velocity = velocity + dt * force / mass;
    }
    
    /** Time integration support: evolve the position
     * according to dx = v * dt.
     *
     * @param dt a double that is the timestep.
     */
    public void leapPosition(double dt) {
        position = position + velocity * dt;
    }

    /** Time integration support: evolve the position
     * according to dx = v * dt + 0.5 * a * dt**2.
     *
     * @param dt a double that is the timestep.
     * @param force a double that is the current force.
     */
    public void leapPosition(double dt, double force) {
        position = position + velocity * dt + 0.5 * force/mass * dt*dt;
    }    
    
    
};
