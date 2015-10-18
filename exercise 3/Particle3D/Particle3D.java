import java.util.Scanner;


public class Particle3D {
	
	private Vector3D position;
	private Vector3D velocity;
	private double mass;
	private String label;
	
	public void setPosition(Vector3D position){
		this.position = position;
	}
	
	public void setVelocty(Vector3D velocity){
		this.velocity = velocity;
	}
	
	public void setMass(double mass) {
		this.mass = mass;
	}
	
	public void setLabel(String label) {
		this.label = label;
	}
	
	public Vector3D getPosition() {
		return position;
	}
	
	public Vector3D getVelocity() {
		return velocity;
	}
	
	public double getMass() {
		return mass;
	}
	
	public String getLabel() {
		return label;
	}
	
	public Particle3D(Vector3D position, Vector3D velocity, double mass, String label) {
		this.position = position;
		this.velocity = velocity;
		this.mass = mass;
		this.label = label;
	}
	
	public Particle3D() {
		this.position = null;
		this.velocity = null;
		this.mass = Double.NaN;
		this.label = null;
	}
	
	public String toString(){
		return "<" + label+"> <"+position.getX()+"> <"+position.getY()+"> <"+position.getZ()+">";
	}
	
	//first position then velocity, expect all components to be input
	public Particle3D(Scanner scanner){
		//expecting exactly 7 doubles - check if true
		double[] vecComp = new double[7];
		for (int i = 0; i < 7; i++) {
			if (scanner.hasNextDouble()) {
				vecComp[i] = scanner.nextDouble();
			} else {
				//not a double here - throw error
			}			
		}
		
		position = new Vector3D(vecComp[0],vecComp[1],vecComp[2]);
		velocity = new Vector3D(vecComp[3],vecComp[4],vecComp[5]);
		mass = vecComp[6];
		label = scanner.next();
	}
	
	public double kineticEnergy() {
		return 0.5 * mass * velocity.magnitudeSquared();
	}
	
	public void jumpVelocity(double dt, Vector3D force) {
		Vector3D forceEffect = new Vector3D(force);
		forceEffect.scalarMultiply(dt/mass);
		velocity = Vector3D.vectorAdd(velocity, forceEffect);		
	}
	
	public void jumpPosition(double dt) {
		Vector3D copyVelocity = new Vector3D(velocity);
		copyVelocity.scalarMultiply(dt);
		position = Vector3D.vectorAdd(position, copyVelocity);
	}
	
	public void jumpPosition(double dt, Vector3D force) {
		Vector3D copyVelocity = new Vector3D(velocity);
		copyVelocity.scalarMultiply(dt);
		Vector3D forceEffect = new Vector3D(force);
		forceEffect.scalarMultiply(dt*dt*0.5/mass);
		velocity = Vector3D.vectorAdd(velocity, copyVelocity);
		velocity = Vector3D.vectorAdd(velocity, forceEffect);
	}
	
	public static Vector3D relativeSeparation(Particle3D p1, Particle3D p2) {
		Vector3D separation = Vector3D.vectorSubtract(p1.getPosition(), p2.getPosition());
		return separation;
	}
	
	//force on particle 1 due to particle 2
	public static Vector3D graviForce(Particle3D p1, Particle3D p2) {
		Vector3D force = Particle3D.relativeSeparation(p2, p1);
		force.scalarDivide(Math.pow(force.magnitude(), 3)*p1.getMass()*p2.getMass());
		return force;
	}
	
	
	//potential at particle p1 due to particle p2
	public static double graviPotential(Particle3D p1, Particle3D p2) {
		return -1 * p1.getMass() * p2.getMass() * Particle3D.relativeSeparation(p2, p1).magnitude();
	}
	

	public static void main(String[] args) {
		/*Vector3D position = new Vector3D(0,0,0);
		Vector3D velocity = new Vector3D(1,1,1);
		Particle3D p = new Particle3D(position, velocity, 1, "hi");
		System.out.println(p);
		p.jumpPosition(0.5);
		System.out.println(p);	*/	
	}
	
	

}
