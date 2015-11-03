import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;


public class Particle3DVerlet {

	public static void main(String[] args) throws IOException {
		
		String outFile = args[0];
        PrintWriter output = new PrintWriter(new FileWriter(outFile));
		
		final double G = 1;
		
		Vector3D particlePos = new Vector3D(1,0,0);
		Vector3D particleVel = new Vector3D(0,1,0);		
		Particle3D particle = new Particle3D(particlePos, particleVel, 0.0001, "orbiting particle");
		
		Vector3D massPos = new Vector3D(0,0,0);
		Vector3D massVel = new Vector3D(0,0,0);
		Particle3D mass = new Particle3D(massPos, massVel, 1, "fixed mass");
		
		Vector3D force = Particle3D.graviForce(particle, mass);
		
		int numstep = 10;
		
		double dt = 0.1;
		
		double t = 0;
		
		output.printf("%10.5f %s\n", t, particle);
		
		for (int i = 0; i < numstep; i++) {
			
			//jump position according to the force and timestep
			particle.jumpPosition(dt, force);
			
			//update force for new position
			Vector3D forceNew = Particle3D.graviForce(particle, mass);
			
			//get average of the two forces
			
			Vector3D averageForce = Vector3D.vectorAdd(force, forceNew);
			averageForce.scalarDivide(2);
			
			//jump velocity by the average force
			particle.jumpVelocity(dt, averageForce);
			
			//set the force to the force at the new position, copy to avoid referencing
			force = new Vector3D(forceNew);
			
			t = t + dt;
			
			output.printf("%10.5f %s\n", t, particle);
			
		}
		
		output.close();

	}

}
