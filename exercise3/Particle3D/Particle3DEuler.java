import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;


public class Particle3DEuler {

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
		
		//force on particle due to mass
		Vector3D force = Particle3D.graviForce(particle, mass);
		
		int numstep = 10;
		
		double dt = 0.1;

		double t = 0;
		
		output.printf("%10.5f %10.5f\n %10.5f %10.5f\n", t, particlePos.getX(),particlePos.getY(),particlePos.getZ());
		
		for (int i = 0; i < numstep; i++) {
			
			particle.jumpPosition(dt);
			
			force = Particle3D.graviForce(particle, mass);
			
			particle.jumpVelocity(dt, force);
			
			t = t + dt;

		        Vector3D newParticlePos = particle.getPosition();
			
			output.printf("%10.5f %10.5f\n %10.5f %10.5f", t, newParticlePos.getX(), newParticlePos.getY(), newParticlePos.getZ());
			
		}
		
		output.close();		
		
	}

}
