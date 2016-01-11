import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Locale;
import java.util.Scanner;


public class NBodyGas {

	public static void main(String[] args) throws IOException {
		//REMOVE LATER
		Locale.setDefault(new Locale("en", "US"));
		
		String paramFile = args[0];
		String energyOutFile = args[1];
		String trajectoryFile = args[2];
		
		
		Scanner simulationData = new Scanner(new File(paramFile));
		PrintWriter energyOutput = new PrintWriter(new FileWriter(energyOutFile));
		PrintWriter trajectoryOutput = new PrintWriter(new FileWriter(trajectoryFile));

		
		//get the parameters from file
		int particleNo = simulationData.nextInt();
		int stepNo = simulationData.nextInt();		
		double timeStep = simulationData.nextDouble();
		//how often to output configuration and energy
		int outputFrequency = simulationData.nextInt();
		double temperature = simulationData.nextDouble();
		double density = simulationData.nextDouble();
		double cutoffDistance = simulationData.nextDouble();
		
		//set up the array to hold the particles, populate, initialise positions and velocities and get the box size
		Particle3D[] particles = new Particle3D[particleNo];
		for (int i = 0; i < particleNo; i++) {
			Particle3D p = new Particle3D();
			p.setMass(1);
			p.setLabel("particle " + i);
			particles[i] = p;
		}
		
		Vector3D boxSize = MDUtilities.setInitialPositions(density, particles);
		MDUtilities.setInitialVelocities(temperature, particles);
		
		
		//data for RDF
		double[] deltaPeaks = new double[150];
		
		//start trajectory file properly - initial position
		
		trajectoryOutput.println(particleNo);
		trajectoryOutput.println("Point = 1");
		for (Particle3D p: particles) {
			p.toVMD(trajectoryOutput);
		}
		
		//start energy file at initial position
		
		//print energies
		double[] energiesInitial = Particle3D.energyArray(particles, boxSize);
		//KE THEN PE THEN TOTAL
		energyOutput.println(energiesInitial[0] + " " + energiesInitial[1] + " " + energiesInitial[2]);
		
		//WHAT SHOULD BE THE VALUE HERE?????
		double dt = timeStep;
		Vector3D[] forces = Particle3D.forceArrayLJ(particles, boxSize);
		
		
		//main verlet loop
		for (int i = 0; i < stepNo; i++) {
			
			Particle3D.jumpPositionArrayVerlet(dt, particles, forces, boxSize);
			
			Vector3D[] new_forces = Particle3D.forceArrayLJ(particles, boxSize);
			
			Vector3D[] temp_forces = new Vector3D[particleNo];
			//sum forces and divide by 2
			for (int j = 0; j < particleNo; j++) {
				temp_forces[j] = Vector3D.vectorAdd(new_forces[j], forces[j]);
				temp_forces[j].scalarDivide(2);
			}
			
			Particle3D.jumpVelocityArray(dt, particles, temp_forces);
			forces = new_forces.clone();
			

			System.out.println("step " + i);
			//only print at specific frequency
			if ((i+1) % outputFrequency == 0) {
				
				
				
				//print trajectory point
				trajectoryOutput.println(particleNo);
				trajectoryOutput.println("Point = " + (i+1));
				for (Particle3D p: particles) {
					p.toVMD(trajectoryOutput);
				}
				
				//print energies
				double[] energies = Particle3D.energyArray(particles, boxSize);
				//KE THEN PE THEN TOTAL
				energyOutput.println(energies[0] + " " + energies[1] + " " + energies[2]);
				
				
				//collect data for RDF plotting
				Particle3D.calculateRDF(particles, boxSize, deltaPeaks);
			}
			
		}
		
		Particle3D.printRDF(deltaPeaks, boxSize, particleNo);
		
		
		//end trajectory file, close all scanners etc.
		simulationData.close();
		energyOutput.close();
		trajectoryOutput.close();
		
		
	}

}
