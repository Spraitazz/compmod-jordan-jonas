import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Scanner;

/**
 * 
 * Computer modelling project B "Liquid & gas simulations" 
 * This is the main method of a Lennard-Jones fluid simulation. It takes
 * an input file with simulation data and produces output files with the
 * energies, trajectories and the radial distribution function. The simulation
 * system is evolved using a velocity verlet algorithm.
 * 
 * @author Jonas Paulavicius
 * @author Louis Headley
 *
 */

public class NBodyGas {

	public static void main(String[] args) throws IOException {		
		
		//parameter and output files as arguments
		String paramFile = args[0];
		String energyOutFile = args[1];
		String trajectoryFile = args[2];
		
		//open file readers and writers
		Scanner simulationData = new Scanner(new File(paramFile));
		PrintWriter energyOutput = new PrintWriter(new FileWriter(energyOutFile));
		PrintWriter trajectoryOutput = new PrintWriter(new FileWriter(trajectoryFile));
		//PrintWriter separations = new PrintWriter(new FileWriter("sep.out"));
		
		//get the parameters from file
		int particleNo = simulationData.nextInt();
		int stepNo = simulationData.nextInt();		
		double dt = simulationData.nextDouble();		
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
		
		//initialize simulation box, setting the temperature, density
		Vector3D boxSize = MDUtilities.setInitialPositions(density, particles);
		MDUtilities.setInitialVelocities(temperature, particles);
		
		//data for RDF
		double[] deltaPeaks = new double[(int) (cutoffDistance*100)];
		
		//start trajectory file properly - initial position		
		trajectoryOutput.println(particleNo);
		trajectoryOutput.println("Point = 1");
		for (Particle3D p: particles) {
			p.toVMD(trajectoryOutput);
		}
		
		//start energy file at initial position
		double[] energiesInitial = Particle3D.energyArray(particles, boxSize, cutoffDistance);
		energyOutput.println(energiesInitial[0] + " " + energiesInitial[1] + " " + (energiesInitial[0]+energiesInitial[1]));
		
		//force, energy arrays
		Vector3D[] forces = Particle3D.forceArrayLJ(particles, boxSize, cutoffDistance);
		Vector3D[] new_forces;
		Vector3D[] temp_forces = new Vector3D[particleNo];
		double[] energies;
		
		//main velocity verlet integration loop
		for (int i = 0; i < stepNo; i++) {
			
			//STEP 1: jump position, calculate force at new position
			Particle3D.jumpPositionArrayVerlet(dt, particles, forces, boxSize);			
			new_forces = Particle3D.forceArrayLJ(particles, boxSize, cutoffDistance);
			
			//STEP 2: calculate intermediate force (old+new)/2
			for (int j = 0; j < particleNo; j++) {
				temp_forces[j] = Vector3D.vectorAdd(new_forces[j], forces[j]);
				temp_forces[j].scalarDivide(2);				
			}
			
			//STEP3: jump velocity, set forces to the new values
			Particle3D.jumpVelocityArray(dt, particles, temp_forces);
			forces = new_forces.clone();			

			//only print at specific frequency
			if ((i+1) % outputFrequency == 0) {
				
				System.out.println("step " + i);
				
				//print trajectory point
				trajectoryOutput.println(particleNo);
				trajectoryOutput.println("Point = " + (i+1));
				for (Particle3D p: particles) {
					p.toVMD(trajectoryOutput);
				}
				
				//print energies
				energies = Particle3D.energyArray(particles, boxSize, cutoffDistance);
				energyOutput.println(energies[0] + " " + energies[1] + " " + (energies[0]+energies[1]));
				
				//separations.println(Particle3D.minSeparation(particles, boxSize));
				
				//collect data for RDF plotting
				Particle3D.calculateRDF(particles, boxSize, deltaPeaks, cutoffDistance);
			}
			
		}
		
		for (int i = 0; i < deltaPeaks.length; i++) {
			//normalize RDF
			deltaPeaks[i] = deltaPeaks[i] / (stepNo/outputFrequency);
			
		}
		
		Particle3D.printRDF(deltaPeaks, boxSize, particleNo, cutoffDistance);
		
		
		//end trajectory file, close all scanners etc.
		simulationData.close();
		energyOutput.close();
		trajectoryOutput.close();
		//separations.close();		
	}
}