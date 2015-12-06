import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Scanner;


public class NBodyGas {

	public static void main(String[] args) throws IOException {
		
		String paramFile = args[0];
		String energyOutFile = args[1];
		String trajectoryFile = args[2];
		
		Scanner simulationData = new Scanner(new File(paramFile));
		PrintWriter energyOutput = new PrintWriter(new FileWriter(energyOutFile));
		PrintWriter trajectoryOutput = new PrintWriter(new FileWriter(trajectoryFile));

		
		//get the parameters from file
		int particleNo = simulationData.nextInt();
		int stepNo = simulationData.nextInt();
		//how often to output configuration and energy
		int outputFrequency = simulationData.nextInt();
		double temperature = simulationData.nextDouble();
		double density = simulationData.nextDouble();
		double cutoffDistance = simulationData.nextDouble();
		
		Particle3D[] particles = new Particle3D[particleNo];
		MDUtilities.setInitialPositions(density, particles);
		MDUtilities.setInitialVelocities(temperature, particles);
		
		//start trajectory file properly
		
		for (int i = 0; i < stepNo; i++) {
			
		}
		
		//end trajectory file, close all scanners etc.
		
	}

}
