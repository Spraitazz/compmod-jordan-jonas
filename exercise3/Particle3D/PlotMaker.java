import java.io.*;
import java.util.*;

/**
 * Part of computer modelling exercise 3, this program takes either the original
 * verlet or euler output file and creates 4 new files with the x-t, y-t, z-t and 
 * x-y data so that this can then be plotted with xmgrace.
 * 
 * @author Louis Headley
 * @author Jonas Paulavicius
 *
 */
public class PlotMaker {


    public static void main(String[] args) throws IOException, FileNotFoundException {
	
    //the filename of the original .out file is the 1st commandline argument
	String fileName = args[0];
	
	//the integer number of steps is the 2nd commandline argument
	int numstep = Integer.valueOf(args[1]);	

	//create scanner to read values	
	Scanner scanner = new Scanner(new File(fileName));



	//output to 4 different files

	String outFile1 = fileName.substring(0, fileName.indexOf(".")) + "tx.out";
	String outFile2 = fileName.substring(0, fileName.indexOf(".")) + "ty.out";
	String outFile3 = fileName.substring(0, fileName.indexOf(".")) + "tz.out";
	String outFile4 = fileName.substring(0, fileName.indexOf(".")) + "xy.out";

	PrintWriter output1 = new PrintWriter(new FileWriter(outFile1));
	PrintWriter output2 = new PrintWriter(new FileWriter(outFile2));
	PrintWriter output3 = new PrintWriter(new FileWriter(outFile3));
	PrintWriter output4 = new PrintWriter(new FileWriter(outFile4));

	double time;
	double x;
	double y;
	double z;
	    
	//step over each line in the original output file
	for(int j = 0; j < numstep; j++){
	    
		while(!scanner.hasNextDouble()) {
			//advance the scanner past the label
			scanner.next();
		}
		
		//collect the data as doubles
	        
	    x = scanner.nextDouble();
	    y = scanner.nextDouble();
	    z = scanner.nextDouble();	
	    time = scanner.nextDouble();	
	    
	    //write into the appropriate file
	    output1.printf("%10.5f %10.5f \n", time , x );
	    output2.printf("%10.5f %10.5f \n", time, y) ;
	    output3.printf("%10.5f %10.5f \n", time, z);
	    output4.printf("%10.5f %10.5f \n", x, y); 
	}
	
	
	//close out files and scanner
	output1.close();
	output2.close();
	output3.close();
	output4.close();
	scanner.close();
    }

}