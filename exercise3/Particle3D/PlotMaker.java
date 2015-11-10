// This is a program which will take a file with 

import java.io.*;
import java.util.*;

public class PlotMaker {


    public static void main(String[] args) throws IOException, FileNotFoundException {
	// Argument is filename
	String fileName = args[0];

	int numstep = 100;



	// Create scanner to read values
	
	Scanner scanner = new Scanner(new File(fileName));



	//program will output to 4 different files

	String outFile1 = "Plot1tx.out";
	String outFile2 = "plot2ty.out";
	String outFile3 = "plot3tz.out";
	String outFile4 = "plot4xy.out";

	PrintWriter Output1 = new PrintWriter(new FileWriter(outFile1));
	PrintWriter Output2 = new PrintWriter(new FileWriter(outFile2));
	PrintWriter Output3 = new PrintWriter(new FileWriter(outFile3));
	PrintWriter Output4 = new PrintWriter(new FileWriter(outFile4));


	// Print values in array to files

	double time;
	double x;
	double y;
	double z;
	    
	
	for(int j=0; j<(numstep); j++){
	    

	    time=scanner.nextDouble();
	    scanner.next();
	    scanner.next();
	    x=scanner.nextDouble();
	    y=scanner.nextDouble();
	    z=scanner.nextDouble();
		
		

	    
	    Output1.printf("%10.5f %10.5f \n", time , x );
	    Output2.printf("%10.5f %10.5f \n", time, y) ;
	    Output3.printf("%10.5f %10.5f \n", time, z);
	    Output4.printf("%10.5f %10.5f \n", x,y); 
	}

	Output1.close();
	Output2.close();
	Output3.close();
	Output4.close();
			   


	
	    


    }

}
