/**
 *
 * Computer Modelling, Exercise 2: Example program to test the Complex class.
 * Creates two Complex numbers from input file, and illustrates the methods
 * implemented in the Complex class.
 *
 * Two methods demonstrate the dangers of assigning Java objects (which is done
 * by reference, not by value) and how to use the copy constructor instead.
 *
 * @author A. Hermann
 * @author A. R. Turner
 * @version "08/2015"
 */

// We need the IO package and Scanner class
import java.io.*;
import java.util.Scanner;

public class ComplexTester {

    /** The main() method, program will run from here. As we are doing file IO,
     *  we need to throw an exception
     */
    public static void main(String[] argv) throws IOException {

        // Create two new complex numbers
	// (i.e. instances of the Complex class)
	Complex A = new Complex();
	Complex B = new Complex();
	
	// Open up the file called "ComplexTester.input"
	// This file should contain four numbers such as "1.0 1.0 2.0 4.0"
	BufferedReader file = new BufferedReader(new FileReader("ComplexTester.input"));
        // Attach a Scanner to the file to parse the numbers
        Scanner scan = new Scanner(file);
	
	// Turn the first four tokens into real and imaginary parts of two complex numbers
        A.setRealPart(scan.nextDouble());
        A.setImagPart(scan.nextDouble());
 
        B.setRealPart(scan.nextDouble());
        B.setImagPart(scan.nextDouble());
	
	// Print out A and B
	System.out.printf("\nA = %s\n", A);
	System.out.printf("B = %s\n", B);
	
	// Test Norm, Conj, Add, Sub, Mult and Divide
	// You should check the results are correct by computing what your input
	// should produce by hand.
	System.out.printf("\nNorm %10.5f %10.5f\n", A.norm(), B.norm());
	System.out.printf("Conj %s, %s\n", A.conj(), B.conj());
	System.out.printf("Add %s\n", Complex.addComplex(A,B));
	System.out.printf("Sub %s\n", Complex.subComplex(A,B));
	System.out.printf("Mult %s\n", Complex.multComplex(A,B));
	System.out.printf("Divide %s\n", Complex.divideComplex(A,B));

	System.out.printf("\nCalling an object assignment test that has a gotcha\n");
	System.out.printf("This illustrates that Java assigns objects by reference not value\n");
	System.out.printf("Explain this to your demonstrator\n");
	GotchaObjectAssignment();

	System.out.printf("\nCalling an object assignment test that is correct\n");
	System.out.printf("This illustrates that Java assigns objects by value if you use a copy constructor\n");
	System.out.printf("Explain this to your demonstrator");
	CorrectObjectAssignment();
	
    };
    
    /** Demonstrates how not to assign object values, because 
     * objects are passed by reference.
     */
    public static void GotchaObjectAssignment() {
	Complex A = new Complex(2,0);	
	Complex B = new Complex(1,0);
	A = B;
	B.setComplex(0,1);
	System.out.printf("\nA = %s\n", A);
	System.out.printf("B = %s\n\n", B);
    };

    /** Demonstrates how to correctly assign the value of one object
     * to another object, by using the copy constructor of the class.
     */ 
    public static void CorrectObjectAssignment() {
	Complex A = new Complex(2,0);	
	Complex B = new Complex(1,0);
	A = new Complex(B);
	B.setComplex(0,1);
	System.out.printf("\nA = %s\n", A);
	System.out.printf("B = %s\n\n", B);
    };

}
