import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;


public class Vector3DTester {
	
	

	public static void main(String[] args) throws FileNotFoundException {
		
		String fileName = args[0];
		Scanner scanner = new Scanner(new File(fileName));
		//first read two vectors and perform simple operations
		double x1 = scanner.nextDouble();
		double y1 = scanner.nextDouble();
		double z1 = scanner.nextDouble();
		double x2 = scanner.nextDouble();
		double y2 = scanner.nextDouble();
		double z2 = scanner.nextDouble();
		
		Vector3D v1 = new Vector3D(x1, y1, z1);
		Vector3D v2 = new Vector3D(x2, y2, z2);
		
		System.out.println("Vector 1: " + v1 + " has magnitude: " + v1.magnitude());
		System.out.println("Vector 2: "+ v2 + " has magnitude: " + v2.magnitude());
		System.out.println("Sum: " + Vector3D.vectorAdd(v1, v2));
		System.out.println("Dot product: " + Vector3D.vectorDot(v1, v2));
		System.out.println("Cross product: " + Vector3D.vectorCross(v1, v2));
		System.out.println(System.lineSeparator()+"Now checking vector identities");
		
		double x3 = scanner.nextDouble();
		double y3 = scanner.nextDouble();
		double z3 = scanner.nextDouble();
		
		Vector3D v3 = new Vector3D(x3, y3, z3);
		Vector3D v4 = Vector3D.vectorCross(v1, v2);
		Vector3D v5 = Vector3D.vectorCross(v2, v1);
		v5.scalarMultiply(-1);
		
		if (Vector3D.equal(v4, v5)){
			System.out.println("First property holds");
		} else {
			System.out.println("Error with showing first property");
		}
		
		Vector3D v6 = Vector3D.vectorCross(v1, Vector3D.vectorAdd(v2, v3));
		Vector3D v7 = Vector3D.vectorAdd(Vector3D.vectorCross(v1, v2),Vector3D.vectorCross(v1, v3));
		
		if (Vector3D.equal(v6, v7)){
			System.out.println("Second property holds");
		} else {
			System.out.println("Error with showing second property");
		}
		
		Vector3D v8 = Vector3D.vectorCross(v1, Vector3D.vectorCross(v2, v3));	
		Vector3D v2copy = new Vector3D(v2);
		Vector3D v3copy = new Vector3D(v3);
		v2copy.scalarMultiply(Vector3D.vectorDot(v1, v3));
		v3copy.scalarMultiply(Vector3D.vectorDot(v1, v2));
		Vector3D v9 = new Vector3D(Vector3D.vectorSubtract(v2copy, v3copy));
		
	    if (Vector3D.equal(v8, v9)){
			System.out.println("Third property holds");
		} else {
			System.out.println("Error with showing third property");
		}
	}

}
