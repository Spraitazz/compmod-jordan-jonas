/**
 * 
 * Computer modelling exercise 2 - a class representing a vector in 3 dimensions
 * with the capabilities to add, subtract, dot and cross two vectors, as well as
 * multiply and divide a vector by a scalar, copy a vector and get its magnitude
 * and magnitude squared
 * 
 * @author Jordan Margetts
 * @author Jonas Paulavicius
 *
 */
public class Vector3D {
	
	private double x,y,z;
	
	public double getX() {
		return x;		
	}
	public double getY() {
		return y;		
	}
	public double getZ() {
		return z;		
	}
	
	public void setX(double x) {
		this.x = x;
	}
	public void setY(double y) {
		this.y = y;
	}
	public void setZ(double z) {
		this.z = z;
	}
	
	/**
	 * Creates a zero-vector if no parameters are specified
	 */
	public Vector3D() {
		x = 0;
		y = 0;
		z = 0;
	}
	
	
	/**
	 * Creates a new vector with the specified parameters
	 * 
	 * @param x double, the magnitude in the x direction
	 * @param y double, the magnitude in the y direction
	 * @param z double, the magnitude in the z direction
	 */
	public Vector3D(double x, double y, double z) {
		this.x = x;
		this.y = y;
		this.z = z;
	}
	
	
	/**
	 * Creates a new Vector3D instance that is a copy of the argument
	 * 
	 * @param vector Vector3D, the vector to make a copy of
	 * @return a new vector that is a copy of the argument instance
	 */
	public Vector3D(Vector3D vector) {
		this.x = vector.getX();
		this.y = vector.getY();
		this.z = vector.getZ();
	}	
	
	
	/**
	 * Returns the length of this vector
	 * 
	 * @return magnitude (length) of this vector
	 */
	public double magnitude() {
		return Math.sqrt(x*x + y*y + z*z);		
	}
	
	
	/**
	 * Returns the length squared of this vector
	 * 
	 * @return magnitude squared of this vector
	 */
	public double magnitudeSquared() {
		return x*x + y*y + z*z;
	}
	
	public String toString() {
		return "("+x+", "+y+", "+z+")";
	}
	
	/**
	 * Multiplies this vector by a given scalar
	 * 
	 * @param scalar double, the scalar to multiply this vector by
	 */
	public void scalarMultiply(double scalar) {
		x = x * scalar;
		y = y * scalar;
		z = z * scalar;
	}
	
	/**
	 * Divides this vector by a given scalar
	 * 
	 * @param scalar double, the scalar to divide this vector by
	 */
	public void scalarDivide(double scalar) {
		x = x / scalar;
		y = y / scalar;
		z = z / scalar;
	}
	
	/**
	 * Static method to add two vectors
	 * 
	 * @param v1 Vector3D, a vector
	 * @param v2 Vector3D, a vector
	 * @return new vector v3 = v1 + v2
	 */
	public static Vector3D vectorAdd(Vector3D v1, Vector3D v2) {
		Vector3D newVector = new Vector3D(v1.getX()+v2.getX(), v1.getY()+v2.getY(), v1.getZ()+v2.getZ());
		return newVector;
	}
	
	/**
	 * Static method to subtract two vectors
	 * 
	 * @param v1 Vector3D, a vector
	 * @param v2 Vector3D, a vector
	 * @return new vector v3 = v1 - v2
	 */
	public static Vector3D vectorSubtract(Vector3D v1, Vector3D v2) {
		Vector3D newVector = new Vector3D(v1.getX()-v2.getX(), v1.getY()-v2.getY(), v1.getZ()-v2.getZ());
		return newVector;
	}
	
	/**
	 * Static method that returns the dot product of the given vectors
	 * 
	 * @param v1 Vector3D, a vector
	 * @param v2 Vector3D, a vector
	 * @return new vector v3 = v1 dot v2
	 */
	public static double vectorDot(Vector3D v1, Vector3D v2) {
		return v1.getX()*v2.getX() + v1.getY()*v2.getY() + v1.getZ()*v2.getZ();		
	}
	
	/**
	 * Static method that returns the cross product (a new vector) of the given vectors
	 * 
	 * @param v1 Vector3D, a vector
	 * @param v2 Vector3D, a vector
	 * @return new vector v3 = v1 cross v2
	 */
	public static Vector3D vectorCross(Vector3D v1, Vector3D v2) {
		double xcomp = v1.getY()*v2.getZ()-v1.getZ()*v2.getY();
		double ycomp = v1.getZ()*v2.getX()-v1.getX()*v2.getZ();
		double zcomp = v1.getX()*v2.getY()-v1.getY()*v2.getX();
		Vector3D newVector = new Vector3D(xcomp, ycomp, zcomp);
		return newVector;
	}
	
	public static boolean equal(Vector3D v1, Vector3D v2) {
		if (v1.getX() == v2.getX() && v1.getY() == v2.getY() && v1.getZ() == v2.getZ()){
			return true;
		} else {
			return false;
		}
	}

	public static void main(String[] args) {
		

	}

}
