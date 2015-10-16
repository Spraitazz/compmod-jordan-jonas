/** 
 * An example class for complex numbers, complete with constructors,
 * setters and getters, instance methods to calculate complex conjugates
 * and norms, and static methods to implement simple arithmetic operations. 
 *
 * @author A. Hermann
 * @author A. R. Turner
 * @version "08/2015"
 *
 */

public class Complex {

    /*
     * Properties
     * Note the "private" keyword used to protect the variables
     */
    private double realPart;
    private double imagPart;

    /*
     *  Constructors
     */

    /** Default constructor. Constructs a new Complex, with uninitialised
     *  real and imaginary parts.
     */
    public Complex () {
        // Set to "not-a-number" to indicate uninitialised
	this.setComplex(Double.NaN ,Double.NaN);
    }

    /** Copy constructor. Constructs a new Complex by copying the real and
     *  imaginary parts of another Complex instance.
     *
     * @param original the Complex to be copied
     */
    public Complex (Complex original) {
	this.setComplex(original.getRealPart(), original.getImagPart());
    }

    /** Explicit constructor. Constructs a new Complex from explicitly
     *  given real and imaginary parts.
     *
     * @param real a double giving the real part of the new Complex
     * @param imag a double giving the imaginary part of the new Complex
     */
    public Complex(double real,double imag) {
	this.setComplex(real, imag);
    }

    /*
     * Setters and getters
     */

    /** Convenient set method to set both real and imaginary part.
     *
     * @param real a double to set the real part
     * @param imag a double to set the imaginary part
     */
    public void setComplex(double real, double imag) {
	this.setRealPart(real);
	this.setImagPart(imag);
    }

    // Setters provide the internal variables

    /** Sets the real part of Complex only.
     *
     * @param real a double to set the real part
     */
    public void setRealPart(double real) { this.realPart = real; }

    /** Sets the imaginary part of Complex only.
     *
     * @param imag a double to set the imaginary part
     */
    public void setImagPart(double imag) { this.imagPart = imag; }
    
    // Getters provide access to private internal variables
 
    /** Gets the real part of a Complex.
     *
     * @return a double instance representing the Complex' real part.
     */
    public double getRealPart() { return this.realPart; }

    /** Gets the imaginary part of a Complex.
     *
     * @return a double instance representing the Complex' imaginary part.
     */
    public double getImagPart() { return this.imagPart; }
    
    
    /** Returns a String representation of the Complex. Methods 
     * called 'toString' are automagically called when an object
     * is printed.<br>
     * If the imaginary part is positive (e.g. for (1.0,3.0)), the output
     * is "1.0 + 3.0i". If the imaginary part is negative (e.g. for (2.0,-4.0)),
     * the output is "2.0 - 4.0i".
     *
     * @return a string representation of the Complex instance
     */
    public String toString() { 
      double real = this.getRealPart();
      double imag = this.getImagPart();
      if (imag >= 0.0) {
        return real + " + " + imag + "i";
      } else {
        return real + " - " + Math.abs(imag) + "i";
      }
    } 

    /*
     * Instance methods
     *
     * These operate on a particular instance of an object
     */

    /** Calculates the square modulus A*A (norm squared) of Complex A.
     * 
     * @return a double representing the Complex' square modulus.
     */ 
    public double normSquared() { return  this.getRealPart() * this.getRealPart() + 
                                          this.getImagPart() * this.getImagPart();  }

    /** Calculates the modulus (norm) |A| of the Complex A.
     *
     * @return a double representing the Complex' modulus.
     */
    public double norm() { return Math.sqrt(this.normSquared());  }
    /** Returns the complex conjugate (p,-q) of Complex A=(p,q).
     *
     * @return a Complex representing the instance's complex conjugate.
     */ 
    public Complex conj() {
	return new Complex(this.getRealPart(),-this.getImagPart());
    }
    
    /* 
     * Static methods
     *
     * The reason we use static methods here that return a Complex
     * is that we can build compound statements such as
     * d = Complex.addComplex(a, Complex.subComplex(b,c)); 
     */
    
    /** Adds two Complex numbers.
     *
     * @param a the first Complex to be added
     * @param b the second Complex to be added
     * @return the complex sum of a and b, a+b.
     */
    public static Complex addComplex(Complex a, Complex b) { 
	return new Complex(a.getRealPart() + b.getRealPart(),
                           a.getImagPart() + b.getImagPart());
    }
    
    /** Subtracts two Complex numbers.
     *
     * @param a the subtrahend
     * @param b the subtractor
     * @return the complex difference between a and b, a-b.
     */
    public static Complex subComplex(Complex a, Complex b) { 
	return new Complex(a.getRealPart() - b.getRealPart(),
                           a.getImagPart() - b.getImagPart());
    }
    
    /** Multiplies two Complex numbers a=(p,q) and b=(r,s) as a*b=(p*r-q*s,p*s+q*r).
     *
     * @param a the first factor
     * @param b the second factor
     * @return the complex product of a and b, a*b.
     */
    public static Complex multComplex(Complex a, Complex b) {
	return new Complex(a.getRealPart()*b.getRealPart() - a.getImagPart()*b.getImagPart(), 
                           a.getRealPart()*b.getImagPart() + a.getImagPart()*b.getRealPart()) ;
    }
    
    // Note that we can overload multiply into two forms:
    // multiply by double, AND multiply by complex
    // Java checks whether you are using double or Complex and calls the right method
    /** Multiplies a Complex a=(p,q) with a double b as a*b = (p*b,q*b).
     *
     * @param a a Complex
     * @param b the scalar factor
     * @return the Complex a scaled by the real number b.
     */
    public static Complex multComplex(Complex a, double b) {
	return new Complex(a.getRealPart()*b, a.getImagPart()*b);
    }

    //
    // Note how we try to build more complicated operations
    // out of the existing and working simple ones instead
    // of relying on new and untested code.
    
    /** Divides Complex a=(p,q) by double b as a/b = (p/b,q/b).
     *
     * @param a a Complex
     * @param b the scalar divisor
     * @return the Complex divided by a scalar
     */
    public static Complex divideComplex(Complex a, double b) {
	return new Complex(a.getRealPart()/b, a.getImagPart()/b);
    }
    
    /** Divides Complex a by Complex b as a/b = a*conj(b)/|b|^2.
     *
     * @param a the Complex dividend
     * @param b the Complex divisor
     * @return the Complex division of a/b
     */
    public static Complex divideComplex(Complex a, Complex b) {
	return divideComplex(multComplex(a,b.conj()), b.normSquared()) ;
    }
    
};

