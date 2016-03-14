/**
 *
 * A class of utility functions for LJ simulations. Used in project B
 * of Computer Modelling.
 * 
 * @author A. Hermann
 * @author A. R. Turner
 * @version "08/2015"
 */

// Import the IO class
import java.io.*;

public class MDUtilities {

    /**
     * Set initial positions on a cubic lattice and returns the simulation box
     * dimensions as a Vector3D object. This routine assumes that you have classes
     * called Vector3D and Particle3D
     * and that you can set the postion of the Particle3D object using a method
     * called 'setPosition' that takes a Vector3D object as its only argument.
     *
     *
     * @param rho density of the system in reduced units
     * @param p array of Particle3D objects
     * @return the box dimensions as Vector3D object
     *
     */
    public static Vector3D setInitialPositions(double rho, Particle3D[] p) {

	// Set the number of atoms
	int nAtom = p.length;
	// Set the box dimensions
	double boxSize = Math.pow(nAtom/rho, 1.0/3.0);
	// Set number of particles in each dimension
	int nDim = (int)Math.pow(nAtom, 1.0/3.0) + 1;
	// Separation between particles
	double delta = boxSize / nDim;

	// Loop over dimensions setting particle positions
	int iAtom = 0;
	for (int ix=0; ix<nDim; ix++) {
	    for (int iy=0; iy<nDim; iy++) {
		for (int iz=0; iz<nDim; iz++) {
		    if (iAtom < nAtom) {
			Vector3D pos = new Vector3D(ix*delta, iy*delta, iz*delta);
			p[iAtom].setPosition(pos);
			iAtom++;
		    }
		}
	    }
	}

	// Print results
	System.out.printf(" %d atoms placed on a cubic lattice\n", iAtom);
	System.out.printf(" Simulation box dimensions = %10.5f %10.5f %10.5f\n", boxSize, boxSize, boxSize);

	// Return the box size
	return new Vector3D(boxSize, boxSize, boxSize);

    }

    /**
     * Set the initial velocities on a Boltzman distribution. This routine 
     * assumes that you have classes called Vector3D and Particle3D
     * and that you can set the velocity of the Particle3D object using a method
     * called 'setVelocity' that takes a Vector3D object as its only argument.
     *
     * @param temp Temperature in reduced units
     * @param p array of Particle3D objects
     *
     */
    public static void setInitialVelocities(double temp, Particle3D[] p) {

        // Set the number of particles
	int nAtom = p.length;

	// Zero the accumulators
        double xv0 = 0.0;
        double yv0 = 0.0;
        double zv0 = 0.0;
	double vsq = 0.0;

	// Loop over particles, setting random velocities
	for (int i=0; i<nAtom; i++) {
	    // Random initial velocity
	    double xvt = Math.random() - 0.5;
	    double yvt = Math.random() - 0.5;
	    double zvt = Math.random() - 0.5;
	    // Set particle velocity
	    p[i].setVelocity(new Vector3D(xvt, yvt, zvt));
	    // Add to total velocity
	    xv0 = xv0 + xvt;
	    yv0 = yv0 + yvt;
	    zv0 = zv0 + zvt;
	    vsq = vsq + xvt*xvt + yvt*yvt + zvt*zvt;
	}

	// Zero the movement of the centre of mass
	xv0 = xv0 / nAtom;
	yv0 = yv0 / nAtom;
	zv0 = zv0 / nAtom;

	// Boltzman factor
	double k = Math.sqrt(3*nAtom*temp/vsq);
    
	// Zero the accumulators
        double xv0Tot = 0.0;
        double yv0Tot = 0.0;
        double zv0Tot = 0.0;
	double v0sq = 0.0;

	// Loop over particles, scaling initial velocities
	for (int i=0; i<nAtom; i++) {
	    // Make sure velocities are scaled correctly
	    Vector3D vtemp = p[i].getVelocity();
	    double xvt = k * (vtemp.getX() - xv0);
	    double yvt = k * (vtemp.getY() - yv0);
	    double zvt = k * (vtemp.getZ() - zv0);
	    // Set particle velocity
	    p[i].setVelocity(new Vector3D(xvt, yvt, zvt));
	    // Add to total velocity
	    xv0Tot = xv0Tot + xvt;
	    yv0Tot = yv0Tot + yvt;
	    zv0Tot = zv0Tot + zvt;
	    v0sq = v0sq + xvt*xvt + yvt*yvt + zvt*zvt;
	}

	// Output the details
	System.out.printf("  Temperature = %10.5f\n", v0sq / (3*nAtom));
	System.out.printf(" COM Velocity = %10.5f %10.5f %10.5f\n", 
			  xv0Tot/nAtom, yv0Tot/nAtom, zv0Tot/nAtom);

	return;
    }

}