//----------------------------------------------------------------------------------------------------------------------
// Copyright (c) 2022 MeteoSwiss, Bruno Zuercher.
// Published under the BSD-3-Clause license.
//----------------------------------------------------------------------------------------------------------------------
package bxr.interpol;

/**
 * The IGridInterpolation interface defines the methods for algorithms that interpolate irregularly spaced
 * observation data to a regular rectangular grid.
 *
 * @author <a href="mailto:bruno.zuercher@meteoswiss.ch">Bruno Z&uuml;rcher</a>
 * @since 2020-12-19
 */
public interface IGridInterpolation {
	
	/**
	 * Specifies the N observation data points that has to be interpolated to a regular rectangular grid.
	 * 
	 * @param pts	    The N x 2 array that specifies the x- and y-coordinates of the observation points.
	 * @param val		The array of length N that specifies the values of the observations.
	 */
	public void setObservations(double pts[][], double val[]);
	
	/**
	 * Specifies the grid, which the observation values are interpolate to.
	 *
     * @param x0        The start coordinate in x-direction.
     * @param y0        The start coordinate in y-direction.
     * @param nX        The number of grid points in x-direction.
     * @param nY        The number of grid points in y-direction.
     * @param step      The step between the grid points (in x- and y-direction).
     */
	public void setGrid(double x0, double y0, int nX, int nY, double step);

	/**
	 * Computes the interpolation of the irregularly spaced observation data specified by
	 * {@link #setObservations(double[][], double[]) to the regular grid given by
	 * {@link #setGrid(double, double, int, int, double)} .
	 * 
	 * @return			The 2-dimensional interpolation data array.
	 */
	public double[][] interpolate();


    /**
     * Returns a short name that describes this IGridInterpolation instance.
     * Useful for naming files.
     * @return          The short name.
     */
	public String getName();
}
