//----------------------------------------------------------------------------------------------------------------------
// Copyright (c) 2022 MeteoSwiss, Bruno Zuercher.
// Published under the BSD-3-Clause license.
//----------------------------------------------------------------------------------------------------------------------
package bxr.interpol.barnesS2;

import bxr.interpol.IGridInterpolation;
import bxr.interpol.barnes.OptConvolBarnesInterpolation;
import bxr.util.projection.LambertConformalProjection;

/**
 * Implements the "optimized convolution Barnes interpolation" algorithm with an algorithmic complexity O(N + W x H)
 * and using the spherical distance measure from S^2.
 * This is done by first mapping the observation points to a Lambert conformal map, applying
 * {@link OptConvolBarnesInterpolation} on this map grid and finally resampling the interpolated data
 * back to the original coordinate system.
 * Thus, we use two approximations to compute Barnes interpolation on S^2, but the delivered results still
 * represents a very good approximation of the exact Barnes interpolation.
 * <br>
 * Following features are embedded:
 * <ul>
 * <li>data range centering (via BaseBarnesInterpolation.setObservations()</li>
 * <li>quantization of result data in interpolate()</li>
 * </ul>
 *
 * @author <a href="mailto:bruno.zuercher@meteoswiss.ch">Bruno Z&uuml;rcher</a>
 * @since 2021-01-23
 */
public class OptConvolBarnesS2Interpolation implements IGridInterpolation {

    /**
     * The used Lambert conformal projection.
     * The scale is chosen such that one unit along the two standard latitudes corresponds to one degree.
     * This allows us to use the same sigma and step as specified for the sphere S^2.
     */
    private final static LambertConformalProjection
        proj = new LambertConformalProjection(11.5, 65.5, 42.5, 180.0/Math.PI, 11.5, 34.5);


    // the statically defined Lambert space grid that completely covers the regular grid used in paper when mapped
    //   with the Lambert conformal projection defined above.
    // TODO: to obtain full generality, this grid definition should be computed dynamically
    /** The x-coordinate of the first grid point of Lambert map grid. */
    private final static double     lam_x0 = -32.0;
    /** The y-coordinate of the first grid point of Lambert map grid. */
    private final static double     lam_y0 = -2.0;
    /** The width of Lambert map grid. */
    private final static double     lam_wX = 64.0;
    /** The height of Lambert map grid. */
    private final static double     lam_wY = 44.0;


    /** The number of points in x-direction on Lambert map grid. */
    private int                     lam_nX;
    /** The number of points in y-direction on Lambert map grid. */
    private int                     lam_nY;


    // specification of regular grid
    /** The start coordinate in x-direction. */
    private double                  reg_x0;
    /** The start coordinate in y-direction. */
    private double                  reg_y0;
    /** The number of grid points in x-direction. */
    private int                     reg_nX;
    /** The number of grid points in y-direction. */
    private int                     reg_nY;

    /** The step between the grid points (in x- and y-direction). Applies also to Lambert grid. */
    private double                  step;

    /** The internally used fast interpolation that is acting on Lambert grid. */
    private OptConvolBarnesInterpolation    optConvolInterpol;

    /** The value for sigma. */
    private double                  sigma;
    /** The number of applied convolutions. */
    private int                     nIter;

    /** The number of observation points. */
    private int                     numPts;


    /**
     * Constructor.
     *
     * @param sigma     The value for sigma.
     * @param nIter     The number of convolutions.
     */
    public OptConvolBarnesS2Interpolation(double sigma, int nIter) {
        this.sigma = sigma;
        this.nIter = nIter;
        this.optConvolInterpol = new OptConvolBarnesInterpolation(sigma, nIter);
    }

    @Override
    public void setObservations(double pts[][], double val[]) {
        this.numPts = pts.length;

        // map first observation points to Lambert map
        for (int k = 0; k < numPts; k++) {
            proj.toMap(pts[k], pts[k], 1);
        }

        // pass these points to internal interpolation instance
        optConvolInterpol.setObservations(pts, val);
    }

    @Override
    public void setGrid(double x0, double y0, int nX, int nY, double step) {
        this.reg_x0 = x0;
        this.reg_y0 = y0;
        this.reg_nX = nX;
        this.reg_nY = nY;

        this.step = step;

        // set grid of internal interpolation instance accordingly
        this.lam_nX = (int)(lam_wX / step);
        this.lam_nY = (int)(lam_wY / step);
        optConvolInterpol.setGrid(lam_x0, lam_y0, lam_nX, lam_nY, step);
    }

    @Override
    public double[][] interpolate() {
        // invoke the two necessary steps
        double      lambertResult[][] = lambertInterpolate();
        return resample(lambertResult);
    }

    /**
     * Invokes the interpolation of the internally hold interpolation instance that acts on the Lambert grid.
     * This basically internal method is publicly exposed in order to allow time measurements.
     *
     * @return      The 2-dimensional interpolation data array.
     */
    public double[][] lambertInterpolate() {
        return optConvolInterpol.interpolate();
    }

    /**
     * Resample data from Lambert conformal grid back to the regular grid.
     * This basically internal method is publicly exposed in order to allow time measurements.
     *
     * @param lambertData   The data on the Lambert grid.
     * @return              The resampled data on the regulas grid.
     */
    public double[][] resample(double lambertData[][]) {
        double              result[][] = new double [reg_nY][reg_nX];

        double              geo[] = new double [2];
        double              map[] = new double [2];
        for (int j = 0; j < reg_nY; j++) {
            // compute y-coordinate in regular grid
            geo[1] = reg_y0 + j*step;
            for (int i = 0; i < reg_nX; i++) {
                // compute x-coordinate in regular grid
                geo[0] = reg_x0 + i*step;
                // lookup corresponding coordinate on Lambert map
                proj.toMap(geo, map, 1);
                // and compute grid point in Lambert map grid
                map[0] = (map[0] - lam_x0) / step;
                map[1] = (map[1] - lam_y0) / step;

                // value in current point of regular grid is bilinear interpolation of 4 neighboring grdi values
                int         ii = (int)(map[0]);
                int         jj = (int)(map[1]);
                double      wx = map[0] - ii;
                double      wy = map[1] - jj;
                result[j][i] = (1.0-wy)*(1.0-wx)*lambertData[jj][ii]
                    + wy*(1.0-wx)*lambertData[jj+1][ii]
                    + wy*wx*lambertData[jj+1][ii+1]
                    + (1.0-wy)*wx*lambertData[jj][ii+1];
            }
        }
        return result;
    }


    @Override
    public String getName() {
        StringBuilder   bld = new StringBuilder();
        bld.append("OptConvolS2_").append(step)
            .append('_').append(sigma).append('_').append(numPts).append('_').append(nIter);

        return bld.toString();
    }


    @Override
    public String toString() {
        StringBuilder       bld = new StringBuilder();
        bld.append(this.getClass().getSimpleName()).append('{');
        bld.append("regGrid=").append(reg_nX).append('x').append(reg_nY).append(',');
        bld.append("lambertGrid=").append(lam_nX).append('x').append(lam_nY).append(',');
        bld.append("step=").append(step).append(',');
        bld.append("sigma=").append(sigma).append(',');
        bld.append("numPts=").append(numPts).append(',');
        bld.append("nIter=").append(nIter).append(',');
        bld.append("rectSize=").append(optConvolInterpol.getRectSize()).append(',');
        bld.append("tailValue=").append(optConvolInterpol.getTailValue(optConvolInterpol.getRectSize())).append(',');
        bld.append("effSigma=").append(optConvolInterpol.getEffectiveSigma()).append('}');
        return bld.toString();
    }
}
