//----------------------------------------------------------------------------------------------------------------------
// Copyright (c) 2022 MeteoSwiss, Bruno Zuercher.
// Published under the BSD-3-Clause license.
//----------------------------------------------------------------------------------------------------------------------
package bxr.interpol.barnes;

/**
 * Implements the "naive Barnes interpolation" algorithm with an algorithmic complexity O(N x W x H).
 * Delivers an exact interpolation result (up to numerical round-off errors).
 *
 * @author <a href="mailto:bruno.zuercher@meteoswiss.ch">Bruno Z&uuml;rcher</a>
 * @since 2020-12-19
 */
public class NaiveBarnesInterpolation extends BaseBarnesInterpolation {

    /** The value for sigma. */
    private double          sigma;


    /**
     * Constructor.
     *
     * @param sigma     The value for sigma.
     */
    public NaiveBarnesInterpolation(double sigma) {
        this.sigma = sigma;
    }


    /**
     * {@inheritDoc}
     * <br>
     * Implements naive Barnes interpolation by a 3-fold nested loop.
     *
     * @return			The 2-dimensional interpolation data array.
     */
    @Override
    public double[][] interpolate() {
        double		            gridVal[][] = new double [nY][nX];

        double                  scale = 2*sigma*sigma;
        // holds the coordinates of the current grid point
        double					c[] = new double [2];
        for (int j = 0; j < nY; j++) {
            // compute y-coordinate
            c[1] = y0 + j*step;
            for (int i = 0; i < nX; i++) {
                // compute x-coordinate
                c[0] = x0 + i*step;

                // loop over all observation points and compute numerator and denominator of equ. (1)
                double			weightTotal = 0.0;
                double			weightedSum = 0.0;
                for (int n = 0; n < numPts; n++) {
                    double		sqrdist = getSqrDist(c, pts[n]);
                    double		weight = Math.exp(-sqrdist/scale);
                    weightedSum += weight*val[n];
                    weightTotal += weight;
                }
                // add offset again to resulting quotient
                gridVal[j][i] = weightedSum / weightTotal + offset;
            }
        }

        return gridVal;
    }

    /**
     * Returns the square of the Euclidean distance between pt0 and pt1.
     *
     * @param pt0       The coordinates of the first point.
     * @param pt1       The coordinates of the second point.
     * @return          The square of the Euclidean distance between the points.
     */
    private static double getSqrDist(double pt0[], double pt1[]) {
        return sqr(pt0[0]-pt1[0]) + sqr(pt0[1] - pt1[1]);
    }

    /**
     * Returns the square of the specified value.
     *
     * @param x     The value.
     * @return      The square of x.
     */
    private static double sqr(double x) {
        return x*x;
    }


    @Override
    public String getName() {
        StringBuilder   bld = new StringBuilder();
        bld.append("Naive_").append(step)
            .append('_').append(sigma).append('_').append(numPts);

        return bld.toString();
    }


    @Override
    public String toString() {
        StringBuilder       bld = new StringBuilder();
        bld.append(this.getClass().getSimpleName()).append('{');
        bld.append("grid=").append(nX).append('x').append(nY).append(',');
        bld.append("step=").append(step).append(',');
        bld.append("sigma=").append(sigma).append(',');
        bld.append("numPts=").append(numPts).append('}');
        return bld.toString();
    }
}
