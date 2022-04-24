//----------------------------------------------------------------------------------------------------------------------
// Copyright (c) 2022 MeteoSwiss, Bruno Zuercher.
// Published under the BSD-3-Clause license.
//----------------------------------------------------------------------------------------------------------------------
package bxr.interpol.barnesS2;

import bxr.interpol.barnes.BaseBarnesInterpolation;

/**
 * Implements the "naive Barnes interpolation" algorithm with an algorithmic complexity O(N x W x H) and using
 * the spherical distance measure from S^2.
 * Delivers an exact interpolation result (up to numerical round-off errors).
 *
 * @author <a href="mailto:bruno.zuercher@meteoswiss.ch">Bruno Z&uuml;rcher</a>
 * @since 2021-01-23
 */
public class NaiveBarnesS2Interpolation extends BaseBarnesInterpolation {

    /** The value for sigma. */
    private double              sigma;


    /**
     * Constructor.
     *
     * @param sigma     The value for sigma.
     */
    public NaiveBarnesS2Interpolation(double sigma) {
        this.sigma = sigma;
    }


    /**
     * {@inheritDoc}
     * <br>
     * Implements naive Barnes interpolation by a 3-fold nested loop and the spherical distance from S^2.
     *
     * @return			The 2-dimensional interpolation data array.
     */
    @Override
    public double[][] interpolate() {
        double		            gridVal[][] = new double [nY][nX];

        double                  scale = 2*sigma*sigma;
        double					lon;
        double                  lat;
        for (int j = 0; j < nY; j++) {
            // compute latitude
            lat = y0 + j*step;
            for (int i = 0; i < nX; i++) {
                // compute longitude
                lon = x0 + i*step;

                // loop over all observation points and compute numerator and denominator of equ. (1)
                double			weightTotal = 0.0;
                double			weightedSum = 0.0;
                for (int n = 0; n < numPts; n++) {
                    double      dist = sphericalDistance(lon, lat, pts[n][0], pts[n][1]);
                    double		weight = Math.exp(-dist*dist/scale);
                    weightedSum += weight*val[n];
                    weightTotal += weight;
                }
                // add offset again to resulting quotient
                gridVal[j][i] = weightedSum / weightTotal + offset;
            }
        }

        return gridVal;
    }

    /** Defines the constant radians per degree. */
    private final static double     RAD_PER_DEGREE = Math.PI / 180.0;

    /**
     * Computes the spherical distance between the specified lon/lat points.
     *
     * @param lon0		The longitude of the first point (in [&deg;]).
     * @param lat0		The latitude of the first point (in [&deg;]).
     * @param lon1		The longitude of the second point (in [&deg;]).
     * @param lat1		The latitude of the second point (in [&deg;]).
     * @return			The spherical distance (in [&deg;]).
     */
    private static double sphericalDistance(double lon0, double lat0, double lon1, double lat1) {
        lon0 *= RAD_PER_DEGREE;
        lat0 *= RAD_PER_DEGREE;
        lon1 *= RAD_PER_DEGREE;
        lat1 *= RAD_PER_DEGREE;

        double		arg = Math.sin(lat0)*Math.sin(lat1) + Math.cos(lat0)*Math.cos(lat1)*Math.cos(lon1-lon0);
        if (arg > 1.0) {
            return 0.0;
        }
        return Math.acos(arg)/RAD_PER_DEGREE;
    }


    @Override
    public String getName() {
        StringBuilder   bld = new StringBuilder();
        bld.append("NaiveS2_").append(step)
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
