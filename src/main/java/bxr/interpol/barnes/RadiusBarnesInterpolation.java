//----------------------------------------------------------------------------------------------------------------------
// Copyright (c) 2022 MeteoSwiss, Bruno Zuercher.
// Published under the BSD-3-Clause license.
//----------------------------------------------------------------------------------------------------------------------
package bxr.interpol.barnes;

import bxr.util.array.DoubleArrayList;
import bxr.util.array.IntegerArrayList;
import bxr.util.kdtree.KdTree;

/**
 * Implements the "radius Barnes interpolation" algorithm with an algorithmic complexity O(N x W x H).
 * This algorithm differs from the naive Barnes interpolation in that it only considers observation points
 * that lie within the radius of influence from the currently considered grid point.
 *
 * @author <a href="mailto:bruno.zuercher@meteoswiss.ch">Bruno Z&uuml;rcher</a>
 * @since 2020-12-19
 */
public class RadiusBarnesInterpolation extends BaseBarnesInterpolation {

    /** The value for sigma. */
    private double              sigma;
    /** The search radius or radius of influence. */
    private double              searchRadius;

    /** The KdTree that holds the locations of the observation points. */
    private KdTree              kdTree;


    /**
     * Constructor.
     *
     * @param sigma             The value for sigma.
     * @param searchRadius      The search radius.
     */
    public RadiusBarnesInterpolation(double sigma, double searchRadius) {
        this.sigma = sigma;
        this.searchRadius = searchRadius;
    }

    /**
     * Computes the weight, above which observation points are taken into consideration.
     * Equal to the Gaussian weight that corresponds to the radius of influence distance.
     *
     * @return      The minimum weight.
     */
    public double getMinimumWeight() {
        double      b = searchRadius / sigma;
        return Math.exp(-b*b/2.0);
    }


    @Override
    public void setObservations(double pts[][], double val[]) {
        super.setObservations(pts, val);

        // build KdTree with location of observation points
        kdTree = new KdTree(pts, true);
    }

    /**
     * {@inheritDoc}
     * <br>
     * Implements radius Barnes interpolation by a 3-fold nested loop, but iterating in the inner-most loop
     * only over those observation points that lie within the radius of influence.
     *
     * @return			The 2-dimensional interpolation data array.
     */
    @Override
    public double[][] interpolate() {
        double		            gridVal[][] = new double [nY][nX];

        KdTree.RadiusSearch	    radiusSearch = kdTree.getRadiusSearch(searchRadius);

        double                  scale = 2*sigma*sigma;
        // holds the coordinates of the current grid point
        double					c[] = new double [2];
        for (int j = 0; j < nY; j++) {
            // compute y-coordinate
            c[1] = y0 + j*step;
            for (int i = 0; i < nX; i++) {
                // compute x-coordinate
                c[0] = x0 + i*step;

                // extract those observation points that lie within radius of influence
                IntegerArrayList    list = radiusSearch.search(c);
                DoubleArrayList     sqrDistList = radiusSearch.getSqrDistances();

                // loop over observation points contained in extracted list and compute numerator and denominator of equ. (1)
                double			weightTotal = 0.0;
                double			weightedSum = 0.0;
                for (int n = list.size()-1; n >= 0 ; n--) {
                    int			index = list.get(n);
                    double		sqrdist = sqrDistList.get(n);
                    double		weight = Math.exp(-sqrdist/scale);
                    weightedSum += weight*val[index];
                    weightTotal += weight;
                }
                // add offset again to resulting quotient
                gridVal[j][i] = weightedSum / weightTotal + offset;
            }
        }

        return gridVal;
    }


    @Override
    public String getName() {
        StringBuilder   bld = new StringBuilder();
        bld.append("Radius_").append(step)
            .append('_').append(sigma).append('_').append(numPts).append('_').append(searchRadius);

        return bld.toString();
    }


    @Override
    public String toString() {
        StringBuilder       bld = new StringBuilder();
        bld.append(this.getClass().getSimpleName()).append('{');
        bld.append("grid=").append(nX).append('x').append(nY).append(',');
        bld.append("step=").append(step).append(',');
        bld.append("sigma=").append(sigma).append(',');
        bld.append("numPts=").append(numPts).append(',');
        bld.append("radius=").append(searchRadius).append('}');
        return bld.toString();
    }
}
