//----------------------------------------------------------------------------------------------------------------------
// Copyright (c) 2022 MeteoSwiss, Bruno Zuercher.
// Published under the BSD-3-Clause license.
//----------------------------------------------------------------------------------------------------------------------
package bxr.interpol.barnes;

import bxr.interpol.IGridInterpolation;

/**
 * Provides the base implementation for Barnes-style interpolations.
 *
 * @author <a href="mailto:bruno.zuercher@meteoswiss.ch">Bruno Z&uuml;rcher</a>
 * @since 2020-12-19
 */
public abstract class BaseBarnesInterpolation implements IGridInterpolation {

    /** The x/y coordinates of the observation points (in [&deg;]). */
    protected double        pts[][];
    /** The observation values. */
    protected double        val[];
    /** The number of observation points. */
    protected int           numPts;

    /** The normalization offset for the set of specified observation values (c.f. section 5.6 Round-Off Error Issues). */
    protected double        offset;

    // grid specification
    /** The start coordinate in x-direction. */
    protected double        x0;
    /** The start coordinate in y-direction. */
    protected double        y0;
    /** The number of grid points in x-direction. */
    protected int           nX;
    /** The number of grid points in y-direction. */
    protected int           nY;
    /** The step between the grid points (in x- and y-direction). */
    protected double        step;


    @Override
    public void setObservations(double pts[][], double val[]) {
        this.pts = pts;
        this.val = val;
        this.numPts = val.length;

        // normalize values
        double      min = val[0];
        double      max = val[0];
        for (int k = 1; k < numPts; k++) {
            if (val[k] < min)  min = val[k];
            else if (val[k] > max)  max = val[k];
        }

        this.offset = (min + max) / 2.0;

        // center range of observation values around 0
        for (int k = 0; k < numPts; k++) {
            val[k] -= offset;
        }
    }

    @Override
    public void setGrid(double x0, double y0, int nX, int nY, double step) {
        this.x0 = x0;
        this.y0 = y0;
        this.nX = nX;
        this.nY = nY;
        this.step = step;
    }
}
