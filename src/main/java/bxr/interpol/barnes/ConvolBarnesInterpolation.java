//----------------------------------------------------------------------------------------------------------------------
// Copyright (c) 2022 MeteoSwiss, Bruno Zuercher.
// Published under the BSD-3-Clause license.
//----------------------------------------------------------------------------------------------------------------------
package bxr.interpol.barnes;

/**
 * Implements the "convolution Barnes interpolation" algorithm with an algorithmic complexity O(N + W x H).
 * The delivered interpolation represents a very good approximation of the exact Barnes interpolation.
 * <br>
 * Following features are embedded:
 * <ul>
 * <li>data range centering (via BaseBarnesInterpolation.setObservations()</li>
 * <li>quantization of result data in interpolate()</li>
 * </ul>
 * Only a rectangular window of fixed size is used, i.e. the actually desired sigma is only approximated.
 *
 * @author <a href="mailto:bruno.zuercher@meteoswiss.ch">Bruno Z&uuml;rcher</a>
 * @since 2020-12-19
 */
public class ConvolBarnesInterpolation extends BaseBarnesInterpolation {

    /** The value for sigma. */
    protected double        sigma;
    /** The number of applied convolutions. */
    protected int           nIter;


    /**
     * Constructor.
     *
     * @param sigma     The value for sigma.
     * @param nIter     The number of convolutions.
     */
    public ConvolBarnesInterpolation(double sigma, int nIter) {
        this.sigma = sigma;
        this.nIter = nIter;
    }

    /**
     * Computes the width 2*T+1 of the rectangular window.
     *
     * @return      The width 2*T+1 of the rectangular window.
     */
    public int getRectSize() {
        double              s = sigma / step;
        int                 halfKernelSize = (int)(Math.sqrt(3.0/nIter)*s + 0.5);
        return 2*halfKernelSize + 1;
    }

    /**
     * Computes the effectively applied sigma, which deviates from the desired sigma, since we rounded the width
     * parameter T to the next integer.
     *
     * @return      The effectively applied sigma.
     */
    public double getEffectiveSigma() {
        int                 rectSize = getRectSize();
        return Math.sqrt(nIter / 12.0 * (rectSize*rectSize - 1)) * step;
    }


    /**
     * {@inheritDoc}
     * <br>
     * Implements algorithm 4 presented in section 4 of paper.
     *
     * @return			The 2-dimensional interpolation data array.
     */
    @Override
    public double[][] interpolate() {
        // the array that stores the numerator field from algorithm 4
        double              valArr[][] = new double[nY][nX];
        // the array that stores the denominator field from algorithm 4
        double              wgtArr[][] = new double[nY][nX];

        // invoke algorithm 1
        injectObservationValues(valArr, wgtArr);

        // the window width 2*T+1
        int                 rectSize = getRectSize();

        // execute algorithm 4
        // convolve rows in x-direction
        double              help[] = new double [nX];
        for (int j = 0; j < nY; j++) {
            // convolve values
            double          res[] = accumulateArray(valArr[j], help, nX, rectSize, nIter);
            if (res != valArr[j]) {
                help = valArr[j];
                valArr[j] = res;
            }

            // convolve weights
            res = accumulateArray(wgtArr[j], help, nX, rectSize, nIter);
            if (res != wgtArr[j]) {
                help = wgtArr[j];
                wgtArr[j] = res;
            }
        }

        // convolve columns in y-direction
        help = new double [nY];
        double              in[] = new double [nY];
        for (int i = 0; i < nX; i++) {
            // convolve values
            // copy data column first from 2-dim array to 1-dim array, which increases the spatial locality of data
            //   and hence increases also the performance
            for (int j = 0; j < nY; j++) {
                in[j] = valArr[j][i];
            }
            double          res[] = accumulateArray(in, help, nY, rectSize, nIter);
            // copy result back to 2-dim array
            for (int j = 0; j < nY; j++) {
                valArr[j][i] = res[j];
            }

            // convolve weights
            // copy data column first from 2-dim array to 1-dim array
            for (int j = 0; j < nY; j++) {
                in[j] = wgtArr[j][i];
            }
            res = accumulateArray(in, help, nY, rectSize, nIter);
            // copy result back to 2-dim array
            for (int j = 0; j < nY; j++) {
                wgtArr[j][i] = res[j];
            }
        }

        // compute limit wgtArr value for which weight > 0.0022, i.e. grid points with greater distance
        //   than 3.5*sigma will evaluate to NaN
        // since we dropped common factors in our computation, we have to revert their cancellation in the
        //   following computation
        double              factor = 1.0;
        for (int k = 1; k < nIter; k++)  factor *= rectSize;
        factor = 12 * factor*factor / 2 / Math.PI / nIter * 0.0022;

        // compute now quotient
        for (int j = 0; j < nY; j++) {
            for (int i = 0; i < nX; i++) {
                if (wgtArr[j][i] > factor) {
                    // consider two things
                    // - add offset again to resulting quotient
                    // - and apply quantization operation:
                    //   here by temporary casting double to float and thus drop around 29 least significant bits
                    valArr[j][i] = (float)(valArr[j][i] / wgtArr[j][i] + offset);

                } else {
                    // grid point more than 3.5*sigma away from nearest observation point
                    valArr[j][i] = Double.NaN;
                }
            }
        }

        return valArr;
    }

    /**
     * Injects the observations values and weights, respectively, into the corresponding fields as described by
     * algorithm 1.
     *
     * @param valArr        The 2-dimensional array that receives the injected obervation values.
     * @param wgtArr        The 2-dimensional array that receives the injected weights.
     */
    protected void injectObservationValues(double valArr[][], double wgtArr[][]) {
        // the grid coordinates of the currently considered sample point
        double			    c[] = new double [2];
        for (int k = 0; k < numPts; k++) {
            // compute grid coordinates of current point
            c[0] = (pts[k][0] - x0) / step;
            c[1] = (pts[k][1] - y0) / step;
            // area factor for discrete Dirac impulse needs not to be taken into account, since it is canceled when dividing fields
            if (c[0] >= 0.0 && c[0] < nX -1 && c[1] >= 0.0 && c[1] < nY -1) {
                // compute weight fractions in x- and y-direction
                int         i = (int)(c[0]);
                double      wx = c[0] - i;
                int         j = (int)(c[1]);
                double      wy = c[1] - j;
                double      w;

                // assign the 4 neighboring grid points their values
                w = (1.0-wx)*(1.0-wy);
                valArr[j][i] += w*val[k];
                wgtArr[j][i] += w;

                w = wx*(1.0-wy);
                valArr[j][i+1] += w*val[k];
                wgtArr[j][i+1] += w;

                w = wx*wy;
                valArr[j+1][i+1] += w*val[k];
                wgtArr[j+1][i+1] += w;

                w = (1.0-wx)*wy;
                valArr[j+1][i] += w*val[k];
                wgtArr[j+1][i] += w;
            }
        }
    }

    /**
     * Computes the nIter-fold convolution of the specified inArr array with a rect-kernel of length rectLen.
     * To obtain the actual convolution with a corresponding uniform distribution, the result would have to be
     * scaled with a factor 1/rectLen^nIter. But this scaling is not implemented, since these factors are canceled
     * when the resulting fields are divided.
     * <br>
     * The contents of the two specified arrays is overwritten.
     *
     * @param inArr         The input data array to be convolved with rect kernel.
     * @param hArr          An auxiliary array of same length.
     * @param arrLen        The length of the two arrays.
     * @param rectLen       The length of the applied rect-kernel, an odd integer number.
     * @param nIter         The number of applied convolutions.
     * @return              The convolved output data array.
     */
    protected double[] accumulateArray(double inArr[], double hArr[], int arrLen, int rectLen, int nIter) {
        // the half window size T
        int             h0 = (rectLen-1) / 2;
        int             h1 = rectLen - h0;
        int             k;
        for (int i = 0; i < nIter; i++) {
            double      accu = 0.0;
            // phase a: window center still outside array
            // accumulate first h0 elements
            for (k = -h0; k < 0; k++) {
                accu += inArr[k + h0];
            }
            // phase b: window center inside array but window does not cover array completely
            // accumulate remaining rectLen elements and write their value into array
            for ( ; k < h1; k++) {
                accu += inArr[k + h0];
                hArr[k] = accu;
            }
            // phase c: window completely contained in array
            // add difference of border elements and write value into array
            for ( ; k < arrLen - h0; k++) {
                accu += (inArr[k + h0] - inArr[k - h1]);
                hArr[k] = accu;
            }
            // phase d (mirroring phase b): window center still inside array but window does not cover array completely
            // de-accumulate elements and write value into array
            for ( ; k < arrLen; k++) {
                accu -= inArr[k - h1];
                hArr[k] = accu;
            }
            // phase e (mirroring phase a): window center left array
            // unnecessary since value is not written

            // hArr contains convolution result of this pass
            // swap arrays and start over next convolution
            double  h[] = inArr;
            inArr = hArr;
            hArr = h;
        }
        return inArr;
    }


    @Override
    public String getName() {
        StringBuilder   bld = new StringBuilder();
        bld.append("Convol_").append(step)
            .append('_').append(sigma).append('_').append(numPts).append('_').append(nIter);

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
        bld.append("nIter=").append(nIter).append(',');
        bld.append("rectSize=").append(getRectSize()).append(',');
        bld.append("effSigma=").append(getEffectiveSigma()).append('}');
        return bld.toString();
    }
}
