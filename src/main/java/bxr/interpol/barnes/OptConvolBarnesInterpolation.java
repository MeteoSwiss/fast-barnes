//----------------------------------------------------------------------------------------------------------------------
// Copyright (c) 2022 MeteoSwiss, Bruno Zuercher.
// Published under the BSD-3-Clause license.
//----------------------------------------------------------------------------------------------------------------------
package bxr.interpol.barnes;

/**
 * Implements the "optimized convolution Barnes interpolation" algorithm with an algorithmic complexity O(N + W x H).
 * The delivered interpolation represents a very good approximation of the exact Barnes interpolation.
 * In contrast to {@link ConvolBarnesInterpolation}, the used rectangular window mask has on both sides a tail value
 * appended, such that the overall variance of the algorithm corresponds exactly to the wanted sigma^2.
 * <br>
 * Following features are embedded:
 * <ul>
 * <li>data range centering (via BaseBarnesInterpolation.setObservations()</li>
 * <li>quantization of result data in interpolate()</li>
 * </ul>
 *
 * @author <a href="mailto:bruno.zuercher@meteoswiss.ch">Bruno Z&uuml;rcher</a>
 * @since 2020-12-29
 */
public class OptConvolBarnesInterpolation extends ConvolBarnesInterpolation {

    /**
     * Constructor.
     *
     * @param sigma     The value for sigma.
     * @param nIter     The number of convolutions.
     */
    public OptConvolBarnesInterpolation(double sigma, int nIter) {
        super(sigma, nIter);
    }

    /**
     * Computes the width 2*T+1 of the rectangular window.
     *
     * @return      The width 2*T+1 of the rectangular window.
     */
    @Override
    public int getRectSize() {
        double              s = sigma / step;
        int                 halfKernelSize = (int)((Math.sqrt(1.0+12*s*s/nIter) - 1) / 2);
        return 2*halfKernelSize + 1;
    }

    /**
     * Computes the window tail value alpha as described by equ. (12) of paper.
     *
     * @param rectSize      The width 2*T+1 of the rectangular window (expected to be odd).
     * @return              The corresponding tail value.
     */
    public double getTailValue(int rectSize) {
        // compute T
        int                 halfKernelSize = (rectSize - 1) / 2;
        // the variance of the pure rectangular window of length 2*T+1
        double              sigmaRectSqr = (halfKernelSize+1)*halfKernelSize/3.0*step*step;
        // slightly rearranged expression from equ. (12)
        return 0.5*rectSize*(sigma*sigma/nIter - sigmaRectSqr) / ((halfKernelSize+1)*(halfKernelSize+1)*step*step - sigma*sigma/nIter);
    }

    /**
     * Computes the effectively applied sigma.
     * In the case of the {@link OptConvolBarnesInterpolation} class, this value corresponds up to some round-off
     * errors to the value of the initially specified sigma.
     *
     * @return      The effectively applied sigma.
     */
    @Override
    public double getEffectiveSigma() {
        int                 rectSize = getRectSize();
        int                 halfKernelSize = (rectSize - 1) / 2;
        double              a = getTailValue(rectSize);
        double              var = nIter * step*step * (halfKernelSize+1) * (2*a*(halfKernelSize+1) + halfKernelSize*rectSize/3.0) /
            (2*(halfKernelSize+a) + 1);
        return Math.sqrt(var);
    }


    /**
     * {@inheritDoc}
     * <br>
     * Implements algorithm 4 presented in section 4 of paper but optimized for a rectangular window with
     * tail value alpha.
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
        // the window tail value alpha
        double              alpha = getTailValue(rectSize);

        // execute algorithm 4
        // convolve rows in x-direction
        double              help[] = new double [nX];
        for (int j = 0; j < nY; j++) {
            // convolve values
            double          res[] = accumulateTailArray(valArr[j], help, nX, rectSize, nIter, alpha);
            if (res != valArr[j]) {
                help = valArr[j];
                valArr[j] = res;
            }

            // convolve weights
            res = accumulateTailArray(wgtArr[j], help, nX, rectSize, nIter, alpha);
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
            double          res[] = accumulateTailArray(in, help, nY, rectSize, nIter, alpha);
            // copy result back to 2-dim array
            for (int j = 0; j < nY; j++) {
                valArr[j][i] = res[j];
            }

            // convolve weights
            // copy data column first from 2-dim array to 1-dim array
            for (int j = 0; j < nY; j++) {
                in[j] = wgtArr[j][i];
            }
            res = accumulateTailArray(in, help, nY, rectSize, nIter, alpha);
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
     * Computes the nIter-fold convolution of the specified inArr array with a rect-kernel of length rectLen and
     * tail values alpha at both ends.
     * To obtain the actual convolution with a corresponding uniform distribution, the result has to be scaled with
     * a factor 1/(rectLen+2*alpha)^nIter. But this scaling is not implemented, since these factors are canceled
     * when the resulting fields are divided.
     * <br>
     * The contents of the two specified arrays is overwritten.
     *
     * @param inArr         The input data array to be convolved with rect kernel.
     * @param hArr          An auxiliary array of same length.
     * @param arrLen        The length of the two arrays.
     * @param rectLen       The length of the applied rect-kernel, an odd integer number.
     * @param nIter         The number of applied convolutions.
     * @param alpha         The tail value of the rect kernel to obtain optimal variance.
     * @return              The convolved output data array.
     */
    private double[] accumulateTailArray(double inArr[], double hArr[], int arrLen, int rectLen, int nIter, double alpha) {
        // the half window size T
        int             h0 = (rectLen-1) / 2;
        int             h0_1 = h0 + 1;
        int             h1 = rectLen - h0;
        int             k;
        for (int i = 0; i < nIter; i++) {
            // to accumulate values under regular part of window (without tails!)
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
                hArr[k] = accu + alpha*inArr[k+h0_1];
            }
            // phase c: window completely contained in array
            // add difference of border elements and write value into array
            for ( ; k < arrLen - h0_1; k++) {
                accu += (inArr[k + h0] - inArr[k - h1]);
                hArr[k] = accu + alpha*(inArr[k-h1]+inArr[k+h0_1]);
            }
            // phase c': very last element
            accu += (inArr[k + h0] - inArr[k - h1]);
            hArr[k] = accu + alpha*inArr[k-h1];
            k++;
            // phase d (mirroring phase b): window center still inside array but window does not cover array completely
            // de-accumulate elements and write value into array
            for ( ; k < arrLen; k++) {
                accu -= inArr[k - h1];
                hArr[k] = accu + alpha*inArr[k-h1];
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
        bld.append("OptConvol_").append(step)
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
        bld.append("tailValue=").append(getTailValue(getRectSize())).append(',');
        bld.append("effSigma=").append(getEffectiveSigma()).append('}');
        return bld.toString();
    }
}
