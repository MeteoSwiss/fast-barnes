//----------------------------------------------------------------------------------------------------------------------
// Copyright (c) 2022 MeteoSwiss, Bruno Zuercher.
// Published under the BSD-3-Clause license.
//----------------------------------------------------------------------------------------------------------------------
package bxr;

import bxr.interpol.IGridInterpolation;
import bxr.interpol.barnesS2.NaiveBarnesS2Interpolation;
import bxr.interpol.barnesS2.OptConvolBarnesS2Interpolation;

/**
 * Main program that allows to invoke the different Barnes interpolation algorithms on S^2 with various setups.
 *
 * @author <a href="mailto:bruno.zuercher@meteoswiss.ch">Bruno Z&uuml;rcher</a>
 * @since 2021-01-23
 */
public class MainS2 {

    /**
     * Defines the main method that can be modified by some parameters.
     *
     * @param args              Not evaluated at all.
     * @throws Exception        In case of errors.
     */
    public static void main(String args[]) throws Exception {

        /// the adjustable parameters //////////////////////////////////////////////////////////////////////////////////

        // one of:  { "NaiveS2", "OptConvolutionS2", }
        String              method = "OptConvolutionS2";

        // one of:  { 54, 218, 872, 3490 }
        int                 numPoints = 3490;

        // one of: { 4.0, 8.0, 16.0, 32.0, 64.0 }
        double              resolution = 32.0;

        // one of:  { 0.25, 0.5, 1.0, 2.0, 4.0 }
        double              sigma = 1.0;

        // one of:  { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 50 }
        // parameter only applies for method "OptConvolutionS2"
        int                 numIter = 4;

        // one of:  { false, true }
        boolean             writeResultFile = false;

        /// do not modify program below this line //////////////////////////////////////////////////////////////////////


        System.out.println("Barnes Interpolation on S^2");
        System.out.println("***************************");


        // read csv file that contains rows with observation data: lat,lon,value
        double              data[][] = Util.readCsvData("input/obs/PressQFF_202007271200_" + numPoints + ".csv");


        // analyze input data
        double              minMaxArr[][] = Util.getColumnMinMax(data);
        System.out.printf("Geographic range: [%.4f - %.4f]\u00b0N x [%.4f - %.4f]\u00b0E\n",
            minMaxArr[0][0], minMaxArr[0][1], minMaxArr[1][0], minMaxArr[1][1]);
        System.out.printf("Value range: [%.1f - %.1f]\n", minMaxArr[2][0], minMaxArr[2][1]);


        // rearrange input data such that it can be used as input for the interpolation algorithms
        // the coordinates of the observation points
        double              obsPts[][] = new double [data.length][2];
        // the observed values
        double              obsVal[] = new double [data.length];
        for (int k = 0; k < obsPts.length; k++) {
            // we need following order: lon,lat
            obsPts[k][0] = data[k][1];
            obsPts[k][1] = data[k][0];

            obsVal[k] = data[k][2];
        }


        // definition of the grid
        double              step = 1.0 / resolution;
        double              x0 = -26.0 + step;
        double              y0 = 34.5;
        double              wX = 75.0;
        double              wY = 37.5;


        // the interpolation method
        IGridInterpolation  interpolMethod = null;
        switch (method) {
            case "NaiveS2":
                interpolMethod = new NaiveBarnesS2Interpolation(sigma);
                break;
            case "OptConvolutionS2":
                interpolMethod = new OptConvolBarnesS2Interpolation(sigma, numIter);
                break;
        }

        interpolMethod.setObservations(obsPts, obsVal);
        interpolMethod.setGrid(x0, y0, (int) (wX / step), (int) (wY / step), step);


        double              result[][];
        if (interpolMethod instanceof OptConvolBarnesS2Interpolation) {
            // invoke the two necessary steps separately in order to measure separate times
            long            start = System.nanoTime();
            double          lambertResult[][] = ((OptConvolBarnesS2Interpolation) interpolMethod).lambertInterpolate();
            long            interpolTime = System.nanoTime();
            result = ((OptConvolBarnesS2Interpolation) interpolMethod).resample(lambertResult);
            long            totalTime = System.nanoTime();

            System.out.println();
            System.out.println("Computing");
            System.out.println(interpolMethod + "\t\t" + ((interpolTime-start)/1000000L) + " ms\t\t"
                + ((totalTime-interpolTime)/1000000L) + " ms\t\t" + ((totalTime-start)/1000000L) + " ms");

        } else {
            long            start = System.nanoTime();
            result = interpolMethod.interpolate();
            long            time = (System.nanoTime() - start) / 1000000L;

            System.out.println();
            System.out.println("Computing");
            System.out.println(interpolMethod + "\t\t" + time + " ms");
        }


        if (writeResultFile) {
            // write interpolated data to file
            String          fileName = "output/grid/" + interpolMethod.getName() + ".gdat";
            Util.writeGridData(fileName, result, x0, y0, step, step);
            System.out.println();
            System.out.printf("Interpolated data written to file \"%s\"\n", fileName);
        }
    }
}
