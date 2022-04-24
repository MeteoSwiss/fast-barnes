//----------------------------------------------------------------------------------------------------------------------
// Copyright (c) 2022 MeteoSwiss, Bruno Zuercher.
// Published under the BSD-3-Clause license.
//----------------------------------------------------------------------------------------------------------------------
package bxr;

import java.util.Arrays;

import bxr.interpol.IGridInterpolation;
import bxr.interpol.barnes.ConvolBarnesInterpolation;
import bxr.interpol.barnes.NaiveBarnesInterpolation;
import bxr.interpol.barnes.RadiusBarnesInterpolation;

/**
 * Performs the time measurements to produce table 3 of paper.
 * Iterates over the different sigma values.
 *
 * @author <a href="mailto:bruno.zuercher@meteoswiss.ch">Bruno Z&uuml;rcher</a>
 * @since 2021-01-09
 */
public class Measurement3 {

    /**
     * Defines the main method.
     *
     * @param args              Not evaluated at all.
     * @throws Exception        In case of errors.
     */
    public static void main(String args[]) throws Exception {

        // the different sigma values
        double              sigmaValues[] = { 0.25, 0.5, 1.0, 2.0, 4.0 };

        // one of:  { 54, 218, 872, 3490 }
        int                 numPoints = 3490;
        // one of: { 4.0, 8.0, 16.0, 32.0, 64.0 }
        double              resolution = 32.0;
        // one of:  { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 50 }
        int                 numIter = 4;

        // the number of measurement repetitions
        int                 repeat = 10;


        System.out.print("Barnes Interpolation Time Measurements 3\n");
        System.out.print("****************************************");


        // give some time to put computer at maximum rest (e.g. install lock screen)
        Thread.sleep(20000L);


        // the loop over the different sigma values
        for (double sigma  : sigmaValues) {
            System.out.println();
            System.out.println();
            System.out.printf("Sigma: %.2f\n", sigma);
            System.out.println();


            // define the interpolation methods to be used
            IGridInterpolation  methods[] = new IGridInterpolation[] {
                new NaiveBarnesInterpolation(sigma),
                new RadiusBarnesInterpolation(sigma, 3.717*sigma),
                new ConvolBarnesInterpolation(sigma, numIter),
            };


            // read csv file that contains rows with observation data: lat,lon,value
            double              data[][] = Util.readCsvData("input/obs/PressQFF_202007271200_" + numPoints + ".csv");

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


            // the test repetition loop
            for (int k = 0; k < repeat; k++) {

                // the loop over the different interpolation methods
                for (IGridInterpolation interpolMethod : methods) {

                    // pass copy of obsVal array since array is modified by interpolation
                    interpolMethod.setObservations(obsPts, Arrays.copyOf(obsVal, obsVal.length));
                    interpolMethod.setGrid(x0, y0, (int) (wX / step), (int) (wY / step), step);

                    long        start = System.nanoTime();
                    double      result[][] = interpolMethod.interpolate();
                    long        time = (System.nanoTime() - start) / 1000000L;

                    System.out.println(interpolMethod + "\t\t" + time + " ms");
                }
            }
        }
    }
}
