//----------------------------------------------------------------------------------------------------------------------
// Copyright (c) 2022 MeteoSwiss, Bruno Zuercher.
// Published under the BSD-3-Clause license.
//----------------------------------------------------------------------------------------------------------------------
package bxr;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;

/**
 * Utility class with different auxiliary methods.
 *
 * @author <a href="mailto:bruno.zuercher@meteoswiss.ch">Bruno Z&uuml;rcher</a>
 * @since 2020-12-21
 */
public class Util {

    /**
     * Reads the N data rows from the specified csv data file.
     * The header line of the file specifies how many data rows and columns follow.
     *
     * @param fileName          The path to the file.
     * @return                  Array of extension rows x columns that contains the read csv data.
     * @throws IOException      In case of an error.
     */
    public static double[][] readCsvData(String fileName) throws IOException {
        BufferedReader          rd = new BufferedReader(new FileReader(fileName));

        // read dimensions from first line of csv file
        String                  lineItems[] = rd.readLine().split(",");
        int                     h = Integer.parseInt(lineItems[0]);
        int                     w = Integer.parseInt(lineItems[1]);

        double                  data[][] = new double [h][w];
        for (int j = 0; j < h; j++) {
            lineItems = rd.readLine().split(",");
            for (int i = 0; i < w; i++) {
                data[j][i] = Double.parseDouble(lineItems[i]);
            }
        }

        return data;
    }


    /**
     * Determines the minimum and maximum of all data columns given a 2-dimensional array.
     * The minimum of column i is returned in result[i][0], the maximum in result[i][1].
     *
     * @param data      The 2-dimensional data array.
     * @return          The corresponding min-max array.
     */
    public static double[][] getColumnMinMax(double data[][]) {
        int             numCol = data[0].length;

        // result array that stores minimum and maximum for each data column
        double          minMaxArr[][] = new double [numCol][2];
        // init array to values of first row
        for (int i = 0; i < numCol; i++) {
            minMaxArr[i][0] = minMaxArr[i][1] = data[0][i];
        }

        for (int j = 1; j < data.length; j++) {
            for (int i = 0; i < numCol; i++) {
                if (data[j][i] < minMaxArr[i][0]) {
                    minMaxArr[i][0] = data[j][i];
                } else if (data[j][i] > minMaxArr[i][1]) {
                    minMaxArr[i][1] = data[j][i];
                }
            }
        }

        return minMaxArr;
    }


    /**
     * Writes the specified 2-dimensional grid data array to a binary file.
     *
     * @param fileName      The file name path to be used.
     * @param data          The 2-dimensional grid data array.
     * @param x0            The x-coordinate of the first grid point.
     * @param y0            The y-coordinate of the first grid point.
     * @param stepX         The step between grid points im x-direction.
     * @param stepY         The step between grid points in y-direction.
     * @throws IOException  In case of an error.
     */
    public static void writeGridData(String fileName, double data[][], double x0, double y0, double stepX, double stepY)
        throws IOException
    {
        File                    file = new File(fileName);
        File                    dir = file.getParentFile();
        if (! dir.isDirectory()) {
            dir.mkdirs();
        }
        DataOutputStream        os = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(file)));

        // write height and width
        int                     h = data.length;
        int                     w = data[0].length;
        os.writeInt(h);
        os.writeInt(w);

        // write grid specification
        os.writeDouble(y0);
        os.writeDouble(x0);
        os.writeDouble(stepY);
        os.writeDouble(stepX);

        // write data array
        for (int j = 0; j < h; j++) {
            for (int i = 0; i < w; i++) {
                os.writeDouble(data[j][i]);
            }
        }

        os.close();
    }
}
