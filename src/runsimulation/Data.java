package runsimulation;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;

// Data class that handles the storage and output of the data generated by the
// simulation
public class Data {
    // example data are poitions as a function of time
    // define the positions and times arrays
    public double[][][] positions;
    public double[][][] velocities;
    public double[] times;
    public String[] objects;
    // constructor takes in the number of values that are to be stored (i.e. iterations of Euler or RK4)
    // and the timestep (to keep time intervals constant)
    public Data(int N, double timestep, String[] names){
        // intiialise arrays with the correct dimensions
        positions = new double[names.length][N][3];
        velocities = new double[names.length][N][3];
        double[] ts = new double[N];
        
        String[] Ns = new String[names.length];
        System.arraycopy(names, 0, Ns, 0, names.length);
        objects = Ns;
        // set all the time values, since they are constant
        for (int i = 0; i < N; i++){
            ts[i] = i*timestep;
        }
        // overwrite the stored array of times
        times = ts;
    }
    
    
    // method to add some position data at some index, where index is the integer
    // multiple of the timestep
    public void addData(double x, double y, double z, double vx, double vy, double vz, int planetIndex, int timeIndex){
        positions[planetIndex][timeIndex][0] = x;
        positions[planetIndex][timeIndex][1] = y;
        positions[planetIndex][timeIndex][2] = z;
        velocities[planetIndex][timeIndex][0] = vx;
        velocities[planetIndex][timeIndex][1] = vy;
        velocities[planetIndex][timeIndex][2] = vz;
    }
    
    // method that dumps all the position data into the console 
    public void output(){
        for (int i = 0; i < positions.length; i++){
            String name = objects[i];
            System.out.println(String.format("Planet: %s",name));
            for (int j = 0; j < positions[0].length; j++){
                double[] pos = positions[i][j];
                double[] vel = velocities[i][j];
                System.out.println(String.format("\tPosition: {%g, %g, %g} m, at t = %g s\n\tVelocity: {%g, %g, %g} m, at t = %g s", pos[0], pos[1], pos[2], times[j], vel[0], vel[1], vel[2], times[j]));
            }
        }
    }
    
    // method to create a CSV file with a desired filename and all of the stored data
    public void writeToCSV(String filename){
        // initialise a stream to write to a file and call it outputFile
        PrintWriter outputFile;
        try {
            // try to create an instance of the stream (creates the CSV file)
            outputFile = new PrintWriter(filename);
        } catch (IOException e) {
            // if there is a problem creating the file, print out an error message and exit
            System.err.println("Failed to open file " + filename + " data was not saved.");
            return;
        }
        
        // First, write out the first line (the column headings for the CSV)
        outputFile.print(",");
        for (String s : objects) {
            outputFile.print(s + ",,,,,,");
        } outputFile.print("\n");
        
        outputFile.print("t [days],");
        for (String s : objects) {
            outputFile.print("x [Au], y [Au], z [Au], vx [Au/day], vy [Au/day], vz [Au/day],");
        } outputFile.print("\n");
                
        //add a loop for x, y and z to include and index, string.fomat
        //start NEW loop
        //outputFile.print("x % [Au], y % [Au], z % [Au]"); // put indexes in at % signs, string.format
        
        // Now make a loop to write the contents of the data to a CSV.
        
        for (int i = 0; i < times.length; i++) {
            outputFile.print(times[i]/86400 + ",");
            for (int n = 0; n < positions.length; n++) {
                outputFile.print(positions[n][i][0]/ (148.28e9) + "," + positions[n][i][1]/ (148.28e9) + "," + positions[n][i][2]/ (148.28e9) + "," + velocities[n][i][0]/ (148.28e9 / (86400)) + "," + velocities[n][i][1]/ (148.28e9 / (86400)) + "," + velocities[n][i][2]/ (148.28e9 / (86400)) + ",");
            } //also do stuff to this with indexes
            outputFile.print("\n");
        }
        
        
        //for (int n = 0; n < positions.length; n++) {
        //    outputFile.println(times[n]/86400 + "," + positions[n][0]/ (148.28e9) + "," + positions[n][1]/ (148.28e9) + "," + positions[n][2]/ (148.28e9));
        //} //also do stuff to this with indexes
        
        //end NEW loop
        
        outputFile.close(); // close the outputFile (stream)
    } 
    
}
