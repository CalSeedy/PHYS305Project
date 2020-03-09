package runsimulation;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import processing.core.*;
import java.util.Random;

// add "extends PApplet" to make our class an extension of PApplet, which lets us call their methods
// and inherit all their properties
public class RunSimulation extends PApplet{
    
    // override settings method in Processing to have the settings we want
    @Override
    public void settings() {
        // sets the window size to 1920 x 1080 pixels (1080p)
        size(1280, 720);
    }
    
    // create an array of file names that we want to open and display
    String[] files = {"Earth_data.csv"};//, "Jupiter_data.csv};

    // create a 3D array of positions for the 10000 steps, for each file
    float[][][] pos = new float[10000][files.length][3];
    // create a 1D array of times
    float[] times = new float[10000];
    
    
    // override the setup method in Processing to provide inital values at the start of the program
    // this method is called only once, at the start of the program
    @Override
    public void setup() {
        
        //frameRate(30);
        
        // load all the data from the files into the premade arrays "pos" and "times"
        // should be noted that the data is stored as floats, so there is lossy conversion
        // (Processing doesnt use doubles if I recall correctly)
        for (int i = 0; i < files.length; i++){
            try {
                BufferedReader csvReader = new BufferedReader(new FileReader(files[i]));
                String row;
                int c = 0;
                while ((row = csvReader.readLine()) != null) {
                    if (c != 0){
                        String[] input = row.split(",");
                        float t = new Float(input[0]);
                        float x = new Float(input[1]);
                        float y = new Float(input[2]);
                        float z = new Float(input[3]);

                        times[c-1] = t;
                        float[] position = {x, y, z};
                        //System.out.println(String.format("%d Pos: (%f, %f, %f)", c-1, x, y, z));
                        pos[c-1][i] = position;
                    }
                    c++;
                }
                csvReader.close();
            } catch (IOException ex) {
                Logger.getLogger(RunSimulation.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        
        
        
        
    }
    
    // create a counter to keep track of the set of data we want to plot
    int a = 0;
    
    // override the draw method to display the objects we want on the screen
    // this method is called once every 1/30th(?) of a second, this rate can be changed with the 
    // "frameRate(X)" method, X is the numebr of frames per second
    @Override
    public void draw() {
        // set the background colour to be white
        background(255);

        // check if we have reached the end of the data
        if (a < pos.length){//< 730){
            // get the current time in the simulation
            float t = times[a];
            // load the Consolas font and set the text font style to be that PFont
            PFont font = createFont("Consolas.ttf", 16);
            textFont(font);
            // set the text size
            textSize(20);
            // add some text at the position x = 4px, y = 30px which shows the time in days to 2d.p.
            text(String.format("Time: %.2f days", t), 4,30);

            // translate the origin of the screen to the centre (rather than the top left)
            translate(width/2.0f, height/2.0f);
            // set the stroke weight (i.e. the outlining) to 0
            strokeWeight(0);
            // set the fill colour to be yellow-ish
            fill(200, 154, 0); // fill(red, green, blue, alpha)
            // create an ellipse at (x = 0, y = 0), with both radii of 10 pixels
            ellipse(0, 0, 10, 10); // ellipse(x, y, horizontal radius, vertical radius)
            
            // for every object (other than the Sun)
            for (int i = 0; i < files.length; i++){
                // find its name by separating the file name into "name" + _ + "data.csv"
                String name = files[i].split("_")[0];
                
                // get its x, y, z values and times them by some scale factor
                float xpos = (pos[a][i][0])*100;
                float ypos = (pos[a][i][1])*100;
                float zpos = pos[a][i][2];
                
                // set stroke weight to be 0 again, to make it affect the current drawn object
                strokeWeight(0);
                // set the fill colour to be black
                fill(0, 0, 0);
                // set the text size
                textSize(14);
                // add some text 10 pixels, and 30 degrees, away from the start of the ellipse
                text(String.format("%s", name), (int)(xpos + 10*Math.cos((float)Math.PI/4.)), (int)(ypos + 10*Math.sin((float)Math.PI/4.)));
                // create an ellipse at x=xpos and y=ypos with a radius of 10 pixels
                ellipse(xpos, ypos, 10, 10);
            }
            // after we display each object, increment the step we are on
            a++;
            
        // if we have reached the end of the data
        } else {
            // set the current position to be 0, looping the data
            a = 0;
        }
    }
    
    public static void main(String[] args) {
        int n = 10000;
        double timestep = (86400.)/2.;
        
        SolarSystem system = new SolarSystem();
        
        double[] E_pos = {1.495978707e11, 0., 0.};
        double[] E_vel = {0., 29.78e3, 0.};
        Body Earth = new Body(E_pos, E_vel, 5.9722e24, 6371e3, "Earth");
        system.addObject(Earth);
        
        /*
        double[] J_pos = {4.*1.495978707e11, 0., 0.};
        double[] J_vel = {0., 2.*29.78e3, 0.};
        Body Jupiter = new Body(J_pos, J_vel, 10.*5.9722e24, 10.*6371e3, "Jupiter");
        system.addObject(Jupiter);
        */
        
        Data storeEarthPos = new Data(n, timestep);
        //Data storeJupiterPos = new Data(n, timestep);
        
        
        for (int i = 0; i < n; i++){
            // REMOVE THE BELOW COMMENT BLOCK TO OUTPUT ALL OBJECT POSITION DATA, AT EACH TIMESTEP.
            // WARNING: SIMULATION BECOMES EXTREMELY SLOW
            /*
            if (i % 1000 == 0){
                System.out.println(i);
            }
            
            for (Body object : system.getObjects()){
                
                
                double[] objPos = object.getPosition();
                double[] objVel = object.getVelocity();
                
                if (i % 1000 == 0){
                
                System.out.println(String.format(
                        "Object: %s | Pos: {%g, %g, %g} m | Vel: {%g, %g, %g} ms^-1", 
                        object.name, objPos[0], objPos[1], objPos[2], objVel[0], objVel[1], objVel[2])    
                );
                }
            }*/
            
            //system.stepEuler(timestep);
            system.stepRK4(timestep);
            
            int E_ind = system.findObjectIndex("Earth");
            double[] E_out = system.getObject(E_ind).getPosition();
            storeEarthPos.addData(E_out[0], E_out[1], E_out[2], i);
            //int J_ind = system.findObjectIndex("Jupiter");
            //double[] J_out = system.getObject(J_ind).getPosition();
            //storeEarthPos.addData(J_out[0], J_out[1], J_out[2], i);
        }
        
        //storeEarthPos.output();
        storeEarthPos.writeToCSV("Earth_data.csv");
        //storeJupiterPos.writeToCSV("Jupiter_data.csv");
        
        // show the applet where our methods are, to help it find the settings, setup and draw methods
        PApplet.main(new String[]{runsimulation.RunSimulation.class.getName()});
        
    }
    
}
