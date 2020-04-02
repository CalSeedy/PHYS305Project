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
public class RunSimulation extends PApplet {
    // override settings method in Processing to have the settings we want
    @Override
    public void settings() {
        // sets the window size to 1920 x 1080 pixels (1080p)
        size(1280, 720);
    }
    
    // create an array of file names that we want to open and display
    String[] files = {"Sun_data.csv","Earth_data.csv","Jupiter_data.csv","Fatty_data.csv"};

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
                if (name.equals("Sun")){
                    fill(200, 154, 0);
                    ellipse(xpos, -ypos, 10, 10);
                } else {
                 // fill(red, green, blue, alpha)
                // set the fill colour to be black
                    fill(0, 0, 0);
                    ellipse(xpos, -ypos, 10, 10);
                }
                // set the text size
                textSize(14);
                // add some text 10 pixels, and 30 degrees, away from the start of the ellipse
                text(String.format("%s", name), (int)(xpos + 10*Math.cos((float)Math.PI/4.)), (int)(-ypos - 10*Math.sin((float)Math.PI/4.)));
                // create an ellipse at x=xpos and y=ypos with a radius of 10 pixels
                
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
        int n = 10000; // number of steps
        double timestep = (86400.)/2.; //1/2 a day is s
        
        SolarSystem system = new SolarSystem(); //initialising the solar system. Ceating a box
        
        Data storeSunPos = new Data(n, timestep);
        Data storeEarthPos = new Data(n, timestep); // where the Earth data will go
        Data storeJupiterPos = new Data(n, timestep); 
        Data storeFattyPos = new Data(n, timestep); //fat planet
        
        /*
        Random r = new Random();
        for (int i = 0; i < 10; i++){
            double[] pos = {r.nextDouble()*(10+10) - 10,
                            r.nextDouble()*(10+10) - 10,
                            r.nextDouble()*(10+10) - 10};
            double[] vel = {r.nextDouble()*(10+10) - 10,
                            r.nextDouble()*(10+10) - 10,
                            r.nextDouble()*(10+10) - 10};
            String name = String.format("Object%d", i);
            Body bod = new Body(pos, vel, 0., 0, name);
            system.addObject(bod);
        }*/
        
        double[] E_pos = {1.495978707e11, 0., 0.}; //position of Earth m = radius or orbit
        double[] E_vel = {0., 29.78e3, 0.}; 
        
        double[] J_pos = {7.78574E11, 0., 0.}; //Jupiter's start position
        double[] J_vel = {0., 13.07e3, 0.}; //Jupiter's velocity
        
        double[] F_pos = {3.495978707e11, 0., 0.}; //position of Fatty m = radius or orbit
        double[] F_vel = {0., 20e3, 0.};        
             
        Body Earth = new Body(E_pos, E_vel, 5.9722e24, 6371e3, "Earth"); //start position and velocity, mass and object radius, "name"
        system.addObject(Earth); //adding eath to the solar system. Creaes a sun in the middle too
        
        Body Jupiter = new Body(J_pos, J_vel, 1.8976E27, 69911e3, "Jupiter");
        system.addObject(Jupiter); //adding Jupiter to the solar system       
        
        Body Fatty = new Body(F_pos, F_vel, 2.*1.98847e30, 69911e3, "Fatty");
        system.addObject(Fatty); //adding Fatty to the solar system
        */
        
        double[] Mars_pos = {227939200.0E3, 0., 0.}; //Mars semi major axis
        double[] Mars_vel = {0., 24.007e3, 0.};
        Body Mars = new Body(Mars_pos, Mars_vel, 6.4171E23, 3389.5e3, "Mars"); //mass, mean radius
        system.addObject(Mars);
        
        double[] V_pos = {108208000E3, 0., 0.}; //Venus semi major axis. Not updated values yet!
        double[] V_vel = {0., 35.02e3, 0.};
        Body Venus = new Body(V_pos, V_vel, 4.8675E24, 6051.8e3, "Venus"); //mass, mean radius
        system.addObject(Venus);
        
        double[] Mercury_pos = {57909050E3, 0., 0.}; //Mercury semi major axis. Not updated values yet!
        double[] Mercury_vel = {0., 47.362e3, 0.};
        Body Mercury = new Body(Mercury_pos, Mercury_vel, 3.3011E23, 2439.7e3, "Mercury"); //mass, mean radius
        system.addObject(Mercury);
        
        double[] S_pos = {1433.53E9, 0., 0.}; //Saturn semi major axis. Not updated values yet!
        double[] S_vel = {0., 9.68e3, 0.};
        Body Saturn = new Body(S_pos, S_vel, 5.6834E26, 58232e3, "Saturn"); //mass, mean radius
        system.addObject(Saturn);
        
        String [] names = {"Sun", "Earth", "Jupiter", "Mars", "Venus"}; //, "Fatty"
        
        for (int i = 0; i < n; i++){
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
            //system.stepEuler(timestep); //step using Euler
            system.stepRK4(timestep); //step using Runge-Kutta            
            //put NEW for loop here
            //for (int j = 0; j < system.getObjects().length-1; j++)){
            //double[] out = system.getObject(j).getPosition(); //for every j value, find the corresponding object
            //then store the data in an appropriate array
            //change Data.java class to put all orbits in one array
            //for the time being, just store them individually and save o individual csvs
            //}
            
            int S_ind = system.findObjectIndex("Sun"); //find where the object is in the list
            double[] S_out = system.getObject(S_ind).getPosition(); //find where it is in the solar system 
            storeSunPos.addData(S_out[0], S_out[1], S_out[2], i); //store the (x,y,z) coordinates
            
            int E_ind = system.findObjectIndex("Earth"); //find where the object is in the list
            double[] E_out = system.getObject(E_ind).getPosition(); //find where it is in the solar system 
            storeEarthPos.addData(E_out[0], E_out[1], E_out[2], i); //store the (x,y,z) coordinates
            
            int J_ind = system.findObjectIndex("Jupiter");
            double[] J_out = system.getObject(J_ind).getPosition();
            storeJupiterPos.addData(J_out[0], J_out[1], J_out[2], i);

            int F_ind = system.findObjectIndex("Fatty");
            double[] F_out = system.getObject(F_ind).getPosition();
            storeFattyPos.addData(F_out[0], F_out[1], F_out[2], i);            
            //end NEW loop
        }
        
        //storeEarthPos.output();
        storeSunPos.writeToCSV("Sun_data.csv");
        storeEarthPos.writeToCSV("Earth_data.csv");
        storeJupiterPos.writeToCSV("Jupiter_data.csv");
        storeFattyPos.writeToCSV("Fatty_data.csv");
        
        PApplet.main(new String[]{runsimulation.RunSimulation.class.getName()});
    }
}
