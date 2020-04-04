package runsimulation;

import java.util.HashSet;
import processing.core.*;
import java.util.Random;
import java.util.Set;

// add "extends PApplet" to make our class an extension of PApplet, which lets us call their methods
// and inherit all their properties
public class RunSimulationRealTime extends PApplet {
    
    SolarSystem sys;
    Hit hits;
    
    // override settings method in Processing to have the settings we want
    @Override
    public void settings() {
        // sets the window size to 1920 x 1080 pixels (1080p)
        size(1280, 720);
    }
    
    // create an array of file names that we want to open and display
    String[] files = {"Sun_data.csv","Earth_data.csv","Jupiter_data.csv","Fatty_data.csv"};
    int n = 2000; // number of steps
    // create a 3D array of positions for the 10000 steps, for each file
    float[][][] pos = new float[n][files.length][3];
    // create a 1D array of times
    float[] times = new float[n];
    
    double timestep = (24*60*60); //f
    // override the setup method in Processing to provide inital values at the start of the program
    // this method is called only once, at the start of the program
    @Override
    public void setup() {
        surface.setTitle("Simulation");
        surface.setResizable(true);
        //frameRate(30);
        /*
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
        */
        
        
        
        
        sys = new SolarSystem(); //initialising the solar sys. Ceating a box
        Data storeSunPos = new Data(n, timestep);
        Data storeEarthPos = new Data(n, timestep); // where the Earth data will go
        Data storeJupiterPos = new Data(n, timestep); 
        //Data storeFattyPos = new Data(n, timestep); //fat planet
        
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
            sys.addObject(bod);
        }*/
        
        
        /*
        double[] E_pos = {1.495978707e11, 0., 0.}; //position of Earth m = radius or orbit
        double[] E_vel = {0., 29.78e3, 0.}; 
        
        double[] J_pos = {7.78574E11, 0., 0.}; //Jupiter's start position
        double[] J_vel = {0., 13.07e3, 0.}; //Jupiter's velocity
              
             
        Body Earth = new Body(E_pos, E_vel, 5.9722e24, 6371e3, "Earth"); //start position and velocity, mass and object radius, "name"
        sys.addObject(Earth); //adding eath to the solar sys. Creaes a sun in the middle too
        
        Body Jupiter = new Body(J_pos, J_vel, 1.8976E27, 69911e3, "Jupiter");
        sys.addObject(Jupiter); //adding Jupiter to the solar system
        */
        
        
        //double[] F_pos = {3.495978707e11, 0., 0.}; //position of Fatty m = radius or orbit
        //double[] F_vel = {0., 20e3, 0.};
        //Body Fatty = new Body(F_pos, F_vel, 2.*1.98847e30, 69911e3, "Fatty");
        //system.addObject(Fatty); //adding Fatty to the solar system
                
        
         double[] Mercury_pos = {2.28E11, 0., 0.}; //Mercury semi major axis. Not updated values yet!
        //double[] Mercury_pos = {2.07E11, 0., 0.}; //Perihelion
        double[] Mercury_vel = {0., 24070., 0.};
        //double[] Mercury_vel = {0., 26500., 0.}; //Max
        Body Mercury = new Body(Mercury_pos, Mercury_vel, 3.3011E+23, 2439.7e+3, "Mercury"); //mass, mean radius
        Mercury.setEccentricity(0.2056);
        sys.addObject(Mercury);
        
        
        double[] V_pos = {1.08E11, 0., 0.}; //Venus semi major axis. Not updated values yet!
        //double[] V_pos = {1.07E11, 0., 0.}; //Perihelion
        double[] V_vel = {0., 35.02e3, 0.};
        //double[] V_vel = {0., 3.53E+04, 0.}; //Max
        Body Venus = new Body(V_pos, V_vel, 4.8675E24, 6051.8e3, "Venus"); //mass, mean radius
        Venus.setEccentricity(0.0068);
        sys.addObject(Venus);
        
        
        double[] E_pos = {1.495978707e11, 0., 0.}; //Average radius of orbit
        //double[] E_pos = {1.52E11, 0., 0.}; //Perihelion
        double[] E_vel = {0., 29780., 0.}; //Mean. Originally set as 29.78e3 m/s
        //double[] E_vel = {0., 30290., 0.}; //Max
        Body Earth = new Body(E_pos, E_vel, 5.9722e24, 6371e3, "Earth"); //start position and velocity, mass and object radius, "name"
        Earth.setEccentricity(0.0167);
        sys.addObject(Earth); //adding eath to the solar system. Creaes a sun in the middle too
        
        double[] Mars_pos = {2.28E11, 0., 0.}; //Mars semi major axis
        //double[] Mars_pos = {2.07E11, 0., 0.}; //Perihelion
        double[] Mars_vel = {0., 24.007e3, 0.}; //Mean
        //double[] Mars_vel = {0., 2.65E04, 0.}; //Max
        Body Mars = new Body(Mars_pos, Mars_vel, 6.4171E+23, 3389.5e+3, "Mars"); //mass, mean radius
        Mars.setEccentricity(0.0934);
        sys.addObject(Mars);
        
        
        double[] J_pos = {7.78574E11, 0., 0.}; //Average
        //double[] J_pos = {7.78574E11, 0., 0.}; //Perihelion
        double[] J_vel = {0., 13.07e3, 0.}; //Average velocity
        //double[] J_vel = {0., 13.07e3, 0.}; //Max
        Body Jupiter = new Body(J_pos, J_vel, 1.8976E27, 69911e3, "Jupiter");
        Jupiter.setEccentricity(0.0484);
        sys.addObject(Jupiter); //adding Jupiter to the solar system       
        
        
        double[] S_pos = {14.3353E11, 0., 0.}; //Saturn semi major axis. Not updated values yet!
        //double[] S_pos = {1.35E12, 0., 0.}; //Perihelion
        double[] S_vel = {0., 9680., 0.}; //Mean
        // double[] S_vel = {0., 10180., 0.}; //Max
        Body Saturn = new Body(S_pos, S_vel, 5.6834E26, 58232E3, "Saturn"); //mass, mean radius
        Saturn.setEccentricity(0.0542);
        sys.addObject(Saturn);
            
        
        double[] U_pos = {28.7246E11, 0., 0.}; //Saturn semi major axis. Not updated values yet!
        //double[] U_pos = {2.74E+12, 0., 0.}; //Perihelion
        double[] U_vel = {0., 6.80E3, 0.};
        //double[] U_vel = {0., 7110.0, 0.}; //Max
        Body Uranus = new Body(U_pos, U_vel, 8.6813E25, 25362E3, "Uranus"); //mass, mean radius
        Uranus.setEccentricity(0.0472);
        sys.addObject(Uranus);
        
        
        double[] N_pos = {44.9506E11, 0., 0.}; //Semi major axis
        //double[] N_pos = {444445E+7, 0., 0.}; //Perihelion
        double[] N_vel = {0., 5430., 0.}; //Mean
        //double[] N_vel = {0., 5500., 0.0}; //Max
        Body Neptune = new Body(N_pos, N_vel, 1.02413E26, 24622000, "Neptune"); //mass, mean radius
        Neptune.setEccentricity(0.0086);
        sys.addObject(Neptune);
        
        double[] P_pos = {59.0638E+11, 0., 0.}; //Semi major axis
        //double[] P_pos = {443682E+7, 0., 0.}; //Perihelion
        double[] P_vel = {0., 4670, 0.}; //Mean
        //double[] P_vel = {0., 6100, 0.}; //Max
        Body Pluto = new Body(P_pos, P_vel, 1.303E22, 1187000, "Pluto"); //mass, mean radius
        Pluto.setEccentricity(0.2488);
        sys.addObject(Pluto);

        Random rand = new Random();
        for (Body obj : sys.getObjects()){
            if (!obj.name.equals("Sun")){
                double randTheta = rand.nextDouble()*2*Math.PI;
                //double randPhi = rand.nextDouble()*Math.PI;

                double[] pos = obj.getPosition();
                double r = Math.sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
                double[] vel = obj.getVelocity();
                double v = Math.sqrt(vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]); 
                
                double[] newPos = {r*Math.cos(randTheta), r*Math.sin(randTheta), 0.}; // {r*Math.cos(randTheta)*Math.sin(randPhi), r*Math.sin(randTheta)*Math.sin(randPhi), r*Math.cos(randPhi)};
                double[] newVel = {v*Math.sin(randTheta), -v*Math.cos(randTheta), 0.}; // {v*Math.sin(randTheta)*Math.sin(randPhi), v*Math.cos(randTheta)*Math.sin(randPhi), v*Math.cos(randPhi)};

                obj.setPosition(newPos[0], newPos[1], newPos[2]);
                obj.setVelocity(newVel[0], newVel[1], newVel[2]);
            }
        }
        
        int thickness = 2;
        double[] pJ = Jupiter.getPosition();
        double rJ = Math.sqrt(pJ[0]*pJ[0] + pJ[1]*pJ[1] + pJ[2]*pJ[2]);
        double[] astLine_pos = {0.,0.,0.};//{2.*rJ + thickness, rJ, 0.};
        
        sys.generateAsteroidLine(astLine_pos[0], astLine_pos[1], astLine_pos[2], thickness, 500, false);
        //Body Fatty = new Body(F_pos, F_vel, 1e10, 69911e3, "Fatty");
        //sys.addObject(Fatty); //adding Fatty to the solar system
        /*
        for (int i = 0; i < n; i++){
            
            if (i % 1000 == 0){
                System.out.println(i);
            }
            for (Body object : sys.getObjects()){
                
                
                double[] objPos = object.getPosition();
                double[] objVel = object.getVelocity();
                
                if (i % 1000 == 0){
                
                System.out.println(String.format(
                        "Object: %s | Pos: {%g, %g, %g} m | Vel: {%g, %g, %g} ms^-1", 
                        object.name, objPos[0], objPos[1], objPos[2], objVel[0], objVel[1], objVel[2])    
                );
                }
            }
            //sys.stepEuler(timestep); //step using Euler
            sys.stepRK4(timestep); //step using Euler
            //put NEW for loop here
            //for (int j = 0; j < sys.getObjects().length-1; j++)){
            //double[] out = sys.getObject(j).getPosition(); //for every j value, find the corresponding object
            //then store the data in an appropriate array
            //change Data.java class to put all orbits in one array
            //for the time being, just store them individually and save o individual csvs
            //}
            
            int S_ind = sys.findObjectIndex("Sun"); //find where the object is in the list
            double[] S_out = sys.getObject(S_ind).getPosition(); //find where it is in the solar system 
            storeSunPos.addData(S_out[0], S_out[1], S_out[2], i); //store the (x,y,z) coordinates
            
            
            int E_ind = sys.findObjectIndex("Earth"); //find where the object is in the list
            double[] E_out = sys.getObject(E_ind).getPosition(); //find where it is in the solar system 
            storeEarthPos.addData(E_out[0], E_out[1], E_out[2], i); //store the (x,y,z) coordinates
            
            int J_ind = sys.findObjectIndex("Jupiter");
            double[] J_out = sys.getObject(J_ind).getPosition();
            storeJupiterPos.addData(J_out[0], J_out[1], J_out[2], i);

            int F_ind = sys.findObjectIndex("Fatty");
            double[] F_out = sys.getObject(F_ind).getPosition();
            storeFattyPos.addData(F_out[0], F_out[1], F_out[2], i);            
            //end NEW loop
        }
        
        //storeEarthPos.output();
        storeSunPos.writeToCSV("Sun_data.csv");
        storeEarthPos.writeToCSV("Earth_data.csv");
        storeJupiterPos.writeToCSV("Jupiter_data.csv");
        storeFattyPos.writeToCSV("Fatty_data.csv");
        
        
        */
        
        sys.genHitArray();
        
    }
    
    // create a counter to keep track of the set of data we want to plot
    int a = 0;
    float scalePlane = 15.f;
    float scalePlot = 0.5f;
    // override the draw method to display the objects we want on the screen
    // this method is called once every 1/30th(?) of a second, this rate can be changed with the 
    // "frameRate(X)" method, X is the numebr of frames per second
    @Override
    public void draw() {
        // set the background colour to be white
        background(255);
        
        // check if we have reached the end of the data
        if (a < 100000){//< 730){
            // get the current time in the simulation
            float t = (float)(timestep*a / 86400);
            // load the Consolas font and set the text font style to be that PFont
            PFont font = createFont("Consolas.ttf", 16);
            textFont(font);
            // set the text size
            textSize(20);
            // add some text at the position x = 4px, y = 30px which shows the time in days to 2d.p.
            text(String.format("Time: %.2f days", t), 4,30);

            // translate the origin of the screen to the centre (rather than the top left)
            translate(width/2.0f, height/2.0f);
            scale(Math.abs(1.5f*(float)Math.sin(a/500.f)) + 1 * scalePlot);
            // set the stroke weight (i.e. the outlining) to 0
            Body[] objs = sys.getObjects();
            // for every object (other than the Sun)
            for (Body obj : objs) {
                // find its name by separating the file name into "name" + _ + "data.csv"
                //String name = files[i].split("_")[0];
                String name = obj.name;
                double[] objPos = obj.getPosition();
                // get its x, y, z values and times them by some scale factor
                float xpos = (float)((objPos[0])*scalePlane / (148.28e9));
                float ypos = (float)((objPos[1])*scalePlane / (148.28e9));
                float zpos = (float)((objPos[2])*scalePlane / (148.28e9));
                // set stroke weight to be 0 again, to make it affect the current drawn object
                strokeWeight(0);
                if (name.equals("Sun")) {
                    fill(250, 254, 76);
                    ellipse(xpos, -ypos, 20*scalePlot, 20*scalePlot);
                } else if (obj.isAsteroid) {
                    // fill(red, green, blue, alpha)
                    // set the fill colour to be black
                    fill(165, 42, 42);
                    ellipse(xpos, -ypos, 3*scalePlot, 3*scalePlot);
                } else if (name.equals("Earth")) {
                    // fill(red, green, blue, alpha)
                    // set the fill colour to be black
                    fill(101, 215, 255);
                    ellipse(xpos, -ypos, 5*scalePlot, 5*scalePlot);
                    
                }else if (name.equals("Jupiter")) {
                    // fill(red, green, blue, alpha)
                    // set the fill colour to be black
                    fill(214, 101, 50);
                    ellipse(xpos, -ypos, 10*scalePlot, 10*scalePlot);
                    
                }else {
                    // fill(red, green, blue, alpha)
                    // set the fill colour to be black
                    fill(0, 0, 0);
                    ellipse(xpos, -ypos, 10*scalePlot, 10*scalePlot);
                }
                // set the text size
                fill(0,0,0);
                textSize(14);
                // add some text 10 pixels, and 30 degrees, away from the start of the ellipse
                if (!obj.isAsteroid) {
                    text(String.format("%s", name), (int)(xpos + 12*Math.cos((float)Math.PI/4.)), (int)(-ypos - 12*Math.sin((float)Math.PI/4.)));
                    // create an ellipse at x=xpos and y=ypos with a radius of 10 pixels
                }
            }
            sys.stepRK4(timestep);
            
            sys.Hits.checkHit(sys);
            // after we display each object, increment the step we are on
            a++;
        // if we have reached the end of the data
        } else {
            // set the current position to be 0, looping the data
            a = 0;
        }
    }
    
    public static void main(String[] args) {
        
        // show the applet where our methods are, to help it find the settings, setup and draw methods
        PApplet.main(new String[]{runsimulation.RunSimulationRealTime.class.getName()});
        
    }
    
}
