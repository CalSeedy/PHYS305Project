package runsimulation;

import processing.core.*;
import java.util.Random;
import processing.event.MouseEvent;

// add "extends PApplet" to make our class an extension of PApplet, which lets us call their methods
// and inherit all their properties
public class RunSimulationRealTime extends PApplet {
    final static double G = 6.67e-11;
    SolarSystem sys;
    Hit hits;
    boolean elliptical = false;
    // override settings method in Processing to have the settings we want
    @Override
    public void settings() {
        // sets the window size to 1920 x 1080 pixels (1080p)
        size(1280, 720);
    }
    
    // create an array of file names that we want to open and display
    String[] names = {"Sun","Mercury","Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"};//"Fatty_data.csv"};
    int n = 2000; // number of steps
    // create a 1D array of times
    float[] times = new float[n];
    
    double timestep = 0.25*(24*60*60); //
    
    Data storeSystem = new Data(n, timestep, names);
    
    // override the setup method in Processing to provide inital values at the start of the program
    // this method is called only once, at the start of the program
    @Override
    public void setup() {
        surface.setTitle("Simulation");
        surface.setResizable(true);
        sys = new SolarSystem(); //initialising the solar sys. Ceating a box
        
        //double[] F_pos = {3.495978707e11, 0., 0.}; //position of Fatty m = radius or orbit
        //double[] F_vel = {0., 20e3, 0.};
        //Body Fatty = new Body(F_pos, F_vel, 2.*1.98847e30, 69911e3, "Fatty");
        //system.addObject(Fatty); //adding Fatty to the solar system
        
        double radScale = 5e2;
        
        double[] Mercury_pos = {5.79E+10, 0., 0.}; //Mercury semi major axis. Not updated values yet!
        //double[] Mercury_pos = {2.07E11, 0., 0.}; //Perihelion
        double[] Mercury_vel = {0., 4.74E+04, 0.};
        //double[] Mercury_vel = {0., 26500., 0.}; //Max
        Body Mercury = new Body(Mercury_pos, Mercury_vel, 3.3011E+23, radScale*2439.7e+3, "Mercury"); //mass, mean radius
        Mercury.setEccentricity(0.2056);
        sys.addObject(Mercury);


        double[] V_pos = {1.08E11, 0., 0.}; //Venus semi major axis. Not updated values yet!
        //double[] V_pos = {1.07E11, 0., 0.}; //Perihelion
        double[] V_vel = {0., 35.02e3, 0.};
        //double[] V_vel = {0., 3.53E+04, 0.}; //Max
        Body Venus = new Body(V_pos, V_vel, 4.8675E24, radScale*6051.8e3, "Venus"); //mass, mean radius
        Venus.setEccentricity(0.0068);
        sys.addObject(Venus);


        double[] E_pos = {1.495978707e11, 0., 0.}; //Average radius of orbit
        //double[] E_pos = {1.52E11, 0., 0.}; //Perihelion
        double[] E_vel = {0., 29780., 0.}; //Mean. Originally set as 29.78e3 m/s
        //double[] E_vel = {0., 30290., 0.}; //Max
        Body Earth = new Body(E_pos, E_vel, 5.9722e24, radScale*6371e3, "Earth"); //start position and velocity, mass and object radius, "name"
        Earth.setEccentricity(0.0167086);
        sys.addObject(Earth); //adding eath to the solar system. Creaes a sun in the middle too

        double[] Mars_pos = {2.28E11, 0., 0.}; //Mars semi major axis
        //double[] Mars_pos = {2.07E11, 0., 0.}; //Perihelion
        double[] Mars_vel = {0., 24.007e3, 0.}; //Mean
        //double[] Mars_vel = {0., 2.65E04, 0.}; //Max
        Body Mars = new Body(Mars_pos, Mars_vel, 6.4171E+23, radScale*3389.5e+3, "Mars"); //mass, mean radius
        Mars.setEccentricity(0.0934);
        sys.addObject(Mars);


        double[] J_pos = {7.78574E11, 0., 0.}; //Average
        //double[] J_pos = {7.78574E11, 0., 0.}; //Perihelion
        double[] J_vel = {0., 13.07e3, 0.}; //Average velocity
        //double[] J_vel = {0., 13.07e3, 0.}; //Max
        Body Jupiter = new Body(J_pos, J_vel, 1.8976E27, radScale*69911e3, "Jupiter");
        Jupiter.setEccentricity(0.0484);
        sys.addObject(Jupiter); //adding Jupiter to the solar system       


        double[] S_pos = {14.3353E11, 0., 0.}; //Saturn semi major axis. Not updated values yet!
        //double[] S_pos = {1.35E12, 0., 0.}; //Perihelion
        double[] S_vel = {0., 9680., 0.}; //Mean
        //double[] S_vel = {0., 10180., 0.}; //Max
        Body Saturn = new Body(S_pos, S_vel, 5.6834E26, radScale*58232E3, "Saturn"); //mass, mean radius
        Saturn.setEccentricity(0.0542);
        sys.addObject(Saturn);


        double[] U_pos = {28.7246E11, 0., 0.}; //Saturn semi major axis. Not updated values yet!
        //double[] U_pos = {2.74E+12, 0., 0.}; //Perihelion
        double[] U_vel = {0., 6.80E3, 0.};
        //double[] U_vel = {0., 7110.0, 0.}; //Max
        Body Uranus = new Body(U_pos, U_vel, 8.6813E25, radScale*25362E3, "Uranus"); //mass, mean radius
        Uranus.setEccentricity(0.0472);
        sys.addObject(Uranus);


        double[] N_pos = {44.9506E11, 0., 0.}; //Semi major axis
        //double[] N_pos = {444445E+7, 0., 0.}; //Perihelion
        double[] N_vel = {0., 5430., 0.}; //Mean
        //double[] N_vel = {0., 5500., 0.0}; //Max
        Body Neptune = new Body(N_pos, N_vel, 1.02413E26, radScale*24622000, "Neptune"); //mass, mean radius
        Neptune.setEccentricity(0.0086);
        sys.addObject(Neptune);

        double[] P_pos = {59.0638E+11, 0., 0.}; //Semi major axis
        //double[] P_pos = {443682E+7, 0., 0.}; //Perihelion
        double[] P_vel = {0., 4670, 0.}; //Mean
        //double[] P_vel = {0., 6100, 0.}; //Max
        Body Pluto = new Body(P_pos, P_vel, 1.303E22, radScale*1187000, "Pluto"); //mass, mean radius
        Pluto.setEccentricity(0.2488);
        sys.addObject(Pluto);
        
        sys.initPeriods();
      
        
        Body[] objs = sys.getObjects();
        int k = 0;
        for (Body o : objs){
            if (!(o.name.contains("Asteroid") || o.name.equals("Sun"))){
                double e = o.getEccentricity();
                double ap = (1.+e)/(1.-e) * o.getPosition()[0];
                o.setAphelion(ap);
                k++;
            }
        }
        
        
        Random rand = new Random();
       
/*        if (elliptical){ // assume position is semi-major axis
            for (Body obj : sys.getObjects()){
               if (!obj.name.equals("Sun")){
                   obj.setPerihelion(obj.getPosition()[0]);
                   double randTheta = rand.nextDouble()*2*Math.PI;
                   //double randPhi = rand.nextDouble()*Math.PI;
                   double sma = (obj.getAphelion() + obj.getPerihelion())/2.;
                   double e = obj.getEccentricity();
                   double r = sma * (1. - e*e) / (1. + e*Math.cos(randTheta));
                   
                   //double[] vel = obj.getVelocity();
                   double M_s = sys.getObjects()[0].getMass();
                   double v = Math.sqrt(G*(M_s+obj.getMass())*((2./r) - (1./sma))); 

                   double[] newPos = {r*Math.cos(randTheta), r*Math.sin(randTheta), 0.}; // {r*Math.cos(randTheta)*Math.sin(randPhi), r*Math.sin(randTheta)*Math.sin(randPhi), r*Math.cos(randPhi)};
                   double[] newVel = {v*Math.sin(randTheta), -v*Math.cos(randTheta), 0.}; // {v*Math.sin(randTheta)*Math.sin(randPhi), v*Math.cos(randTheta)*Math.sin(randPhi), v*Math.cos(randPhi)};

                   obj.setPosition(newPos[0], newPos[1], newPos[2]);
                   obj.setVelocity(newVel[0], newVel[1], newVel[2]);
               }
           }
            
            
        } else {
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
        } */
        
        int thickness = 200;
        double[] pJ = Pluto.getPosition();
        double rJ = Math.sqrt(pJ[0]*pJ[0] + pJ[1]*pJ[1] + pJ[2]*pJ[2]);
        //double[] astLine_pos = {0.,0.,0.};//{2.*rJ + thickness, rJ, 0.};

        //sys.generateAsteroidLine(astLine_pos[0], astLine_pos[1], astLine_pos[2], thickness, 100, true);
        sys.generateAsteroidCircle(0.,0.,0., 1.5*rJ, 1000, false);
        
        sys.genHitArray();
        translate(width/2.0f, height/2.0f);
        pushMatrix();
    }
    
    // create a counter to keep track of the set of data we want to plot
    int a = 0;
    float scalePlane = 50.f;
    float scalePlot = 1f;
    // override the draw method to display the objects we want on the screen
    // this method is called once every 1/30th(?) of a second, this rate can be changed with the 
    // "frameRate(X)" method, X is the numebr of frames per second
    @Override
    public void draw() {
        // translate the origin of the screen to the centre (rather than the top left)
        
        
        // set the background colour to be white
        background(255);
        
        // check if we have reached the end of the data
        if (a < 2000){//< 730){
            // get the current time in the simulation
            float t = (float)(timestep*a / 86400);
            // load the Consolas font and set the text font style to be that PFont
            PFont font = createFont("Consolas.ttf", 16);
            textFont(font);
            // set the text size
            textSize(20);
            // add some text at the position x = 4px, y = 30px which shows the time in days to 2d.p.
            text(String.format("Time: %.2f days", t), 4,30);
            //scale(Math.abs(1.5f*(float)Math.sin(a/500.f)) + 1 * scalePlot);
            // set the stroke weight (i.e. the outlining) to 0
            popMatrix();
            Body[] objs = sys.getObjects();
            // for every object (other than the Sun)
            int p = 0;
            for (Body obj : objs) {
                // find its name by separating the file name into "name" + _ + "data.csv"
                //String name = files[i].split("_")[0];
                String name = obj.name;
                double[] objPos = obj.getPosition();
                double[] objVel = obj.getVelocity();
                if(!obj.isAsteroid){
                    storeSystem.addPosData(objPos[0], objPos[1], objPos[2], objVel[0], objVel[1], objVel[2], p, a);
                    p++;
                }
                
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
               }else if (name.equals("Earth")) {
                 // fill(red, green, blue, alpha)
                // set the fill colour to be black
                    fill(101, 215, 255);
                    ellipse(xpos, -ypos, 5, 5);
                }else if (name.equals("Mercury")) {
                 // fill(red, green, blue, alpha)
                // set the fill colour to be black
                    fill(126, 125, 108);
                    ellipse(xpos, -ypos, 2, 2);
                }else if (name.equals("Venus")) {
                 // fill(red, green, blue, alpha)
                // set the fill colour to be black
                    fill(251, 189, 195);
                    ellipse(xpos, -ypos, 5, 5);
                }else if (name.equals("Mars")) {
                 // fill(red, green, blue, alpha)
                // set the fill colour to be black
                    fill(247, 0, 0);
                    ellipse(xpos, -ypos, 5, 5);
                }else if (name.equals("Jupiter")) {
                 // fill(red, green, blue, alpha)
                // set the fill colour to be black
                    fill(214, 101, 50);
                    ellipse(xpos, -ypos, 10, 10);
                }else if (name.equals("Saturn")) {
                 // fill(red, green, blue, alpha)
                // set the fill colour to be black
                    fill(249, 212, 80);
                    ellipse(xpos, -ypos, 8, 8);
                }else if (name.equals("Uranus")) {
                 // fill(red, green, blue, alpha)
                // set the fill colour to be black
                    fill(147, 212, 180);
                    ellipse(xpos, -ypos, 7, 7);
                }else if (name.equals("Neptune")) {
                 // fill(red, green, blue, alpha)
                // set the fill colour to be black
                    fill(120, 114, 180);
                    ellipse(xpos, -ypos, 7, 7);
                }else if (name.equals("Pluto")) {
                 // fill(red, green, blue, alpha)
                // set the fill colour to be black
                    fill(127, 114, 42);
                    ellipse(xpos, -ypos, 3, 3);
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
            sys.cleanAsteroids();
            sys.updatePeriods();
            sys.Hits.checkHit(sys);
            // after we display each object, increment the step we are on
            a++;
        // if we have reached the end of the data
        } else {
            storeSystem.writePosToCSV("Simulation.csv");
            // set the current position to be 0, looping the data
            a = 0;
        }
        pushMatrix();
    }
    
    @Override
    public void mouseWheel(MouseEvent event){
        float e = event.getCount();
        if (scalePlane - e > 0.0f){
            scalePlane -= e;
        } else {
            scalePlane -= e/10.0f;
        }
        
        
    }
    
    @Override
    public void mouseDragged(){
        popMatrix();
        translate(mouseX-pmouseX, mouseY-pmouseY);
        pushMatrix();
    }
    
    public static void main(String[] args) {
        
        // show the applet where our methods are, to help it find the settings, setup and draw methods
        PApplet.main(new String[]{runsimulation.RunSimulationRealTime.class.getName()});
        
    }
    
}
