package runsimulation;

import java.util.Random;

// SolarSystem class manages all the objects that are inside the system,
// handles the addition/removal of any objects and has the responsibility
// of stepping all objects forward by 1 timestep with either Euler or RK4.
public class SolarSystem {
    
    final static double G = 6.67e-11;   // Create a const. for G in SI [m]^3 [kg]^-1 [s]^-2 
    
    private Body[] objects;             // store a array of all bodies in the system
    private long astNum;                // Numbers of asteroids
    
    
    // blank constructor to initialise the system with a star ("Sun")
    // Sun has the same properties (im SI) as the Sun, from Google
    public SolarSystem(){
        // when creating a new solar system, create a "Sun" at the centre
        Body[] new_objects = new Body[1];   // initialise a new array of Bodies of length = 1
        double[] S_pos = {0., 0., 0.};      // define the position vector (origin for Sun)
        double[] S_vel = {0., 0., 0.};      // define the velocity vector (no motion for Sun)
        
        // create a new body that has the above vectors and the Google'd values for Solar mass and radius
        Body Sun = new Body(S_pos, S_vel, 1.98847e30, 696342e3, "Sun");
        
        // set the first Body in the array to be this "Sun"
        new_objects[0] = Sun;
        
        // overwrite the null array
        objects = new_objects;
        
        astNum = 0;
    }
    
    
    // method to add a new object to the current system: takes a Body as input,
    // creates a new array, copies the old set of Bodies, adds the new one,
    // then finally overwrites the array of Bodies that is stored
    public void addObject(Body new_Body){
        // copy over the current array of objects (Bodies)
        Body[] current = objects;
        
        // create a new array of Bodies that has 1 extra space
        Body[] new_objects = new Body[current.length+1];
        
        // copy over the current Bodies to the new array
        for(int i = 0; i < current.length; i++){
            new_objects[i] = current[i];           
        }
        // add the new Body
        new_objects[current.length] = new_Body;
        
        // overwrite the stored array
        objects = new_objects;
    }
    
    
    // method to remove an object in the current system: takes an int as input,
    // creates a new array, copies the old set of Bodies (except the one we dont 
    // want), then finally overwrites the array of Bodies that is stored 
    public void removeObject(int index){
        // copy over the current array of objects (Bodies)
        Body[] current = objects;
        
        // create a new array of Bodies that has 1 less space
        Body[] new_objects = new Body[current.length-1];
        
        // create a counter to keep track of the index in the new array
        int k = 0;
        
        // loop over every Body in the current array
        for(int i = 0; i < current.length; i++){
            // and so long as the Body isn't the one we want to eject
            if (i != index){
                // copy that Body over to the new array at the index k
                new_objects[k] = current[i];
                
                // increment the index, k, by 1
                k++;
            }
        }
        // overwrite the stored array of Bodies
        objects = new_objects;      
    }
    
    // getter for the array of objects
    public Body[] getObjects(){
        return objects;
    }
    
    // getter for a specific Body at some index
    public Body getObject(int index){
        return objects[index];
    }
    
    // private method to return the magnitude of a 3D vector
    private double magnitude(double[] vector){
        // R = {x, y, z}, |R| = sqrt(x^2 + y^2 + z^2), |R|^2 = R.R
        // --> |R| = sqrt(R.R)
        return Math.sqrt(dot(vector, vector));
    }
    
    // private method to return the dot product of two 3D vectors
    private double dot(double[] vector1, double[] vector2){
        // A.B = (A_x * B_x) + (A_y * B_y) + (A_z * B_z)
        return (vector1[0]*vector2[0] + vector1[1]*vector2[1] + vector1[2]*vector2[2]);
    }
    
    
    // method to calculate the angle between two vectors (Bodies)
    // takes two Body instances as input, returns the angle theta
    public double angleBetween2Bodies(Body b1, Body b2){
        // define the values that are calculated later (reserve some space)
        double mag1, mag2, theta;
        
        // get the position of Body 1 and Body 2
        double[] pos1 = b1.getPosition();
        double[] pos2 = b2.getPosition();
        
        // calculate the respective magnitudes of Body 1 and 2
        mag1 = magnitude(pos1);
        mag2 = magnitude(pos2);
        
        // A.B = |A||B|cos(theta)
        // --> theta = cos^-1 (A.B / (|A||B|))
        theta = Math.acos((dot(pos1, pos2))/(mag1*mag2));
        
        return theta;
    }
    
    // method to find the index of a Body in bodies through the chosen name
    // returns first occurance of matching Body and -1 if not found.
    public int findObjectIndex(String identifier){
        // create a counter to keep track of current position in objects
        int c = 0;
        
        // while we havent reached the end of the array of Bodies
        while (c < objects.length){
            // for each body, b, in objects (one at a time)
            for (Body b : objects){
                
                // check if the Body has the same name as the inputted "identifier"
                if (b.name.equals(identifier)){
                    // if it does, return the counter as the index of the Body
                    return c;
                } else {
                    // otherwise, increment the counter by 1
                    c++;
                }
            }
        }
        // if we reach the end, the "identifier" wasn't able to be matched with
        // any of the names of the Bodies, we should then return some easily handled
        // value for the index. -1 will throw an ArrayIndexOutOfBounds exception.
        return -1;
    }
    
    
    // method for the implementation of the Euler integration method, taking in
    // some timestep as dt
    public void stepEuler(double timestep){
        // keep track of current body we are referencing
        int j = 0;
              
        // for every body in our system/ array
        for (Body b : objects){
            // first initialise the acceleration, velocity and position vectors
            // A, V, X respectively
            double[] A = {0.,0.,0.};
            double[] V = {0.,0.,0.};
            double[] X = {0.,0.,0.};
            
            // loop over all Bodies in the system/ array
            for (int i = 0; i < objects.length; i++){
                
                // make sure we aren't going to include the current Body
                if (i != j){
                    //System.out.println(String.format("Position: {%g, %g, %g}", X[0], X[1], X[2]));
                    
                    // get the other Body's mass
                    double mass = objects[i].getMass();
                    
                    // calculate the distance vector from the current Body to the other Body
                    double[] distVec = b.vectorToBody(objects[i]);
                    
                    // calculate the magnitude of this distance
                    double d = magnitude(distVec);
                    
                    // create a factor, from the gravitational force equation
                    // --> F(vector) = m_1*a = ((-G*m_1*m_2)/(|r|^3)) * r(vector)
                    // --> factor = (-G*m_2)/(|r|^3) --> a(vector) = factor * r(vector)
                    double factor = (-G*mass)/(d*d*d);
                    
                    
                    //Uncomment Block below for Start position value
                    /*
                    System.out.println(String.format("Start: %s", b.name));
                    System.out.println(String.format("Position: {%g, %g, %g}", objPos[0], objPos[1], objPos[2]));
                    */
                    
                    
                    // add the acceleration to its respective components
                    A[0] += factor * distVec[0];
                    A[1] += factor * distVec[1];
                    A[2] += factor * distVec[2];
                    
                    // Euler integration at step n:
                    // dr(vector)/dt|_(n+1) = dr(vector)/dt|_n + (d^2 r(vector)/dt^2) * dt
                    // r(vector)|_(n+1) = r(vector)|_n + (dr(vector)/dt)|_(n+1) * dt
                    // or --> V_new = V_current + A * dt, X_new = X_current + V_new *dt

                    // get the other Body's position and velocity
                    double[] objVel = b.getVelocity();
                    double[] objPos = b.getPosition();

                    // calculate new velocity components
                    V[0] += objVel[0] + A[0]*timestep;
                    V[1] += objVel[1] + A[1]*timestep;
                    V[2] += objVel[2] + A[2]*timestep;

                    // calculate new position components
                    X[0] += objPos[0] + V[0]*timestep;
                    X[1] += objPos[1] + V[1]*timestep;
                    X[2] += objPos[2] + V[2]*timestep;


                    //Uncomment Block below for end position values
                    /*
                    System.out.println(String.format("End: %s", b.name));
                    System.out.println(String.format("Position: {%g, %g, %g}", X[0], X[1], X[2]));
                    */

                    // update the object's velocity
                    b.setVelocity(V[0], V[1], V[2]);

                    // update the object's position
                    b.setPosition(X[0], X[1], X[2]);
                    
                }
                
                
                
                
                
                
                
                
            }
            // increment the counter, keeping in sync with the current Body
            j++;
        }
    }
    
    private double[] fvector(double[][] fullstate){
        double[] fvec = new double[objects.length*8];
        for (int i = 0; i < objects.length; i++){
            fvec[i*8 + 0] = fullstate[i][3];
            fvec[i*8 + 1] = fullstate[i][4];
            fvec[i*8 + 2] = fullstate[i][5];
            for (int j = 0; j < objects.length; j++){
                if (i != j){
                    //double m1 = fullstate[i][6];
                    double m2 = fullstate[j][6];
                    
                    double[] r21 = {fullstate[i][0] - fullstate[j][0], fullstate[i][1] - fullstate[j][1], fullstate[i][2] - fullstate[j][2]};
                    double d = Math.sqrt(r21[0]*r21[0] + r21[1]*r21[1] + r21[2]*r21[2]);
        
                    double factor = -G*m2 / (d*d*d);
                    
                    // a12 = -G*m1*m2 / (|r21|)^3 * r21
                    double Ax = factor * r21[0];
                    double Ay = factor * r21[1];
                    double Az = factor * r21[2];
                    
                    
                    fvec[i*8 + 3] = Ax;
                    fvec[i*8 + 4] = Ay;
                    fvec[i*8 + 5] = Az;
                }
            }
        }
        return fvec;
    }
    
    /*
    private double[] Gravity(double[] state1, double[] state2){
        double m1 = state1[6];
        double m2 = state2[6];
        double[] r1 = {state1[0], state1[1], state1[2]};
        double[] r2 = {state2[0], state2[1], state2[2]};
        
        double[] r21 = {state2[0]-state1[0], state2[1]-state1[1],state2[2]-state1[2]};
        double d = Math.sqrt(r21[0]*r21[0] + r21[1]*r21[1] + r21[2]*r21[2]);
        
        double factor = -G*m1*m2 / (d*d*d);
        double Fx = factor * r21[0];
        double Fy = factor * r21[1];
        double Fz = factor * r21[2];
        double[] F = {Fx, Fy, Fz};
        return F;
    }
    */
    
    // method for the implementation of the 4th Order Runge-Kutt integration 
    // method, taking in some timestep as dt
    public void stepRK4(double timestep){
        // TODO: ADD 4th Order Runge-Kutta method for propogation
        //double[] a = {0., 0.5, 0.5, 1.};
        //double[] b = {1./6., 1./3., 1./3., 1./6.};
        double[][] state_initial = new double[objects.length][8];
        double[][] state_second, state_third, state_fourth, state_final;
        state_second = state_third = state_fourth = state_final = state_initial;
        int count = 0;
        for (Body obj : objects){
            double[] state = obj.getState();
            state_initial[count] = state;
            count++;
        }
        
        
        // calculating k1 and second state
        double[] fvec1 = fvector(state_initial);
        double[] k1 = new double[fvec1.length];
        count = 0;
        for (double val : fvec1){
            k1[count] = val*timestep;
            count++;
        }
        
        count = 0;
        for (double[] state : state_initial){
            for (int i = 0; i < state.length; i++){
                if (i != 6 || i != 7){
                    state_second[count][i] = state[i] + 0.5*k1[count*8 + i]; 
                } else {
                    state_second[count][i] = state_initial[count][i];
                }
            }
            count++;
        }
        
        
        // calculating k2 and the 3rd state
        double[] fvec2 = fvector(state_second);
        double[] k2 = new double[fvec2.length];
        count = 0;
        for (double val : fvec2){
            k2[count] = val*timestep;
            count++;
        }
        
        count = 0;
        for (double[] state : state_initial){
            for (int i = 0; i < state.length; i++){
                if (i != 6 || i != 7){
                    state_third[count][i] = state[i] + 0.5*k2[count*8 + i]; 
                } else {
                    state_third[count][i] = state_initial[count][i];
                } 
            }
            count++;
        }
        
        
        
        // calculating k3 and 4th state
        double[] fvec3 = fvector(state_third);
        double[] k3 = new double[fvec3.length];
        count = 0;
        for (double val : fvec3){
            k3[count] = val*timestep;
            count++;
        }
        
        count = 0;
        for (double[] state : state_initial){
            for (int i = 0; i < state.length; i++){
                if (i != 6 || i != 7){
                    state_fourth[count][i] = state[i] + k3[count*8 + i]; 
                } else {
                    state_fourth[count][i] = state_initial[count][i];
                }  
            }
            count++;
        }
        
        
        // calculating k4 and final state
        double[] fvec4 = fvector(state_fourth);
        double[] k4 = new double[fvec4.length];
        count = 0;
        for (double val : fvec4){
            k4[count] = val*timestep;
            count++;
        }
        
        count = 0;
        for (double[] state : state_initial){
            for (int i = 0; i < state.length; i++){
                if (i != 6 || i != 7){
                    state_final[count][i] = state[i] + (1./6.)*k1[count*8 + i] + (1./3.)*k2[count*8 + i] + (1./3.)*k3[count*8 + i] + (1./6.)*k4[count*8 + i]; 
                } else {
                    state_final[count][i] = state_initial[count][i];
                } 
            }
            count++;
        }
        
        count = 0;
        for (Body obj : objects){
            double[] new_state = state_final[count];
            obj.setState(new_state);
            count++;
        }
        
        
        
    }
    
    
    // method that creates a new instance of an asteroid (which is just a Body)
    // the asteroid is generated with a random position on a unit sphere and always
    // has a velocity towards the Earth
    public void generateAsteroid(){
        // Find the object that is the furthest away (that isnt an asteroid)
        double furthest = 0.;
        for (Body b : objects){
            double[] pos = b.getPosition();
            double distance = magnitude(pos);
            if ((distance > furthest) && !(b.name.contains("Asteroid"))){
                furthest = distance;
            }
        }
        
        Random rand = new Random();
        
        // Use Spherical Co-ordinates r, theta, phi
        double phi = rand.nextDouble()*(2.*Math.PI);    // generate random angle
        double theta = rand.nextDouble()*Math.PI;       // generate another random angle
        double x, y, z;
        // set the absolute distance from centre of a sphere to be twice the distance as the furthest object
        double r = furthest*2.;                         // set abs. distance to be twice as far away as furthest object in SolarSystem
        
        // 3-D generation
        /*
        x = r*Math.sin(theta)*Math.cos(phi);
        y = r*Math.sin(theta)*Math.sin(phi);
        z = r*Math.cos(theta);
        */
        
        // 2D - generation
        x = r*Math.cos(phi);
        y = r*Math.sin(phi);
        z = 0.;
        
        // set the asteroid position
        double[] a_pos = {x, y, z};

        
        // Find out where the Earth is
        int ind = findObjectIndex("Earth");
        Body E = getObject(ind);
        double[] E_pos = E.getPosition();
        
        // create a velocity such that the asteroid heads towards Earth
        double vx, vy, vz;
        vx = ((x - E_pos[0])/1e7) * rand.nextGaussian() + ((x - E_pos[0])/1e6);
        vy = ((y - E_pos[1])/1e7) * rand.nextGaussian() + ((x - E_pos[1])/1e6);
        vz = ((z - E_pos[2])/1e7) * rand.nextGaussian() + ((x - E_pos[2])/1e6);
        double[] a_vel = {vx , vy, vz};
        
        // set the asteroid's radius such that it is a Gaussian about (50 +/- 30)m
        double a_radius = rand.nextGaussian()*33. + 50.;
        
        // use Density * Volume to calculate mass; assume Density of asteroids = 5000 kg/m^3
        // Volume of a sphere = 4/3 * PI * r^3
        double a_mass = 5000. * (4./3.) * Math.PI * a_radius * a_radius * a_radius;
        
        // create the asteroid with the above calculated values with the name AsteroidX, where X
        // is the current number of asteroids in the system
        Body ast = new Body(a_pos, a_vel, a_mass, a_radius, String.format("Asteroid%d", astNum));
        
        // add the asteroid to the system and increment the number of asteroids
        addObject(ast);
        astNum++; 
    }
    
}
