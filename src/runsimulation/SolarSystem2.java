
package runsimulation;
import java.util.ArrayList;
// SolarSystem class manages all the objects that are inside the system,
// handles the addition/removal of any objects and has the responsibility
// of stepping all objects forward by 1 timestep with either Euler or RK4.
public class SolarSystem2 {
    
    final static double G = 6.67e-11; // Create a const. for G in SI [m]^3 [kg]^-1 [s]^-2 
    
    private Body[] objects;         // store a array of all bodies in the system
    
    // blank constructor to initialise the system with a star ("Sun")
    // Sun has the same properties (im SI) as the Sun, from Google
    public SolarSystem2(){
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
            new_objects[i]= current[i];           
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
        return objects; //returns whole list
    }
    
    // getter for a specific Body at some index
    public Body getObject(int index){
        return objects[index]; //retun the single body
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
        return -1; //so you know it hasn't worked
    }
    
    
    // method for the implementation of the Euler integration method, taking in
    // some timestep as dt
    public void stepEuler(double timestep){
        
        /*for every body:
        loop over every OTHER body
            add the foces acting from all other bodies / the accelelerations due to the other bodies
            enact the foce on the body (a each time interval)
            find the total acceleration and directeion acing on every body
            create a list of the accelertion for each body (x,y,z)
        update the velocities and then positions
        */
        
        double[][] accelerations = new double[objects.length][3]; //inituialising an array of accelerations for the bodies
        double[][] position1 = new double[objects.length][3]; //array of start positions
        double[][] positions2 = new double[objects.length][3]; //array of final positions
        double[][] v1 = new double[objects.length][3]; //array of start velocities
        double[][] v2 = new double[objects.length][3]; //array of final velocities      

        //loops over every body and every other body:
            for (int j = 0; j < objects.length; j++){ //j will be m1
                for (int i = 0; i < objects.length; i++){ //i will be m2
                    if (i != j){
                        double[] positionVec = objects[i].vectorToBody(objects[j]); //the vector between the bodies (from m2 to m1)
                        double d = objects[i].distanceToBody(objects[j]); //distance betweem m2 and m1
                        double factor = G*objects[i].getMass()/ (Math.pow(d, 3)); //the acceleration factor
                        accelerations[j][0] += factor * positionVec[0]; //the acceleration on m1 due to m2 in the x-direction is added on in every loop over j
                        accelerations[j][1] += factor * positionVec[1];
                        accelerations[j][2] += factor * positionVec[2];
                        
                        v1[j] = objects[j].getVelocity(); //saving the initial velocity of m1
                        position1[j] = objects[j].getPosition(); //saving the initial position of m1
                    }
                }                
            }
           
    for (int i = 0; i < accelerations.length; i++){
        v2[i][0] = v1[i][0] + accelerations[i][0]*timestep; //the final acceleration of m1 in the x-direction
        v2[i][1] = v1[i][1] + accelerations[i][1]*timestep;
        v2[i][2] = v1[i][2] + accelerations[i][2]*timestep;  
        
        positions2[i][0] = position1[i][0] + v2[i][0]*timestep; //the final velocity of m1 in the x-direction
        positions2[i][1] = position1[i][1] + v2[i][1]*timestep;                        
        positions2[i][2] = position1[i][2] + v2[i][2]*timestep;
        
        objects[i].setPosition(positions2[i][0], positions2[i][1], positions2[i][2]); //put the final position into the array 
        objects[i].setVelocity(v2[i][0], v2[i][1], v2[i][2]); //put the final velocity into the array
        }
    }    
      
    
    public void Hit(){
        int [] hits = new int[objects.length]; //will store the number of hits on each object
        
        for (int j = 0; j < objects.length; j++){
            for (int i = 0; i < objects.length; i++){
                double d = objects[i].distanceToBody(objects[j]);
                double radii = objects[i].getRadius() + objects[j].getRadius();
                if ( d <= radii){
                    
                    double m1 = objects[i].getMass();
                    double m2 = objects[j].getMass();
                    
                    if (m1 / m2 > Math.pow(10,6)) { //if m2 is much smaller than m1
                        //m2 will merge with m1 and transfer its momentum to m1

                        double[] momentum1 = objects[i].getMomentum();
                        double[] momentum2 = objects[j].getMomentum();
                        double[] new_momentum =  {0.,0.,0.};
                        new_momentum[0] = momentum1[0]+momentum2[0];
                        new_momentum[1] = momentum1[1]+momentum2[1];
                        new_momentum[2] = momentum1[2]+momentum2[2];
                    
                        double new_mass = objects[i].getMass()+objects[j].getMass();
                    
                        String name1 = objects[j].getName(); //the names of the objects which have collided
                        String name2 = objects[i].getName();
                        objects[j].updateMass(objects[i].getMass());
                        objects[j].updateVelocity(new_momentum);
                        removeObject(i);

                        hits[j] ++;
                        hits[i] ++;
                    
                    //update the mass and momentum x of the first body to include the mass of the second
                    //update the speed of the first body
                    //delete the seocnd body
                    //record what body hit what
                    //or should I just delete the body with 'asteroid' in its name??
                    }
                }
            }
        }
    }
    
    // method for the implementation of the 4th Order Runge-Kutt integration 
    // method, taking in some timestep as dt
    public void stepRK4(double timestep){
        // TODO: ADD 4th Order Runge-Kutta method for propogation 
    }
    
    
}
