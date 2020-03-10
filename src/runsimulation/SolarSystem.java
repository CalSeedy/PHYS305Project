package runsimulation;

// SolarSystem class manages all the objects that are inside the system,
// handles the addition/removal of any objects and has the responsibility
// of stepping all objects forward by 1 timestep with either Euler or RK4.
public class SolarSystem {
    
    final static double G = 6.67e-11; // Create a const. for G in SI [m]^3 [kg]^-1 [s]^-2 
    
    private Body[] objects;         // store a array of all bodies in the system
    
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
        // keep track of current body we are referencing
        int j = 0;
        
        // for every body in our system/ array
        for (Body b : objects){
            
            
            // loop over all Bodies in the system/ array
            for (int i = 0; i < objects.length; i++){
                // first initialise the acceleration, velocity and position vectors
                // A, V, X respectively
                double[] A = {0.,0.,0.}; //acceleration
                double[] V = {0.,0.,0.}; //velocity
                double[] X = {0.,0.,0.}; //position
                
                // make sure we aren't going to include the current Body ie it doesn't interact with itself
                if (i != j){
                    
                    //System.out.println(String.format("Position: {%g, %g, %g}", X[0], X[1], X[2]));
                    
                    double mass = objects[i].getMass(); //get the other Body's mass
                    
                    // calculate the distance vector from the current Body to the other Body
                    double[] distVec = b.vectorToBody(objects[i]);
                    
                    // calculate the magnitude of this distance
                    double d = magnitude(distVec);
                    
                    // create a factor, from the gravitational force equation
                    // --> F(vector) = m_1*a = ((-G*m_1*m_2)/(|r|^3)) * r(vector)
                    // --> factor = -G*m_1*m_2)/(|r|^3) --> a(vector) = factor * r(vector)
                    double factor = (-G*mass)/(d*d*d); //acceleration factor
                    
                    // get the other Body's position and velocity

                    
                    //Uncomment Block below for Start position value
                    /*
                    System.out.println(String.format("Start: %s", b.name));
                    System.out.println(String.format("Position: {%g, %g, %g}", objPos[0], objPos[1], objPos[2]));
                    */
                    
                    
                    // add the acceleration to its respective components
                    A[0] += factor * distVec[0];
                    A[1] += factor * distVec[1];
                    A[2] += factor * distVec[2];
                    
                }
                    double[] objVel = b.getVelocity();
                    double[] objPos = b.getPosition(); 
                    
                    // Euler integration at step n:
                    // dr(vector)/dt|_(n+1) = dr(vector)/dt|_n + (d^2 r(vector)/dt^2) * dt
                    // r(vector)|_(n+1) = r(vector)|_n + (dr(vector)/dt)|_(n+1) * dt
                    // or --> V_new = V_current + A * dt, X_new = X_current + V_new *dt
                    
                    //put this part outside of the loop
                    // calculate new velocity components
                    V[0] += objVel[0] + A[0]*timestep;
                    V[1] += objVel[1] + A[1]*timestep;
                    V[2] += objVel[2] + A[2]*timestep;
                    
                    // calculate new position components
                    X[0] += objPos[0] + V[0]*timestep;
                    X[1] += objPos[1] + V[1]*timestep;
                    X[2] += objPos[2] + V[2]*timestep;
                    
                    // update the object's position
                    b.setPosition(X[0], X[1], X[2]);
                    
                    // update the object's velocity
                    b.setVelocity(V[0], V[1], V[2]);
                    
                    //Uncomment Block below for end position values
                    /*
                    System.out.println(String.format("End: %s", b.name));
                    System.out.println(String.format("Position: {%g, %g, %g}", X[0], X[1], X[2]));
                    */
                
            }
            // increment the counter, keeping in sync with the current Body
            j++;
        }
    }
    
    /*for every body:
        loop over every OTHER body
            add the foces acting from all other bodies / the accelelerations due to the other bodies
            enact the foce on the body (a each time interval)
            find the total acceleration and directeion acing on every body
            create a list of the accelertion for each body (x,y,z)
    update the velocities and then positions
*/
            
    
    // method for the implementation of the 4th Order Runge-Kutt integration 
    // method, taking in some timestep as dt
    public void stepRK4(double timestep){
        // TODO: ADD 4th Order Runge-Kutta method for propogation 
    }
    
    
}
