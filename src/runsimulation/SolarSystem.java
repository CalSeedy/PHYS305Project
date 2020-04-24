
package runsimulation;

import java.util.Random;

// SolarSystem class manages all the objects that are inside the system,
// handles the addition/removal of any objects and has the responsibility
// of stepping all objects forward by 1 timestep with either Euler or RK4.
public class SolarSystem {
    
    final static double G = 6.67e-11; // Create a const. for G in SI [m]^3 [kg]^-1 [s]^-2 
    
    private Body[] objects;         // store a array of all bodies in the system
    private int astNum = 0;         // store number of asteroids in the system
    public double time = 0.; 
    public Hit Hits;
    
    // blank constructor to initialise the system with a star ("Sun")
    // Sun has the same properties (im SI) as the Sun, from Google
    public SolarSystem() {
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
    
    public void genHitArray(){
        Hit hit = new Hit(this);
        Hits = hit;
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
    public double magnitude(double[] vector){
        // R = {x, y, z}, |R| = sqrt(x^2 + y^2 + z^2), |R|^2 = R.R
        // --> |R| = sqrt(R.R)
        return Math.sqrt(dot(vector, vector));
    }
    
    // private method to return the dot product of two 3D vectors
    public double dot(double[] vector1, double[] vector2){
        // A.B = (A_x * B_x) + (A_y * B_y) + (A_z * B_z)
        return (vector1[0]*vector2[0] + vector1[1]*vector2[1] + vector1[2]*vector2[2]);
    }
    
    public double[] cross(double[] vec1, double[] vec2){
        double[] out = new double[3];
        if (vec1.length == vec2.length && vec1.length == 3){
            out[0] = (vec1[1]*vec2[2])-(vec1[2]*vec2[1]);
            out[1] = -(vec1[0]*vec2[2])-(vec1[2]*vec2[0]);
            out[2] = (vec1[0]*vec2[1])-(vec1[1]*vec2[0]);
            
        }
        return out;
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
        
        time += timestep;
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
    
    private double[] fvector(double[][] fullstate){
        
        double[][] accelerations = new double[objects.length][3]; //initialising an array of accelerations for the bodies
        
        double[] fvec = new double[objects.length*8];
        for (int i = 0; i < objects.length; i++){
            fvec[i*8 + 0] = fullstate[i][3];
            fvec[i*8 + 1] = fullstate[i][4];
            fvec[i*8 + 2] = fullstate[i][5];
            
            for (int j = 0; j < objects.length; j++){
                if (i != j) {
                    if ((objects[i].isAsteroid) && (objects[j].isAsteroid)){
                        accelerations[i][0] += 0.; //the acceleration on m1 due to m2 in the x-direction is added on in every loop over j
                        accelerations[i][1] += 0.;
                        accelerations[i][2] += 0.;
                    } else {
                        //double m1 = fullstate[i][6];
                        double m2 = fullstate[j][6];

                        double[] r21 = {fullstate[i][0] - fullstate[j][0], fullstate[i][1] - fullstate[j][1], fullstate[i][2] - fullstate[j][2]};
                        double d = Math.sqrt(r21[0]*r21[0] + r21[1]*r21[1] + r21[2]*r21[2]);

                        double factor = -G*m2 / (d*d*d);

                        // a12 = -G*m2 / (|r21|)^3 * r21

                        accelerations[i][0] += factor * r21[0]; //the acceleration on m1 due to m2 in the x-direction is added on in every loop over j
                        accelerations[i][1] += factor * r21[1];
                        accelerations[i][2] += factor * r21[2];

                    }
                }
            }
        }
        
        for (int i = 0; i < accelerations.length; i++){
            fvec[i*8 + 3] = accelerations[i][0];
            fvec[i*8 + 4] = accelerations[i][1];
            fvec[i*8 + 5] = accelerations[i][2];
        }
        return fvec;
    }
    
    private double[] calculateK(double[][] state, double dt){
        double[] fvec = fvector(state);
        double[] k = new double[fvec.length];
        int count = 0;
        for (double val : fvec){
            k[count] = val*dt/2;
            count++;
        }
        
        return k;
        
    }
    
    private double[][] nextState(double[][] initial, double[][] other, double[] k, int step){
        double multiple;
        switch(step){
            
            case 0:
                multiple = 0.0;
            case 1:
                multiple = 0.5;
            case 2:
                multiple = 0.5;
            case 3:
                multiple = 1.0;
                
            default:
                multiple = 0.0;
                
        }
        
        int count = 0;
        double[][] output = new double[initial.length][initial[0].length];
        for (double[] state : other){
            for (int i = 0; i < state.length; i++){
                if (i != 6 || i != 7){
                    output[count][i] = state[i] + multiple*k[count*8 + i]; 
                } else {
                    output[count][i] = initial[count][i];
                }
            }
            count++;
        }
        
        return output;
    }
    
    
    // method for the implementation of the 4th Order Runge-Kutt integration 
    // method, taking in some timestep as dt
    public void stepRK4(double timestep){
        // TODO: ADD 4th Order Runge-Kutta method for propogation
        double[][] state_initial = new double[objects.length][8];
        //double[][] state_second, state_third, state_fourth, state_final;
        //state_second = state_third = state_fourth = state_final = state_initial;
        int count = 0;
        for (Body obj : objects){
            double[] state = obj.getState();
            state_initial[count] = state;
            count++;
        }
        
        
        // calculating k1 and second state
        double[] k1 = calculateK(state_initial, timestep);
        double[][] state_second = nextState(state_initial, state_initial, k1, 1);

        
        // calculating k2 and the 3rd state
        double[] k2 = calculateK(state_second, timestep);
        double[][] state_third = nextState(state_initial, state_second, k2, 2);
        
        
        // calculating k3 and 4th state
        double[] k3 = calculateK(state_third, timestep);
        double[][] state_fourth = nextState(state_initial, state_third, k2, 3);
        
        
        // calculating k4 and final state
        double[] k4 = calculateK(state_fourth, timestep);
        
        double[][] state_final = new double[state_initial.length][state_initial[0].length];
        count = 0;
        for (double[] state : state_initial){
            for (int i = 0; i < state.length; i++){
                if (i != 6 || i != 7){
                    state_final[count][i] = state[i] + (1./6.)*k1[count * state.length + i] + (1./3.)*k2[count*state.length + i] + (1./3.)*k3[count*state.length + i] + (1./6.)*k4[count*state.length + i]; 
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
        
        
        time += timestep;
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
            if ((distance > furthest) && !(b.isAsteroid)){
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
        
        // set the asteroid's radius such that it is a Gaussian about (1000 +/- 100)m
        double a_radius = Math.abs(rand.nextGaussian()*100. + 1000.);
        
        // use Density * Volume to calculate mass; assume Density of asteroids = 5000 kg/m^3
        // Volume of a sphere = 4/3 * PI * r^3
        double a_mass = 5000 * (4./3.) * Math.PI * a_radius * a_radius * a_radius;
        
        // create the asteroid with the above calculated values with the name AsteroidX, where X
        // is the current number of asteroids in the system
        Body ast = new Body(a_pos, a_vel, a_mass, a_radius, String.format("Asteroid%d", astNum));
        
        // add the asteroid to the system and increment the number of asteroids
        addObject(ast);
        astNum++; 
    }
    
    // method that creates a new instance of an asteroid (which is just a Body)
    // the asteroid is generated with a random position on a unit sphere and always
    // has a velocity towards the Earth
    public void generateAsteroidLine(double x_in, double y_in, double z_in, int thickness, int N, boolean target){
        /*
        // Find the object that is the furthest away (that isnt an asteroid)
        double furthest = 0.;
        for (Body b : objects){
            double[] pos = b.getPosition();
            double distance = magnitude(pos);
            if ((distance > furthest) && !(b.name.contains("Asteroid"))){
                furthest = distance;
            }
        }
        */
        Random rand = new Random();
        for (int i = 0; i < N; i++){
            double x, y, z;
            x = x_in - rand.nextDouble()*(thickness*148.28e7)*0.99;
            y = y_in - rand.nextDouble()*(2.*7.78574E11);
            z = z_in; // can do same as above but with some "width/height"
            
            double[] a_pos = {x, y, z};
            double[] a_vel = {0., 0., 0.};
            double vx, vy, vz;

            if (target){
                // Find out where the Earth is
                int ind = findObjectIndex("Earth");
                Body E = getObject(ind);
                double[] E_pos = E.getPosition();
                // create a velocity such that the asteroid heads towards Earth
                vx = (a_pos[0] - E_pos[0])/ 1e8;
                vy = -(a_pos[1] - E_pos[1])/ 1e8; //+ (0.5*(x - E_pos[1]))/7.78574E11)*1e3 *(rand.nextDouble());
                vz = 0.;//(a_pos[2] - E_pos[2])/ 1e7;
                a_vel[0] = vx;
                a_vel[1] = vy;
                a_vel[2] = vz;
            } else {
                vx = -1.7e4;
                a_vel[0] = vx;
            }
            
            // set the asteroid's radius such that it is a Gaussian about (50 +/- 33)m
            double a_radius = Math.abs(rand.nextGaussian()*100. + 1000.);
            
            // use Density * Volume to calculate mass; assume Density of asteroids = 5000 kg/m^3
            // Volume of a sphere = 4/3 * PI * r^3
            double a_mass = 5e3 *(4./3.) * Math.PI * a_radius * a_radius * a_radius;
            Body ast = new Body(a_pos, a_vel, a_mass, a_radius, String.format("Asteroid%d", astNum));
            
            // add the asteroid to the system and increment the number of asteroids
            addObject(ast);
            astNum++; 
            
        }
    }
    
    
    public void generateAsteroidCircle(double x_in, double y_in, double z_in, double radius, int N, boolean target){
        Random rand = new Random();
        for (int i = 0; i < N; i++){
            double x, y, z, dx, dy, dz;
            double theta = rand.nextDouble() * Math.PI * 2;
            double phi = rand.nextDouble() * Math.PI;
            double r = rand.nextDouble()*(N * 148.28e7 / (Math.PI * 2)) + radius ;
            
            x = x_in + r*Math.cos(theta);//*Math.sin(phi);
            y = y_in + r*Math.sin(theta);//*Math.sin(phi);
            z = z_in;// + r*Math.cos(phi);
            
            double[] a_pos = {x, y, z};
            double[] a_vel = {0., 0., 0.};
            double vx, vy, vz;

            if (target){
                // Find out where the Earth is
                int ind = findObjectIndex("Earth");
                Body E = getObject(ind);
                double[] E_pos = E.getPosition();
                // create a velocity such that the asteroid heads towards Earth
                vx = (a_pos[0] - E_pos[0])/ 1e8;
                vy = -(a_pos[1] - E_pos[1])/ 1e8; //+ (0.5*(x - E_pos[1]))/7.78574E11)*1e3 *(rand.nextDouble());
                vz = 0.;//(a_pos[2] - E_pos[2])/ 1e7;
                a_vel[0] = vx/10.;
                a_vel[1] = vy/10.;
                a_vel[2] = vz/10.;
            } else {
                vx = -x / 5e7;
                vy = -y / 5e7;
                vz = -z / 5e7;
                a_vel[0] = vx;
                a_vel[1] = vy;
                a_vel[2] = vz;
                
            }
            
            // set the asteroid's radius such that it is a Gaussian about (50 +/- 33)m
            double a_radius = Math.abs(rand.nextGaussian()*100. + 1000.);
            
            // use Density * Volume to calculate mass; assume Density of asteroids = 5000 kg/m^3
            // Volume of a sphere = 4/3 * PI * r^3
            double a_mass = 5e3 *(4./3.) * Math.PI * a_radius * a_radius * a_radius;
            Body ast = new Body(a_pos, a_vel, a_mass, a_radius, String.format("Asteroid%d", astNum));
            
            // add the asteroid to the system and increment the number of asteroids
            addObject(ast);
            astNum++; 
            
        }
    }
    
    
    
    public void cleanAsteroids(){
        for(Body b : objects){
            if (b.isAsteroid){
                double[] pos = b.getPosition();
                double d = magnitude(pos);
                if (d >= 59.0638E+11 *2.5){ // 2.5 * pluto orbit radius
                    removeObject(findObjectIndex(b.name));
                }
            }
        }
    }
    
    
    public Hit getHits(){
        return Hits;
    }
}
