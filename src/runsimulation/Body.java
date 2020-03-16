package runsimulation;

public class Body {
    private double[] position;              // x,y,z position [m]
    private double[] velocity;              // velocity in x, y, z [m][s]^-1
    private double mass;                    // object's masss [kg]
    private double radius;                  // object's radius [m]
    public String name;                     // name to refer to object
    
    // empty contructor
    public Body() {
        // initialise the class variables
        // easy to handle errors by setting mass as a -ve value
        mass = -1.;     
        radius = 0.;
        double[] pos = {0., 0., 0.};
        double[] vel = pos;
        position = pos;
        velocity = vel;
        name = "UnspecifiedObject";
    }
    
    // constructor that takes in the intial position, velocity, mass and object radius in SI
    public Body(double[] pos_in, double[] vel_in, double mass_in, double radius_in) {
        // use the arguments to intialise the class variables
        mass = mass_in;     
        radius = radius_in;
        position = pos_in;
        velocity = vel_in;
        name = "UnspecifiedObject";
    }
    
    // same as Body(double[], double[], double, double), 
    // but the object can be asigned a chosen name.
    public Body(double[] pos_in, double[] vel_in, double mass_in, double radius_in, String name_in) {
        mass = mass_in;     
        radius = radius_in;
        position = pos_in;
        velocity = vel_in;
        name = name_in;
    }
    
    
    // method to return (x, y, z) vector from the object to the Body b
    public double[] vectorToBody(Body b){
        double[] vec = {position[0] - b.position[0], position[1] - b.position[1], position[2] - b.position[2]};
        return vec;
    }
    
    // method to calculate the absolute distance to the Body b
    public double distanceToBody(Body b){
        // get the displacement vector, vec, to the Body b
        double[] vec = vectorToBody(b);
        // return the magnitude of the displacement vector
        return Math.sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    }
    
    // getter mehtod to return all values in a 1D array, execpt for the name (arbitrarily) 
    public double[] getState(){
        double[] state = {position[0], position[1], position[2], velocity[0], velocity[1], velocity[2], mass, radius};
        return state;
    }
    
    
    // setter mehtod to overwrite all values in the Body, execpt for the name (arbitrarily) 
    public void setState(double[] state){
        position[0] = state[0];
        position[1] = state[1];
        position[2] = state[2];
        velocity[0] = state[3];
        velocity[1] = state[4];
        velocity[2] = state[5];
        mass = state[6];
        radius = state[7];        
    }
    
    // getter for Body position
    public double[] getPosition(){
        return position;
    }
    
    // setter for position, takes in compoents as arguments/ parameters
    public void setPosition(double x, double y, double z){
        //System.out.println(String.format("Position: {%g, %g, %g}", x, y, z));
        position[0] = x;       
        position[1] = y;
        position[2] = z;
    }
    
    // getter for Body velocity
    public double[] getVelocity(){
        return velocity;
    }
    
    // setter for velocity, takes in compoents as arguments/ parameters
    public void setVelocity(double vx, double vy, double vz){
        velocity[0] = vx;       
        velocity[1] = vy;
        velocity[2] = vz;
    }
    
    // getter for the Body's name
    public String getName(){
        return name;
    }
    
    // setter for Body's name
    public void setName(String name_in){
        name = name_in;
    }
    
    // getter for the Body's mass
    public double getMass(){
        return mass;
    }
    
    public double getRadius(){ //getter for the Body's radius
        return radius;
    }
    
    public double [] getMomentum(){  //getter for the Body's momentum
        double[] p = {0.,0.,0.,};
        double [] velocity = getVelocity();
        double mass = getMass(); //do I have to call the method to use 'mass'?
        p[0] = mass * velocity[0];
        p[1] = mass * velocity[1];
        p[2] = mass * velocity[2];
        return p;
    }
    
    public void updateMass(double m){ //add a mass to th body's mass. this must be called before updateVelocity when two bodies collide
        mass += m;
    }
    
    //update the velocity after the momentum has been changed
    public void updateVelocity(double [] momentum){ 
        double vx = momentum[0] / getMass();
        double vy = momentum[1] / getMass();
        double vz = momentum[2] / getMass();
        setVelocity(vx, vy, vz);
    }
}
