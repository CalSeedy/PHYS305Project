package runsimulation;

import java.util.Random;

public class RunSimulation {
    public static void main(String[] args) {
        int n = 10000; // number of steps
        double timestep = (86400.)/2.; //1/2 a day is s
        
        SolarSystem system = new SolarSystem(); //initialising the solar system. Ceating a box
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
             
        Body Earth = new Body(E_pos, E_vel, 5.9722e24, 6371e3, "Earth"); //star position and velocity, mass and object radius, "name"
        system.addObject(Earth); //adding eath to the solar system. Creaes a sun in the middle too
        
        Body Jupiter = new Body(J_pos, J_vel, 1.8976E27, 69911, "Jupiter");
        system.addObject(Jupiter); //adding Jupiter to the solar system       
        
        Body Fatty = new Body(F_pos, F_vel, 2E29, 69911, "Fatty");
        system.addObject(Fatty); //adding Fatty to the solar system
        
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
            system.stepEuler(timestep); //step using Euler
            
            //put NEW for loop here
            //for (int j = 0; j < system.getObjects().length-1; j++)){
            //double[] out = system.getObject(j).getPosition(); //for every j value, find the corresponding object
            //then store the data in an appropriate array
            //change Data.java class to put all orbits in one array
            //for the time being, just store them individually and save o individual csvs
            //}
        
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
        storeEarthPos.writeToCSV("Earth_data.csv");
        storeJupiterPos.writeToCSV("Jupiter_data.csv");
        storeFattyPos.writeToCSV("Fatty_data.csv");
    }
}
