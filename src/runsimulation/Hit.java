
package runsimulation;


class Hit {

    int [] hits;   //will store the number of hits on each object
    String[] planetNames;
    
    int[] planetIndexes;
    
    public Hit(SolarSystem sys){
        Body[] objs = sys.getObjects();
        int c = 0;
        for (Body b : objs){
            if (b.isAsteroid){
                c++;
            }
        }
        
        String[] names = new String[objs.length - c];
        int[] idxs = new int[objs.length - c];
        //maybe create an aray of planet names
        hits = new int [names.length]; //objects.length
        //for (int i : hits){
        //    i = 0;
        //}
        
        planetNames = names;
        planetIndexes = idxs;
        
        index(sys);
        
    }
    
    private void index(SolarSystem sys){
        Body[] objs = sys.getObjects();
        String[] names = new String[planetNames.length];
        int[] idxs = new int[planetIndexes.length];
        int k = 0;
        for (Body b : objs){
            if (!b.isAsteroid){
                names[k] = b.name;
                idxs[k] = sys.findObjectIndex(b.name);
            }
        }
    }
    
    private int getIndex(Body[] arr, String identifier){
        // create a counter to keep track of current position in objects
        int c = 0;
        
        // while we havent reached the end of the array of Bodies
        while (c < arr.length){
            // for each body, b, in objects (one at a time)
            for (Body b : arr){
                
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
        return -1;
    }
    
    
    public void checkHit(SolarSystem sys){
        Body[] objects = sys.getObjects();
        
        for (int j = 0; j < objects.length; j++){
            for (int i = j; i < objects.length; i++){
                Body b1 = objects[i];
                Body b2 = objects[j];
                if (i != j){
                    double d = b1.distanceToBody(b2);
                    double radii = b1.getRadius() + b2.getRadius();

                    if (d <= radii){
                        //System.out.println(String.format("%e: %s collided with %s", sys.time, b1.name, b2.name));
                        if (b1.isAsteroid && !b2.isAsteroid){
                            //System.out.println(String.format("\t %s was removed!", b1.name));
                            sys.removeObject(sys.findObjectIndex(b1.name));
                            int idx = sys.findObjectIndex(b2.name); 
                            hits[idx] ++;
                        } else if (!b1.isAsteroid && b2.isAsteroid){
                            //System.out.println(String.format("\t %s was removed!", b2.name));
                            sys.removeObject(sys.findObjectIndex(b2.name));
                            int idx = sys.findObjectIndex(b1.name);
                            hits[idx] ++;
                            
                        } else {
                            // do stuff for planet-planet collisions
                        }
                    } else if (d <= 1.5*radii){
                        if (b1.isAsteroid){
                            //System.out.println(String.format("%e: %s barely passed by %s", sys.time, b1.name, b2.name));
                        } else {
                            //System.out.println(String.format("%e: %s barely passed by %s", sys.time, b2.name, b1.name));
                        }
                    }
                }
            }       
        }
        System.out.println("\n\n");
        for (int i = 0; i < hits.length; i++){
            
            System.out.println(String.format("%s :  %d", planetNames[i], hits[i]));
        }
    }
    
    /*
    public void Hit(){
        
        
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
                        
                        //Scalett is working on the below
                        ///*if (name1.indexOf("Asteroid") == -1){ //if the string "Asteroid" is not in the name
                            //if name1 is in the list 'names':
                            if (names.contains(name1)){ //google isn't helping me!!! I want to see whether this list contains name1
                                //find the index
                                //add 1 to the same index in the hit list
                            }
                            //else:
                                //add the name to the list 'names'
                                
                                //find the index of that name
                                //add 1 to the same index in the hit list
                                //OR add the same to the end of the 'names' list and add 1 to the end of the 'hits' list
                                names.add(name1); //add the name to the end of the 'names' list
                                hits.add(1); //add one to the end of the 'hits' list
                        }
                        */
                        /*
                        removeObject(i);

                        //hits[j] ++;
                        //hits[i] ++;
                    
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
    */
}            
    

