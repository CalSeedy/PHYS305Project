
package runsimulation;


class Hit {

    int [] hits;   //will store the number of hits on each object
    int [] nearMisses;   //will store the number of nearMisses on each object
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
        nearMisses = new int [names.length]; //objects.length
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
        int[] misses = new int[planetIndexes.length];
        int k = 0;
        for (Body b : objs){
            if (!b.isAsteroid){
                names[k] = b.name;
                idxs[k] = sys.findObjectIndex(b.name);
                misses[k] = nearMisses[idxs[k]];
                k++;
            }
        }
        
        planetNames = names;
        planetIndexes = idxs;
        nearMisses = misses;
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
                    double d = Math.abs(b1.distanceToBody(b2));
                    double radii = b1.getRadius() + b2.getRadius();

                    if (d <= radii){
                        //System.out.println(String.format("%e: %s collided with %s", sys.time, b1.name, b2.name));
                        if (b1.isAsteroid && !b2.isAsteroid){
                            //System.out.println(String.format("\t %s was removed!", b1.name));
                            
                            int idx = sys.findObjectIndex(b2.name); 
                            hits[idx] ++;
                            
                            sys.removeObject(sys.findObjectIndex(b1.name));
                        } else if (!b1.isAsteroid && b2.isAsteroid){
                            //System.out.println(String.format("\t %s was removed!", b2.name));
                            
                            int idx = sys.findObjectIndex(b1.name);
                            hits[idx] ++;
                            sys.removeObject(sys.findObjectIndex(b2.name));
                        } else {
                            // do stuff for planet-planet collisions
                        }
                    } else if (d <= 1.5*radii){
                        if (b1.isAsteroid && !b2.isAsteroid){
                            //System.out.println(String.format("%e: %s barely passed by %s", sys.time, b1.name, b2.name));
                            int idx = sys.findObjectIndex(b2.name); 
                            nearMisses[idx] ++;
                        } else if (!b1.isAsteroid && b2.isAsteroid){
                            //System.out.println(String.format("%e: %s barely passed by %s", sys.time, b2.name, b1.name));
                            int idx = sys.findObjectIndex(b1.name);
                            nearMisses[idx] ++;
                            
                        } else {
                            // do stuff for planet-planet collisions
                        }
                    }
                }
            }       
        }
    }
    
    public int[] getHits() {
        return hits;
    }
    
    public int[] getMisses() {
        return nearMisses;
    }
    
    public String[] getNames() {
        return planetNames;
    }
    
    
    public void display(){
        System.out.println("\n\n");
        for (int i = 0; i < hits.length; i++){
            System.out.println(String.format("%s :  %d (%d)", planetNames[i], hits[i], nearMisses[i]));
        }
    }
}        
    

