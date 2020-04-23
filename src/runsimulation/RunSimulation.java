package runsimulation;

import java.util.Random;

public class RunSimulation {
    final static int ITERATIONS = 10;
    final static int STEPS = 5000;
    
    final static double G = 6.67e-11;
    
    public static void main(String[] args) {
        boolean elliptical = false;
        double timestep = 1*(24*60*60);
        String[] names = {"Sun","Mercury","Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"};//"Fatty_data.csv"};
        
        int[] totalHits = new int[names.length];
        for (int item : totalHits){
            item = 0;
        }
        
        int[] totalMisses = new int[names.length];
        for (int item : totalMisses){
            item = 0;
        }
        
        for (int k = 0; k < ITERATIONS; k++){
            SolarSystem sys = new SolarSystem();
            Hit hits;
            double[] Mercury_pos = {5.79E+10, 0., 0.}; //Mercury semi major axis. Not updated values yet!
            //double[] Mercury_pos = {2.07E11, 0., 0.}; //Perihelion
            double[] Mercury_vel = {0., 4.74E+04, 0.};
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
            Earth.setEccentricity(0.0167086);
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
            //double[] S_vel = {0., 10180., 0.}; //Max
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

            Body[] objs = sys.getObjects();
            int c = 0;
            for (Body o : objs){
                if (!(o.name.contains("Asteroid") || o.name.equals("Sun"))){
                    double e = o.getEccentricity();
                    double ap = (1.+e)/(1.-e) * o.getPosition()[0];
                    o.setAphelion(ap);
                    c++;
                }
            }


            Random rand = new Random();
            if (elliptical){ // assume position is semi-major axis
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
            }

            int thickness = 200;
            double[] pJ = Pluto.getPosition();
            double rJ = Math.sqrt(pJ[0]*pJ[0] + pJ[1]*pJ[1] + pJ[2]*pJ[2]);
            double[] astLine_pos = {0.,0.,0.};//{2.*rJ + thickness, rJ, 0.}
            sys.generateAsteroidCircle(0.,0.,0., rJ*2, 1000, false);
            sys.genHitArray();
            //Data storeSystem = new Data(n, timestep, names);
            
            int a = 0;
            System.out.println(String.format("Iteration: %d",k));
            while (a < STEPS){
                sys.Hits.checkHit(sys);
                sys.stepRK4(timestep);
                sys.cleanAsteroids();
                if (a % 1000 == 0){
                    System.out.println(String.format("\tStep: %d", a));
                }
                a++;
            }
            //sys.Hits.display();
            
            int[] currentHits = sys.Hits.getHits();
            int[] currentMisses = sys.Hits.getMisses();
            String[] currentNames = sys.Hits.getNames();

            for (int i = 0; i< currentHits.length; i++){
                if (currentNames[i].equals(names[i])){
                    totalHits[i] += currentHits[i];
                    totalMisses[i] += currentMisses[i];    
                } else {
                    String s1 = names[i];
                    for (int j = 0; j < currentNames.length; j++){
                        String s2 = currentNames[j];
                        if (s1.equals(s2)){
                            totalHits[i] += currentHits[j];
                            totalMisses[i] += currentMisses[j];
                        }    
                    }
                    
                    
                }
            }
            /*
            for (int i = 0; i < totalHits.length; i++){
                System.out.println(String.format("%s :  %d (%d)", names[i], totalHits[i], totalMisses[i]));
            }
            */
        }
    }
}
    
    
    
